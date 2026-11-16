#!/usr/bin/env python3
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""GSLC map-coordinate adaptation of Xin Wang's run7 LA/LB ionosphere workflow.

Gomba et al. (2016) frequency combination; original workflow by Yajun Zhang /
Xin Wang. See UPSTREAM_PROVENANCE.md and LICENSE (GPL-3.0).
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import numpy as np
from pyproj import CRS
from scipy.interpolate import RegularGridInterpolator

import make_intf_nsr_gslc as core
from gslc_config import read_iono_config
from nearest_grid_c_exact import reproduce_c_nearest
from plot_gslc_gmt import plot_grid, plot_dependencies, read_plot_grid

ROOT = Path(__file__).resolve().parent
UPSTREAM = "https://github.com/Xin-Wang-520/gmtsar-nisar-dinsar-workflow"


def sample_grid(data, x, y, target_x, target_y, method="linear", edge_pixels=0):
    """Sample true map coordinates. Never relabel or stretch grid extents.

    edge_pixels=None allows nearest edge extrapolation for reconstructing a
    coarse filter grid inside its original extent. Finite margins are used
    when transferring the final B correction back to A.
    """
    dx, dy = x[1]-x[0], y[0]-y[1]
    interp = RegularGridInterpolator((y[::-1], x), np.asarray(data)[::-1],
                                     method=method, bounds_error=False, fill_value=np.nan)
    result = np.full((len(target_y), len(target_x)), np.nan, dtype="f4")
    margin_x = np.inf if edge_pixels is None else edge_pixels*dx
    margin_y = np.inf if edge_pixels is None else edge_pixels*dy
    for row in range(0, len(target_y), 256):
        yy, xx = np.meshgrid(target_y[row:row+256], target_x, indexing="ij")
        allowed = ((xx >= x[0]-margin_x) & (xx <= x[-1]+margin_x) &
                   (yy >= y[-1]-margin_y) & (yy <= y[0]+margin_y))
        pts = np.column_stack([np.clip(yy.ravel(), y[-1], y[0]), np.clip(xx.ravel(), x[0], x[-1])])
        values = interp(pts).reshape(xx.shape)
        result[row:row+len(yy)] = np.where(allowed, values, np.nan)
    return result


def nearest_fill(data, radius=0):
    # The upstream C-exact tie order expects a south-to-north NetCDF array.
    # Core GSLC arrays are north-to-south; normalize before calling the helper.
    array = np.where(np.isfinite(data), data, np.nan).astype("f4")
    filled, _, _ = reproduce_c_nearest(array[::-1].copy(), radius, 512)
    return filled[::-1].copy()


def combine_ionosphere(phi_a, phi_b, f_a, f_b):
    if not (np.isfinite(f_a) and np.isfinite(f_b) and 0 < f_a < f_b):
        raise ValueError("Expected finite frequencies f_B > f_A > 0")
    r = f_a/f_b
    # Algebraically identical to fb*(fb*phi_A-fa*phi_B)/(fb**2-fa**2).
    return (np.asarray(phi_a, dtype="f8")-r*np.asarray(phi_b, dtype="f8"))/(1-r*r)


def grid_bounds(x, y):
    dx, dy = (x[1]-x[0])/2, (y[0]-y[1])/2
    return [float(x[0]-dx), float(x[-1]+dx), float(y[-1]-dy), float(y[0]+dy)]


def read_checked(path, x=None, y=None):
    values, xx, yy = core.read_grid(path)
    if x is not None and (not np.array_equal(x, xx) or not np.array_equal(y, yy)):
        raise ValueError(f"Grid does not match its band's coordinates: {path}")
    return values, xx, yy


def estimate(la_dir, lb_dir, output, config, *, gmt=None, pdf=True, lonlat=True,
             mask_polygon=None, mask_mode='exclude'):
    la_dir, lb_dir, out = Path(la_dir).resolve(), Path(lb_dir).resolve(), Path(output).resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"Output must be a new or empty directory: {out}")
    if lonlat:
        from reproject_gslc import reproject_dependencies
        reproject_dependencies(pdf=pdf, gmt=gmt)
    reports = [json.loads((d/"run.json").read_text()) for d in (la_dir, lb_dir)]
    ra, rb = reports
    for band, r in zip("AB", reports):
        if r['status'] != 'complete' or r['arguments']['frequency'] != band or 'unwrapping' not in r:
            raise ValueError(f"Band {band} requires a complete unwrapped result")
    for key in ('reference_product', 'secondary_product'):
        if ra[key]['zeroDopplerStartTime'] != rb[key]['zeroDopplerStartTime']:
            raise ValueError("A/B are not the same acquisition pair")
    if ra['arguments']['polarization'] != rb['arguments']['polarization']:
        raise ValueError("A/B polarizations differ")
    if ra['phase_convention'] != rb['phase_convention'] or ra['phase_convention'] != 'reference * conjugate(secondary)':
        raise ValueError("A/B phase convention is incompatible")
    crs = CRS.from_epsg(ra['epsg'])
    if ra['epsg'] != rb['epsg'] or not crs.is_projected or any(a.unit_name != 'metre' for a in crs.axis_info):
        raise ValueError("Ionosphere filters require a common projected CRS with meter axes")
    a, xa, ya = read_checked(la_dir/'unwrap.grd')
    b, xb, yb = read_checked(lb_dir/'unwrap.grd')
    f_a, f_b = [r['reference_product']['center_frequency_hz'] for r in reports]
    combine_ionosphere(0., 0., f_a, f_b)  # Frequency validation before writing.
    ca, _, _ = read_checked(la_dir/'corr.grd', xa, ya)
    cb, _, _ = read_checked(lb_dir/'corr.grd', xb, yb)
    cca, _, _ = read_checked(la_dir/'conncomp.grd', xa, ya)
    ccb, _, _ = read_checked(lb_dir/'conncomp.grd', xb, yb)
    a = np.where(cca > 0, a, np.nan)
    b = np.where(ccb > 0, b, np.nan)
    spatial_mask = None
    roi_a, roi_b = None, None
    if mask_polygon is not None:
        from preview_iono_mask import polygon_selection
        vertices = np.asarray(mask_polygon, dtype='f8')
        if vertices.size < 6 or vertices.size % 2 or not np.isfinite(vertices).all():
            raise ValueError('Mask polygon needs at least three finite E N pairs in boundary order')
        vertices = vertices.reshape(-1, 2)
        if np.array_equal(vertices[0], vertices[-1]): vertices = vertices[:-1]
        area = abs(np.sum(vertices[:,0]*np.roll(vertices[:,1],-1)-vertices[:,1]*np.roll(vertices[:,0],-1)))/2
        if len(vertices)<3 or area<1 or len(np.unique(vertices,axis=0))!=len(vertices):
            raise ValueError('Mask polygon is degenerate; order its vertices around the boundary')
        if mask_mode not in ('keep','exclude'):
            raise ValueError('Mask mode must be keep or exclude')
        roi_a, roi_b = polygon_selection(xa, ya, vertices), polygon_selection(xb, yb, vertices)
        if mask_mode == 'exclude': roi_a, roi_b = ~roi_a, ~roi_b
        a, b = np.where(roi_a,a,np.nan), np.where(roi_b,b,np.nan)
        ca, cb = np.where(roi_a,ca,np.nan), np.where(roi_b,cb,np.nan)
        if min(np.isfinite(a).sum(),np.isfinite(b).sum()) < 16:
            raise ValueError('Polygon leaves fewer than 16 valid samples in A or B')
        spatial_mask = dict(kind='polygon',mode=mask_mode,vertices_xy_m=vertices.tolist(),epsg=ra['epsg'],
                            retained_valid_A=int(np.isfinite(a).sum()),retained_valid_B=int(np.isfinite(b).sum()))
    dependencies = plot_dependencies(gmt) if pdf else {'gmt': core.executable('gmt', gmt)}
    out.mkdir(parents=True, exist_ok=True)
    report = {
        'status':'running', 'version':core.VERSION, 'operation':'GSLC A/B ionospheric estimate at frequency A',
        'method_source':UPSTREAM, 'la_directory':str(la_dir), 'lb_directory':str(lb_dir),
        'created_utc':datetime.now(timezone.utc).isoformat(), 'config':config,
        'frequencies_hz':{'A':f_a,'B':f_b}, 'epsg':ra['epsg'],
        'phase_convention':ra['phase_convention'], 'amplification':float(f_b*f_b/(f_b*f_b-f_a*f_a)),
        'phase_reference':'Direct input unwrapped phases, no automatic 2pi/constant datum alignment',
        'grid_transfer':'True map-coordinate interpolation; no extent relabeling',
        'warnings':list(dict.fromkeys(ra.get('warnings',[])+rb.get('warnings',[]))),
        'spatial_mask':spatial_mask,
    }
    report['warnings'].append('A/B integer-cycle errors and spatial filtering differences can enter the ionosphere estimate; validate the correction for each scene.')
    start = time.monotonic()
    def write(name, values, x=xb, y=yb, title='', unit='radian'):
        core.write_grid(out/(name+'.grd'), values, x, y, crs, title or name, unit)
    try:
        if spatial_mask is not None:
            write('roi_mask_A',np.where(np.isfinite(a),1,np.nan),xa,ya,'Polygon selection intersected with valid A unwrap','1')
            write('roi_mask_B',np.where(np.isfinite(b),1,np.nan),xb,yb,'Polygon selection intersected with valid B unwrap','1')
        with tempfile.TemporaryDirectory(prefix='.iono_',dir=out) as temporary:
            work = Path(temporary)
            gmtbin = dependencies['gmt']
            region = '-R'+'/'.join(f'{v:.14g}' for v in grid_bounds(xb,yb))
            # All GMT operations already use UTM/polar meters, -D0 Cartesian.
            def filter_grid(values, kind, wavelength, coarse=False):
                core.write_grid(work/'input.grd',values,xb,yb,crs)
                cmd = [gmtbin,'grdfilter','input.grd','-D0',f'-F{kind}{wavelength:g}/{wavelength:g}',
                       '-Ni','-Gcoarse.grd']
                if coarse:
                    inc = wavelength/config['coarse_samples_per_wavelength']
                    cmd += [f'-I{inc:g}/{inc:g}']
                core.run_tool(cmd,work,out/'iono_gmt.log')
                z,x,y = core.read_grid(work/'coarse.grd')
                if np.count_nonzero(np.isfinite(z)) < 4 or min(z.shape)<2:
                    raise ValueError('Ionosphere coarse filter has too few samples; reduce wavelength or use a larger scene')
                return sample_grid(z,x,y,xb,yb,'nearest',None)
            a_on_b = sample_grid(a,xa,ya,xb,yb,'linear')
            if roi_b is not None:
                a_on_b = np.where(roi_b,a_on_b,np.nan)
            write('unwrap_A_onB_raw',a_on_b)
            if config['prefilter_a'] > 0:
                a_on_b = filter_grid(a_on_b,'g',config['prefilter_a'])
            if roi_b is not None:
                a_on_b = np.where(roi_b,a_on_b,np.nan)
            write('unwrap_A_onB',a_on_b)
            corr_a = sample_grid(ca,xa,ya,xb,yb,'nearest')
            if roi_b is not None:
                corr_a = np.where(roi_b,corr_a,np.nan)
            mean_corr = (corr_a+cb)/2
            write('corr_A_onB',corr_a,unit='1')
            write('corr_mean',mean_corr,unit='1')
            raw = combine_ionosphere(a_on_b,b,f_a,f_b)
            write('iono_correction',raw)
            good = np.isfinite(raw) & np.isfinite(mean_corr) & (mean_corr>=config['coherence_threshold'])
            if np.count_nonzero(good)<16:
                raise ValueError('Fewer than 16 common reliable A/B ionosphere samples')
            median = float(np.median(raw[good]))
            halfwidth = math.pi*f_a*f_b/(f_b*f_b-f_a*f_a)
            good &= abs(raw-median)<=halfwidth
            quality = np.where(good,raw,np.nan).astype('f4')
            write('iono_quality',quality)
            write('iono_quality_mask',np.where(good,1,np.nan),unit='1')
            report.update(quality_pixels=int(good.sum()),quality_fraction=float(good.mean()),
                          outlier_center_rad=median,outlier_halfwidth_rad=halfwidth)
            working = nearest_fill(quality)
            wavelength = config['filter_wavelength']
            for kind in ('m','b','m'):
                core.LOG.info('Ionosphere: %g m %s filter',wavelength,kind)
                smoothed = filter_grid(working,kind,wavelength,True)
                working = nearest_fill(np.where(good,quality,smoothed))
            boxcar = filter_grid(working,'b',wavelength,True)
            filtered_b = filter_grid(boxcar,'g',wavelength,True)
            # Upstream's final Gaussian coarse grid is reconstructed with surface.
            xyz = work/'coarse.xyz'
            with xyz.open('wb') as stream:
                subprocess.run([gmtbin,'grd2xyz',str(work/'coarse.grd'),'-s'],stdout=stream,check=True)
            core.run_tool([gmtbin,'surface',str(xyz),region,
                           f'-I{xb[1]-xb[0]:.14g}/{yb[0]-yb[1]:.14g}',
                           '-r',f'-T{config["surface_tension"]:g}','-Gsurface.grd'],work,out/'iono_gmt.log')
            z,x,y = core.read_grid(work/'surface.grd')
            filtered_b = sample_grid(z,x,y,xb,yb,'linear',None)
            if np.any(~np.isfinite(filtered_b)):
                filtered_b = nearest_fill(filtered_b)
            write('iono_correction_filtered_LB',filtered_b if roi_b is None else np.where(roi_b,filtered_b,np.nan))
            iono_a = sample_grid(filtered_b,xb,yb,xa,ya,'linear',1.5)
            iono_a = np.where(np.isfinite(a),iono_a,np.nan)
            write('iono_correction_filt',iono_a,xa,ya,'Ionospheric correction at frequency A on original A grid')
            corrected = a-iono_a
            write('unwrap_corrected',corrected,xa,ya,'Input A unwrapped phase minus ionospheric estimate')
            wrap = lambda z: np.angle(np.exp(1j*np.asarray(z,dtype='f8')))
            for name,z,x,y in [('phase_LA_before',a,xa,ya),('phase_LB_before',b,xb,yb),
                               ('phase_iono_correction',iono_a,xa,ya),('phase_LA_after',corrected,xa,ya)]:
                write(name,wrap(z),x,y)
            report['edge_extrapolation_note']='Only up to 1.5 B pixel spacings allowed when returning the filtered field to A; unsupported A pixels stay NaN.'
            report['corrected_a_valid_fraction']=float(np.isfinite(corrected).mean())
            report['filter_sequence']=['median','boxcar','median','boxcar','Gaussian','surface']
        wavelength_a = core.C/f_a
        report['los_iono'] = core.phase_to_los(out/'iono_correction_filt.grd', out/'los_iono.grd', wavelength_a,
                                               title='Ionospheric correction at A expressed as LOS')
        report['los_corrected'] = core.phase_to_los(out/'unwrap_corrected.grd', out/'los_corrected.grd', wavelength_a)
        if pdf:
            report['plots']=plot_ionosphere(out,la_dir,lb_dir,dependencies,frequencies=report['frequencies_hz'])
        if lonlat:
            from reproject_gslc import reproject_ionosphere
            report['lonlat'] = reproject_ionosphere(out,la_dir,lb_dir,report,pdf=pdf,gmt=gmt)
        report['status']='complete'
    except BaseException as exc:
        report['status'],report['error']='failed',str(exc)
        raise
    finally:
        report['elapsed_seconds']=time.monotonic()-start
        (out/'iono_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    return report


def plot_ionosphere(out, la, lb, dependencies, *, geographic=False, frequencies=None):
    """Unshaded GMT 2x2 figures with a colorbar below each map panel."""
    out, la, lb = Path(out), Path(la), Path(lb)
    suffix = '_ll' if geographic else ''
    if frequencies is None:
        frequencies = {band: json.loads((directory/'run.json').read_text())['reference_product']['center_frequency_hz']
                       for band, directory in (('A', la), ('B', lb))}
    factors = {band: -core.C/frequency/(4*np.pi)*1000 for band, frequency in frequencies.items()}
    before = ([out/f'unwrap_LA_before{suffix}.grd', out/f'unwrap_LB_before{suffix}.grd'] if geographic
              else [la/'unwrap.grd', lb/'unwrap.grd'])
    unwrapped = before + [out/f'iono_correction_filt{suffix}.grd', out/f'unwrap_corrected{suffix}.grd']
    wrapped = [out/f'{name}{suffix}.grd' for name in
               ('phase_LA_before', 'phase_LB_before', 'phase_iono_correction', 'phase_LA_after')]
    los = before + [out/f'los_iono{suffix}.grd', out/f'los_corrected{suffix}.grd']
    pages = []
    with tempfile.TemporaryDirectory(prefix='.iono_plot_', dir=out) as temp:
        work = Path(temp)
        (work/'user').mkdir()
        env = dict(os.environ, GMT_USERDIR=str(work/'user'), GMT_TMPDIR=str(work))
        gmt, gs = dependencies['gmt'], dependencies['gs']
        for mode, sources in (('phase', unwrapped), ('wrapped', wrapped), ('los', los)):
            stem = {'phase': 'iono_correction_2x2', 'wrapped': 'iono_wrapped_phase_2x2',
                    'los': 'iono_los_2x2'}[mode] + suffix
            titles = (['LA before', 'LB before', 'Ionospheric estimate at A', 'LA after subtraction'] if mode != 'los'
                      else ['LA LOS before', 'LB LOS before', 'Ionospheric correction as A LOS', 'LA LOS after subtraction'])
            multipliers = [factors['A'], factors['B'], 1000, 1000] if mode == 'los' else [1]*4
            first = True
            ps = work/(stem+'.ps')
            with (out/'iono_plot.log').open('a') as log:
                subprocess.run([gmt, 'set', 'GMT_THEME', 'classic', 'GMT_HISTORY', 'false',
                                'COLOR_NAN', '255/255/254', 'MAP_FRAME_TYPE', 'plain',
                                'FORMAT_GEO_MAP', 'ddd.xxF'], cwd=work, env=env,
                               stdout=log, stderr=log, check=True)
                def run(args, text=None, append=True):
                    with ps.open('ab' if append else 'wb') as stream:
                        subprocess.run(list(map(str, args)), cwd=work, env=env,
                                       input=text.encode() if text else None, stdout=stream, stderr=log, check=True)
                _, xa, ya = read_plot_grid(sources[0], geographic)
                bounds = np.array(grid_bounds(xa, ya))/(1 if geographic else 1000)
                region = '-R'+'/'.join(map(str, bounds))
                width = 11.0
                height = width*(bounds[3]-bounds[2])/(bounds[1]-bounds[0])
                if geographic:
                    height /= math.cos(math.radians((bounds[3]+bounds[2])/2))
                projection = f'-JM{width}c' if geographic else f'-JX{width}c/{height}c'
                axes = ['-Bxa', '-Bya'] if geographic else ['-Bxa+lEasting (km)', '-Bya+lNorthing (km)']
                for i, (source, title, multiplier) in enumerate(zip(sources, titles, multipliers)):
                    z, x, y = read_plot_grid(source, geographic)
                    z = z*multiplier
                    finite = z[np.isfinite(z)]
                    if not finite.size:
                        raise ValueError(f'No finite values available for plotting {source}')
                    grid, palette = work/f'panel_{i}.grd', work/f'panel_{i}.cpt'
                    plot_grid(grid, z, x, y, geographic)
                    if mode == 'wrapped':
                        scale = ['-Crainbow', '-T-3.15/3.15/0.05']
                    else:
                        lo, hi = np.quantile(finite, [.02, .98])
                        if hi-lo < 1e-5:
                            lo, hi = float(lo)-1, float(hi)+1
                        scheme = ['-Cpolar'] if mode == 'los' else ['-Cseis', '-I']
                        scale = [*scheme, f'-T{lo:.12g}/{hi:.12g}/{(hi-lo)/100:.12g}']
                    with palette.open('wb') as cp:
                        subprocess.run([gmt, 'makecpt', *scale, '-Z', '-D'], cwd=work, env=env,
                                       stdout=cp, stderr=log, check=True)
                    px, py = 2+(i%2)*14, 4+(1-i//2)*(height+4.5)
                    run([gmt, 'grdimage', grid, region, projection, f'-C{palette}', '-nn+a',
                         *axes, '-BWSen', f'-Xa{px}c', f'-Ya{py}c', '--FONT_TITLE=11p',
                         '--FONT_LABEL=10p', '--FONT_ANNOT_PRIMARY=9p', '--PS_MEDIA=A1',
                         *(['-P'] if first else ['-O']), '-K'], append=not first)
                    first = False
                    label = '-Bxa+lLOS (mm; positive toward satellite)' if mode == 'los' else '-Bxa+lPhase (rad)'
                    run([gmt, 'psscale', f'-C{palette}', f'-Dx{px+width/2}c/{py-1.65}c+w8c/0.25c+jTC+h',
                         label, '-Xa0c', '-Ya0c', '--FONT_LABEL=10p', '--FONT_ANNOT_PRIMARY=9p',
                         '--MAP_LABEL_OFFSET=0.15c', '-O', '-K'])
                    run([gmt, 'pstext', '-R0/30/0/40', '-JX30c/40c', '-Xa0c', '-Ya0c',
                         '-F+f11p,Helvetica-Bold+jBC', '-N', '-O', '-K'],
                        f'{px+width/2} {py+height+0.75} {title}\n')
                header = {'wrapped': 'Wrapped phases', 'phase': 'Unwrapped phases; individual color scales',
                          'los': 'LOS conversion; individual color scales'}[mode]
                coordinates = ' | WGS 84 lon/lat' if geographic else ''
                run([gmt, 'pstext', '-R0/30/0/40', '-JX30c/40c', '-Xa0c', '-Ya0c',
                     '-F+f13p,Helvetica-Bold+jBC', '-N', '-O', '-K'],
                    f'15 {2*height+10.6} NISAR GSLC A/B ionosphere | {header}{coordinates}\n')
                footer = ('LOS = -wavelength/(4 pi) x phase; each band uses its own wavelength. No offset subtraction.' if mode == 'los'
                          else 'Input SNAPHU phases used directly; validate DEM/processing differences and integer cycles.')
                run([gmt, 'pstext', '-R0/30/0/40', '-JX30c/40c', '-Xa0c', '-Ya0c',
                     '-F+f8p,Helvetica+jBC', '-N', '-O'], f'15 0.15 {footer}\n')
                subprocess.run([gmt, 'psconvert', str(ps), '-Tf', '-P', '-A+m0.2c', f'-D{out}',
                                f'-F{stem}', f'-G{gs}'], cwd=work, env=env, stdout=log, stderr=log, check=True)
            pages.append(stem+'.pdf')
    return pages


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('la',type=Path)
    parser.add_argument('lb',type=Path)
    parser.add_argument('-o','--output',type=Path,required=True)
    parser.add_argument('--config',type=Path,default=ROOT/'configs/config_IONO.txt')
    parser.add_argument('--gmt')
    parser.add_argument('--no-pdf',action='store_true')
    parser.add_argument('--mask-polygon',nargs='+',type=float,metavar='COORD')
    parser.add_argument('--mask-mode',choices=('exclude','keep'),default='exclude')
    parser.add_argument('--lonlat',action=argparse.BooleanOptionalAction,default=True,
                        help='Export longitude/latitude GRD and figures; enabled by default')
    args=parser.parse_args()
    import logging
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    try:
        r=estimate(args.la,args.lb,args.output,read_iono_config(args.config),gmt=args.gmt,
                   pdf=not args.no_pdf,lonlat=args.lonlat,
                   mask_polygon=args.mask_polygon,mask_mode=args.mask_mode)
        print(json.dumps({'status':r['status'],'output':str(args.output),'elapsed_seconds':r['elapsed_seconds']}))
    except (ValueError,OSError,RuntimeError,KeyError,ImportError,subprocess.CalledProcessError) as exc:
        parser.exit(1,f'{exc}\n')


if __name__=='__main__':
    main()
