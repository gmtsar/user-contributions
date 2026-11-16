#!/usr/bin/env python3
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""Preview an oblique north-side mask on existing A/B LOS grids; inputs stay read-only."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import netCDF4
import numpy as np
from pyproj import CRS

from make_intf_nsr_gslc import read_grid, write_grid
from plot_gslc_gmt import plot_grid, plot_dependencies


def normalize_polygon(values):
    vertices = np.asarray(values, dtype='f8')
    if vertices.size < 6 or vertices.size % 2 or not np.isfinite(vertices).all():
        raise ValueError('Polygon needs at least three finite E N pairs in boundary order')
    vertices = vertices.reshape(-1, 2)
    if np.array_equal(vertices[0], vertices[-1]): vertices = vertices[:-1]
    area = abs(np.sum(vertices[:,0]*np.roll(vertices[:,1],-1)-vertices[:,1]*np.roll(vertices[:,0],-1)))/2
    if len(vertices)<3 or area<1 or len(np.unique(vertices,axis=0))!=len(vertices):
        raise ValueError('Polygon is degenerate; order vertices around its boundary')
    return vertices


def polygon_selection(x, y, vertices):
    """Ray crossing with boundary inclusion; compute in 256-row stripes."""
    selected = np.zeros((len(y), len(x)), dtype=bool)
    xx = x[None, :]
    for row in range(0, len(y), 256):
        yy = y[row:row+256, None]
        inside = np.zeros((len(yy), len(x)), dtype=bool)
        edge = np.zeros_like(inside)
        for a, b in zip(vertices, np.roll(vertices, -1, axis=0)):
            dx, dy = b-a
            if dy != 0:
                inside ^= ((a[1]>yy)!=(b[1]>yy)) & (xx < a[0]+(yy-a[1])*dx/dy)
            cross = (xx-a[0])*dy-(yy-a[1])*dx
            edge |= ((np.abs(cross)<=1e-6*np.hypot(dx,dy)) &
                     (xx>=min(a[0],b[0])-1e-6) & (xx<=max(a[0],b[0])+1e-6) &
                     (yy>=min(a[1],b[1])-1e-6) & (yy<=max(a[1],b[1])+1e-6))
        selected[row:row+len(yy)] = inside | edge
    return selected


def preview_band(source, output, p1, p2, binaries, *, polygon=None, mode="keep"):
    report = json.loads((source/'run.json').read_text())
    if report.get('status') != 'complete' or not report.get('unwrapping'):
        raise ValueError(f'Run4 result is not recorded complete: {source}')
    data, x, y = read_grid(source/'los.grd')
    with netCDF4.Dataset(source/'los.grd') as ds:
        if ds['z'].units != 'm':
            raise ValueError('Expected LOS in meters')
        crs = CRS.from_wkt(ds['crs'].spatial_ref)
    if crs.to_epsg() != 32646:
        raise ValueError('These proposed points are for EPSG:32646; choose matching coordinates for other CRSs')
    slope = (p2[1]-p1[1])/(p2[0]-p1[0])
    boundary = lambda e: p1[1]+slope*(e-p1[0])
    selected = (y[:, None] >= boundary(x)[None, :]) if polygon is None else polygon_selection(x,y,polygon)
    keep = np.isfinite(data) & (selected if mode=='keep' else ~selected)
    selection = (('north side' if mode=='keep' else 'south side') if polygon is None else
                 ('polygon interior' if mode=='keep' else 'outside polygon'))
    if not keep.any():
        raise ValueError(f'The proposed line excludes all valid pixels: {source}')
    output.mkdir(parents=True)
    after = np.where(keep, data, np.nan).astype('f4')
    write_grid(output/'mask.grd', np.where(keep, 1, np.nan), x, y, crs,
               f'PREVIEW: retain {selection}; 1=retain, NaN=exclude', '1')
    write_grid(output/'los_masked.grd', after, x, y, crs, f'PREVIEW ONLY: LOS retained {selection}', 'm')
    gmt, gs = binaries['gmt'], binaries['gs']
    dx, dy = x[1]-x[0], y[0]-y[1]
    bounds = np.array([x[0]-dx/2, x[-1]+dx/2, y[-1]-dy/2, y[0]+dy/2])/1000
    region = '-R'+'/'.join(map(str, bounds))
    width = 12.0
    height = width*(bounds[3]-bounds[2])/(bounds[1]-bounds[0])
    projection = f'-JX{width}c/{height}c'
    with tempfile.TemporaryDirectory(prefix='.plot_', dir=output) as temporary:
        work = Path(temporary)
        (work/'user').mkdir()
        env = dict(os.environ, GMT_USERDIR=str(work/'user'), GMT_TMPDIR=str(work))
        for name, z in (('before', data), ('after', after)):
            plot_grid(work/f'{name}.grd', z*1000, x, y)
        palette = work/'los.cpt'
        ps = work/'mask_preview.ps'
        with (output/'gmt_preview.log').open('w') as log:
            def run(args, *, text=None, first=False):
                with ps.open('wb' if first else 'ab') as stream:
                    subprocess.run(list(map(str, args)), cwd=work, env=env,
                        input=None if text is None else text.encode(), stdout=stream, stderr=log, check=True)
            subprocess.run([gmt,'set','GMT_THEME','classic','GMT_HISTORY','false',
                'COLOR_NAN','255/255/254','MAP_FRAME_TYPE','plain'], cwd=work, env=env,
                stdout=log, stderr=log, check=True)
            if (source/'los.cpt').is_file():
                shutil.copyfile(source/'los.cpt', palette)
            else:
                finite = data[np.isfinite(data)]*1000
                lo, hi = float(finite.min()), float(finite.max())
                if hi-lo < 1: lo, hi = lo-1, hi+1
                with palette.open('wb') as stream:
                    subprocess.run([gmt,'makecpt','-Cpolar',f'-T{lo}/{hi}/{(hi-lo)/100}','-Z','-D'],
                        cwd=work, env=env, stdout=stream, stderr=log, check=True)
            for i, name in enumerate(('before','after')):
                px, py = 2.2+i*16.0, 4.0
                run([gmt,'grdimage',work/f'{name}.grd',region,projection,f'-C{palette}','-nn+a',
                     '-Bxa+lEasting (km)','-Bya+lNorthing (km)','-BWSen',
                     f'-Xa{px}c',f'-Ya{py}c','--FONT_LABEL=10p','--FONT_ANNOT_PRIMARY=9p',
                     '--PS_MEDIA=A1',*(['-P'] if i==0 else ['-O']),'-K'],first=i==0)
                outline = (np.array([[e,boundary(e*1000)/1000] for e in bounds[:2]]) if polygon is None else
                           np.vstack([polygon,polygon[0]])/1000)
                line = ''.join(f'{e} {n}\n' for e,n in outline)
                run([gmt,'psxy',region,projection,f'-Xa{px}c',f'-Ya{py}c','-W1.2p,black,--','-O','-K'],text=line)
                if i==0:
                    vertices = (p1,p2) if polygon is None else polygon
                    points = ''.join(f'{e/1000} {n/1000}\n' for e,n in vertices)
                    run([gmt,'psxy',region,projection,f'-Xa{px}c',f'-Ya{py}c','-Sc0.15c','-Gblack','-O','-K'],text=points)
                    labels = ''.join(f'{e/1000} {n/1000} P{j}\n' for j,(e,n) in enumerate(vertices,1))
                    run([gmt,'pstext',region,projection,f'-Xa{px}c',f'-Ya{py}c','-F+f9p,Helvetica-Bold+jBL','-D0.1c/0.1c','-N','-O','-K'],text=labels)
                run([gmt,'psscale',f'-C{palette}',f'-Dx{px+width/2}c/{py-1.65}c+w9c/0.25c+jTC+h+e',
                     '-Bxa+lLOS (mm; positive toward satellite)','-Xa0c','-Ya0c',
                     '--FONT_LABEL=9p','--FONT_ANNOT_PRIMARY=9p','-O','-K'])
                title = 'Original LOS + proposed boundary' if i==0 else f'Masked LOS: retain {selection}'
                run([gmt,'pstext','-R0/40/0/50','-JX40c/50c','-Xa0c','-Ya0c',
                     '-F+f11p,Helvetica-Bold+jBC','-N','-O','-K'],text=f'{px+width/2} {py+height+0.7} {title}\n')
            band = report['arguments']['frequency']
            run([gmt,'pstext','-R0/40/0/50','-JX40c/50c','-Xa0c','-Ya0c',
                 '-F+f13p,Helvetica-Bold+jBC','-N','-O','-K'],
                text=f'16.2 {height+5.9} NISAR GSLC {band} | mask preview | EPSG:32646\n')
            note = (f'P1=({p1[0]:g}, {p1[1]:g}) m; P2=({p2[0]:g}, {p2[1]:g}) m' if polygon is None else
                    f'{len(polygon)}-vertex polygon; mode={mode}; vertices saved in preview.json')
            run([gmt,'pstext','-R0/40/0/50','-JX40c/50c','-Xa0c','-Ya0c',
                 '-F+f9p,Helvetica+jBC','-N','-O'],
                text=f'16.2 0.2 PREVIEW ONLY | {note} | shared color scale\n')
            subprocess.run([gmt,'psconvert',ps,'-Tf','-P','-A+m0.25c',f'-D{output}',
                            '-Fmask_preview',f'-G{gs}'],cwd=work,env=env,stdout=log,stderr=log,check=True)
    result = dict(source=str(source), p1_m=p1, p2_m=p2, epsg=32646, keep=selection, mode=mode,
                  polygon_m=None if polygon is None else polygon.tolist(),
                  valid_before=int(np.isfinite(data).sum()), retained=int(keep.sum()),
                  removed=int((np.isfinite(data)&~keep).sum()), preview_only=True,
                  pdf=str(output/'mask_preview.pdf'))
    (output/'preview.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,default=Path.cwd())
    parser.add_argument('-o','--output',type=Path,default=Path('iono_mask_preview'))
    parser.add_argument('--frequency',choices=('A','B','AB'),default='AB')
    parser.add_argument('--p1',nargs=2,type=float,default=[520000.,3510000.],metavar=('E','N'))
    parser.add_argument('--p2',nargs=2,type=float,default=[730000.,3570000.],metavar=('E','N'))
    parser.add_argument('--polygon',nargs='+',type=float,metavar='COORD',
                        help='Polygon vertices E1 N1 E2 N2 ... in UTM meters, ordered around the boundary')
    parser.add_argument('--mode',choices=('keep','exclude'),default='keep',
                        help='Keep or exclude polygon interior; for a line, apply to its north side')
    parser.add_argument('--gmt')
    args = parser.parse_args()
    if not np.isfinite(args.p1+args.p2).all() or args.p1[0] >= args.p2[0]:
        parser.error('Use finite UTM meter coordinates with P2 east of P1')
    polygon = None
    if args.polygon is not None:
        if len(args.polygon)<6 or len(args.polygon)%2 or not np.isfinite(args.polygon).all():
            parser.error('--polygon requires at least three finite E N pairs in boundary order')
        polygon = np.array(args.polygon).reshape(-1,2)
        if np.array_equal(polygon[0],polygon[-1]): polygon=polygon[:-1]
        area = abs(np.sum(polygon[:,0]*np.roll(polygon[:,1],-1)-polygon[:,1]*np.roll(polygon[:,0],-1)))/2
        if len(polygon)<3 or area<1 or len(np.unique(polygon,axis=0))!=len(polygon):
            parser.error('Polygon is degenerate or has repeated vertices; order vertices around the boundary')
    project = args.project.resolve()
    output = (args.output if args.output.is_absolute() else project/args.output).resolve()
    if output.exists(): parser.error(f'Preview output already exists; choose a new -o directory: {output}')
    try:
        pair = json.loads((project/'gslc_project.json').read_text())['pair']
        binaries = plot_dependencies(args.gmt)
        results = {}
        for band in args.frequency:
            source = project/('LA' if band=='A' else 'LB')/'intf'/pair
            if not (source/'los.grd').is_file():
                if args.frequency=='AB':
                    print(f'SKIP {band}: LOS not available yet: {source}',flush=True)
                    continue
                raise ValueError(f'Missing completed LOS: {source}')
            state = json.loads((source/'run.json').read_text())
            if args.frequency=='AB' and (state.get('status')!='complete' or not state.get('unwrapping')):
                print(f'SKIP {band}: native unwrapping is not recorded complete yet',flush=True)
                continue
            target = output/('LA' if band=='A' else 'LB')
            results[band] = preview_band(source,target,args.p1,args.p2,binaries,polygon=polygon,mode=args.mode)
            print(json.dumps(results[band],indent=2),flush=True)
        if not results: raise ValueError('No completed A/B LOS results available')
        (output/'preview.json').write_text(json.dumps(results,indent=2)+'\n')
    except (ValueError,OSError,KeyError,subprocess.CalledProcessError) as exc:
        parser.exit(1,f'{exc}\n')


if __name__=='__main__':
    main()
