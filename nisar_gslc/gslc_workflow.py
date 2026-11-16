#!/usr/bin/env python3
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""GSLC LA/LB workflow with six consecutive stages, run1 through run6."""
import argparse
import copy
import fcntl
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import h5py
import netCDF4
import numpy as np
from pyproj import CRS

import make_intf_nsr_gslc as core
import gslc_pair as entry
from gslc_config import load_text_pair, read_iono_config
from iono_gslc import estimate, nearest_fill

ROOT = Path(__file__).resolve().parent
CONFIG_NAMES = {'A':'config_GSLC_A.txt','B':'config_GSLC_B.txt'}


def json_read(path):
    return json.loads(Path(path).read_text())


def project_state(project):
    project = Path(project).resolve()
    state = json_read(project/'gslc_project.json')
    paths = [project/p for p in state['inputs']]
    return project, state, paths


def band_dir(project, state, band):
    return project/('LA' if band=='A' else 'LB')/'intf'/state['pair']


def config_paths(project):
    return (project/'LA'/CONFIG_NAMES['A'], project/'LB'/CONFIG_NAMES['B'])


def options(args, project, state, inputs, band, inspect=False):
    pa,pb = config_paths(project)
    cfg = load_text_pair(pa,pb,band)
    argv = ['full',*map(str,inputs),'--frequency',band,'--no-pdf']
    for key in ('phasefilt','snaphu','gmt'):
        if getattr(args,key,None):
            argv += ['--'+key,getattr(args,key)]
    cli = entry.parser().parse_args(argv)
    opt = entry.band_options(cli,cfg,band,inspect)
    opt.output = band_dir(project,state,band)
    opt.bbox = state.get('bbox')
    opt.pdf = not args.no_pdf
    if getattr(args,'threshold',None) is not None:
        opt.correlation_threshold = args.threshold
    if getattr(args,'reference_bbox',None) is not None:
        opt.reference_bbox = args.reference_bbox
    return opt


def prepare(args):
    project = args.project.resolve()
    raw = project/'raw'
    if bool(args.reference) != bool(args.secondary):
        raise ValueError('Specify both --reference and --secondary, or place exactly two GSLC .h5 files in raw/')
    files = [args.reference.resolve(),args.secondary.resolve()] if args.reference else sorted(raw.glob('*.h5'))
    if len(files)!=2:
        raise ValueError(f'Place exactly two NISAR GSLC .h5 files in {raw}; found {len(files)}')
    if files[0].resolve()==files[1].resolve():
        raise ValueError('The two raw/ inputs must resolve to different GSLC HDF5 files')
    if files[0].name == files[1].name:
        raise ValueError('Input filenames must differ so raw/ links remain unambiguous')
    from reproject_gslc import reproject_dependencies
    reproject_dependencies(pdf=not args.no_pdf, gmt=args.gmt)
    acquisitions=[]
    for file in files:
        with h5py.File(file,'r') as f:
            p = core.Product(f,file,'A','HH')
            acquisitions.append((p.meta['zeroDopplerStartTime'],file.resolve()))
    acquisitions.sort()
    files = [p for _,p in acquisitions]
    days = [datetime.fromisoformat(t[:19]).strftime('%Y%j') for t,_ in acquisitions]
    if days[0] == days[1]:
        raise ValueError('This workflow requires two different acquisition dates')
    state = {'version':core.VERSION,'pair':'_'.join(days),'bbox':args.bbox,
             'input_times':[t for t,_ in acquisitions], 'inputs':[f'raw/{p.name}' for p in files]}
    manifest = project/'gslc_project.json'
    if manifest.exists():
        previous=json_read(manifest)
        if any(previous.get(k)!=state[k] for k in ('pair','bbox','input_times','inputs')):
            raise ValueError('Existing project belongs to another pair/ROI; choose a new project directory')
        core.LOG.info('Project already prepared; existing configurations retained: %s',project)
        return previous
    for band in 'AB':
        conf = project/('LA' if band=='A' else 'LB')/CONFIG_NAMES[band]
        if conf.exists():
            raise ValueError(f'Unmanaged configuration already exists: {conf}')
    # Validate both bands using templates before creating directories.
    cfg=load_text_pair(ROOT/'configs/config_GSLC_A.txt',ROOT/'configs/config_GSLC_B.txt')
    cli=entry.parser().parse_args(['inspect',*map(str,files),'--no-pdf'])
    for band in 'AB':
        opt=entry.band_options(cli,cfg,band,True)
        opt.bbox=args.bbox
        opt.stage='filter'  # A small diagnostic ROI may need fewer unwrap tiles, edited after preparation.
        core.process_one(opt,emit_inspection=False)
    project.mkdir(parents=True,exist_ok=True)
    raw.mkdir(exist_ok=True)
    for file in files:
        dest=raw/file.name
        if dest.exists() or dest.is_symlink():
            if dest.resolve()!=file:
                raise ValueError(f'raw/ contains a conflicting filename: {dest}')
        else:
            dest.symlink_to(file)
    for band in 'AB':
        channel=project/('LA' if band=='A' else 'LB')
        (channel/'raw').mkdir(parents=True,exist_ok=True)
        (channel/'intf').mkdir(exist_ok=True)
        shutil.copy2(ROOT/'configs'/CONFIG_NAMES[band],channel/CONFIG_NAMES[band])
        for file in files:
            (channel/'raw'/file.name).symlink_to(Path('../../raw')/file.name)
    shutil.copy2(ROOT/'configs/config_IONO.txt',project/'config_IONO.txt')
    entry.write_report(manifest,state)
    core.LOG.info('Prepared GSLC %s: %s',state['pair'],project)
    return state


def process(args):
    project,state,inputs=project_state(args.project)
    opts=[]
    for band in args.frequency:
        opt=options(args,project,state,inputs,band,True)
        opt.stage='filter'
        core.process_one(opt,emit_inspection=False)
        if opt.output.exists() and any(opt.output.iterdir()):
            raise ValueError(f'Run2 output exists; use a new project for filter/sampling comparisons: {opt.output}')
        opt.inspect=False
        opts.append(opt)
    result={'status':'running','step':'run2','bands':{}}
    try:
        for opt in opts:
            result['bands'][opt.frequency]=core.process_one(opt)
            source=project/('LA' if opt.frequency=='A' else 'LB')/CONFIG_NAMES[opt.frequency]
            shutil.copy2(source,opt.output/f'config_GSLC_{opt.frequency}_used.txt')
        result['status']='complete'
    except BaseException as exc:
        result['status'],result['error']='failed',str(exc)
        raise
    finally:
        entry.write_report(project/'run2_report.json',result)
    return result


def quality_mask(directory,threshold):
    phase,x,y=core.read_grid(directory/'phasefilt.grd')
    corr,cx,cy=core.read_grid(directory/'corr.grd')
    mask,mx,my=core.read_grid(directory/'mask.grd')
    if not (np.array_equal(x,cx) and np.array_equal(x,mx) and np.array_equal(y,cy) and np.array_equal(y,my)):
        raise ValueError('Phase, correlation and mask grids do not share coordinates')
    good=np.isfinite(phase) & np.isfinite(corr) & (corr>=threshold) & np.isfinite(mask) & (mask>0)
    if (directory/'landmask.grd').exists():
        from iono_gslc import read_checked
        land,_,_=read_checked(directory/'landmask.grd',x,y)
        good &= np.isfinite(land) & (land>0)
    return phase,corr,good,x,y


def preview(args):
    project,state,inputs=project_state(args.project)
    result={}
    for band in args.frequency:
        directory=band_dir(project,state,band)
        report=json_read(directory/'run.json')
        opt=options(args,project,state,inputs,band)
        phase,corr,good,x,y=quality_mask(directory,opt.correlation_threshold)
        out=directory/args.preview_name
        if out.exists() and any(out.iterdir()):
            raise ValueError(f'Preview exists; use --preview-name with a new name: {out}')
        out.mkdir(parents=True,exist_ok=True)
        crs=CRS.from_epsg(report['epsg'])
        core.write_grid(out/'mask.grd',np.where(good,1,np.nan),x,y,crs,'Pre-SNAPHU valid/correlation/optional land mask')
        core.write_grid(out/'phasefilt.grd',np.where(good,phase,np.nan),x,y,crs,'Pre-SNAPHU phase','radian')
        shutil.copy2(directory/'corr.grd',out/'corr.grd')
        report=copy.deepcopy(report)
        report['operation']='pre-SNAPHU preview; no new phase filtering'
        report['quality_threshold']=opt.correlation_threshold
        report['quality_pixels']=int(good.sum())
        if not args.no_pdf:
            from plot_gslc_gmt import plot_products
            report['plotting']=plot_products(out,report,gmt_dependencies(args))
        entry.write_report(out/'run.json',report)
        result[band]=str(out)
    return result


def gmt_dependencies(args):
    from plot_gslc_gmt import plot_dependencies
    return plot_dependencies(args.gmt)


def verify_filter_configuration(opt,report):
    old=report['arguments']
    expected=[opt.output_spacing[0],-opt.output_spacing[1]]
    if not np.allclose(expected,report['output_spacing_xy_m'],rtol=0,atol=1e-8):
        raise ValueError('Output spacing changed after run2; rerun run2 in a new project')
    for key in ('filter_wavelength','goldstein_mode','filter_patch','min_valid_fraction'):
        if getattr(opt,key) != old[key]:
            raise ValueError(f'{key} changed after run2; rerun run2 in a new project')
    if opt.goldstein_mode=='fixed' and (opt.filter_alpha if opt.filter_alpha is not None else .5)!=old['filter_alpha']:
        raise ValueError('Fixed Goldstein alpha changed after run2')


def reference_outputs(directory,opt,report):
    if opt.reference_bbox is None:
        report['los'] = core.phase_to_los(directory/'unwrap.grd', directory/'los.grd', report['wavelength_m'])
        return
    phase,x,y=core.read_grid(directory/'unwrap.grd')
    cc,_,_=core.read_grid(directory/'conncomp.grd')
    xmin,xmax,ymin,ymax=opt.reference_bbox
    if not (xmin<xmax and ymin<ymax):
        raise ValueError('Invalid stable reference bbox')
    box=(y[:,None]>=ymin)&(y[:,None]<=ymax)&(x[None,:]>=xmin)&(x[None,:]<=xmax)
    good=box&np.isfinite(phase)&(cc>0)
    ids=np.unique(cc[good])
    if len(ids)!=1:
        raise ValueError('Stable reference must contain reliable pixels in exactly one connected-component label')
    zero=float(np.median(phase[good]))
    relative=np.where(cc==ids[0],phase-zero,np.nan)
    crs=CRS.from_epsg(report['epsg'])
    core.write_grid(directory/'unwrap_ref.grd',relative,x,y,crs,'Phase relative to user-provided stable reference','radian')
    core.write_grid(directory/'los.grd',-report['wavelength_m']/(4*np.pi)*relative,x,y,crs,
                    'Relative LOS, positive toward satellite','m')
    report['reference']={'phase_rad':zero,'component':int(ids[0]),'bbox':opt.reference_bbox,
                         'selection_is_user_supplied_not_verified_stable':True}


def unwrap_existing(directory,opt,interpolation,*,pdf=True,gmt=None):
    """Use saved phase/correlation; preserve run2 scientific grids exactly."""
    directory=Path(directory)
    report=json_read(directory/'run.json')
    if report['status']!='complete' or not (directory/'phasefilt.grd').exists():
        raise ValueError(f'Run2/phase input is incomplete: {directory}')
    if any((directory/p).exists() for p in ('unwrap.grd','conncomp.grd','snaphu_work')):
        raise ValueError(f'Unwrap output already exists: {directory}')
    if opt.correlation_threshold==0:
        raise ValueError('Threshold 0 means skip unwrapping; select a positive threshold for run4/run6')
    verify_filter_configuration(opt,report)
    validation=copy.copy(opt)
    validation.inspect=True
    validation.stage='unwrap'
    core.process_one(validation,emit_inspection=False)
    phase,corr,good,x,y=quality_mask(directory,opt.correlation_threshold)
    opt.looks=report['arguments']['looks']
    plan=core.unwrap_plan(opt,phase.shape)
    core.enforce_unwrap_guard(opt,plan)
    binary=core.executable('snaphu',opt.snaphu)
    crs=CRS.from_epsg(report['epsg'])
    previous=copy.deepcopy(report)
    result=copy.deepcopy(report)
    solver_options=copy.copy(opt)
    solver_options.reference_bbox=None
    solver_options.reference_pixel=None
    start=time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix='.unwrap_',dir=directory) as temp:
            work=Path(temp)
            solver_phase=np.where(good,phase,np.nan)
            if interpolation['interpolate_phase']:
                solver_phase=nearest_fill(solver_phase,interpolation['nearest_grid_radius'])
                solver_options.correlation_threshold=0
            core.write_grid(work/'phasefilt.grd',solver_phase,x,y,crs,'SNAPHU input; interpolated only when enabled','radian')
            core.write_grid(work/'corr.grd',np.where(good,corr,0),x,y,crs,'Solver correlation','1')
            core.unwrap_interferogram(solver_options,work,x,y,crs,report['wavelength_m'],result,binary,emit_los=False)
            # Remove interpolation-only pixels from scientific outputs.
            labels,_,_=core.read_grid(work/'conncomp.grd')
            unwrapped,_,_=core.read_grid(work/'unwrap.grd')
            accepted=good & np.isfinite(unwrapped) & (labels>0)
            core.write_grid(work/'unwrap.grd',np.where(accepted,unwrapped,np.nan),x,y,crs,
                            'Unwrapped phase; independent ambiguities may remain','radian')
            core.write_grid(work/'conncomp.grd',np.where(accepted,labels,0),x,y,crs,'Reliable connected-component label; tile-local when tiled')
            ids,counts=np.unique(np.where(accepted,labels,0),return_counts=True)
            result['unwrapping']['component_counts']={str(int(i)):int(n) for i,n in zip(ids,counts)}
            result['unwrapping']['accepted_pixels_after_original_mask']=int(accepted.sum())
            result['unwrapping']['interpolation_enabled']=bool(interpolation['interpolate_phase'])
            result['unwrapping']['nearest_grid_radius_pixels']=interpolation['nearest_grid_radius']
            if interpolation['interpolate_phase']:
                core.write_grid(work/'phasefilt_interp.grd',solver_phase,x,y,crs,'Solver-only nearest filled phase','radian')
            for name in ('unwrap.grd','conncomp.grd','snaphu_work','snaphu.log','phasefilt_interp.grd'):
                if (work/name).exists():
                    shutil.move(str(work/name),str(directory/name))
        result['arguments'].update(stage='unwrap',requested_stage='unwrap',
                                   correlation_threshold=opt.correlation_threshold,
                                   unwrap_tiles=opt.unwrap_tiles,unwrap_overlap=opt.unwrap_overlap,
                                   unwrap_nproc=opt.unwrap_nproc,effective_looks=opt.effective_looks,
                                   unwrap_mode=opt.unwrap_mode,defomax=opt.defomax)
        result['unwrapping_plan']=plan
        reference_outputs(directory,opt,result)
        if pdf:
            from plot_gslc_gmt import plot_products
            unchanged=() if result.get('ionospheric_correction') else ('display_amp','phase','corr','phasefilt')
            result['plotting']=plot_products(directory,result,gmt_dependencies(argparse.Namespace(gmt=gmt)),
                                            overwrite=True,reuse_layers=unchanged)
        result['status']='complete'
        result['unwrap_step_elapsed_seconds']=time.monotonic()-start
        entry.write_report(directory/'filter_run.json',previous)
        entry.write_report(directory/'run.json',result)
        (directory/'unwrap_config_used.json').write_text(json.dumps({
            'arguments':result['arguments'],'interpolation':interpolation},indent=2)+'\n')
        return result
    except BaseException as exc:
        entry.write_report(directory/'unwrap_failed.json',{'status':'failed','error':str(exc)})
        raise


def reuse_completed_unwrap(directory, opt, interpolation):
    """Reuse only a recorded complete unwrap with matching processing settings."""
    directory = Path(directory)
    report = json_read(directory/'run.json')
    if report.get('status') != 'complete' or not report.get('unwrapping'):
        raise ValueError(f'Existing unwrap is not recorded complete; retained for diagnosis: {directory}')
    verify_filter_configuration(opt, report)
    for key in ('correlation_threshold', 'unwrap_tiles', 'unwrap_overlap', 'unwrap_mode', 'defomax', 'effective_looks'):
        if report['arguments'].get(key) != getattr(opt, key):
            raise ValueError(f'{key} differs from completed unwrap; refusing to reuse {directory}')
    old = report['unwrapping']
    if bool(old.get('interpolation_enabled')) != bool(interpolation['interpolate_phase']):
        raise ValueError('Interpolation settings changed since completed unwrap')
    if interpolation['interpolate_phase'] and old.get('nearest_grid_radius_pixels') != interpolation['nearest_grid_radius']:
        raise ValueError('Interpolation radius changed since completed unwrap')
    if report.get('reference', {}).get('bbox') != opt.reference_bbox:
        raise ValueError('Reference area differs from completed unwrap')
    expected = tuple(report['output_shape'])
    coordinates = None
    for name in ('unwrap', 'conncomp', 'los'):
        path = directory/f'{name}.grd'
        if not path.is_file():
            raise ValueError(f'Cannot reuse incomplete unwrap; missing {path}')
        with netCDF4.Dataset(path) as ds:
            if ds['z'].shape != expected or CRS.from_wkt(ds['crs'].spatial_ref).to_epsg() != report['epsg']:
                raise ValueError(f'Existing unwrap geometry conflicts with its report: {path}')
            xy = (np.asarray(ds['x'][:]), np.asarray(ds['y'][:]))
            if coordinates is not None and not all(np.array_equal(a, b) for a, b in zip(coordinates, xy)):
                raise ValueError(f'Existing unwrap products have inconsistent grids: {path}')
            coordinates = xy
    core.LOG.info('Reusing completed %s unwrap and LOS; SNAPHU is not rerun: %s', opt.frequency, directory)
    return report


def unwrap(args):
    if args.lonlat:
        from reproject_gslc import reproject_dependencies
        reproject_dependencies(pdf=not args.no_pdf, gmt=args.gmt)
    project,state,inputs=project_state(args.project)
    interpolation=read_iono_config(project/'config_IONO.txt')
    result={}
    for band in args.frequency:
        opt=options(args,project,state,inputs,band)
        directory = band_dir(project,state,band)
        existing = any((directory/name).exists() for name in ('unwrap.grd', 'conncomp.grd', 'snaphu_work'))
        if args.resume and existing:
            result[band] = reuse_completed_unwrap(directory,opt,interpolation)
        else:
            result[band]=unwrap_existing(directory,opt,interpolation,pdf=not args.no_pdf,gmt=args.gmt)
        if args.lonlat:
            from reproject_gslc import reproject
            result[band]['lonlat'] = reproject(directory,pdf=not args.no_pdf,gmt=args.gmt,resume=args.resume)
    entry.write_report(project/'run4_report.json',{'status':'complete','bands':result})
    return {'status':'complete','pair':state['pair']}


def ionosphere(args):
    project,state,_=project_state(args.project)
    cfg=read_iono_config(args.iono_config or project/'config_IONO.txt')
    output=args.output or project/('iono_correction_mask' if args.mask_polygon else 'iono_correction')
    return estimate(band_dir(project,state,'A'),band_dir(project,state,'B'),output,cfg,
                    gmt=args.gmt,pdf=not args.no_pdf,lonlat=args.lonlat,
                    mask_polygon=args.mask_polygon,mask_mode=args.mask_mode)


def mask_source_identity(project, state, bands='AB'):
    identities = {}
    for band in bands:
        for name in ('unwrap.grd','los.grd','corr.grd','conncomp.grd'):
            path = band_dir(project,state,band)/name
            stat = path.stat()
            identities[str(path)] = {'bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns}
    return identities


def mask_workflow(args):
    """Mode 1 publishes a preview/selection; mode 2 uses that exact selection."""
    from preview_iono_mask import normalize_polygon, preview_band
    project,state,_ = project_state(args.project)
    destination = args.output or Path('iono_correction_mask')
    root = (destination if destination.is_absolute() else project/destination).resolve()
    if root == project or root in project.parents or any(
            root == project/name or project/name in root.parents for name in ('raw','LA','LB')):
        raise ValueError('Choose a separate mask workflow output directory')
    root.mkdir(parents=True,exist_ok=True)
    with (root/'.run5mask.lock').open('a+') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(f'Another run5mask task is already using {root}') from exc
        lock.seek(0);lock.truncate();lock.write(str(os.getpid()));lock.flush()
        selection_path = root/'mask_selection.json'
        saved = json_read(selection_path) if selection_path.is_file() else None
        if saved and (saved.get('workflow')!='run5mask_two_stage' or saved.get('project')!=str(project)
                      or saved.get('pair')!=state['pair']):
            raise ValueError('Existing mask selection belongs to another project/pair')
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        def archive(path):
            if path.exists():
                history = root/'history';history.mkdir(exist_ok=True)
                destination = history/(path.name+'_'+stamp)
                path.rename(destination)
                core.LOG.info('Previous %s retained at %s',path.name,destination)
        if args.mode == '1':
            if args.no_pdf:
                raise ValueError('Mode 1 is a visual preview; remove --no-pdf')
            points = args.mask_polygon if args.mask_polygon is not None else (saved or {}).get('vertices_xy_m')
            if points is None:
                raise ValueError('First preview requires --mask-polygon E1 N1 E2 N2 E3 N3 E4 N4')
            vertices = normalize_polygon(points)
            action = args.mask_mode or (saved or {}).get('mask_mode','exclude')
            binaries = gmt_dependencies(args)
            results = {}
            with tempfile.TemporaryDirectory(prefix='.preview_',dir=root) as temporary:
                prepared = Path(temporary)/'preview'
                for band in 'AB':
                    source = band_dir(project,state,band)
                    report = json_read(source/'run.json')
                    if report.get('status')!='complete' or not report.get('unwrapping') or not (source/'los.grd').is_file():
                        core.LOG.info('Skipping %s preview: native unwrap/LOS not complete yet',band)
                        continue
                    results[band] = preview_band(source,prepared/('LA' if band=='A' else 'LB'),
                        [520000.,3510000.],[730000.,3570000.],binaries,polygon=vertices,mode=action)
                if not results:
                    raise ValueError('No complete native A/B LOS results available for preview')
                selection = dict(workflow='run5mask_two_stage',version=core.VERSION,project=str(project),
                    pair=state['pair'],epsg=32646,mask_mode=action,vertices_xy_m=vertices.tolist(),
                    created_utc=datetime.now(timezone.utc).isoformat(),
                    source_identity=mask_source_identity(project,state,results))
                (prepared/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
                archive(root/'preview')
                prepared.rename(root/'preview')
            for band,result in results.items():
                directory = root/'preview'/('LA' if band=='A' else 'LB')
                result['pdf'] = str(directory/'mask_preview.pdf')
                entry.write_report(directory/'preview.json',result)
            entry.write_report(root/'preview/preview.json',results)
            temporary_selection = root/'.mask_selection.tmp'
            entry.write_report(temporary_selection,selection)
            temporary_selection.replace(selection_path)
            return dict(status='complete',mode=1,output=str(root),selection=str(selection_path),
                        previews={band:r['pdf'] for band,r in results.items()},
                        next='Run the same run5mask entry with mode 2 after reviewing the preview')
        if saved is None or not (root/'preview/selection.json').is_file():
            raise ValueError('Run mode 1 and review its preview before mode 2')
        if json_read(root/'preview/selection.json') != saved:
            raise ValueError('Preview and saved selection differ; rerun mode 1')
        if args.mask_polygon is not None and not np.array_equal(normalize_polygon(args.mask_polygon),saved['vertices_xy_m']):
            raise ValueError('Polygon differs from the preview; update mode 1 before mode 2')
        if args.mask_mode is not None and args.mask_mode != saved['mask_mode']:
            raise ValueError('Mask action differs from the preview; update mode 1 before mode 2')
        if mask_source_identity(project,state,saved_bands(saved,project,state)) != saved['source_identity']:
            raise ValueError('Native A/B inputs changed since preview; rerun mode 1')
        cfg = read_iono_config(args.iono_config or project/'config_IONO.txt')
        result_dir = root/'result'
        archive(result_dir)
        entry.write_report(root/'run5mask_state.json',dict(status='running',mode=2,pid=os.getpid(),result=str(result_dir)))
        try:
            result = estimate(band_dir(project,state,'A'),band_dir(project,state,'B'),result_dir,cfg,
                gmt=args.gmt,pdf=not args.no_pdf,lonlat=args.lonlat,
                mask_polygon=saved['vertices_xy_m'],mask_mode=saved['mask_mode'])
            result['mask_preview_selection'] = saved
            entry.write_report(result_dir/'iono_report.json',result)
        except BaseException as exc:
            entry.write_report(root/'run5mask_state.json',dict(status='failed',mode=2,error=str(exc)))
            raise
        entry.write_report(root/'run5mask_state.json',dict(status='complete',mode=2,result=str(result_dir)))
        return dict(status='complete',mode=2,output=str(result_dir),selection=str(selection_path))


def saved_bands(selection, project, state):
    return ''.join(band for band in 'AB' if str(band_dir(project,state,band)/'los.grd') in selection['source_identity'])


def correct(args):
    if args.lonlat:
        from reproject_gslc import reproject_dependencies
        reproject_dependencies(pdf=not args.no_pdf, gmt=args.gmt)
    project,state,inputs=project_state(args.project)
    iono_dir=(args.iono_dir or project/'iono_correction').resolve()
    selection_path=iono_dir/'mask_selection.json'
    selection=json_read(selection_path) if selection_path.is_file() else None
    if selection and selection.get('workflow')=='run5mask_two_stage':
        iono_dir=iono_dir/'result'
        if not (iono_dir/'iono_report.json').is_file():
            raise ValueError('Mask preview exists but mode 2 result is missing; run run5mask mode 2 first')
    iono_report=json_read(iono_dir/'iono_report.json')
    if iono_report['status']!='complete':
        raise ValueError('Run5 ionosphere estimate is incomplete')
    if selection:
        region=iono_report.get('spatial_mask') or {}
        if region.get('mode')!=selection['mask_mode'] or region.get('vertices_xy_m')!=selection['vertices_xy_m']:
            raise ValueError('Current mask preview differs from the correction; rerun mode 2 before run6')
        recorded=iono_report.get('mask_preview_selection')
        if recorded and any(recorded.get(key)!=selection.get(key) for key in
                ('project','pair','epsg','mask_mode','vertices_xy_m','source_identity')):
            raise ValueError('Preview inputs differ from the correction; rerun mode 2 before run6')
    recorded=iono_report.get('mask_preview_selection')
    if recorded and mask_source_identity(project,state,saved_bands(recorded,project,state))!=recorded['source_identity']:
        raise ValueError('Native inputs changed after the mask correction; regenerate preview and mode 2')
    la=band_dir(project,state,'A')
    if Path(iono_report['la_directory'])!=la:
        raise ValueError('Run5 belongs to a different LA result')
    report=copy.deepcopy(json_read(la/'filter_run.json' if (la/'filter_run.json').exists() else la/'run.json'))
    out=(args.output or project/'LA_iono_corrected').resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Corrected output must be new or empty: {out}')
    out.mkdir(parents=True,exist_ok=True)
    for name in ('phase.grd','phasefilt.grd'):
        shutil.copy2(iono_dir/'phase_LA_after.grd',out/name)
    for name in ('corr.grd','mask.grd','landmask.grd'):
        if (la/name).exists():
            shutil.copy2(la/name,out/name)
    if iono_report.get('spatial_mask') is not None:
        roi,x,y=core.read_grid(iono_dir/'roi_mask_A.grd')
        corr,cx,cy=core.read_grid(out/'corr.grd')
        if not (np.array_equal(x,cx) and np.array_equal(y,cy)):
            raise ValueError('Run5 ROI mask does not match the original A grid')
        keep=np.isfinite(roi)&(roi>0)
        crs=CRS.from_epsg(report['epsg'])
        core.write_grid(out/'mask.grd',np.where(keep,1,np.nan),x,y,crs,'Run5 polygon selection on valid A unwrap','1')
        core.write_grid(out/'corr.grd',np.where(keep,corr,np.nan),x,y,crs,'Original A coherence within run5 polygon selection','1')
    report['operation']='Ionosphere-corrected LA phase; re-unwrapping without reapplying Gaussian/Goldstein'
    report['ionospheric_correction']={'source':str(iono_dir),'report':str(iono_dir/'iono_report.json')}
    if iono_report.get('spatial_mask') is not None:
        report['ionospheric_correction']['spatial_mask']=iono_report['spatial_mask']
    missing=report.get('processing_recipe',{}).get('unimplemented_legacy_options',[])
    if 'ionospheric correction' in missing:
        missing.remove('ionospheric correction')
    report.pop('unwrapping',None)
    report['status']='complete'
    entry.write_report(out/'run.json',report)
    opt=options(args,project,state,inputs,'A')
    opt.output=out
    result=unwrap_existing(out,opt,read_iono_config(project/'config_IONO.txt'),pdf=not args.no_pdf,gmt=args.gmt)
    if args.lonlat:
        from reproject_gslc import reproject
        result['lonlat']=reproject(out,pdf=not args.no_pdf,gmt=args.gmt)
    return {'status':'complete','output':str(out),'los':str(out/'los.grd'),
            'los_lonlat':str(out/'lonlat/los_ll.grd') if args.lonlat else None}


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='step',required=True)
    for name in ('prepare','process','preview','unwrap','iono','correct','mask'):
        s=sub.add_parser(name)
        s.add_argument('--project',type=Path,default=Path.cwd(),help='Data project directory; default current directory')
        s.add_argument('--no-pdf',action='store_true')
        s.add_argument('--gmt')
        if name=='prepare':
            s.add_argument('--reference',type=Path,help='Optional input override; normally put two GSLC .h5 files in project/raw/')
            s.add_argument('--secondary',type=Path,help='Optional second input override; set together with --reference')
            s.add_argument('--bbox',nargs=4,type=float,help='Optional diagnostic ROI in product CRS meters; default full coverage')
        if name in ('process','preview','unwrap'):
            s.add_argument('--frequency',choices=('A','B','AB'),default='AB')
        if name in ('process','unwrap','correct'):
            s.add_argument('--phasefilt')
            s.add_argument('--snaphu')
        if name in ('preview','unwrap','correct'):
            s.add_argument('--threshold',type=float)
        if name=='preview':
            s.add_argument('--preview-name',default='presnaphu')
        if name=='unwrap':
            s.add_argument('--resume',action='store_true',
                           help='Reuse completed matching native unwraps and retry recorded geographic exports')
        if name in ('unwrap','iono','correct','mask'):
            s.add_argument('--lonlat',action=argparse.BooleanOptionalAction,default=True,
                           help='Export longitude/latitude GRD and figures; enabled by default')
        if name in ('unwrap','correct'):
            s.add_argument('--reference-bbox',nargs=4,type=float)
        if name in ('iono','correct','mask'):
            s.add_argument('-o','--output',type=Path)
        if name=='mask':
            s.add_argument('mode',choices=('1','2'),help='1: preview and save selection; 2: estimate using saved selection')
        if name in ('iono','mask'):
            s.add_argument('--iono-config',type=Path)
            s.add_argument('--mask-polygon',nargs='+',type=float,metavar='COORD',
                           help='Polygon E1 N1 E2 N2 ... in product CRS meters, ordered around its boundary')
            s.add_argument('--mask-mode',choices=('exclude','keep'),default='exclude' if name=='iono' else None,
                           help='Exclude or retain polygon interior for A/B ionosphere estimation')
        if name=='correct':
            s.add_argument('--iono-dir',type=Path)
    return p


def main(argv=None):
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    p=parser()
    args=p.parse_args(argv)
    if getattr(args,'threshold',None) is not None and not 0<=args.threshold<=1:
        p.error('threshold must lie in [0,1]')
    functions=dict(prepare=prepare,process=process,preview=preview,unwrap=unwrap,iono=ionosphere,correct=correct,mask=mask_workflow)
    try:
        result=functions[args.step](args)
        print(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False))
    except (ValueError,OSError,RuntimeError,KeyError,ImportError,subprocess.CalledProcessError) as exc:
        core.LOG.error('%s',exc)
        return 1
    except KeyboardInterrupt:
        core.LOG.error('Interrupted; partial output directories are retained for diagnosis')
        return 130
    return 0


if __name__=='__main__':
    raise SystemExit(main())
