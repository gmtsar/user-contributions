#!/usr/bin/env python3
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""Export native GSLC phase/LOS and ionosphere products as EPSG:4326 GRD."""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import logging
import tempfile
import shutil
import subprocess

import netCDF4
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window
from pyproj import CRS

from make_intf_nsr_gslc import GridWriter, VERSION

from plot_gslc_gmt import plot_dependencies, plot_products

LAYERS = ("display_amp", "phase", "corr", "phasefilt", "unwrap", "mask", "conncomp", "unwrap_ref", "los")


_RASTERIO_CHECKED = False
LOG = logging.getLogger('nisar_gslc')


def reproject_dependencies(*, pdf=False, gmt=None, gdalwarp=None):
    """Check the Python warp path before expensive processing starts."""
    global _RASTERIO_CHECKED
    result = {'backend': 'rasterio + netCDF4'}
    if gdalwarp is not None:
        warp = shutil.which(str(gdalwarp))
        if not warp:
            raise ValueError(f'Cannot find requested gdalwarp: {gdalwarp}')
        result.update(gdalwarp=warp, backend='external gdalwarp')
    if not _RASTERIO_CHECKED:
        # Exercise the actual NC4 reader and coordinate transformation, not
        # merely module imports or the external gdalwarp executable's presence.
        with tempfile.TemporaryDirectory(prefix='gslc_crs_check_') as temporary:
            path = Path(temporary)/'probe.grd'
            x = 574260 + np.arange(8)*40.0
            y = 3625340 - np.arange(8)*40.0
            writer = GridWriter(path, x, y, CRS.from_epsg(32646), 'CRS dependency check')
            try:
                writer[:] = np.ones((8, 8), dtype='f4')
            finally:
                writer.close()
            try:
                with rasterio.open(path) as src, WarpedVRT(src, crs='EPSG:4326', nodata=np.nan) as vrt:
                    if not np.isfinite(vrt.read(1)).any():
                        raise ValueError('CRS dependency check produced no valid pixels')
            except (rasterio.errors.RasterioError, ValueError) as exc:
                raise ValueError(f'rasterio cannot read/reproject NC4 GRD: {exc}') from exc
        _RASTERIO_CHECKED = True
    if pdf:
        result.update(plot_dependencies(gmt))
    return result


def source_identity(inputs):
    return {name: {'path': str(Path(path).resolve()), 'bytes': Path(path).stat().st_size,
                   'mtime_ns': Path(path).stat().st_mtime_ns} for name, path in inputs.items()}


def prepare_export(outdir, source, operation, inputs, *, resume=False, pdf=False, geotiff=False):
    """Reuse verified complete exports or archive a failed export before retry."""
    if not outdir.exists():
        return None
    if not outdir.is_dir():
        raise ValueError(f'Output is not a directory: {outdir}')
    if not any(outdir.iterdir()):
        return None
    if not resume:
        raise ValueError(f'Output must be new/empty; use --resume to retry a recorded export: {outdir}')
    manifest = outdir/'run.json'
    if not manifest.is_file():
        raise ValueError(f'Cannot resume an unrecognized export directory: {outdir}')
    old = json.loads(manifest.read_text())
    if old.get('operation') != operation or Path(old.get('source_directory', '')).resolve() != source:
        raise ValueError(f'Existing export belongs to another source/operation: {outdir}')
    data_ok = (old.get('status') == 'complete' and old.get('export_version') == VERSION
               and old.get('source_identity') == source_identity(inputs)
               and set(old.get('layers', {})) == set(inputs)
               and all((outdir/f'{name}_ll.grd').is_file() for name in inputs)
               and bool(old.get('geotiff')) == bool(geotiff))
    plots = old.get('plotting', {})
    pages = plots.get('pdfs', old.get('plots', []))
    if plots.get('multipage_pdf'):
        pages = pages + [plots['multipage_pdf']]
    pdf_ok = not pdf or (bool(pages) and all((outdir/name).is_file() for name in pages))
    tif_ok = not geotiff or all((outdir/f'{name}_ll.tif').is_file() for name in inputs)
    if data_ok and pdf_ok and tif_ok:
        LOG.info('Reusing complete geographic export: %s', outdir)
        return old
    if old.get('status') not in ('failed', 'running', 'complete'):
        raise ValueError(f'Cannot resume export with unknown status: {outdir}')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    backup = outdir.with_name(outdir.name + '_retry_' + stamp)
    outdir.rename(backup)
    LOG.info('Previous geographic export retained at %s; retrying %s', backup, outdir)
    return None


def warp_rasterio_layers(inputs, outdir, epsg, *, geotiff=False):
    """Stream a virtual nearest-neighbor warp into CF/GMT netCDF4 grids.

    This path does not call external GDAL tools or require GDAL's NC writer.
    Each read is limited to 256 output rows, including for full-scene grids.
    """
    groups, layers, commands = {}, {}, []
    with (outdir/'reprojection.log').open('w', encoding='utf-8') as log:
        log.write(f'Backend: rasterio {rasterio.__version__}, GDAL {rasterio.__gdal_version__}; netCDF4 writer\n')
        log.flush()
        with rasterio.Env(GDAL_CACHEMAX=64*1024*1024):
            for name, path in inputs.items():
                path = Path(path)
                target = outdir/f'{name}_ll.grd'
                log.write(f'Reprojecting {path} -> {target}\n'); log.flush()
                with netCDF4.Dataset(path) as ds:
                    var = next(v for v in ds.variables.values() if v.ndim == 2)
                    unit, title = getattr(var, 'units', '1'), getattr(var, 'long_name', name)
                with ExitStack() as stack:
                    src = stack.enter_context(rasterio.open(path))
                    if src.count != 1 or not src.crs or src.crs.to_epsg() != epsg:
                        raise ValueError(f'{path}: input CRS/band count conflicts with its report')
                    key = (src.shape, tuple(src.transform), src.crs.to_wkt())
                    options = dict(crs='EPSG:4326', nodata=np.nan, dtype='float32',
                                   resampling=Resampling.nearest, warp_mem_limit=64)
                    if key in groups:
                        shape, transform = groups[key]
                        options.update(height=shape[0], width=shape[1], transform=transform)
                    vrt = stack.enter_context(WarpedVRT(src, **options))
                    if vrt.crs.to_epsg() != 4326 or vrt.transform.a <= 0 or vrt.transform.e >= 0:
                        raise ValueError('Expected a north-up longitude/latitude output grid')
                    groups[key] = (vrt.shape, vrt.transform)
                    t = vrt.transform
                    x = t.c + (np.arange(vrt.width)+0.5)*t.a
                    y = t.f + (np.arange(vrt.height)+0.5)*t.e
                    if min(len(x), len(y)) < 2:
                        raise ValueError('Geographic output has fewer than two samples on an axis')
                    writer = GridWriter(target, x, y, CRS.from_epsg(4326), title, unit)
                    stack.callback(writer.close)
                    writer.file.reprojection_backend = 'rasterio WarpedVRT + netCDF4'
                    writer.file.resampling = 'nearest_neighbor'
                    writer.file.source_grid = str(path)
                    dst = None
                    if geotiff:
                        dst = stack.enter_context(rasterio.open(outdir/f'{name}_ll.tif', 'w', driver='GTiff',
                            width=vrt.width, height=vrt.height, count=1, dtype='float32',
                            crs=vrt.crs, transform=t, nodata=np.nan, compress='deflate', predictor=3, tiled=True))
                    count = 0
                    for row in range(0, vrt.height, 256):
                        window = Window(0, row, vrt.width, min(256, vrt.height-row))
                        values = vrt.read(1, window=window, masked=True).filled(np.nan)
                        writer[row:row+values.shape[0]] = values
                        count += int(np.isfinite(values).sum())
                        if dst is not None:
                            dst.write(values, 1, window=window)
                    layers[name] = dict(source=str(path), grid=target.name, units=unit,
                        shape=list(vrt.shape), bounds_wsen=list(vrt.bounds), spacing_degrees=[t.a, t.e],
                        finite_pixels=count, geotiff=f'{name}_ll.tif' if geotiff else None)
                    commands.append(dict(backend='rasterio WarpedVRT + netCDF4', source=str(path), target=str(target),
                                         source_epsg=epsg, target_epsg=4326, resampling='nearest_neighbor', rows_per_read=256))
                log.write(f'Completed {name}: {count} finite pixels\n'); log.flush()
                LOG.info('Geographic export completed: %s', target.name)
    (outdir/'reprojection_commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    last = next(reversed(layers.values()))
    return dict(layers=layers, source_grid_groups=len(groups), output_shape=last['shape'],
                output_bounds_wsen=last['bounds_wsen'], output_spacing_xy_degrees=last['spacing_degrees'],
                export_backend='rasterio WarpedVRT + netCDF4')


def warp_layers(inputs, outdir, epsg, *, gdalwarp=None, geotiff=False):
    """Warp each native grid; matching input grids must retain matching outputs."""
    if not inputs:
        raise ValueError("No input grids to reproject")
    outdir = Path(outdir)
    if gdalwarp is None:
        return warp_rasterio_layers(inputs, outdir, epsg, geotiff=geotiff)
    env = dict(os.environ)
    for key in ("PROJ_LIB", "PROJ_DATA"):
        env.pop(key, None)
    commands, groups, layers = [], {}, {}
    with (outdir/"reprojection.log").open("w", encoding="utf-8") as log:
        for name, path in inputs.items():
            path = Path(path)
            with rasterio.open(path) as src:
                if src.count != 1 or not src.crs or src.crs.to_epsg() != epsg:
                    raise ValueError(f"{path}: input CRS/band count conflicts with its report")
                source_key = (src.shape, tuple(src.transform), src.crs.to_wkt())
            with netCDF4.Dataset(path) as src:
                variable = next(v for v in src.variables.values() if v.ndim == 2)
                unit = getattr(variable, 'units', '1')
                title = getattr(variable, 'long_name', name)
            target = outdir/f"{name}_ll.grd"
            command = [gdalwarp, "-s_srs", f"EPSG:{epsg}", "-t_srs", "EPSG:4326",
                       "-r", "near", "-of", "netCDF", "-ot", "Float32", "-dstnodata", "nan",
                       "-co", "FORMAT=NC4C", "-co", "COMPRESS=DEFLATE", str(path), str(target)]
            commands.append(command)
            log.write(json.dumps(command)+"\n"); log.flush()
            if subprocess.run(command, env=env, stdout=log, stderr=log).returncode:
                log.flush()
                detail = (outdir/"reprojection.log").read_text(errors="replace")[-4000:]
                raise RuntimeError(f"gdalwarp failed for {name}; see {outdir/'reprojection.log'}\n{detail}")
            # Retain physical units: phase is radians; LOS is meters.
            with netCDF4.Dataset(target, 'a') as converted:
                for v in converted.variables.values():
                    if v.ndim == 2:
                        v.units, v.long_name = unit, title
            with ExitStack() as stack:
                ds = stack.enter_context(rasterio.open(target))
                grid = (ds.shape, ds.transform, ds.crs)
                if not ds.crs or ds.crs.to_epsg() != 4326:
                    raise ValueError(f"{target}: expected EPSG:4326")
                if source_key in groups and groups[source_key] != grid:
                    raise ValueError("Matching native grids produced inconsistent geographic grids")
                groups[source_key] = grid
                dst = None
                if geotiff:
                    profile = dict(driver="GTiff", width=ds.width, height=ds.height, count=1,
                                   dtype="float32", crs=ds.crs, transform=ds.transform,
                                   nodata=np.nan, compress="deflate", predictor=3, tiled=True)
                    dst = stack.enter_context(rasterio.open(outdir/f"{name}_ll.tif", "w", **profile))
                    dst.update_tags(source_grid=str(path), resampling="nearest_neighbor")
                count = 0
                for _, window in ds.block_windows(1):
                    values = ds.read(1, window=window, masked=True).filled(np.nan)
                    count += int(np.isfinite(values).sum())
                    if dst is not None:
                        dst.write(values, 1, window=window)
                layers[name] = {"source": str(path), "grid": target.name, "units": unit,
                                "shape": list(ds.shape), "bounds_wsen": list(ds.bounds),
                                "spacing_degrees": [ds.transform.a, ds.transform.e],
                                "finite_pixels": count,
                                "geotiff": f"{name}_ll.tif" if geotiff else None}
    (outdir/"reprojection_commands.json").write_text(json.dumps(commands, indent=2)+"\n")
    last = next(reversed(layers.values()))
    return {"layers": layers, "source_grid_groups": len(groups), "output_shape": last['shape'],
            "output_bounds_wsen": last['bounds_wsen'], "output_spacing_xy_degrees": last['spacing_degrees']}


def preview_pdfs(outdir, names):
    renderer = shutil.which("pdftoppm")
    if not renderer:
        return []  # GMT PDFs do not require Poppler; previews are optional.
    outputs = []
    for name in names:
        path = Path(outdir)/name
        subprocess.run([renderer, "-scale-to", "1600", "-singlefile", "-png", str(path),
                        str(path.with_suffix(""))], check=True, capture_output=True)
        outputs.append(path.with_suffix('.png').name)
    return outputs


def reproject(directory, output=None, *, pdf=False, geotiff=False, gmt=None, gdalwarp=None, resume=False):
    source = Path(directory).resolve()
    outdir = Path(output).resolve() if output else source/"lonlat"
    original = json.loads((source/"run.json").read_text())
    dependencies = reproject_dependencies(pdf=pdf, gmt=gmt, gdalwarp=gdalwarp)
    inputs = {name: source/f"{name}.grd" for name in LAYERS if (source/f"{name}.grd").is_file()}
    if not inputs:
        raise ValueError("No supported input layers found")
    old = prepare_export(outdir, source, "reproject_existing_results", inputs, resume=resume, pdf=pdf, geotiff=geotiff)
    if old is not None:
        return old
    outdir.mkdir(parents=True, exist_ok=True)
    report = dict(original)
    report.update(status="running", operation="reproject_existing_results", epsg=4326,
                  created_utc=datetime.now(timezone.utc).isoformat(), source_directory=str(source),
                  source_run_json=str(source/"run.json"), source_epsg=original['epsg'],
                  source_output_spacing_xy_m=original['output_spacing_xy_m'],
                  resampling="nearest_neighbor", geotiff=geotiff,
                  export_version=VERSION, source_identity=source_identity(inputs))
    report.pop('output_spacing_xy_m', None)
    try:
        report.update(warp_layers(inputs, outdir, original['epsg'],
                                  gdalwarp=dependencies.get('gdalwarp'), geotiff=geotiff))
        if report['source_grid_groups'] != 1:
            raise ValueError("Single-frequency products must share one native grid")
        if pdf:
            report['plotting'] = plot_products(outdir, report, dependencies, geographic=True)
            report['png_previews'] = preview_pdfs(outdir, report['plotting']['pdfs'])
        report['status'] = 'complete'
    except BaseException as exc:
        report['status'], report['error'] = 'failed', str(exc)
        raise
    finally:
        (outdir/'run.json').write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    return report


def reproject_ionosphere(directory, la, lb, original, *, pdf=True, gmt=None, output=None, gdalwarp=None, resume=False):
    source = Path(directory).resolve()
    outdir = Path(output).resolve() if output else source/'lonlat'
    dependencies = reproject_dependencies(pdf=pdf, gmt=gmt, gdalwarp=gdalwarp)
    inputs = {p.stem: p for p in sorted(source.glob('*.grd'))}
    inputs.update(unwrap_LA_before=Path(la)/'unwrap.grd', unwrap_LB_before=Path(lb)/'unwrap.grd')
    old = prepare_export(outdir, source, 'reproject_ionosphere_products', inputs, resume=resume, pdf=pdf)
    if old is not None:
        return old
    outdir.mkdir(parents=True, exist_ok=True)
    report = dict(original)
    report.update(status='running', operation='reproject_ionosphere_products', epsg=4326,
                  source_epsg=original['epsg'], source_directory=str(source),
                  resampling='nearest_neighbor', geotiff=False,
                  export_version=VERSION, source_identity=source_identity(inputs))
    try:
        report.update(warp_layers(inputs, outdir, original['epsg'], gdalwarp=dependencies.get('gdalwarp')))
        if pdf:
            from iono_gslc import plot_ionosphere
            report['plots'] = plot_ionosphere(outdir, la, lb, dependencies, geographic=True,
                                              frequencies=original['frequencies_hz'])
            report['png_previews'] = preview_pdfs(outdir, report['plots'])
        report['status'] = 'complete'
    except BaseException as exc:
        report['status'], report['error'] = 'failed', str(exc)
        raise
    finally:
        (outdir/'run.json').write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('-o', '--output', type=Path)
    parser.add_argument('--pdf', action='store_true')
    parser.add_argument('--geotiff', action='store_true')
    parser.add_argument('--gmt')
    parser.add_argument('--gdalwarp', help='Optional external backend; default uses rasterio and netCDF4')
    parser.add_argument('--resume', action='store_true', help='Reuse verified exports or retain a recorded failed export and retry')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    try:
        if (args.directory/'iono_report.json').is_file():
            original = json.loads((args.directory/'iono_report.json').read_text())
            report = reproject_ionosphere(args.directory, original['la_directory'], original['lb_directory'],
                                          original, pdf=args.pdf, gmt=args.gmt, output=args.output,
                                          gdalwarp=args.gdalwarp, resume=args.resume)
        else:
            report = reproject(args.directory, args.output, pdf=args.pdf, geotiff=args.geotiff,
                               gmt=args.gmt, gdalwarp=args.gdalwarp, resume=args.resume)
        print(json.dumps({k: report[k] for k in ('status', 'source_directory', 'output_shape', 'output_bounds_wsen')}, indent=2))
    except (ValueError, OSError, RuntimeError, KeyError, rasterio.errors.RasterioError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'{exc}\n')


if __name__ == '__main__':
    main()
