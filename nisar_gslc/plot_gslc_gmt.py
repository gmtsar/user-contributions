#!/usr/bin/env python3
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""Draw geocoded GSLC products using the local GMTSAR p2p GMT style.

UTM plots use temporary km coordinates. The --geographic mode draws previously
reprojected *_ll.grd grids with longitude/latitude axes. No radar FLIPUD is used.
Every layer uses unshaded colors with a bottom colorbar; no illumination is added.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import netCDF4
import numpy as np
from pyproj import CRS

from make_intf_nsr_gslc import read_grid, write_grid


def plot_dependencies(gmt=None):
    gmt = gmt or shutil.which("gmt")
    if not gmt or not Path(gmt).is_file():
        raise ValueError("GMT is required for --pdf; put gmt on PATH or specify --gmt")
    gs = shutil.which("gs")
    if not gs:
        candidate = Path(gmt).resolve().parent/"gs"
        gs = str(candidate) if candidate.is_file() else None
    if not gs:
        raise ValueError("Ghostscript (gs) is required for GMT PDF conversion")
    return {"gmt": str(Path(gmt).resolve()), "gs": str(Path(gs).resolve())}


def read_plot_grid(path, geographic=False):
    if not geographic:
        return read_grid(path)
    import rasterio
    with rasterio.open(path) as ds:
        if ds.count != 1 or not ds.crs or ds.crs.to_epsg() != 4326:
            raise ValueError(f"{path}: geographic plotting requires a single-band EPSG:4326 grid")
        t = ds.transform
        if t.b != 0 or t.d != 0 or t.a <= 0 or t.e >= 0:
            raise ValueError(f"{path}: expected a north-up unrotated geographic grid")
        data = ds.read(1, masked=True).filled(np.nan)
        return data, t.c+(np.arange(ds.width)+0.5)*t.a, t.f+(np.arange(ds.height)+0.5)*t.e


def plot_grid(path, data, x, y, geographic=False):
    """Temporary GMT grid in km or correctly declared longitude/latitude."""
    factor = 1 if geographic else 1000
    with netCDF4.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        ds.Conventions, ds.node_offset = "CF-1.8", np.int32(1)
        for name, axis in (("x", x/factor), ("y", y/factor)):
            ds.createDimension(name, len(axis))
            var = ds.createVariable(name, "f8", (name,))
            var[:] = axis
            var.units = ("degrees_east" if name == "x" else "degrees_north") if geographic else "km"
            var.axis = name.upper()
            if geographic:
                var.standard_name = "longitude" if name == "x" else "latitude"
            half = abs(axis[1]-axis[0])/2
            var.actual_range = [min(axis)-half, max(axis)+half]
        z = ds.createVariable("z", "f4", ("y", "x"), fill_value=np.float32(np.nan))
        z[:] = data
        finite = data[np.isfinite(data)]
        z.actual_range = [float(finite.min()), float(finite.max())]


def plot_products(outdir, report, binaries=None, overwrite=False, geographic=False, *, reuse_layers=()):
    """Plot available layers; callers may reuse named unchanged layers and rebuild the bundle."""
    outdir = Path(outdir).resolve()
    binaries = binaries or plot_dependencies()
    band = report["arguments"]["frequency"]
    pol = report["arguments"]["polarization"]
    names = ["display_amp", "phase", "corr", "phasefilt", "unwrap", "los"]
    suffix = "_ll" if geographic else ""
    bundle = outdir/f"gslc_{band}_{pol}{suffix}.pdf"
    targets = [outdir/f"{name}{suffix}.pdf" for name in names] + [bundle]
    if not overwrite and any(p.exists() for p in targets):
        raise ValueError("PDF output already exists; standalone plotter requires --overwrite to replace it")
    previous = {}
    if reuse_layers and (outdir/"plots.json").is_file():
        try:
            previous = json.loads((outdir/"plots.json").read_text(encoding="utf-8"))
        except (ValueError, OSError):
            pass  # Missing plot metadata is repaired by drawing the layers again.
    commands, pages, skipped = [], [], []
    reused, generated = [], []
    with tempfile.TemporaryDirectory(prefix=".gmt_", dir=outdir) as temporary:
        work = Path(temporary)
        (work/"user").mkdir()
        env = dict(os.environ, GMT_USERDIR=str(work/"user"), GMT_TMPDIR=str(work))
        with (outdir/"gmt_plot.log").open("w", encoding="utf-8") as log:
            def run(argv, output=None, append=False, text=None):
                argv = [str(a) for a in argv]
                commands.append(argv)
                log.write(json.dumps(argv)+"\n"); log.flush()
                if output:
                    with Path(output).open("ab" if append else "wb") as stream:
                        result = subprocess.run(argv, cwd=work, env=env, input=text.encode() if text else None,
                                                stdout=stream, stderr=log)
                else:
                    result = subprocess.run(argv, cwd=work, env=env, input=text.encode() if text else None,
                                            stdout=log, stderr=log)
                if result.returncode:
                    raise RuntimeError(f"GMT plotting command failed ({result.returncode}); see {outdir/'gmt_plot.log'}")

            gmt, gs = binaries["gmt"], binaries["gs"]
            run([gmt, "set", "GMT_THEME", "classic", "GMT_HISTORY", "false",
                 "FONT_ANNOT_PRIMARY", "10p,Helvetica", "FONT_LABEL", "11p,Helvetica",
                 "FONT_TITLE", "13p,Helvetica-Bold", "MAP_FRAME_TYPE", "plain",
                 "MAP_FRAME_PEN", "1p", "MAP_TICK_LENGTH_PRIMARY", "0.08i",
                 "MAP_LABEL_OFFSET", "0.10i", "MAP_ANNOT_OFFSET_PRIMARY", "0.06i",
                 "COLOR_NAN", "255/255/254", "COLOR_BACKGROUND", "black",
                 "COLOR_FOREGROUND", "white", "PS_MEDIA", "A3", "PS_CHAR_ENCODING", "ISOLatin1+"])
            if geographic:
                run([gmt, "set", "FORMAT_GEO_MAP", "ddd.xxF", "MAP_DEGREE_SYMBOL", "degree"])
            # Same display transform as filter.csh: sqrt(abs(complex interferogram)).
            complex_names = ("realfilt", "imagfilt") if (outdir/"realfilt.grd").exists() else ("real", "imag")
            if (not geographic and not (outdir/"display_amp.grd").exists()
                    and all((outdir/f"{name}.grd").exists() for name in complex_names)):
                real, x, y = read_grid(outdir/f"{complex_names[0]}.grd")
                imag, _, _ = read_grid(outdir/f"{complex_names[1]}.grd")
                amplitude = np.sqrt(np.hypot(real, imag))
                write_grid(outdir/"display_amp.grd", amplitude, x, y, CRS.from_epsg(report["epsg"]),
                           "sqrt(abs(complex interferogram used for phase)); p2p amplitude display")
                del real, imag, amplitude
            start = report["reference_product"]["zeroDopplerStartTime"][:10]
            end = report["secondary_product"]["zeroDopplerStartTime"][:10]
            heading = f"NISAR GSLC | {start} to {end} | {band} / {pol}"
            dx, dy = report.get("output_spacing_xy_m", report.get("source_output_spacing_xy_m", (None, None)))
            titles = {"display_amp": "Interferometric amplitude", "phase": "Wrapped phase",
                      "corr": "Coherence", "phasefilt": "Goldstein filtered phase",
                      "unwrap": "Unwrapped phase", "los": "LOS displacement"}
            corrected = bool(report.get("ionospheric_correction"))
            if corrected:
                titles.update(phase="Ionosphere-corrected wrapped phase",
                              phasefilt="Corrected phase supplied to SNAPHU",
                              unwrap="Re-unwrapped corrected phase",
                              los="Ionosphere-corrected LOS displacement")
            elif report.get("operation", "").startswith("pre-SNAPHU"):
                titles["phasefilt"] = "Pre-SNAPHU masked phase"
            mixed = report["reference_product"].get("demSource") != report["secondary_product"].get("demSource")
            limits = {}
            for name in names:
                stem = name+suffix
                path = outdir/f"{stem}.grd"
                if not path.exists():
                    continue
                pdf, palette = outdir/f"{stem}.pdf", outdir/f"{stem}.cpt"
                if (name in reuse_layers and previous.get("shading") == "none"
                        and previous.get("colorbar_position") == "bottom"
                        and pdf.name in previous.get("pdfs", [])
                        and name in previous.get("color_limits", {})
                        and pdf.is_file() and pdf.stat().st_size > 100
                        and palette.is_file() and palette.stat().st_size > 0
                        and pdf.stat().st_mtime_ns >= max(path.stat().st_mtime_ns, palette.stat().st_mtime_ns)):
                    pages.append(pdf)
                    reused.append(pdf.name)
                    limits[name] = previous["color_limits"][name]
                    log.write(f"Reusing unchanged layer: {pdf.name}\n"); log.flush()
                    continue
                data, x, y = read_plot_grid(path, geographic)
                finite = data[np.isfinite(data)]
                if not finite.size:
                    skipped.append({"layer": name, "reason": "no finite accepted pixels"})
                    continue
                if name == "los":
                    data, finite = data*1000, finite*1000
                grid = work/f"{stem}.nc"
                plot_grid(grid, data, x, y, geographic)
                factor = 1 if geographic else 1000
                xmin, xmax = (x[0]-(x[1]-x[0])/2)/factor, (x[-1]+(x[1]-x[0])/2)/factor
                ymin, ymax = (y[-1]+(y[1]-y[0])/2)/factor, (y[0]-(y[1]-y[0])/2)/factor
                aspect = (ymax-ymin)/(xmax-xmin)
                if geographic:
                    aspect /= math.cos(math.radians((ymax+ymin)/2))
                width = min(6.5, 7/aspect)
                height = width*aspect
                region = f"-R{xmin:.12g}/{xmax:.12g}/{ymin:.12g}/{ymax:.12g}"
                projection = f"-JM{width:.6g}i" if geographic else f"-JX{width:.6g}i/{height:.6g}i"
                colorbar_width = min(5, width*0.8)
                if name == "display_amp":
                    upper = max(float(3*np.std(finite)), float(np.quantile(finite, 0.5)), 1e-12)
                    run([gmt, "grd2cpt", grid, "-Cgray", "-Z", "-D", f"-L0/{upper:.12g}"], output=palette)
                    bar = ["-Bx0+lAmplitude (histogram equalized)"]
                    limits[name] = [0, upper]
                elif name == "corr":
                    run([gmt, "makecpt", "-Cgray", "-T0/0.8/0.1", "-Z", "-D"], output=palette)
                    bar = ["-Bxa0.2+lCorrelation"]
                    limits[name] = [0, 0.8]
                elif name in ("phase", "phasefilt"):
                    run([gmt, "makecpt", "-Crainbow", "-T-3.15/3.15/0.1", "-Z", "-D"], output=palette)
                    bar = ["-Bxa1.57+lPhase", "-By+lrad"]
                    limits[name] = [-3.15, 3.15]
                else:
                    mean, std = float(np.mean(finite)), float(np.std(finite))
                    lo, hi = math.floor(mean-2*std), math.ceil(mean+2*std)
                    if lo == hi:
                        lo, hi = lo-1, hi+1
                    step = max(1, math.ceil((hi-lo)/100))
                    scheme = ["-Cseis", "-I"] if name == "unwrap" else ["-Cpolar"]
                    run([gmt, "makecpt", *scheme, f"-T{lo}/{hi}/{step}", "-Z", "-D"], output=palette)
                    bar = ["-Bxa+lUnwrapped phase", "-By+lrad"] if name == "unwrap" else [
                        "-Bxa+lLOS displacement (range decrease)", "-By+lmm"]
                    limits[name] = [lo, hi]
                # Reuse the source p2p color scale when changing only map coordinates.
                if geographic and report.get("source_directory"):
                    original = Path(report["source_directory"])/f"{name}.cpt"
                    source_plots = Path(report["source_directory"])/"plots.json"
                    if original.is_file():
                        shutil.copyfile(original, palette)
                        if source_plots.is_file():
                            limits[name] = json.loads(source_plots.read_text())["color_limits"][name]
                ps = work/f"{stem}.ps"
                if geographic:
                    def tick_step(span, count):
                        base = 10**math.floor(math.log10(span/count))
                        return next(v*base for v in (1, 2, 5, 10) if v*base >= span/count)
                    xtick, ytick = tick_step(xmax-xmin, 5), tick_step(ymax-ymin, 4)
                    precision = max(2, -math.floor(math.log10(min(xtick, ytick))))
                    run([gmt, "set", "FORMAT_GEO_MAP", "ddd."+"x"*precision+"F"])
                    axes = [f"-Bxa{xtick:g}f{xtick/5:g}", f"-Bya{ytick:g}f{ytick/5:g}"]
                else:
                    axes = ["-Bxaf+lEasting (km)", "-Byaf+lNorthing (km)"]
                run([gmt, "grdimage", grid, region, projection, f"-C{palette}",
                     *axes, "-BWSen", *(["-nn+a"] if geographic else []),
                     "-X1.3i", "-Y2.0i", "-P", "-K"], output=ps)
                extend = "+ef" if name in ("corr", "display_amp") else "+e" if name in ("unwrap", "los") else ""
                run([gmt, "psscale", region, projection, f"-C{palette}",
                     f"-Dx{width/2:g}i/-0.65i+w{colorbar_width:g}i/0.18i+jTC+h{extend}",
                     *bar, "-O", "-K"], output=ps, append=True)
                def text_line(value, location, offset, font, justify, last=False):
                    run([gmt, "pstext", region, projection, f"-F+f{font}+j{justify}",
                         f"-D0i/{offset}i", "-N", "-O", *([] if last else ["-K"])], output=ps, append=True,
                        text=f"{(xmin+xmax)/2:.12g} {location:.12g} {value}\n")
                text_line(heading, ymax, 0.70, "13p,Helvetica-Bold", "BC")
                subtitle = (f"{titles[name]} | WGS 84 lon/lat | EPSG:4326" if geographic else
                            f"{titles[name]} | {dx:g} x {abs(dy):g} m | EPSG:{report['epsg']}")
                text_line(subtitle, ymax, 0.44, "10p,Helvetica", "BC")
                footer = "Diagnostic subset: mixed processing and DEM versions." if mixed else "Diagnostic interferogram in the GSLC map grid."
                notes = {"unwrap": "Arbitrary phase datum; connected components have independent ambiguities.",
                         "los": ("Positive toward satellite; zero is the user-supplied reference area."
                                 if report.get('reference') else
                                 "LOS = -wavelength / (4 pi) x unwrapped phase; positive toward satellite."),
                         "corr": "P2P grayscale: values at or above 0.8 are white; invalid pixels are pale cream.",
                         "display_amp": "Display transform: sqrt(abs(complex interferogram used for phase)).",
                         "phase": "Phase convention: reference * conjugate(secondary).",
                         "phasefilt": ("Modified Goldstein: alpha = 1 - patch coherence."
                                       if report['arguments'].get('goldstein_mode') == 'coherence'
                                       else f"Goldstein filter alpha = {report['arguments']['filter_alpha']:g}.")}
                gaussian = report.get("gaussian_filter", {})
                if gaussian.get("enabled") and name in ("phase", "phasefilt", "corr"):
                    notes[name] += f" Gaussian wavelength: {gaussian['wavelength_m']:g} m."
                if corrected and name in ("phase", "phasefilt"):
                    notes[name] = "Ionospheric estimate subtracted; Gaussian/Goldstein not repeated. Original A grid retained."
                text_line(footer, ymin, -1.35, "8p,Helvetica,gray30", "TC")
                text_line(notes[name], ymin, -1.52, "8p,Helvetica,gray30", "TC", last=True)
                run([gmt, "psconvert", ps, "-Tf", "-P", "-A+m0.10i", f"-D{outdir}", f"-F{stem}", f"-G{gs}"])
                if not pdf.is_file() or pdf.stat().st_size < 100:
                    raise RuntimeError(f"GMT did not produce {pdf}")
                pages.append(pdf)
                generated.append(pdf.name)
            if not pages:
                raise ValueError("No finite products available to plot")
            run([gs, "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite", "-dAutoRotatePages=/None",
                 "-dAutoFilterColorImages=false", "-dColorImageFilter=/FlateEncode",
                 "-dAutoFilterGrayImages=false", "-dGrayImageFilter=/FlateEncode",
                 "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false",
                 f"-sOutputFile={bundle}", *pages])
    result = {"style": "GMTSAR p2p palettes and annotations; unshaded geocoded grids",
              "shading": "none",
              "colorbar_position": "bottom",
              "pdfs": [p.name for p in pages], "multipage_pdf": bundle.name,
              "reused_pdfs": reused, "generated_pdfs": generated,
              "color_limits": limits, "skipped": skipped,
              "coordinate_note": ("EPSG:4326 geographic grids; Mercator map with longitude/latitude axes; nearest-neighbor image sampling."
                                  if geographic else "Temporary plotting axes in km; scientific grids remain in EPSG CRS meters.")}
    (outdir/"gmt_commands.json").write_text(json.dumps(commands, indent=2), encoding="utf-8")
    (outdir/"plots.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="A single-frequency output directory containing run.json")
    parser.add_argument("--gmt")
    parser.add_argument("--geographic", action="store_true", help="Plot already-reprojected *_ll.grd grids in EPSG:4326")
    parser.add_argument("--overwrite", action="store_true", help="Replace plotting products only; does not rerun the interferogram")
    args = parser.parse_args()
    try:
        report = json.loads((args.directory/"run.json").read_text())
        result = plot_products(args.directory, report, plot_dependencies(args.gmt), args.overwrite, args.geographic)
        print(json.dumps(result, indent=2))
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
