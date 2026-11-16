#!/usr/bin/env python3
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""Portable A/B pair helper, with explicit input paths and per-band settings."""
import argparse
from datetime import datetime, timezone
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import time

import make_intf_nsr_gslc as core
from gslc_config import load_text_pair

ROOT = Path(__file__).resolve().parent
COMMON = {"polarization", "stage", "goldstein_mode", "filter_patch", "filter_alpha",
          "correlation_threshold", "unwrap_mode", "defomax", "unwrap_overlap", "unwrap_nproc",
          "min_valid_fraction", "block_size", "filter_block_size", "max_unwrap_pixels"}
BAND = {"looks", "output_spacing", "filter_wavelength", "effective_looks", "unwrap_tiles"}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("inspect", "full"), help="inspect is read-only; full processes all common coverage")
    p.add_argument("reference", type=Path, help="Earlier GSLC HDF5")
    p.add_argument("secondary", type=Path, help="Later GSLC HDF5")
    p.add_argument("-o", "--output", type=Path, help="New or empty output root; required for full")
    p.add_argument("--config", type=Path, help="Legacy JSON configuration; exclusive with --config-a/--config-b")
    p.add_argument("--config-a", type=Path, help="LA key = value configuration; default configs/config_GSLC_A.txt")
    p.add_argument("--config-b", type=Path, help="LB key = value configuration; default configs/config_GSLC_B.txt")
    p.add_argument("--frequency", choices=("A", "B", "AB"), default="AB")
    p.add_argument("--stage", choices=("intf", "filter", "unwrap"), help="Override configured last processing stage")
    p.add_argument("--unwrap-nproc", type=int, help="Override simultaneous SNAPHU processes (within each band)")
    p.add_argument("--no-pdf", action="store_true", help="Skip GMT figures; still output GRD")
    p.add_argument("--geotiff", action="store_true", help="Also output GeoTIFF; off by default")
    p.add_argument("--keep-snaphu-files", action="store_true")
    p.add_argument("--reference-bbox", nargs=4, type=float, metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
                   help="User-verified stable region in product CRS meters; tiled labels limit LOS coverage")
    for name in ("phasefilt", "snaphu", "gmt"):
        p.add_argument("--"+name, help="Override installed executable")
    p.add_argument("--json", action="store_true", help="Emit full preflight metadata in inspect mode")
    return p


def load_config(path):
    cfg = json.loads(path.read_text(encoding="utf-8"))
    if set(cfg)-{"description", "notes", "common", "bands", "source_files"}:
        raise ValueError("Unknown top-level configuration key")
    if set(cfg["common"])-COMMON or set(cfg["bands"]) != {"A", "B"}:
        raise ValueError("Configuration requires common options and bands A/B; unknown options are rejected")
    for band, options in cfg["bands"].items():
        if (set(options)-(BAND|COMMON) or not {"filter_wavelength", "unwrap_tiles"} <= set(options)
                or ("looks" in options) == ("output_spacing" in options)):
            raise ValueError(f"Invalid band {band} options")
    cfg.pop("source_files",None)  # Saved effective values replay without requiring old source-file paths.
    return cfg


def configured_pair(args):
    if args.config is not None:
        if args.config_a is not None or args.config_b is not None:
            raise ValueError("Use --config JSON or --config-a/--config-b text files, not both")
        return load_config(args.config)
    return load_text_pair(args.config_a or ROOT/"configs/config_GSLC_A.txt",
                          args.config_b or ROOT/"configs/config_GSLC_B.txt", args.frequency)


def band_options(args, cfg, band, inspect):
    values = {**cfg["common"], **cfg["bands"][band], "frequency": band}
    for name in ("stage", "unwrap_nproc", "reference_bbox", "phasefilt", "snaphu", "gmt"):
        if getattr(args, name) is not None:
            values[name] = getattr(args, name)
    argv = [str(args.reference), str(args.secondary)]
    for name, value in values.items():
        argv.append("--"+name.replace("_", "-"))
        argv.extend(map(str, value if isinstance(value, list) else [value]))
    if args.output:
        argv += ["--output", str(args.output/f"frequency{band}_{values.get('polarization', 'HH')}")]
    if inspect:
        argv.append("--inspect")
    if not args.no_pdf:
        argv.append("--pdf")
    if args.geotiff:
        argv.append("--geotiff")
    if args.keep_snaphu_files:
        argv.append("--keep-snaphu-files")
    # Deliberately no bbox or native window: this entry always uses the full common grid.
    return core.parser().parse_args(argv)


def preflight(args, cfg):
    bands = list(args.frequency) if args.frequency == "AB" else [args.frequency]
    planned = {}
    for band in bands:
        options = band_options(args, cfg, band, True)
        report = core.process_one(options, emit_inspection=False)
        if "unwrapping_plan" in report:
            core.enforce_unwrap_guard(options, report["unwrapping_plan"])
        planned[band] = report
    # Resolve all requested external tools before any output is created.
    checked_snaphu = set()
    for report in planned.values():
        stage = report["arguments"]["stage"]
        if stage in ("filter", "unwrap"):
            core.executable("phasefilt", args.phasefilt)
        if stage == "unwrap":
            binary = core.executable("snaphu", args.snaphu)
            if binary not in checked_snaphu and max(report["unwrapping_plan"]["tiles_yx"]) > 1:
                result = subprocess.run([binary, "--help"], capture_output=True, text=True, timeout=15)
                help_text = result.stdout+result.stderr
                if not all(option in help_text for option in ("--tile", "--nproc", "--tiledir")):
                    raise ValueError("SNAPHU lacks the required native tiling options; use a compatible 2.x build")
                checked_snaphu.add(binary)
    if not args.no_pdf:
        from plot_gslc_gmt import plot_dependencies
        plot_dependencies(args.gmt)
    if args.geotiff:
        import rasterio  # noqa: F401
    return planned


def print_plan(planned):
    print("Full common coverage; A/B processed sequentially; all dimensions are rows x columns.", flush=True)
    for band, report in planned.items():
        ny, nx = report["output_shape"]
        dx, dy = report["output_spacing_xy_m"]
        wave = report["arguments"]["filter_wavelength"]
        print(f"{band}: {ny} x {nx} ({ny*nx:,} pixels), X/Y {dx:g}/{abs(dy):g} m; Gaussian {wave:g} m", flush=True)
        if "unwrapping_plan" in report:
            p = report["unwrapping_plan"]
            print(f"   SNAPHU tiles Y/X={p['tiles_yx']}, overlap={p['overlap_yx_pixels']}, "
                  f"processes={p['processes']}, largest tile={p['max_tile_shape_yx']} "
                  f"({p['max_tile_pixels']:,} pixels)", flush=True)


def write_report(path, report):
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    temporary.replace(path)


def run(args):
    if args.mode == "full" and args.output is None:
        raise ValueError("full requires --output")
    cfg = configured_pair(args)
    planned = preflight(args, cfg)
    if args.mode == "inspect":
        if args.json:
            print(json.dumps(planned, indent=2, ensure_ascii=False, allow_nan=False))
        else:
            print_plan(planned)
            for note in cfg.get("notes", []):
                print("Note: "+note)
        return planned
    out = args.output.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"Output must be a new or empty directory: {out}")
    print_plan(planned)
    for note in cfg.get("notes", []):
        core.LOG.warning("%s", note)
    out.mkdir(parents=True, exist_ok=True)
    for band, source in cfg.get("source_files", {}).items():
        shutil.copy2(source, out/f"config_GSLC_{band}_used.txt")
    (out/"config_used.json").write_text(json.dumps(cfg, indent=2)+"\n", encoding="utf-8")
    report = {"version": core.VERSION, "status": "running", "coverage": "full common grid",
              "started_utc": datetime.now(timezone.utc).isoformat(), "config": cfg,
              "bands": {}, "preflight": planned,
              "note": "A/B independent interferograms. No automatic ionospheric/DEM correction or stable reference."}
    start = time.monotonic()
    write_report(out/"batch.json", report)
    try:
        for band in planned:
            report["active_band"] = band
            write_report(out/"batch.json", report)
            report["bands"][band] = core.process_one(band_options(args, cfg, band, False))
            write_report(out/"batch.json", report)
        report["status"] = "complete"
        report.pop("active_band", None)
    except BaseException as exc:
        report["status"], report["error"] = "failed", str(exc)
        raise
    finally:
        report["elapsed_seconds"] = time.monotonic()-start
        write_report(out/"batch.json", report)
    return report


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parser().parse_args(argv)
    try:
        run(args)
    except (ValueError, OSError, RuntimeError, KeyError, TypeError, subprocess.TimeoutExpired) as exc:
        core.LOG.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        core.LOG.error("Interrupted; inspect batch.json and per-band run.json for partial outputs")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
