# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""Strict readers for commentable GMTSAR-style GSLC configuration files."""
from pathlib import Path
import math

TEXT_TYPES = {
    "frequency": str, "polarization": str, "stage": str,
    "filter_wavelength": float, "output_spacing_x": float, "output_spacing_y": float,
    "goldstein_mode": str, "filter_alpha": float, "filter_patch": int,
    "correlation_threshold": float, "unwrap_mode": str, "defomax": float,
    "effective_looks": float, "unwrap_tiles_y": int, "unwrap_tiles_x": int,
    "unwrap_overlap_y": int, "unwrap_overlap_x": int, "unwrap_nproc": int,
    "min_valid_fraction": float, "block_size": int, "filter_block_size": int,
    "max_unwrap_pixels": int,
}

IONO_TYPES = {
    "filter_wavelength": float, "prefilter_a": float, "coherence_threshold": float,
    "surface_tension": float, "coarse_samples_per_wavelength": int,
    "interpolate_phase": int, "nearest_grid_radius": int,
}


def read_band_config(path, band):
    path = Path(path)
    values = {}
    for number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.partition("#")[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{number}: expected key = value")
        key, value = (v.strip() for v in line.split("=", 1))
        if key not in TEXT_TYPES:
            raise ValueError(f"{path}:{number}: unsupported parameter {key!r}; use GSLC map-grid parameters")
        if key in values:
            raise ValueError(f"{path}:{number}: duplicate parameter {key!r}")
        try:
            values[key] = TEXT_TYPES[key](value)
        except ValueError as exc:
            raise ValueError(f"{path}:{number}: invalid value for {key}: {value!r}") from exc
    if values.pop("frequency", band) != band:
        raise ValueError(f"{path}: expected frequency = {band}")
    for key in ("filter_wavelength", "output_spacing_x", "output_spacing_y"):
        if key not in values:
            raise ValueError(f"{path}: missing required parameter {key}")
    for target, first, second in (
            ("output_spacing", "output_spacing_x", "output_spacing_y"),
            ("unwrap_tiles", "unwrap_tiles_y", "unwrap_tiles_x"),
            ("unwrap_overlap", "unwrap_overlap_y", "unwrap_overlap_x")):
        if (first in values) != (second in values):
            raise ValueError(f"{path}: set both {first} and {second}")
        if first in values:
            values[target] = [values.pop(first), values.pop(second)]
    return values


def load_text_pair(path_a, path_b, bands="AB"):
    paths = {"A": Path(path_a), "B": Path(path_b)}
    options = {b: read_band_config(paths[b], b) for b in bands}
    return {
        "description": "GMTSAR-style key = value GSLC configuration; full common coverage",
        "notes": ["Output X/Y spacing and Gaussian wavelength are separate parameters.",
                  "SNAPHU effective_looks values are provisional, not measured independent looks.",
                  "Native tiled connected-component labels remain tile-local."],
        "common": {}, "bands": options,
        "source_files": {b: str(paths[b].resolve()) for b in bands},
    }


def read_iono_config(path):
    values = {}
    for n, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.partition("#")[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{n}: expected key = value")
        key, raw = [v.strip() for v in line.split("=", 1)]
        if key not in IONO_TYPES or key in values:
            raise ValueError(f"{path}:{n}: unknown or duplicate parameter {key!r}")
        try:
            values[key] = IONO_TYPES[key](raw)
        except ValueError as exc:
            raise ValueError(f"{path}:{n}: invalid {key}: {raw}") from exc
    defaults = dict(filter_wavelength=20000., prefilter_a=1000., coherence_threshold=.1,
                    surface_tension=.1, coarse_samples_per_wavelength=8,
                    interpolate_phase=0, nearest_grid_radius=300)
    defaults.update(values)
    v = defaults
    if any(not math.isfinite(float(x)) for x in v.values()):
        raise ValueError("Ionosphere parameters must be finite")
    if (v['filter_wavelength'] <= 0 or v['prefilter_a'] < 0 or not 0 <= v['coherence_threshold'] <= 1
            or not 0 <= v['surface_tension'] <= 1 or v['coarse_samples_per_wavelength'] < 2
            or v['interpolate_phase'] not in (0, 1) or v['nearest_grid_radius'] < 1):
        raise ValueError("Invalid ionosphere/interpolation configuration")
    return defaults
