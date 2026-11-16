# GSLC workflow provenance

The staged LA/LB layout, original RSLC run7 Gomba frequency combination (GSLC run5), mean-coherence mask,
frequency-dependent outlier interval, median/boxcar/median/boxcar/Gaussian
sequence and final GMT surface reconstruction follow Xin Wang's
[gmtsar-nisar-dinsar-workflow](https://github.com/Xin-Wang-520/gmtsar-nisar-dinsar-workflow)
at revision `f106c873b8dde30db5823fc6783e2367cd9a2f66` (retrieved 2026-11-15).

- Method: Gomba et al. (2016), as implemented by that workflow.
- Subsequent script development: Yajun Zhang.
- NISAR LA/LB workflow integration: Xin Wang, USTC, Hefei, China.
- `nearest_grid_c_exact.py` is reproduced from that revision, with its author
  header retained. This distribution uses GNU GPL version 3; see `LICENSE`.
- Other numerical programs (GMTSAR, GMT, SNAPHU, NumPy, SciPy, h5py, netCDF4,
  pyproj and GDAL) keep their respective authorship and licenses.

## Changes for GSLC

Inputs are geocoded complex HDF5 products, not RSLC radar-coordinate files.
Center frequencies, CRS and actual map spacings come from HDF5 / recorded
`run.json` metadata. The code does not fabricate PRMs or relabel radar ranges.
A is sampled at B's actual map coordinates, using linear interpolation of
unwrapped phase and nearest-neighbor sampling of coherence. Filter widths are
in the existing map CRS meters. Reconstruction returns to A's actual centers.
Only a narrow (1.5 B pixel) margin permits edge extrapolation; larger
unsupported areas remain NaN. The final A valid-data footprint is restored.

Gaussian/Goldstein processing before GSLC run4 keeps the selected GSLC recipe.
Solver-only nearest interpolation is optional, disabled by default, and uses
the upstream C-exact helper after normalizing array row orientation. The
scientific output retains the original correlation/valid-data mask.

GSLC is already geocoded. Run4-run6 export geographic GRD by default; this is
a map reprojection, not a radar-to-map geocoding step. `--no-lonlat` disables
that export on the individual staged entry points. Native grids retain their CRS.
LOS is converted directly as `-wavelength/(4*pi) * unwrapped_phase`, positive
toward the satellite, without subtracting a mean or selecting a reference area.
The default retains the inherited phase datum and is not an absolute displacement
reference. An explicitly supplied reference area remains optional. Tiled
connected-component labels remain local to the SNAPHU tiles. Full-scene seams,
cycle ambiguities and phase references need evaluation.

The original workflow and local GMTSAR installation are not edited. Its
run9 manual isolated-area masking is not part of this adaptation.
