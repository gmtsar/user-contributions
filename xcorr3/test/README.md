# Launcher integration tests

Portable tests (PowerShell: use `python`):

```sh
python3 -B test/test_xcorr3.py
```

For real Linux tests, install the launcher and requested backends first. Generate
small synthetic PRM/SLC/grid fixtures with the existing MT functional suite:

```sh
python3 backends/mt/test/test_modes.py --original "$(command -v xcorr)" \
  --candidate "$HOME/.local/bin/xcorr_mt" --work /tmp/xcorr3-fixtures-new
python3 test/test_integration.py --prefix "$HOME/.local" \
  --original "$(command -v xcorr)" --fixtures /tmp/xcorr3-fixtures-new/fixtures \
  --work /tmp/xcorr3-integration-new
python3 test/test_geocode.py --prefix "$HOME/.local" \
  --fixtures /tmp/xcorr3-fixtures-new/fixtures --work /tmp/xcorr3-geocode-new
```

Use a new work directory for each run. These Linux suites require NumPy, GMT,
C shell, the GMTSAR projection script, both built backends, and an NVIDIA GPU.
Inputs/results stay in the test directories. Integration cases compare native
output with direct backend/original results across frequency/time/real/grid
modes, interpolation settings, CPU worker counts, config and csh routing.
Stock `-nointerp` has known undefined offset fractions: compare coordinates and
correlation, and use the independently known (3,2) translation for MT/CC offsets.

Failure cases preserve previous output and test exact file-size limits, missing
inputs, ignored child exit codes, failed/missing projection products, and empty
filters. The geocoding suite uses a synthetic radar/geographic lookup and compares
all four grid arrays, including NaNs, to the direct CC result. It does not measure
earthquake geocoding accuracy or validate a full TOPS/ESD processing chain.

See [validation results](../VALIDATION.md) for the real Dingri data checks and
the distinction between fresh small-grid CPU checks and full-grid CUDA checks.
