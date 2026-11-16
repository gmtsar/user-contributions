# Server Use and Recovery

**English** | [简体中文](SERVER_RUN.zh-CN.md) · [Home](README.md)

## Install the repository copy

Transfer the `nisar_gslc/` package to a server directory, then enter the
package directory. Input data and processing results are separate from the source
files and are excluded by `.gitignore`.

Activate the existing GMTSAR/SAR environment. Install only missing Python
dependencies if they are not already available:

```bash
python -m pip install -r requirements.txt
python -c 'import numpy,scipy,h5py,netCDF4,pyproj,rasterio; print(pyproj.CRS.from_epsg(32646)); print(rasterio.crs.CRS.from_epsg(4326))'
command -v csh phasefilt snaphu gmt gs
export NISAR_GSLC_PYTHON="$(command -v python)"
chmod +x *.csh
```

Keep executables on `PATH`. Individual process/unwrap entry points also accept
`--phasefilt`, `--snaphu` or `--gmt` overrides where applicable. The all-stage
controller accepts only its documented controller options; it does not forward
arbitrary numerical options.

The default geographic export uses rasterio and netCDF4 internally. It does
not require external `gdalwarp` HDF5/netCDF drivers. Optional `pdftoppm` enables
PNG previews but is not required for GRD or PDF output.

## Run in the repository directory

Put exactly two GSLC `.h5` files or symbolic links in `raw/`. Review the
`configs/` templates for a fresh project, then execute:

```bash
nohup ./run1-6_gslc_all.csh > run1-6_gslc_all.log 2>&1 &
tail -f run1-6_gslc_all.log
```

For project-specific configuration before run2:

```bash
./run1_prepare_GSLC_LA_LB.csh
# Edit LA/config_GSLC_A.txt, LB/config_GSLC_B.txt and config_IONO.txt.
nohup ./run1-6_gslc_all.csh --from-step 2 > run2-6.log 2>&1 &
```

## Keep the data project outside the repository

Replace the example path below and put the two acquisitions in that project's
`raw/`. Run commands from the code directory:

```bash
export GSLC_PROJECT=/path/to/data_project
mkdir -p "$GSLC_PROJECT/raw"
# Place the two GSLC files or links in "$GSLC_PROJECT/raw" first.
./run1_prepare_GSLC_LA_LB.csh --project "$GSLC_PROJECT"
# Edit configurations under "$GSLC_PROJECT".
nohup ./run1-6_gslc_all.csh --project "$GSLC_PROJECT" \
  --from-step 2 > run2-6.log 2>&1 &
```

Every individual stage, including a direct Python call, must use that same
`--project` value. For example:

```bash
python gslc_workflow.py process --project "$GSLC_PROJECT"
```

Do not run the individual processing command concurrently with the all-stage
controller for the same project.

## Logs and expected outputs

The outer log shows start/finish/failure for each stage. It prints the detailed
log directory, `PROJECT/logs/all_UTC_TIMESTAMP_PID/`. Inspect `run2.log`,
`run4.log`, etc. there. The controller stops after a failing stage.

Default final outputs are:

```text
PROJECT/LA_iono_corrected/unwrap.grd
PROJECT/LA_iono_corrected/los.grd
PROJECT/LA_iono_corrected/unwrap.pdf
PROJECT/LA_iono_corrected/los.pdf
PROJECT/LA_iono_corrected/lonlat/unwrap_ll.grd
PROJECT/LA_iono_corrected/lonlat/los_ll.grd
```

A/B are processed sequentially. Within a band, the configured SNAPHU tiles run
with `unwrap_nproc` simultaneous processes. No fixed full-scene runtime is
guaranteed; coverage, coherence, grid spacing and server resources all matter.

## Continue after a failure

The controller supports `--from-step N` for N=1…6. It still runs through step 6.
Existing outputs are protected; the option is not a general overwrite switch.

If native unwrapping completed but geographic export failed, use:

```bash
nohup ./run1-6_gslc_all.csh --from-step 4 --resume \
  > continue_run4-6.log 2>&1 &
```

`--resume` applies to run4. It checks the saved completion report and matching
settings, reuses complete unwrap/LOS grids, preserves a recorded failed
geographic directory, and retries export. It then continues incomplete bands
and run5/run6. It does not resume a partly solved SNAPHU tile job or silently
replace existing run5/run6 outputs.

To retry only run4 or only a geographic export:

```bash
./run4_unwrap_GSLC_LA_LB.csh --resume
python reproject_gslc.py LA/intf/<pair> --resume --pdf
```

Replace `<pair>` with the actual directory. Use `--project` or explicit paths
as appropriate for an external data project.

For a new ionospheric-filter experiment, edit `config_IONO.txt` and choose new
output directories:

```bash
./run5_iono_GSLC_LA_LB.csh -o iono_correction_trial
./run6_corrected_GSLC_LA.csh \
  --iono-dir iono_correction_trial -o LA_iono_corrected_trial
```

The optional polygon branch uses a saved preview and mode 2 result instead;
see the [workflow guide](WORKFLOW_GSLC.md#5-optional-polygon-mask-preview-first).
Direct Python alternatives are in the [README](README.md#alternative-call-python-directly).

## Scope of checks

The copied source corresponds to v0.6.5. Local diagnostic runs checked the
A/B processing route, masked correction, geographic export recovery and plot
reuse on a subset. A/B effective looks remain provisional. Runtime and
scientific validity on a new full scene must be assessed from its own results;
see the interpretation notes in the [README](README.md#interpretation-and-validation).
