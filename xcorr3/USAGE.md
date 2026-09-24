# More about using xcorr3

[Back to README](README.md)

Skip the MT or CUDA commands if that backend is not needed. Add PATH to your
shell startup file. MT defaults target Ubuntu x86-64; override its Makefile
include/library settings for other installations. CUDA overrides include
`CUDA_HOME` and `CUDA_ARCH` (e.g. `make cc CUDA_ARCH=sm_75`).
Geocoding additionally needs GMT, GMTSAR `proj_ra2ll.csh`, C shell and a matching
`trans.dat`; it is enabled only with the CC option `-geocode`.

## Usage

```sh
xcorr3 master.PRM secondary.PRM [options]
xcorr3 master.PRM secondary.PRM --backend mt -nproc 6 -nx 20 -ny 50
xcorr3 --backend cc master.PRM secondary.PRM -xsearch 128 -ysearch 128
xcorr3 --help
```

Both input files are required; missing inputs print usage/examples and exit 2.
SLC paths in the PRMs must be accessible from the working directory.

| Option | Meaning |
|---|---|
| `--backend original\|mt\|cc` | Choose the executable; no automatic fallback |
| `-nproc N` / `--nproc N` | MT processes, or `auto`; CLI overrides config |
| `-nx N -ny N` | Range/azimuth sample counts |
| `-freq` | Frequency correlation; mapped to the CC default |
| `-xsearch N -ysearch N` | Search radii |
| `-range_interp N -interp N` | Interpolation factors |
| `-noshift -norange -nointerp` | Ignore coarse shifts / disable interpolation |
| `--config FILE` | Read backend/process settings from an existing GMTSAR config |
| `--dry-run` | Print the selection and arguments without processing |

Append `xcorr_backend = mt` and `xcorr_nproc = 6` once to your existing config
(see [example](xcorr-options.config)). CLI selection wins; absent settings use
`original`/`auto`. Process settings in the config apply only to MT.

```sh
xcorr3 --config config.txt master.PRM secondary.PRM -nx 20 -ny 50
# Optional: route bare xcorr calls in a foreground, non-TOPS workflow.
xcorr3 --config config.txt --run p2p_processing.csh ALOS a b config.txt
```

`--run` uses a temporary child PATH. Absolute xcorr paths or scripts that reset
PATH bypass it. Stock TOPS geometric/ESD alignment does not call xcorr and is
not accelerated by this launcher. Intercepted xcorr failures remain nonzero
even if an outer script ignores the error and subsequently returns success.

Native output formats are preserved: normally `freq_xcorr.dat`, five columns
for original/MT and six for CC (additional `peak_snr`). CC does not support
time-domain or real/grid input modes; its results can differ from CPU results.
Use one job per output directory and check exit status before using results.

[Tests and validation](VALIDATION.md) · [MT details](backends/mt/README.md) ·
[CC details](backends/cc/README.md). Python help/routing tests run in PowerShell;
the scientific executables and `--run` require Linux/POSIX.
