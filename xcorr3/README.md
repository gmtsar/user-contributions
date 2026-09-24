# xcorr3

[中文](README.zh-CN.md)

If you already use GMTSAR's `xcorr`, xcorr3 lets you run the correlation with
multiple CPU processes (`mt`) or an NVIDIA GPU (`cc`). You still provide the
two PRM files and the usual correlation options. You can also select the
original `xcorr` to compare results; it remains the default.

## Installation

The commands below run on Linux, including WSL2. You need Python 3 and a
working GMTSAR installation. For the CPU version, you also need GCC/Make and
the GMT, LAPACK, BLAS and libtiff development libraries used by your GMTSAR build.

From this directory, set `GMTSAR_HOME` to your built GMTSAR source tree:

```sh
make install-mt GMTSAR_HOME=/usr/local/GMTSAR PREFIX="$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

This installs `xcorr3` and `xcorr_mt`. Add the PATH line to your shell startup
file to keep it for future sessions. Your existing `xcorr` is left in place.

For the GPU version, you need an NVIDIA driver, CUDA Toolkit, a compatible
C/C++ compiler, Make, pkg-config and GLib development headers. Then run:

```sh
make cc                              # CUDA in /usr/local/cuda
make install-cc PREFIX="$HOME/.local"
```

Use `make cc-system` instead of `make cc` for system CUDA in `/usr/bin`.
If you only want the launcher, use `make install PREFIX="$HOME/.local"`.

## Try it

In your data directory, replace the example names with your two PRM files.
The SLC paths they contain must be accessible from that directory.

```sh
xcorr3 master.PRM secondary.PRM                          # original GMTSAR
xcorr3 master.PRM secondary.PRM --backend mt -nproc 6     # six CPU processes
xcorr3 master.PRM secondary.PRM --backend cc             # NVIDIA GPU
```

Add `-nx 20 -ny 50` to sample 20 range positions and 50 azimuth positions,
or `-xsearch 128 -ysearch 128` to set the search radii. Run `xcorr3 --help`
for more options; running it without both inputs also prints the usage.

The usual output is `freq_xcorr.dat`. Original/MT write five columns; CC adds
`peak_snr` as a sixth. CC supports frequency-domain SLC correlation and can
give different results from the CPU versions. We compared them using the
Dingri earthquake data; see the [results and figures](VALIDATION.md).

To choose a backend in your existing GMTSAR config or use it from a processing
script, see [the usage notes](USAGE.md). This does not accelerate GMTSAR's
TOPS geometric/ESD alignment, which does not call `xcorr`.

The launcher and MT use [GPL-3.0-or-later](LICENSE); MT is based on GMTSAR
`xcorr`. CC retains its [MIT license and CUI Hao/Jazz-0626 notices](backends/cc/LICENSE).
