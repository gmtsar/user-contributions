# xcorr_mt backend

CPU multiprocess implementation of GMTSAR xcorr using POSIX fork(). It preserves
location order and links the same GMTSAR/GMT numerical libraries as stock xcorr.

Install and select this backend through [xcorr3](../../README.md):

```sh
# From the xcorr3 package root
make install-mt GMTSAR_HOME=/usr/local/GMTSAR PREFIX="$HOME/.local"
xcorr3 master.PRM secondary.PRM --backend mt -nproc 6 -nx 20 -ny 50
```

The binary can also be called directly as `xcorr_mt`. Worker precedence is an
explicit `-nproc`, then `OMP_NUM_THREADS`, then online CPUs minus two (at least
one). Workers are processes, not OpenMP threads, and are capped by CPUs/points.
The GMTSAR parser still uses fixed-size filename buffers; use short input paths.

Native GMTSAR arguments are passed through, including time-domain and real/grid
modes. `-time` writes `time_xcorr.dat`; normal output is `freq_xcorr.dat`.
See [Dingri results](VALIDATION.md) and [methods](VALIDATION_METHODS.md).

[GPL-3.0-or-later](LICENSE); derived from GMTSAR xcorr.
