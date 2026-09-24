# Validation methods and compatibility boundaries

Current results and acceptance status: [validation summary](VALIDATION.md).

## Reference and environment

The local reference is the unmodified installed GMTSAR `xcorr` executable,
SHA256 `4158d29678bbf79e69083cf511eba1639a862fe164b6a3b101c418146538f33c`.
The installed GMTSAR source checkout is
`90503c37b427e9f8c80157a4653367ec52df7f14`; its generated configuration reports
6.2.0. The candidate links that installation's existing `libgmtsar.a` and GMT
6.5.0. The test environment is Ubuntu 24.04 in WSL2, with 8 visible logical CPUs
and approximately 23 GiB RAM. Automatic worker selection uses 6 processes.
`OMP_NUM_THREADS` is unset and `OPENBLAS_NUM_THREADS=1` for both programs.
The linked `libgmtsar.a` has SHA256
`7b49161eb5dece78e873be9ecc3d6ee132baf57c907dfeafc1fcd43e90be3489`.
An unrelated pre-existing local edit in `update_PRM_sub.c` was left unchanged;
the reference executable itself is preserved byte-for-byte.

These measurements concern this environment, not every GMTSAR release or CPU.
The old source project report used a different computer and is not reused as
evidence for this validation.

## Real inputs and processing boundary

The before/after original Sentinel-1 SAFE products are:

- `S1A_IW_SLC__1SDV_20250105T121424_20250105T121452_057309_070D2E_EF04.SAFE`
- `S1A_IW_SLC__1SDV_20250117T121423_20250117T121451_057484_07141C_8E02.SAFE`

All three VV subswaths are regenerated using `make_s1a_tops` mode 1 without
shift grids. No existing resampled secondary image is reused. `calc_dop_orb`
and `SAT_baseline` supply identical initial geometry to both correlation
programs. The full PRM-declared scene is used; no data-dependent crop or
point selection is made to improve agreement.

GMTSAR rounds the last zero-based output line index down to a multiple of four
when writing the PRM. Its burst writer may retain 1–4 additional complete lines
in the SLC. These extra lines are recorded; both executables use the same PRM
extent. Every requested correlation patch is checked against both image bounds.

## Matrix and acceptance

| Configuration | Grid | Search parameters |
|---|---|---|
| Standard | 20 x 50 (1000 points) | x=128, y=128 |
| Wide | 20 x 50 (1000 points) | x=128, y=256 |
| Dense | 40 x 100 (4000 points) | x=128, y=128 |

All three configurations cover IW1, IW2 and IW3. The representative IW2 standard
and wide configurations use three paired repetitions with alternating execution
order. IW2 standard additionally checks 1, 2, 3 and 4 workers.

Checks require complete finite five-column outputs and identical point ordering
and coordinates. Exact SHA256 equality is reported separately. At points with
`corr >= 18` in either output, acceptance masks must agree, offset differences
must not exceed 1/16 pixel and correlation differences must not exceed 0.01.
The maximum difference between the two fitted offset models over the image
must be at most 0.01 pixel. All lower-quality differences are retained for
inspection; they are never silently dropped from difference counts.

Wall-clock timing includes process startup, computation and output publication.
Runs are sequential, without competing benchmark jobs. The representative
standard configuration must achieve at least 1.5x median speedup across three
paired repetitions. Other configurations and memory costs are also reported.
This measures `xcorr`, not total InSAR workflow runtime.

The saved `/usr/bin/time` records include peak RSS (`%M`, in KiB). This is not
a measurement of the simultaneous total memory of the parent and all workers;
do not interpret it as the complete parallel job's memory requirement. The
per-worker allocation estimate in the README is useful when selecting a worker
count, but also excludes some library workspaces and runtime overhead.

## Why byte equality is not a universal promise

In the tested original, `print_results` reads fractional offsets even when
`-nointerp` disables the function that assigns them. A known integer-translation
fixture exposed extremely large erroneous offsets in both the original and the
initial parallel implementation. The candidate now initializes those fractions
to zero for every point. The regression oracle in this mode is the known
translation `(3, 2)` pixels plus reference coordinate/correlation agreement,
rather than reproduction of uninitialized memory.

Separately, inherited `highres_corr.c` can leave interpolation cells unwritten
near correlation-window boundaries. Numerical comparison and downstream fitting
are therefore necessary even when most cases are byte-identical. That library
has not been changed for this validation.

## Repeating checks with other data

Use the same built GMTSAR installation for both executables. Keep separate
working directories, containing identical PRM files and read-only input links.
Run `test/test_modes.py` first, then use representative complete scenes and
`test/compare_outputs.py` to inspect differences. For example:

```sh
# Execute each command in its own otherwise identical input directory.
/usr/bin/time -v xcorr master.PRM aligned.PRM -nx 20 -ny 50 -xsearch 128 -ysearch 128
/usr/bin/time -v xcorr_mt master.PRM aligned.PRM -nx 20 -ny 50 -xsearch 128 -ysearch 128
fitoffset.csh 3 3 freq_xcorr.dat 18 > fitoffset.txt
```

Preserve the complete output files, exact commands, input/executable hashes,
exit codes and all repetitions. A zero exit code alone is not a correctness or
performance result. This validation does not certify TOPS ESD or phase accuracy.
