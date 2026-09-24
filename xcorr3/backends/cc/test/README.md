# Regression tests

These tests compile the actual command-line and PRM parsers into a small host
harness. They do not need CUDA, a GPU, or SLC imagery, and do not allocate image
or search-window buffers.

Requirements: Python 3, GCC, `pkg-config`, and GLib development headers (Ubuntu:
`build-essential pkg-config libglib2.0-dev`). Run from the repository root:

```sh
python3 test/test_inputs.py --work /tmp/xcorr-cc-input-tests \
  --sanitize address,undefined
```

Use `--source /path/to/checkout` to test another checkout; `--sanitize none` or
`--sanitize undefined` can be used if AddressSanitizer is unavailable. Each run
copies the production input sources into the work directory and records their
SHA-256 hashes in `results.json` alongside every case, exit code, and diagnostic.

Coverage includes preserved defaults and option overrides; missing, unknown,
empty, malformed, out-of-range and overflowed parameters; PRM line boundaries,
long values and paths with spaces/equals; required-key handling; strict numeric
conversion; safe derived dimensions; and deterministic malformed numeric cases.
Every invalid case must fail with a diagnostic identifying its option or key,
without a signal, timeout, or sanitizer report. Successful cases free returned
paths, allowing leak checking of the complete input layer.

These checks cover input handling only. They do not establish numerical CUDA
equivalence or validate output-file error handling.

`test_gpu_inputs.py` additionally exercises the production executable with GPU
visibility disabled for malformed inputs, checks preservation of an existing
output, and validates real SLC filenames containing spaces and `=`. The positive
cases require a GPU and the 1536×1024 known-shift `(3,2)` complex-SLC fixtures:

```sh
python3 test/test_gpu_inputs.py --binary ./xcorr_cc \
  --fixtures /path/to/known-shift/fixtures --work /tmp/xcorr-cc-production-inputs
```

## Output and geocoding failures

Run the actual host output/postprocessing code under ASan, UBSan and leak
checking (also requires G++). Use a fresh work directory for each run:

```sh
python3 test/test_outputs.py --work /tmp/xcorr-cc-output-tests --real-tools
```

Fault injection covers write/flush/sync/close errors, each external command,
missing or all-NaN grids, invalid output targets, and publication rollback with
and without previous outputs. Every handled failure must return nonzero, retain
previous results, and clean its private workspace. Mock commands test failure
propagation; `--real-tools` adds a synthetic mapping using installed GMT and
GMTSAR (`proj_ra2ll.csh`). Omit that flag when these tools are unavailable.
The synthetic case does not validate earthquake displacement accuracy.

Exercise the production CUDA executable's output handling with the same GPU
fixtures used above:

```sh
python3 test/test_output_binary.py --binary ./xcorr_cc \
  --fixtures /path/to/known-shift/fixtures --work /tmp/xcorr-cc-production-outputs
```

This checks the known shift, file-size limit, missing transformation, invalid
target types, and an unwritable directory. Run as a normal user, not root.
