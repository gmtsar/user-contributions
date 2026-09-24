#!/usr/bin/env python3
"""Compile the real input layer and test it without CUDA or SLC data.

Run on Linux/WSL with GCC, pkg-config and GLib development headers. Each case
runs in a fresh subprocess, so a crash cannot be mistaken for a rejected input.
The build copies the production files into --work for a reproducible snapshot.
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shlex
import shutil
import subprocess
import sys


MASTER = {
    "SLC_file": "master.SLC", "num_rng_bins": "4096", "num_patches": "2",
    "num_valid_az": "4096", "PRF": "1000",
}
SECONDARY = {**MASTER, "SLC_file": "secondary.SLC", "PRF": "1001",
             "rshift": "-16", "ashift": "3"}
SANITIZER = re.compile(r"AddressSanitizer|UndefinedBehaviorSanitizer|LeakSanitizer|runtime error:")


def prm(values: dict[str, str], newline: str = "\n", final: bool = True) -> bytes:
    return (newline.join(f"{key} = {value}" for key, value in values.items())
            + (newline if final else "")).encode()


@dataclasses.dataclass
class Case:
    name: str
    options: list[str] = dataclasses.field(default_factory=list)
    master: bytes = dataclasses.field(default_factory=lambda: prm(MASTER))
    secondary: bytes = dataclasses.field(default_factory=lambda: prm(SECONDARY))
    error_token: str | None = None
    expected: dict = dataclasses.field(default_factory=dict)
    argv: list[str] | None = None
    help_only: bool = False


def cases() -> list[Case]:
    result = [
        Case("defaults", expected={
            "m_nx": 4096, "m_ny": 8192, "s_nx": 4096, "s_ny": 8192,
            "x_offset": -16, "y_offset": 3, "xsearch": 64, "ysearch": 64,
            "nxl": 16, "nyl": 32, "ri": 2, "interp_factor": 16,
            "nyquist_split": 0, "snr_thr": 10, "psnr_thr": 5,
            "do_geocode": 0, "do_blockmedian": 1, "throttle_ms": 2,
            "throttle_every": 64, "device": 0,
            "m_path": "master.SLC", "s_path": "secondary.SLC",
            "astretcha": .001}),
        Case("explicit_zero_filters", ["-snr", "0", "-psnr", "0", "-throttle", "0"],
             expected={"snr_thr": 0, "psnr_thr": 0, "throttle_ms": 0}),
        Case("signed_zero_filters", ["-snr", "-0", "-psnr", "+0", "-throttle", "-0"],
             expected={"snr_thr": 0, "psnr_thr": 0, "throttle_ms": 0}),
        Case("positive_sign", ["-nx", "+20"], expected={"nxl": 20}),
        Case("arbitrary_interp_factor", ["-interp", "3"], expected={"interp_factor": 3}),
        Case("nointerp_ignores_derived_size", ["-interp", "2147483647", "-nointerp"],
             expected={"interp_factor": 1}),
        Case("explicit_controls", ["-nx", "20", "-ny", "50", "-xsearch", "128",
             "-ysearch", "256", "-range_interp", "4", "-interp", "8",
             "-snr", "18", "-psnr", "6", "-throttle", "5", "-throttle_every", "10"],
             expected={"nxl": 20, "nyl": 50, "xsearch": 128, "ysearch": 256,
                       "ri": 4, "interp_factor": 8, "snr_thr": 18, "psnr_thr": 6,
                       "throttle_ms": 5, "throttle_every": 10}),
        Case("override_flags", ["-range_interp", "4", "-interp", "8", "-norange",
             "-nointerp", "-nyquist_split", "-geocode", "-no_blockmedian"],
             expected={"ri": 1, "interp_factor": 1, "nyquist_split": 1,
                       "do_geocode": 1, "do_blockmedian": 0}),
        Case("no_geocode_wins", ["-geocode", "-no_geocode"], expected={"do_geocode": 0}),
        Case("noshift_omitted_keys", ["-noshift"], secondary=prm(MASTER),
             expected={"x_offset": 0, "y_offset": 0}),
        Case("no_trailing_newline", master=prm(MASTER, final=False),
             secondary=prm(SECONDARY, final=False), expected={"y_offset": 3}),
        Case("crlf", master=prm(MASTER, "\r\n"), secondary=prm(SECONDARY, "\r\n")),
        Case("blank_comments", master=b" \t\n# ignored = value\n" + prm(MASTER) + b"\n \t\n"),
        Case("duplicate_last_assignment", master=prm(MASTER) + b"num_rng_bins = 8192\n",
             expected={"m_nx": 8192}),
        Case("empty_unrelated", master=prm(MASTER) + b"unrelated = \n"),
        Case("help", argv=["-help"], help_only=True),
        Case("help_double_dash", argv=["--help"], help_only=True),
        Case("no_args_help", argv=[], help_only=True),
        Case("double_dash_paths", argv=["--", "-master.PRM", "-secondary.PRM"]),
        Case("extra_positional", argv=["master.PRM", "secondary.PRM", "extra.PRM"],
             error_token="PRM"),
        Case("missing_secondary", argv=["master.PRM"], error_token="PRM"),
        Case("unknown_option", ["-bogus"], error_token="bogus"),
        Case("unknown_option_value", ["-bogus=2"], error_token="bogus"),
        Case("flag_with_argument", ["-noshift=1"], error_token="noshift"),
        Case("invalid_device", ["-af", "metal"], error_token="af"),
        Case("missing_device_value", ["-af"], error_token="af"),
        Case("missing_file", argv=["absent.PRM", "secondary.PRM"], error_token="absent.PRM"),
        Case("master_read_error", argv=[".", "secondary.PRM"], error_token="PRM"),
        Case("secondary_read_error_cleanup", argv=["master.PRM", "."], error_token="PRM"),
        Case("empty_key", master=b"= value\n" + prm(MASTER), error_token="PRM"),
        Case("whitespace_key", master=b"  \t = value\n" + prm(MASTER), error_token="PRM"),
        Case("malformed_line", master=prm(MASTER) + b"not an assignment\n", error_token="PRM"),
        Case("embedded_nul", master=prm(MASTER) + b"junk = a\x00b\n", error_token="PRM"),
        Case("empty_file", master=b"", error_token="SLC_file"),
        Case("dimension_product_overflow", master=prm({**MASTER, "num_patches": "2147483647"}),
             error_token="PRM"),
        Case("sampling_product_overflow", ["-nx", "2147483647", "-ny", "2"],
             error_token="nx"),
        Case("ny_sampling_overflow", ["-ny", "2147483647"], error_token="ny"),
        Case("interpolation_size_overflow", ["-interp", "2147483647"], error_token="interp"),
        Case("range_interpolation_size_overflow", ["-range_interp", "1073741824"],
             error_token="range_interp"),
        Case("xsearch_size_overflow", ["-xsearch", "1073741824"], error_token="xsearch"),
        Case("ysearch_size_overflow", ["-ysearch", "1073741824"], error_token="ysearch"),
        Case("prf_ratio_overflow", master=prm({**MASTER, "PRF": "1e-200"}),
             secondary=prm({**SECONDARY, "PRF": "1e200"}), error_token="PRF"),
        Case("prf_ratio_cancellation", master=prm({**MASTER, "PRF": "1e200"}),
             secondary=prm({**SECONDARY, "PRF": "1e-200"}), error_token="PRF"),
        Case("large_streamable_image", master=prm({**MASTER, "num_valid_az": "10000000"}),
             expected={"m_ny": 20000000}),
        Case("master_row_cache_overflow", master=prm({**MASTER, "num_rng_bins": "10000000"}),
             error_token="num_rng_bins"),
        Case("secondary_row_cache_overflow", secondary=prm({**SECONDARY, "num_rng_bins": "10000000"}),
             error_token="num_rng_bins"),
    ]
    # These values are parsed but never used to read or allocate any SLC data.
    long_path = "directory with spaces/" + "x" * 400 + " = scene.SLC"
    result.append(Case("long_path_spaces_equals", master=prm({**MASTER, "SLC_file": long_path}),
                       expected={"m_path": long_path}))
    result.append(Case("long_unknown_value", master=prm(MASTER) + b"metadata = " + b"a" * 4096 + b"\n"))
    for device, value in (("cuda", 1), ("opencl", 2), ("cpu", 3)):
        result.append(Case(f"device_{device}", ["-af", device], expected={"device": value}))
    numeric = ("nx", "ny", "xsearch", "ysearch", "range_interp", "interp", "snr", "psnr",
               "throttle", "throttle_every")
    positive = set(numeric) - {"snr", "psnr", "throttle"}
    invalid = ("", "abc", "1.5", "1x", "+", "+-1", "-1", " 1", "1 ", "2147483648", "4294967296",
               "999999999999999999999999999999999")
    for option in numeric:
        result.append(Case(f"missing_value_{option}", [f"-{option}"], error_token=option))
        for index, value in enumerate(invalid):
            result.append(Case(f"numeric_{option}_{index}", [f"-{option}", value], error_token=option))
        if option in positive:
            result.append(Case(f"zero_{option}", [f"-{option}", "0"], error_token=option))
    for option in ("xsearch", "ysearch", "range_interp"):
        result.append(Case(f"non_power_two_{option}", [f"-{option}", "3"], error_token=option))
    for option in ("xsearch", "ysearch"):
        for value in (1, 2, 4):
            result.append(Case(f"small_subpixel_{option}_{value}", [f"-{option}", str(value)],
                               error_token=option))
            result.append(Case(f"small_integer_{option}_{value}", [f"-{option}", str(value), "-nointerp"],
                               expected={option: value, "interp_factor": 1}))
        result.append(Case(f"minimum_subpixel_{option}", [f"-{option}", "8"], expected={option: 8}))
    for label, source in (("master", MASTER), ("secondary", SECONDARY)):
        for key in source:
            changed = dict(source)
            del changed[key]
            values = {label: prm(changed)}
            result.append(Case(f"missing_{label}_{key}", error_token=key, **values))
            result.append(Case(f"empty_{label}_{key}", error_token=key,
                               **{label: prm({**source, key: " \t "})}))
    for key in ("num_rng_bins", "num_patches", "num_valid_az", "rshift", "ashift"):
        for index, value in enumerate(("abc", "1junk", "1.2", "2147483648", "-2147483649", "1e2")):
            result.append(Case(f"integer_{key}_{index}", secondary=prm({**SECONDARY, key: value}),
                               error_token=key))
        if key.startswith("num_"):
            for value in ("0", "-1"):
                result.append(Case(f"nonpositive_{key}_{value}", secondary=prm({**SECONDARY, key: value}),
                                   error_token=key))
    for index, value in enumerate(("abc", "1x", "nan", "NaN", "inf", "-inf", "1e9999", "1e-9999", "0", "-1")):
        result.append(Case(f"float_PRF_{index}", master=prm({**MASTER, "PRF": value}), error_token="PRF"))
    # Deterministic bounded malformed corpus: strict conversions must reject the
    # entire token; each generated token has a digit prefix plus a nonnumeric suffix.
    rng = random.Random(20260924)
    for index in range(32):
        option = rng.choice(numeric)
        value = str(rng.randrange(1, 1000000)) + rng.choice(["junk", "_", "x42", "..", "=3"])
        result.append(Case(f"corpus_cli_{index}", [f"-{option}", value], error_token=option))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--work", type=Path, required=True,
                        help="Isolated build and fixture directory; existing unrelated files are not deleted")
    parser.add_argument("--sanitize", choices=("none", "undefined", "address,undefined"), default="undefined")
    parser.add_argument("--cc", default="gcc")
    args = parser.parse_args()
    source, work = args.source.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    snapshot = work / "source"
    snapshot.mkdir(exist_ok=True)
    hashes = {}
    for name in ("xcorr2_args.c", "xcorr2_args.h", "prm_helper.c", "xcorr2.h"):
        data = (source / name).read_bytes()
        hashes[name] = hashlib.sha256(data).hexdigest()
        (snapshot / name).write_bytes(data)
    shutil.copyfile(Path(__file__).with_name("parser_harness.c"), snapshot / "parser_harness.c")
    flags = shlex.split(subprocess.check_output(["pkg-config", "--cflags", "--libs", "glib-2.0"], text=True))
    binary = work / "parser_harness"
    command = [args.cc, "-std=gnu11", "-O1", "-g", "-Wall", "-Wextra", "-Werror", "-fno-omit-frame-pointer"]
    if args.sanitize != "none":
        command += [f"-fsanitize={args.sanitize}", "-fno-sanitize-recover=all"]
    # Non-PIE avoids ASan shadow-map startup collisions on some WSL hosts.
    if "address" in args.sanitize:
        command += ["-no-pie"]
    command += [str(snapshot / name) for name in ("parser_harness.c", "xcorr2_args.c", "prm_helper.c")]
    command += ["-o", str(binary), *flags, "-lm"]
    subprocess.run(command, check=True)
    environment = {**os.environ, "LC_ALL": "C.UTF-8", "ASAN_OPTIONS": "detect_leaks=1:halt_on_error=1",
                   "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1"}
    records = []
    for case in cases():
        case_dir = work / "cases" / case.name
        case_dir.mkdir(parents=True, exist_ok=True)
        for name, data in (("master.PRM", case.master), ("secondary.PRM", case.secondary),
                           ("-master.PRM", case.master), ("-secondary.PRM", case.secondary)):
            (case_dir / name).write_bytes(data)
        argv = case.argv if case.argv is not None else ["master.PRM", "secondary.PRM", *case.options]
        failures = []
        try:
            completed = subprocess.run([str(binary), *argv], cwd=case_dir, env=environment,
                                       text=True, errors="replace", capture_output=True, timeout=5)
            if completed.returncode < 0:
                failures.append(f"terminated by signal {-completed.returncode}")
            if SANITIZER.search(completed.stderr):
                failures.append("sanitizer diagnostic")
            if case.error_token is not None:
                if completed.returncode == 0:
                    failures.append("invalid input returned success")
                if case.error_token.casefold() not in completed.stderr.casefold():
                    failures.append(f"diagnostic does not identify {case.error_token!r}")
            else:
                if completed.returncode != 0:
                    failures.append(f"valid input exited {completed.returncode}")
                elif case.help_only:
                    if "xcorr_cc" not in completed.stdout:
                        failures.append("missing help output")
                else:
                    try:
                        normalized = json.loads(completed.stdout)
                        for key, value in case.expected.items():
                            if normalized.get(key) != value:
                                failures.append(f"{key}: expected {value!r}, got {normalized.get(key)!r}")
                    except (ValueError, TypeError) as error:
                        failures.append(f"invalid normalized JSON: {error}")
            records.append({"name": case.name, "pass": not failures, "failures": failures,
                            "returncode": completed.returncode, "stderr": completed.stderr,
                            "stdout": completed.stdout, "argv": argv})
        except subprocess.TimeoutExpired:
            records.append({"name": case.name, "pass": False, "failures": ["5-second timeout"], "argv": argv})
    failed = [record for record in records if not record["pass"]]
    report = {"pass": not failed, "total": len(records), "passed": len(records) - len(failed),
              "failed": len(failed), "sanitize": args.sanitize, "source_sha256": hashes,
              "compile_command": command, "cases": records}
    (work / "results.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"{report['passed']}/{len(records)} passed; {len(failed)} failed ({args.sanitize})")
    for record in failed:
        print(f"FAIL {record['name']}: {'; '.join(record['failures'])}")
    print(f"Detailed results: {work / 'results.json'}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
