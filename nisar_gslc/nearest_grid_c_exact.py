#!/usr/bin/env python3
# =============================================================================
# Author       : Xin Wang
# Affiliation  : University of Science and Technology of China (USTC)
# Location     : Hefei, China
# Email        : xinw11@mail.ustc.edu.cn
# GitHub       : https://github.com/Xin-Wang-520
# Uploaded     : 2026-08-01
# =============================================================================

"""Fast numerical reproduction of GMTSAR nearest_grid.

Usage:
    nearest_grid_c_exact.py input.grd output.grd [search_radius]

The output grid reproduces the original C program's numerical behavior:

* Existing non-NaN samples are unchanged.
* NaNs are filled from the Euclidean-nearest non-NaN sample.
* Equidistant candidates use the same directional order as find_nearest()
  in GMTSAR's nearest_grid.c.
* A positive radius reproduces the C loop's boundary behavior, which checks
  the first squared-distance ring beyond radius**2.  On an integer grid that
  ring is radius**2 + 1.

The NetCDF container is copied from the input, so coordinate variables and
grid metadata are retained.  File bytes and history metadata need not match
the GMT-written C output, but z values and NaN locations are intended to.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError as exc:
    raise SystemExit(
        "ERROR: netCDF4 is required. Install it with: "
        "python3 -m pip install netCDF4"
    ) from exc

try:
    from scipy.ndimage import distance_transform_edt
except ImportError as exc:
    raise SystemExit(
        "ERROR: SciPy is required. Install it with: "
        "python3 -m pip install scipy"
    ) from exc


class FillError(RuntimeError):
    """Invalid input or a failed grid operation."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "快速复现 GMTSAR nearest_grid 的最近邻和等距离选择顺序。"
        )
    )
    parser.add_argument("input", type=Path, help="输入 GMT/NetCDF 网格")
    parser.add_argument("output", type=Path, help="输出网格")
    parser.add_argument(
        "search_radius",
        nargs="?",
        type=int,
        default=0,
        help="搜索半径（像元）；0与原C程序相同，表示覆盖整幅网格",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="允许覆盖已经存在的输出文件",
    )
    parser.add_argument(
        "--block-rows",
        type=int,
        default=512,
        help="等距离选择的分块行数，默认512",
    )
    return parser.parse_args()


def find_grid_variable(dataset: Dataset) -> str:
    if "z" in dataset.variables and dataset.variables["z"].ndim == 2:
        return "z"

    candidates = [
        name
        for name, variable in dataset.variables.items()
        if variable.ndim == 2
        and np.issubdtype(variable.dtype, np.number)
    ]
    if len(candidates) != 1:
        found = ", ".join(candidates) if candidates else "无"
        raise FillError(f"无法唯一确定二维网格变量；找到：{found}")
    return candidates[0]


def c_ring_offsets(distance_squared: int) -> tuple[tuple[int, int], ...]:
    """Return one squared-distance ring in nearest_grid.c candidate order.

    GMT's in-memory grid row 0 is the northern row, while the CF/NetCDF
    variable read by netCDF4 is ordered from the lower y coordinate upward.
    Consequently, the C program's positive ``i`` offset is a negative NumPy
    row offset.  Column signs are unchanged.
    """

    if distance_squared <= 0:
        return ()

    offsets: list[tuple[int, int]] = []
    maximum_x = math.isqrt(distance_squared)
    minimum_x = math.isqrt(distance_squared // 2)

    for x_offset in range(minimum_x, maximum_x + 1):
        y_squared = distance_squared - x_offset * x_offset
        if y_squared < 0:
            continue

        y_offset = math.isqrt(y_squared)
        if y_offset * y_offset != y_squared:
            continue
        if y_offset > x_offset:
            continue

        if y_offset == 0:
            offsets.extend(
                (
                    (-x_offset, 0),
                    (x_offset, 0),
                    (0, x_offset),
                    (0, -x_offset),
                )
            )
        elif x_offset != y_offset:
            offsets.extend(
                (
                    (-x_offset, y_offset),
                    (x_offset, y_offset),
                    (-x_offset, -y_offset),
                    (x_offset, -y_offset),
                    (-y_offset, x_offset),
                    (-y_offset, -x_offset),
                    (y_offset, x_offset),
                    (y_offset, -x_offset),
                )
            )
        else:
            offsets.extend(
                (
                    (-x_offset, x_offset),
                    (x_offset, x_offset),
                    (-x_offset, -x_offset),
                    (x_offset, -x_offset),
                )
            )

    if not offsets:
        raise FillError(
            f"内部错误：平方距离 {distance_squared} 没有整数网格偏移"
        )
    return tuple(offsets)


def fill_one_distance_group(
    source: np.ndarray,
    output: np.ndarray,
    rows: np.ndarray,
    columns: np.ndarray,
    distance_squared: int,
    offsets: tuple[tuple[int, int], ...],
) -> int:
    """Fill a group using exactly the C candidate order."""

    unresolved = np.ones(rows.size, dtype=bool)
    n_rows, n_columns = source.shape
    filled = 0

    for row_offset, column_offset in offsets:
        if not np.any(unresolved):
            break

        candidate_rows = rows + row_offset
        candidate_columns = columns + column_offset
        inside = (
            unresolved
            & (candidate_rows >= 0)
            & (candidate_rows < n_rows)
            & (candidate_columns >= 0)
            & (candidate_columns < n_columns)
        )
        if not np.any(inside):
            continue

        inside_indices = np.flatnonzero(inside)
        valid_indices = inside_indices[
            ~np.isnan(
                source[
                    candidate_rows[inside_indices],
                    candidate_columns[inside_indices],
                ]
            )
        ]
        if valid_indices.size == 0:
            continue

        output[rows[valid_indices], columns[valid_indices]] = source[
            candidate_rows[valid_indices],
            candidate_columns[valid_indices],
        ]
        unresolved[valid_indices] = False
        filled += int(valid_indices.size)

    if np.any(unresolved):
        raise FillError(
            "内部错误：距离变换给出的最近距离没有找到有效候选点；"
            f"distance_squared={distance_squared}, "
            f"unresolved={int(np.count_nonzero(unresolved))}"
        )

    return filled


def reproduce_c_nearest(
    data: np.ndarray,
    search_radius: int,
    block_rows: int,
) -> tuple[np.ndarray, int, int]:
    if data.ndim != 2:
        raise FillError(f"只支持二维网格，当前维数：{data.ndim}")
    if search_radius < 0:
        raise FillError("search_radius不能小于0")
    if block_rows <= 0:
        raise FillError("block_rows必须大于0")

    # The C source tests isnan(), not isfinite().  Infinity is therefore
    # treated as an existing sample here as well.
    invalid = np.isnan(data)
    invalid_count = int(np.count_nonzero(invalid))
    if invalid_count == data.size:
        raise FillError("输入网格没有任何非NaN像元")
    if invalid_count == 0:
        return data, 0, 0

    # EDT supplies only the exact minimum Euclidean distance.  It does not
    # choose the value, because SciPy and nearest_grid.c break ties differently.
    distances = distance_transform_edt(invalid)

    if search_radius == 0:
        maximum_distance_squared: int | None = None
    else:
        # nearest_grid.c enters find_nearest while the previous ring is <= R^2,
        # then also evaluates the next representable integer-grid ring.
        maximum_distance_squared = search_radius * search_radius + 1

    # Keep source and output separate.  The C program always searches its
    # original input array and never uses a value filled earlier in the scan.
    output = data.copy()
    filled_count = 0
    ring_cache: dict[int, tuple[tuple[int, int], ...]] = {}
    n_rows = data.shape[0]

    for row0 in range(0, n_rows, block_rows):
        row1 = min(row0 + block_rows, n_rows)
        block_invalid = invalid[row0:row1]
        if not np.any(block_invalid):
            continue

        block_distance_squared = np.rint(
            distances[row0:row1] * distances[row0:row1]
        ).astype(np.int64)

        fillable = block_invalid.copy()
        if maximum_distance_squared is not None:
            fillable &= block_distance_squared <= maximum_distance_squared
        if not np.any(fillable):
            continue

        local_rows, columns = np.nonzero(fillable)
        rows = local_rows + row0
        keys = block_distance_squared[local_rows, columns]

        order = np.argsort(keys, kind="stable")
        rows = rows[order]
        columns = columns[order]
        keys = keys[order]

        starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1]) + 1]
        ends = np.r_[starts[1:], keys.size]

        for start, end in zip(starts, ends):
            distance_squared = int(keys[start])
            offsets = ring_cache.get(distance_squared)
            if offsets is None:
                offsets = c_ring_offsets(distance_squared)
                ring_cache[distance_squared] = offsets

            filled_count += fill_one_distance_group(
                data,
                output,
                rows[start:end],
                columns[start:end],
                distance_squared,
                offsets,
            )

    return output, invalid_count, filled_count


def fill_grid(
    input_path: Path,
    output_path: Path,
    search_radius: int,
    *,
    overwrite: bool,
    block_rows: int,
) -> None:
    input_path = input_path.resolve()
    output_path = output_path.resolve()

    if not input_path.is_file():
        raise FillError(f"输入网格不存在：{input_path}")
    if input_path == output_path:
        raise FillError("输入和输出不能是同一个文件")
    if output_path.exists() and not overwrite:
        raise FillError(f"输出已存在，拒绝覆盖：{output_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(
        f".{output_path.name}.tmp-{os.getpid()}"
    )

    total_started = time.perf_counter()
    with Dataset(input_path, "r") as source_dataset:
        grid_name = find_grid_variable(source_dataset)
        raw = source_dataset.variables[grid_name][:]
        if np.ma.isMaskedArray(raw):
            data = np.asarray(raw.filled(np.nan))
        else:
            data = np.asarray(raw).copy()

    calculation_started = time.perf_counter()
    output, invalid_count, filled_count = reproduce_c_nearest(
        data,
        search_radius,
        block_rows,
    )
    calculation_elapsed = time.perf_counter() - calculation_started

    try:
        if temporary_path.exists():
            temporary_path.unlink()
        shutil.copy2(input_path, temporary_path)
        with Dataset(temporary_path, "r+") as target_dataset:
            target_dataset.variables[grid_name][:] = output
        os.replace(temporary_path, output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    total_elapsed = time.perf_counter() - total_started
    print(
        "C-exact nearest-neighbour grid written: "
        f"{output_path} "
        f"(radius={search_radius}, input_nan={invalid_count}, "
        f"filled={filled_count}, calculation={calculation_elapsed:.2f}s, "
        f"total={total_elapsed:.2f}s)",
        flush=True,
    )


def main() -> int:
    args = parse_args()
    fill_grid(
        args.input,
        args.output,
        args.search_radius,
        overwrite=args.overwrite,
        block_rows=args.block_rows,
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except FillError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
