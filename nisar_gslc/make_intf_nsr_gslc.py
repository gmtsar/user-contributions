#!/usr/bin/env python3
# GSLC workflow adaptation and maintenance: Xin Wang
# Affiliation: University of Science and Technology of China (USTC), Hefei, China
# Contact: xinw11@mail.ustc.edu.cn | https://github.com/Xin-Wang-520
# Attribution and license: AUTHORS.md, UPSTREAM_PROVENANCE.md, LICENSE
"""Independent NISAR GSLC pair processing in the product's map coordinates.

No GMTSAR installation files are changed. The phase convention is
reference * conjugate(secondary). See README.md for scientific limitations.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import logging
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import h5py
import netCDF4
import numpy as np
from pyproj import CRS
from scipy.ndimage import convolve1d
from scipy.signal import fftconvolve

LOG = logging.getLogger("nisar_gslc")
C = 299792458.0
VERSION = "0.6.5"
POL_BITS = {p: 1 << (8 + i) for i, p in enumerate(
    ("HH", "HV", "VH", "VV", "LH", "LV", "RH", "RV"))}


def scalar(value):
    """Decode small HDF5 metadata without evaluating its contents."""
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace").rstrip("\x00")
    if isinstance(value, np.ndarray):
        return scalar(value.item()) if value.size == 1 else [scalar(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def metadata_value(h5, path, default=None):
    return scalar(h5[path][()]) if path in h5 else default


def bool_value(value):
    if value is None:
        return None
    if str(value).strip().lower() in ("true", "1"):
        return True
    if str(value).strip().lower() in ("false", "0"):
        return False
    return None


def warn(report, message):
    if message not in report["warnings"]:
        report["warnings"].append(message)
        LOG.warning(message)


def regular_axis(values, name):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)):
        raise ValueError(f"{name}: expected at least two finite pixel-center coordinates")
    step = float(values[1] - values[0])
    if step == 0 or not np.allclose(np.diff(values), step, rtol=1e-8, atol=abs(step)*1e-6):
        raise ValueError(f"{name}: grid must be monotonic and regularly spaced")
    return values, step


def complex_fields(dataset):
    typ = dataset.id.get_type()
    if typ.get_class() != h5py.h5t.COMPOUND or typ.get_nmembers() != 2:
        raise ValueError(f"{dataset.name}: expected complex samples, not detected power")
    names = [typ.get_member_name(i).decode() for i in range(2)]
    for real, imag in (("r", "i"), ("real", "imag"), ("real", "imaginary")):
        if set(names) == {real, imag}:
            for i in range(2):
                if typ.get_member_type(i).get_class() != h5py.h5t.FLOAT:
                    raise ValueError(f"{dataset.name}: complex components must be floating point")
            return real, imag
    raise ValueError(f"{dataset.name}: unknown complex component names {names}")


def read_complex(dataset, start, shape, fields):
    """HDF5 converts each component to f32; also handles NISAR complex32/f16."""
    dtype = np.dtype([(fields[0], "f4"), (fields[1], "f4")])
    buf = np.empty(shape, dtype=dtype)
    space = dataset.id.get_space()
    space.select_hyperslab(start, shape)
    memspace = h5py.h5s.create_simple(shape)
    dataset.id.read(memspace, space, buf, mtype=h5py.h5t.py_create(dtype))
    return buf[fields[0]] + np.complex64(1j) * buf[fields[1]]


class Product:
    def __init__(self, h5, path, frequency, polarization):
        self.h5, self.path = h5, Path(path).resolve()
        self.polarization = polarization
        self.root = "/science/LSAR/GSLC"
        self.group_path = f"{self.root}/grids/frequency{frequency}"
        if self.group_path not in h5 or polarization not in h5[self.group_path]:
            raise ValueError(f"{path}: missing {self.group_path}/{polarization}")
        self.group = h5[self.group_path]
        self.data = self.group[polarization]
        self.fields = complex_fields(self.data)
        self.x, sx = regular_axis(self.group["xCoordinates"][:], "xCoordinates")
        self.y, sy = regular_axis(self.group["yCoordinates"][:], "yCoordinates")
        if self.data.shape != (len(self.y), len(self.x)):
            raise ValueError(f"{path}: sample shape does not match coordinate vectors")
        self.shape = self.data.shape
        self.flip_x, self.flip_y = sx < 0, sy > 0
        if self.flip_x:
            self.x = self.x[::-1]
        if self.flip_y:
            self.y = self.y[::-1]
        self.dx, self.dy = abs(sx), -abs(sy)
        for key, step in (("xCoordinateSpacing", sx), ("yCoordinateSpacing", sy)):
            if key in self.group and not np.isclose(float(self.group[key][()]), step, rtol=1e-6):
                raise ValueError(f"{path}: {key} conflicts with pixel-center coordinates")
        self.epsg = int(self.group["projection"][()])
        self.crs = CRS.from_epsg(self.epsg)
        if not self.crs.is_projected or any(abs(a.unit_conversion_factor-1) > 1e-10 for a in self.crs.axis_info[:2]):
            raise ValueError("First release requires a projected grid in meters (UTM or polar stereographic)")
        self.frequency_hz = float(self.group["centerFrequency"][()])
        if not np.isfinite(self.frequency_hz) or self.frequency_hz <= 0:
            raise ValueError("Invalid centerFrequency")
        self.mask = self.group.get("mask")
        self.exception = self.group.get("inputDataExceptionMask")
        if self.mask is not None and self.mask.shape != self.shape:
            raise ValueError(f"{path}: mask shape {self.mask.shape} differs from SLC {self.shape}")
        self.new_exception = False
        if self.exception is not None:
            desc = scalar(self.exception.attrs.get("description", b""))
            self.new_exception = "validity mask bits" in desc.lower()
            if self.exception.dtype.itemsize > 1 and not self.new_exception:
                raise ValueError("Unknown multi-byte inputDataExceptionMask encoding")
        self.exception_mapping = None
        if self.exception is not None and self.exception.shape != self.shape:
            self.exception_mapping = self.legacy_exception_mapping(frequency)
        identification = "/science/LSAR/identification"
        self.meta = {k: metadata_value(h5, f"{identification}/{k}") for k in (
            "productType", "isGeocoded", "trackNumber", "frameNumber", "orbitPassDirection",
            "lookDirection", "zeroDopplerStartTime", "zeroDopplerEndTime", "granuleId",
            "compositeReleaseId", "productVersion", "productSpecificationVersion", "isMixedMode")}
        params = f"{self.root}/metadata/processingInformation/parameters"
        self.meta["corrections"] = {k: metadata_value(h5, f"{params}/{k}") for k in (
            "ellipsoidalFlatteningApplied", "topographicFlatteningApplied",
            "rangeIonosphericGeolocationCorrectionApplied", "azimuthIonosphericGeolocationCorrectionApplied",
            "dryTroposphericGeolocationCorrectionApplied", "wetTroposphericGeolocationCorrectionApplied")}
        self.meta["demSource"] = metadata_value(h5, f"{self.root}/metadata/processingInformation/inputs/demSource")
        self.meta.update(path=str(self.path), dataset=self.data.name, shape=list(self.shape), epsg=self.epsg,
                         spacing_m=[self.dx, self.dy], center_frequency_hz=self.frequency_hz,
                         exception_mask_encoding="polarization_validity_bits" if self.new_exception else "legacy_anomaly_codes",
                         exception_mask_grid_mapping=self.exception_mapping)

    def legacy_exception_mapping(self, frequency):
        """Recover the known old B LUT on A's geogrid, preserving all anomaly bits.

        ISCE3 v0.25.16 GslcWriter omitted frequency= in geocode_lut; its
        default selects A's geogrid. This is a quality-mask operation only,
        never an interpolation of complex SLC samples.
        """
        message = (f"{self.path}: inputDataExceptionMask shape {self.exception.shape} "
                   f"differs from SLC {self.shape}; not a verified legacy B-on-A grid")
        a = self.h5.get(f"{self.root}/grids/frequencyA")
        if frequency != "B" or self.new_exception or self.exception.dtype != np.dtype("u1") or a is None:
            raise ValueError(message)
        ax, dx = regular_axis(a["xCoordinates"][:], "legacy mask xCoordinates")
        ay, dy = regular_axis(a["yCoordinates"][:], "legacy mask yCoordinates")
        if self.exception.shape != (len(ay), len(ax)) or int(a["projection"][()]) != self.epsg:
            raise ValueError(message)
        fx, fy = dx < 0, dy > 0
        ax, ay = ax[::-1] if fx else ax, ay[::-1] if fy else ay
        dx, dy = abs(dx), abs(dy)
        ratio = np.array([-self.dy/dy, self.dx/dx])
        if np.any(ratio < 1) or not np.allclose(ratio, np.rint(ratio), rtol=0, atol=1e-8):
            raise ValueError(message + "; non-integer footprint ratio")
        ry, rx = (int(v) for v in np.rint(ratio))
        # Require matching full pixel-edge extents, not just matching dimensions.
        edges_a = [ax[0]-dx/2, ax[-1]+dx/2, ay[-1]-dy/2, ay[0]+dy/2]
        edges_b = [self.x[0]-self.dx/2, self.x[-1]+self.dx/2,
                   self.y[-1]+self.dy/2, self.y[0]-self.dy/2]
        if ((len(ay), len(ax)) != (self.shape[0]*ry, self.shape[1]*rx)
                or not np.allclose(edges_a, edges_b, rtol=0, atol=1e-5)):
            raise ValueError(message + "; pixel-edge extents do not match")
        return {"source_grid": "frequencyA", "target_grid": "frequencyB",
                "method": "bitwise_OR_of_all_fine_pixels_in_each_B_footprint",
                "factors_yx": [ry, rx], "source_flip_yx": [bool(fy), bool(fx)],
                "source_shape": list(self.exception.shape)}

    def read_mapped_exception(self, row, col, height, width):
        mapping = self.exception_mapping
        ry, rx = mapping["factors_yx"]
        fy, fx = mapping["source_flip_yx"]
        h, w = height*ry, width*rx
        r, c = row*ry, col*rx
        if fy:
            r = self.exception.shape[0]-r-h
        if fx:
            c = self.exception.shape[1]-c-w
        e = self.exception[r:r+h, c:c+w]
        if fy:
            e = e[::-1]
        if fx:
            e = e[:, ::-1]
        e = e.reshape(height, ry, width, rx)
        return np.bitwise_or.reduce(np.bitwise_or.reduce(e, axis=3), axis=1)

    def read(self, row, col, height, width):
        r = self.shape[0]-row-height if self.flip_y else row
        c = self.shape[1]-col-width if self.flip_x else col
        data = read_complex(self.data, (r, c), (height, width), self.fields)
        valid = np.isfinite(data) & (np.abs(data) > 0)
        slices = np.s_[r:r+height, c:c+width]
        if self.mask is not None:
            m = self.mask[slices]
            valid &= (m > 0) & (m < 255)
        if self.exception is not None and self.exception_mapping is None:
            e = self.exception[slices].astype(np.uint32)
            if self.new_exception:
                valid &= ((e & 255) == 0) & ((e & POL_BITS[self.polarization]) != 0)
            else:
                valid &= e == 0
        if self.flip_y:
            data, valid = data[::-1], valid[::-1]
        if self.flip_x:
            data, valid = data[:, ::-1], valid[:, ::-1]
        if self.exception_mapping is not None:
            valid &= self.read_mapped_exception(row, col, height, width) == 0
        return data, valid


@dataclass
class Region:
    row: int
    col: int
    height: int
    width: int
    secondary_row_offset: int
    secondary_col_offset: int


def pair_region(ref, sec, args, report):
    if ref.path == sec.path:
        raise ValueError("Reference and secondary must be different files")
    for product in (ref, sec):
        if product.meta["productType"] != "GSLC" or bool_value(product.meta["isGeocoded"]) is not True:
            raise ValueError(f"{product.path}: input must be an identified geocoded GSLC product")
        for key in ("ellipsoidalFlatteningApplied", "topographicFlatteningApplied"):
            value = bool_value(product.meta["corrections"][key])
            if value is False or (value is None and not args.assume_flattened):
                raise ValueError(f"{product.path}: {key} is {value}; this workflow requires flattened GSLCs. "
                                 "--assume-flattened applies only to missing flags, never explicit False.")
            if value is None:
                warn(report, f"User assumed missing {key}=True for {product.path.name}")
        if product.mask is None:
            warn(report, f"No subswath mask in {product.path.name}; using finite, nonzero samples only")
        if product.exception_mapping is not None:
            warn(report, f"Legacy B exception mask in {product.path.name} uses A's geogrid; "
                 "all anomaly bits within each B pixel footprint are conservatively OR-combined. "
                 "Recovery is inferred from ISCE3 writer behavior and verified matching grid extents.")
    for key in ("trackNumber", "frameNumber", "orbitPassDirection", "lookDirection"):
        if ref.meta[key] is None or sec.meta[key] is None:
            warn(report, f"Cannot verify {key}: metadata missing")
        elif ref.meta[key] != sec.meta[key]:
            raise ValueError(f"Incompatible {key}: {ref.meta[key]} vs {sec.meta[key]}")
    for key in ("compositeReleaseId", "productSpecificationVersion", "productVersion"):
        if ref.meta[key] != sec.meta[key]:
            warn(report, f"Different {key}: {ref.meta[key]} vs {sec.meta[key]}")
    if not ref.meta["demSource"] or not sec.meta["demSource"]:
        warn(report, "Flattening DEM provenance missing; a common DEM cannot be verified")
    elif ref.meta["demSource"] != sec.meta["demSource"]:
        warn(report, "Different flattening DEM provenance: residual topographic/geolocation effects may remain. "
             "Outputs are diagnostic until the DEM difference is evaluated or harmonized.")
    if ref.meta["corrections"] != sec.meta["corrections"]:
        warn(report, "The two products have different processing correction flags")
    if ref.meta["zeroDopplerStartTime"] and sec.meta["zeroDopplerStartTime"]:
        if ref.meta["zeroDopplerStartTime"] >= sec.meta["zeroDopplerStartTime"]:
            warn(report, "Reference is not earlier than secondary; phase/LOS signs follow the supplied order")
    if ref.epsg != sec.epsg or not np.allclose([ref.dx, ref.dy], [sec.dx, sec.dy], rtol=1e-8, atol=1e-8):
        raise ValueError("Different projections or pixel spacings; phase-preserving resampling is required")
    if not np.isclose(ref.frequency_hz, sec.frequency_hz, rtol=1e-8):
        raise ValueError("Different center frequencies; common-band processing is required")
    offsets = np.array([(ref.y[0]-sec.y[0])/ref.dy, (ref.x[0]-sec.x[0])/ref.dx])
    if not np.allclose(offsets, np.rint(offsets), atol=1e-5, rtol=0):
        raise ValueError("Fractional-pixel grid offset: no interpolation is performed in this release")
    ro, co = (int(v) for v in np.rint(offsets))
    r0, c0 = max(0, -ro), max(0, -co)
    r1, c1 = min(len(ref.y), len(sec.y)-ro), min(len(ref.x), len(sec.x)-co)
    if args.window:
        r, c, h, w = args.window
        if h <= 0 or w <= 0 or r < r0 or c < c0 or r+h > r1 or c+w > c1:
            raise ValueError("--window must lie entirely within the common coverage (north-up reference indices)")
        r0, c0, r1, c1 = r, c, r+h, c+w
    if args.bbox:
        xmin, xmax, ymin, ymax = args.bbox
        if not (xmin < xmax and ymin < ymax):
            raise ValueError("--bbox order is XMIN XMAX YMIN YMAX in product CRS meters")
        c0 = max(c0, int(np.searchsorted(ref.x, xmin)))
        c1 = min(c1, int(np.searchsorted(ref.x, xmax, side="right")))
        r0 = max(r0, int(np.searchsorted(-ref.y, -ymax)))
        r1 = min(r1, int(np.searchsorted(-ref.y, -ymin, side="right")))
    ly, lx = args.looks
    h, w = (r1-r0)//ly*ly, (c1-c0)//lx*lx
    if h < ly*2 or w < lx*2:
        raise ValueError("Common region is empty or smaller than two output pixels per axis")
    if h != r1-r0 or w != c1-c0:
        warn(report, "Incomplete look blocks at the south/east edge were dropped")
    return Region(r0, c0, h, w, int(ro), int(co))


def registration_check(ref, sec, region, report):
    """Amplitude NCC diagnostic only; does not move or resample either SLC."""
    size = min(256, region.height, region.width)
    if size < 64:
        warn(report, "ROI too small for amplitude registration diagnostics")
        return {"status": "insufficient_area", "samples": []}
    samples = []
    for fy, fx in ((0.2, 0.2), (0.2, 0.8), (0.5, 0.5), (0.8, 0.2), (0.8, 0.8)):
        row = region.row + int((region.height-size)*fy)
        col = region.col + int((region.width-size)*fx)
        a, va = ref.read(row, col, size, size)
        b, vb = sec.read(row+region.secondary_row_offset, col+region.secondary_col_offset, size, size)
        valid = va & vb
        if valid.mean() < 0.9:
            continue
        window = np.outer(np.hanning(size), np.hanning(size))
        processed = []
        for z in (a, b):
            amp = np.log(np.maximum(abs(z), np.finfo(np.float32).tiny))
            amp = np.where(valid, amp-np.mean(amp[valid]), 0)*window
            processed.append(amp)
        denom = np.linalg.norm(processed[0])*np.linalg.norm(processed[1])
        if denom < 1e-15:
            continue
        corr = fftconvolve(processed[0], processed[1][::-1, ::-1], mode="full")/denom
        radius, center = 8, size-1
        box = corr[center-radius:center+radius+1, center-radius:center+radius+1]
        y, x = np.unravel_index(np.argmax(box), box.shape)
        if y in (0, 2*radius) or x in (0, 2*radius):
            continue
        def peak_offset(v):
            denominator = v[0]-2*v[1]+v[2]
            return 0.0 if abs(denominator) < 1e-12 else float(np.clip(0.5*(v[0]-v[2])/denominator, -0.5, 0.5))
        samples.append({"reference_row": row, "reference_col": col, "ncc": float(box[y, x]),
                        "secondary_shift_y_px": float(y-radius)+peak_offset(box[y-1:y+2, x]),
                        "secondary_shift_x_px": float(x-radius)+peak_offset(box[y, x-1:x+2])})
    reliable = [s for s in samples if s["ncc"] >= 0.15]
    result = {"status": "diagnostic_only", "samples": samples, "reliable_samples": len(reliable),
              "convention": "estimated shift to apply to secondary; native north-up pixel units; not applied"}
    if len(reliable) < 3:
        warn(report, "Amplitude registration check has fewer than three reliable patches; alignment is unverified")
    else:
        shift = [float(np.median([s[k] for s in reliable])) for k in ("secondary_shift_y_px", "secondary_shift_x_px")]
        result["median_secondary_shift_yx_px"] = shift
        if max(abs(v) for v in shift) > 0.25:
            warn(report, "Amplitude check indicates >0.25 native pixel residual shift; precision coregistration may be needed")
    return result


class GridWriter:
    def __init__(self, path, x, y, crs, description, units="1"):
        self.file = netCDF4.Dataset(path, "w", format="NETCDF4_CLASSIC")
        d = self.file
        d.Conventions, d.node_offset = "CF-1.8", np.int32(1)
        d.title, d.history = description, f"NISAR_GSLC {VERSION}; {datetime.now(timezone.utc).isoformat()}"
        d.phase_convention = "reference * conjugate(secondary)"
        d.createDimension("x", len(x)); d.createDimension("y", len(y))
        for name, values in (("x", x), ("y", y)):
            var = d.createVariable(name, "f8", (name,))
            var[:] = values
            var.standard_name = ("longitude" if name == "x" else "latitude") if crs.is_geographic else f"projection_{name}_coordinate"
            var.units = ("degrees_east" if name == "x" else "degrees_north") if crs.is_geographic else "m"
            var.axis = name.upper()
            # GMT/GDAL use actual_range as outer bounds when node_offset=1.
            half_pixel = abs(float(values[1]-values[0]))/2
            var.actual_range = np.array([min(values)-half_pixel, max(values)+half_pixel])
        projection = d.createVariable("crs", "i4")
        projection.assignValue(0)
        projection.setncatts(crs.to_cf())
        projection.spatial_ref = crs.to_wkt()
        projection.epsg_code = f"EPSG:{crs.to_epsg()}"
        dx, dy = float(x[1]-x[0]), float(y[1]-y[0])
        projection.GeoTransform = f"{x[0]-dx/2:.12g} {dx:.12g} 0 {y[0]-dy/2:.12g} 0 {dy:.12g}"
        self.z = d.createVariable("z", "f4", ("y", "x"), fill_value=np.float32(np.nan),
                                  zlib=True, complevel=2, chunksizes=(min(256, len(y)), min(256, len(x))))
        self.z.long_name, self.z.units, self.z.grid_mapping = description, units, "crs"
        self.minimum, self.maximum = math.inf, -math.inf

    def __setitem__(self, key, values):
        self.z[key] = values
        finite = np.asarray(values)[np.isfinite(values)]
        if finite.size:
            self.minimum = min(self.minimum, float(finite.min()))
            self.maximum = max(self.maximum, float(finite.max()))

    def close(self):
        self.z.actual_range = np.array([self.minimum, self.maximum] if self.minimum <= self.maximum
                                       else [np.nan, np.nan], dtype="f4")
        self.file.close()


def write_grid(path, data, x, y, crs, description="", units="1"):
    writer = GridWriter(path, x, y, crs, description, units)
    try:
        writer[:] = data
    finally:
        writer.close()


def read_grid(path):
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        data = np.asarray(ds["z"][:])
        x, y = np.asarray(ds["x"][:]), np.asarray(ds["y"][:])
        if y[0] < y[-1]:
            data, y = data[::-1], y[::-1]
        if x[0] > x[-1]:
            data, x = data[:, ::-1], x[::-1]
        return data, x, y


def blocks(shape, size):
    for r in range(0, shape[0], size):
        for c in range(0, shape[1], size):
            yield r, c, min(size, shape[0]-r), min(size, shape[1]-c)


def multilook(a, b, valid, looks, min_fraction):
    ly, lx = looks
    shape = (a.shape[0]//ly, ly, a.shape[1]//lx, lx)
    a, b = np.where(valid, a, 0), np.where(valid, b, 0)
    count = valid.reshape(shape).sum(axis=(1, 3))
    total = lambda z: z.reshape(shape).sum(axis=(1, 3), dtype=np.float64)
    cross = a * np.conj(b)
    cross_sum = total(cross.real) + 1j*total(cross.imag)
    p1, p2 = total(abs(a)**2), total(abs(b)**2)
    denominator = np.maximum(count, 1)
    good = (count >= math.ceil(ly*lx*min_fraction)) & (count > 0) & (p1 > 0) & (p2 > 0)
    cross_mean = cross_sum/denominator
    corr = np.divide(abs(cross_sum), np.sqrt(p1*p2), out=np.zeros(count.shape), where=(p1*p2)>0)
    out = {"real": cross_mean.real, "imag": cross_mean.imag, "phase": np.angle(cross_mean),
           "corr": np.clip(corr, 0, 1), "amp1": p1/denominator, "amp2": p2/denominator}
    out = {k: np.where(good, v, np.nan).astype("f4") for k, v in out.items()}
    out.update(mask=np.where(good, 1.0, np.nan).astype("f4"), nlooks=count.astype("f4"),
               valid_fraction=(count/(ly*lx)).astype("f4"))
    return out


def gaussian_axis(wavelength, spacing, stride):
    """GMTSAR's sigma=w/5.3 and approximately +/-2 sigma support, in map pixels.

    Even strides put output centers halfway between native pixel centers.
    Use a symmetric even kernel there; odd strides reproduce the odd kernel
    in make_gaussian_filter.c. No phase interpolation is involved.
    """
    sigma = wavelength / 5.3 / abs(spacing)
    n = int(4 * sigma)
    n += 1 if n % 2 == 0 else 0
    if stride % 2 == 0:
        n += 1
    offsets = np.arange(n, dtype=float) - (n-1)/2
    exponent = -0.5*(offsets/sigma)**2
    weights = np.exp(exponent-exponent.max())
    return weights/weights.sum()


def gaussian_plan(args, ref):
    enabled = args.stage in ("filter", "unwrap") and args.filter_wavelength > 0
    result = {"enabled": enabled, "wavelength_m": args.filter_wavelength,
              "definition": "GMTSAR make_gaussian_filter: sigma=wavelength/5.3; approximately +/-2 sigma",
              "operation": "joint-valid Gaussian weighting of native complex cross-products and powers, then sample output centers",
              "output_spacing_is_separate": True}
    if enabled:
        kernels = [gaussian_axis(args.filter_wavelength, d, s)
                   for d, s in zip((-ref.dy, ref.dx), args.looks)]
        result.update(sigma_m=args.filter_wavelength/5.3,
                      native_kernel_shape_yx=[len(k) for k in kernels],
                      weights_y=kernels[0].tolist(), weights_x=kernels[1].tolist(),
                      nominal_weighted_sample_count=float(1/np.prod([np.sum(k*k) for k in kernels])),
                      sample_count_is_not_measured_independent_looks=True,
                      correlation_postfilter="GMTSAR fill.3x3: outer product of [1,2,1]/4; original mask retained",
                      roi_edges="read native halo outside selected ROI, limited to common product coverage")
    return result


def sample_separable(data, ky, kx, anchors_y, anchors_x, looks):
    """Convolve only required rows/columns; strided views avoid full kernel cubes."""
    windows = np.lib.stride_tricks.sliding_window_view(data, len(kx), axis=1)
    columns = slice(int(anchors_x[0])-len(kx)//2,
                    int(anchors_x[-1])-len(kx)//2+1, looks[1])
    reduced_x = np.einsum("ijk,k->ij", windows[:, columns, :], kx, dtype=np.float64)
    windows = np.lib.stride_tricks.sliding_window_view(reduced_x, len(ky), axis=0)
    rows = slice(int(anchors_y[0])-len(ky)//2,
                 int(anchors_y[-1])-len(ky)//2+1, looks[0])
    return np.einsum("ijk,k->ij", windows[rows, :, :], ky, dtype=np.float64)


def gaussian_tile(ref, sec, region, row, col, h, w, args, kernels):
    """Read one halo once, form both raw block means and Gaussian products."""
    ky, kx = kernels
    ly, lx = args.looks
    anchors_y = row + np.arange(h)*ly + ly//2
    anchors_x = col + np.arange(w)*lx + lx//2
    r0 = min(row, int(anchors_y[0])-len(ky)//2)
    c0 = min(col, int(anchors_x[0])-len(kx)//2)
    r1 = max(row+h*ly, int(anchors_y[-1])+(len(ky)-1)//2+1)
    c1 = max(col+w*lx, int(anchors_x[-1])+(len(kx)-1)//2+1)
    ro, co = region.secondary_row_offset, region.secondary_col_offset
    rr0, cc0 = max(r0, 0, -ro), max(c0, 0, -co)
    rr1 = min(r1, ref.shape[0], sec.shape[0]-ro)
    cc1 = min(c1, ref.shape[1], sec.shape[1]-co)
    a = np.zeros((r1-r0, c1-c0), dtype="complex64")
    b = np.zeros_like(a)
    valid = np.zeros(a.shape, dtype=bool)
    aa, va = ref.read(rr0, cc0, rr1-rr0, cc1-cc0)
    bb, vb = sec.read(rr0+ro, cc0+co, rr1-rr0, cc1-cc0)
    dst = np.s_[rr0-r0:rr1-r0, cc0-c0:cc1-c0]
    joint = va & vb
    a[dst], b[dst], valid[dst] = np.where(joint, aa, 0), np.where(joint, bb, 0), joint
    del aa, bb, va, vb, joint
    crop = np.s_[row-r0:row-r0+h*ly, col-c0:col-c0+w*lx]
    raw = multilook(a[crop], b[crop], valid[crop], args.looks, args.min_valid_fraction)
    ay, ax = anchors_y-r0, anchors_x-c0
    sample = lambda v: sample_separable(v, ky, kx, ay, ax, args.looks)
    support = np.clip(sample(valid), 0, 1)
    denom = np.maximum(support, np.finfo(float).tiny)
    cross = a*np.conj(b)
    z = (sample(cross.real)+1j*sample(cross.imag))/denom
    p1, p2 = sample(abs(a)**2)/denom, sample(abs(b)**2)/denom
    good = (np.isfinite(raw["mask"]) & (support >= args.min_valid_fraction-1e-12)
            & (p1 > 0) & (p2 > 0))
    coherence = np.divide(abs(z), np.sqrt(p1*p2), out=np.zeros(z.shape), where=(p1*p2)>0)
    filtered = {"realfilt": z.real, "imagfilt": z.imag, "phase": np.angle(z),
                "corr": np.clip(coherence, 0, 1), "amp1": p1, "amp2": p2,
                "mask": np.ones(z.shape)}
    filtered = {k: np.where(good, v, np.nan).astype("f4") for k, v in filtered.items()}
    filtered["filter_support"] = support.astype("f4")
    raw.update(filtered)
    return raw


def smooth_correlation(outdir, x, y, crs, report):
    """Match fill.3x3 weighting, renormalizing NaNs without filling masked holes."""
    source = outdir/"corr_gaussian.grd"
    (outdir/"corr.grd").rename(source)
    hist, count = np.zeros(100, dtype=np.int64), 0
    kernel = np.array([1., 2., 1.])/4
    with netCDF4.Dataset(source) as ds, ExitStack() as stack:
        ds.set_auto_mask(False)
        writer = GridWriter(outdir/"corr.grd", x, y, crs, "Gaussian coherence followed by GMTSAR fill.3x3 weighting", "1")
        stack.callback(writer.close)
        for r in range(0, len(y), 256):
            r0, r1 = max(0, r-1), min(len(y), r+257)
            data = ds["z"][r0:r1]
            valid = np.isfinite(data)
            def smooth(v):
                return convolve1d(convolve1d(v, kernel, axis=1, mode="constant"),
                                  kernel, axis=0, mode="constant")
            numerator, support = smooth(np.where(valid, data, 0).astype(float)), smooth(valid.astype(float))
            result = np.divide(numerator, support, out=np.full(data.shape, np.nan), where=valid & (support > 0))
            result = result[r-r0:r-r0+min(256, len(y)-r)].astype("f4")
            writer[r:r+len(result)] = result
            finite = result[np.isfinite(result)]
            hist += np.histogram(finite, bins=100, range=(0, 1))[0]
            count += len(finite)
    report["coherence_before_fill3x3"] = report["coherence"]
    report["coherence"] = {"valid_pixels": count, "valid_fraction": count/(len(x)*len(y)),
                           "histogram_0_to_1": hist.tolist(), "histogram_bin_width": 0.01,
                           "median_approx": float((np.searchsorted(np.cumsum(hist), count/2)+0.5)/100)}


def make_interferogram(ref, sec, region, args, outdir, x, y, report):
    names = {"real": ("Multilooked complex interferogram: real", "1"),
             "imag": ("Multilooked complex interferogram: imaginary", "1"),
             "phase": ("Wrapped reference-minus-secondary phase", "radian"),
             "corr": ("Sample coherence before Goldstein filtering", "1"),
             "amp1": ("Mean reference power (not amplitude)", "1"),
             "amp2": ("Mean secondary power (not amplitude)", "1"),
             "mask": ("Joint valid mask: 1 valid, NaN invalid", "1"),
             "nlooks": ("Number of contributing valid samples, not independent looks", "1"),
             "valid_fraction": ("Joint valid sample fraction", "1")}
    gaussian = report["gaussian_filter"]["enabled"]
    if gaussian:
        names.update(realfilt=("Gaussian weighted complex interferogram: real", "1"),
                     imagfilt=("Gaussian weighted complex interferogram: imaginary", "1"),
                     filter_support=("Fraction of Gaussian kernel weight supported by jointly valid native samples", "1"),
                     amp1=("Gaussian weighted reference power (not amplitude)", "1"),
                     amp2=("Gaussian weighted secondary power (not amplitude)", "1"))
        kernels = [np.array(report["gaussian_filter"][f"weights_{axis}"]) for axis in ("y", "x")]
    ly, lx = args.looks
    hist = np.zeros(100, dtype=np.int64)
    valid_count = 0
    last_message = time.monotonic()
    with ExitStack() as stack:
        writers = {}
        for name, (description, units) in names.items():
            writer = GridWriter(outdir/f"{name}.grd", x, y, ref.crs, description, units)
            stack.callback(writer.close)
            writers[name] = writer
        for r, c, h, w in blocks((len(y), len(x)), args.block_size):
            row, col = region.row+r*ly, region.col+c*lx
            if gaussian:
                result = gaussian_tile(ref, sec, region, row, col, h, w, args, kernels)
            else:
                a, va = ref.read(row, col, h*ly, w*lx)
                b, vb = sec.read(row+region.secondary_row_offset, col+region.secondary_col_offset, h*ly, w*lx)
                result = multilook(a, b, va & vb, args.looks, args.min_valid_fraction)
            for name, values in result.items():
                writers[name][r:r+h, c:c+w] = values
            coherence = result["corr"]
            finite = coherence[np.isfinite(coherence)]
            hist += np.histogram(finite, bins=100, range=(0, 1))[0]
            valid_count += len(finite)
            if time.monotonic()-last_message > 20:
                LOG.info("Interferogram: output row %d/%d", r+h, len(y))
                last_message = time.monotonic()
    if not valid_count:
        raise ValueError("No jointly valid pixels in selected region")
    report["coherence"] = {"valid_pixels": valid_count, "valid_fraction": valid_count/(len(y)*len(x)),
                           "histogram_0_to_1": hist.tolist(), "histogram_bin_width": 0.01,
                           "median_approx": float((np.searchsorted(np.cumsum(hist), valid_count/2)+0.5)/100)}
    if gaussian:
        smooth_correlation(outdir, x, y, ref.crs, report)
    elif ly*lx == 1:
        warn(report, "Single-sample coherence is identically 1; use multiple looks for meaningful coherence")


def executable(name, override):
    result = shutil.which(override or name)
    if result is None:
        raise ValueError(f"Cannot find executable {override or name}")
    return str(Path(result).resolve())


def run_tool(command, cwd, log_path):
    with open(log_path, "a", encoding="utf-8") as log:
        log.write("\nCOMMAND " + json.dumps(command) + "\n"); log.flush()
        completed = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f"{Path(command[0]).name} failed ({completed.returncode}); see {log_path}")


def phasefilt_selfcheck(binary, outdir, crs):
    """Protect against incompatible binaries, orientation errors and edge artifacts."""
    n = 96
    x = np.arange(n, dtype=float)+0.5
    y = np.arange(n, dtype=float)[::-1]+0.5
    expected = 0.23 + 0.021*np.arange(n)[:, None] + 0.013*np.arange(n)[None, :]
    with tempfile.TemporaryDirectory(prefix="phasefilt-check-", dir=outdir) as tmp:
        tmp = Path(tmp)
        write_grid(tmp/"real.grd", np.cos(expected), x, y, crs)
        write_grid(tmp/"imag.grd", np.sin(expected), x, y, crs)
        run_tool([binary, "-real", "real.grd", "-imag", "imag.grd", "-alpha", "0", "-psize", "32"],
                 tmp, outdir/"phasefilt.log")
        phase, _, _ = read_grid(tmp/"filtphase.grd")
        error = np.max(abs(np.angle(np.exp(1j*(phase[32:-32, 32:-32]-expected[32:-32, 32:-32])))))
        if not np.isfinite(error) or error > 1e-4:
            raise RuntimeError(f"Installed phasefilt failed alpha=0 identity check (error {error}); installation unchanged")
    return float(error)


def filter_interferogram(args, outdir, x, y, crs, report, binary):
    report["phasefilt_identity_max_error_rad"] = phasefilt_selfcheck(binary, outdir, crs)
    pad, size = args.filter_patch, args.filter_block_size
    adaptive = args.goldstein_mode == "coherence"
    gaussian = report["gaussian_filter"]["enabled"]
    source_names = {"real": "realfilt" if gaussian else "real",
                    "imag": "imagfilt" if gaussian else "imag", "mask": "mask"}
    if adaptive:
        source_names.update(amp1="amp1", amp2="amp2")
    report["goldstein_filter"] = {
        "mode": args.goldstein_mode, "patch_pixels": pad,
        "alpha": "1 - patch_weighted_coherence" if adaptive else args.filter_alpha,
        "input_complex": [source_names["real"]+".grd", source_names["imag"]+".grd"],
        "matches_local_filter_csh_amplitude_inputs": adaptive}
    with ExitStack() as stack:
        inputs = {}
        for name, source_name in source_names.items():
            ds = stack.enter_context(netCDF4.Dataset(outdir/f"{source_name}.grd"))
            ds.set_auto_mask(False)
            inputs[name] = ds["z"]
        writer = GridWriter(outdir/"phasefilt.grd", x, y, crs, "Goldstein filtered phase", "radian")
        stack.callback(writer.close)
        for r, c, h, w in blocks((len(y), len(x)), size):
            r0, r1, c0, c1 = max(0, r-pad), min(len(y), r+h+pad), max(0, c-pad), min(len(x), c+w+pad)
            pr0, pc0 = max(0, pad-r), max(0, pad-c)
            pr1, pc1 = max(0, r+h+pad-len(y)), max(0, c+w+pad-len(x))
            # Align the FFT patch lattice at tile boundaries and protect outer edges.
            pr1 += (-h) % pad
            pc1 += (-w) % pad
            with tempfile.TemporaryDirectory(prefix="phasefilt-tile-", dir=outdir) as tmp:
                tmp = Path(tmp)
                for name in ("real", "imag") + (("amp1", "amp2") if adaptive else ()):
                    values = np.nan_to_num(inputs[name][r0:r1, c0:c1], nan=0.0)
                    values = np.pad(values, ((pr0, pr1), (pc0, pc1)), mode="edge")
                    tx = x[0]+(c-pad+np.arange(values.shape[1]))*(x[1]-x[0])
                    ty = y[0]+(r-pad+np.arange(values.shape[0]))*(y[1]-y[0])
                    write_grid(tmp/f"{name}.grd", values, tx, ty, crs)
                command = [binary, "-v", "-real", "real.grd", "-imag", "imag.grd",
                           "-alpha", str(args.filter_alpha), "-psize", str(pad)]
                if adaptive:
                    # Both power files override -alpha in phasefilt.c, as in filter.csh.
                    command += ["-amp1", "amp1.grd", "-amp2", "amp2.grd"]
                run_tool(command, tmp, outdir/"phasefilt.log")
                filtered, _, _ = read_grid(tmp/"filtphase.grd")
                core = filtered[pad:pad+h, pad:pad+w]
                valid = np.isfinite(inputs["mask"][r:r+h, c:c+w])
                writer[r:r+h, c:c+w] = np.where(valid, core, np.nan)
            LOG.info("Phase filter: output tile row=%d col=%d", r, c)


def unwrap_plan(args, shape):
    """Use SNAPHU's tile geometry; the size guard applies to each solver tile."""
    tiles = [int(v) for v in args.unwrap_tiles]
    overlap = [int(v) if n > 1 else 0 for v, n in zip(args.unwrap_overlap, tiles)]
    if min(tiles) < 1 or min(args.unwrap_overlap) < 0 or args.unwrap_nproc < 1:
        raise ValueError("Unwrap tile counts/processes must be positive and overlaps nonnegative")
    tile_shape = [math.ceil((size+(n-1)*ov)/n) for size, n, ov in zip(shape, tiles, overlap)]
    last_shape = [size-(n-1)*(tile-ov) for size, n, tile, ov in zip(shape, tiles, tile_shape, overlap)]
    if min(last_shape) < 2 or any(n > 1 and ov >= tile for n, ov, tile in zip(tiles, overlap, tile_shape)):
        raise ValueError("Unwrap tiles/overlaps do not fit the output grid")
    pixels = math.prod(tile_shape)
    return {"tiles_yx": tiles, "overlap_yx_pixels": overlap,
            "processes": min(args.unwrap_nproc, math.prod(tiles)), "max_tile_shape_yx": tile_shape,
            "max_tile_pixels": pixels, "total_pixels": math.prod(shape),
            "within_size_guard": pixels <= args.max_unwrap_pixels,
            "assembly": "SNAPHU native region assembly" if math.prod(tiles) > 1 else "single tile"}


def enforce_unwrap_guard(args, plan):
    if not plan["within_size_guard"]:
        raise ValueError(f"SNAPHU solver tile has {plan['max_tile_pixels']:,} pixels, above the size guard "
                         f"({args.max_unwrap_pixels:,}). Increase --unwrap-tiles, select a smaller region, "
                         "or explicitly raise --max-unwrap-pixels after checking available memory. "
                         "--stage filter does not run SNAPHU.")


def phase_to_los(source, destination, wavelength, *, title='LOS from unwrapped phase; positive toward satellite'):
    """Convert phase to LOS in meters by direct scaling, without any offset."""
    if not np.isfinite(wavelength) or wavelength <= 0:
        raise ValueError("LOS conversion requires a positive finite wavelength")
    factor = -float(wavelength)/(4*np.pi)
    with netCDF4.Dataset(source) as ds:
        ds.set_auto_mask(False)
        x, y = np.asarray(ds['x'][:]), np.asarray(ds['y'][:])
        crs = CRS.from_wkt(ds['crs'].spatial_ref)
        writer = GridWriter(destination, x, y, crs, title, "m")
        try:
            writer.file.los_scale_m_per_radian = factor
            writer.file.radar_wavelength_m = float(wavelength)
            writer.file.phase_offset_subtracted = 0.0
            for row in range(0, len(y), 256):
                phase = np.asarray(ds['z'][row:row+256], dtype='f8')
                writer[row:row+256] = np.where(np.isfinite(phase), factor*phase, np.nan)
        finally:
            writer.close()
    return {'formula':'-wavelength/(4*pi)*phase', 'units':'m',
            'wavelength_m':float(wavelength), 'scale_m_per_radian':factor,
            'phase_offset_subtracted_rad':0.0, 'positive':'toward satellite'}


def unwrap_interferogram(args, outdir, x, y, crs, wavelength, report, binary, *, emit_los=True):
    shape = (len(y), len(x))
    plan = unwrap_plan(args, shape)
    enforce_unwrap_guard(args, plan)
    work = outdir/"snaphu_work"
    work.mkdir()
    valid_count = 0
    with netCDF4.Dataset(outdir/"phasefilt.grd") as ph, netCDF4.Dataset(outdir/"corr.grd") as co:
        ph.set_auto_mask(False); co.set_auto_mask(False)
        with open(work/"phase.bin", "wb") as pf, open(work/"corr.bin", "wb") as cf, open(work/"mask.bin", "wb") as mf:
            for r in range(0, len(y), 256):
                phase, corr = ph["z"][r:r+256], co["z"][r:r+256]
                valid = np.isfinite(phase) & np.isfinite(corr) & (corr >= args.correlation_threshold)
                np.where(valid, phase, 0).astype("=f4").tofile(pf)
                np.where(valid, corr, 0).astype("=f4").tofile(cf)
                valid.astype("i1").tofile(mf)
                valid_count += int(valid.sum())
    if valid_count < 16:
        raise ValueError("Fewer than 16 valid pixels above the unwrapping coherence threshold")
    effective = args.effective_looks or int(np.prod(args.looks))
    if args.effective_looks is None:
        warn(report, f"SNAPHU NCORRLOOKS={effective} is a provisional nominal output-block sample count, "
             "not measured independent looks; GSLC resampling and Gaussian weighting change the statistics. "
             "Use --effective-looks after estimating correlation statistics.")
    config = ["INFILEFORMAT FLOAT_DATA", "OUTFILEFORMAT FLOAT_DATA", "CORRFILEFORMAT FLOAT_DATA",
              "BYTEMASKFILE mask.bin", "CONNCOMPOUTTYPE UINT", "MINCONNCOMPFRAC 0.001", "MAXNCOMPS 1000",
              f"LAMBDA {wavelength:.16g}", f"NCORRLOOKS {effective:.8g}", f"DEFOMAX_CYCLE {args.defomax:.8g}",
              "# SMOOTH/DEFO are phase costs; do not interpret map axes as slant range/azimuth.",
              "# No radar geometry, baseline, or legacy C-band configuration is imported."]
    (work/"snaphu.conf").write_text("\n".join(config)+"\n", encoding="utf-8")
    selected_mode = args.unwrap_mode
    if selected_mode == "auto":
        selected_mode = "smooth" if args.defomax == 0 else "deformation"
    mode = "-s" if selected_mode == "smooth" else "-d"
    LOG.info("SNAPHU: %s, %d valid pixels; log: snaphu.log", selected_mode, valid_count)
    command = [binary, "phase.bin", str(len(x)), "-f", "snaphu.conf", "-c", "corr.bin", "-o", "unwrap.bin",
               "-g", "conncomp.bin", mode, "-v"]
    if math.prod(plan["tiles_yx"]) > 1:
        command += ["--tile", *map(str, plan["tiles_yx"]), *map(str, plan["overlap_yx_pixels"]),
                    "--nproc", str(plan["processes"]), "--tiledir", "tiles"]
        if args.keep_snaphu_files:
            command += ["-k"]
        warn(report, "SNAPHU assembles overlapping tiles into one phase grid. Its connected-component "
             "labels remain tile-local: different labels need not mean a phase discontinuity. "
             "Direct LOS conversion retains all accepted pixels; optional reference-area products retain the selected label.")
    solver_start = time.monotonic()
    run_tool(command, work, outdir/"snaphu.log")
    solver_seconds = time.monotonic()-solver_start
    LOG.info("SNAPHU completed in %.1f s", solver_seconds)
    for name in ("unwrap.bin", "conncomp.bin"):
        if (work/name).stat().st_size != np.prod(shape)*4:
            raise RuntimeError(f"Unexpected SNAPHU output size: {name}")
    phase = np.memmap(work/"unwrap.bin", dtype="=f4", mode="r", shape=shape)
    labels = np.memmap(work/"conncomp.bin", dtype="=u4", mode="r", shape=shape)
    mask = np.memmap(work/"mask.bin", dtype="i1", mode="r", shape=shape)
    wrapped = np.memmap(work/"phase.bin", dtype="=f4", mode="r", shape=shape)
    counts = {}
    rewrap_max = 0.0
    with ExitStack() as stack:
        uw = GridWriter(outdir/"unwrap.grd", x, y, crs, "Unwrapped phase; each connected component has its own ambiguity", "radian")
        cc = GridWriter(outdir/"conncomp.grd", x, y, crs, "SNAPHU connected component; 0 unreliable or invalid")
        stack.callback(uw.close); stack.callback(cc.close)
        for r in range(0, len(y), 256):
            good = (mask[r:r+256] != 0) & (labels[r:r+256] > 0) & np.isfinite(phase[r:r+256])
            if np.any(good):
                residual = np.angle(np.exp(1j*(phase[r:r+256][good].astype('f8')-wrapped[r:r+256][good])))
                rewrap_max = max(rewrap_max, float(np.max(np.abs(residual))))
            uw[r:r+256] = np.where(good, phase[r:r+256], np.nan)
            chunk_labels = np.where(good, labels[r:r+256], 0)
            cc[r:r+256] = chunk_labels
            ids, ns = np.unique(chunk_labels, return_counts=True)
            for label, count in zip(ids, ns):
                counts[int(label)] = counts.get(int(label), 0)+int(count)
    report["unwrapping"] = {"mode": selected_mode, "valid_input_pixels": valid_count,
                            "component_counts": counts, "effective_looks_used": effective,
                            "tile_plan": plan, "rewrap_max_error_rad": rewrap_max,
                            "solver_elapsed_seconds": solver_seconds,
                            "component_label_scope": "tile-local" if math.prod(plan["tiles_yx"]) > 1 else "whole input grid",
                            "rewrap_check_note": "Modulo-phase consistency does not prove integer-cycle correctness."}
    if not any(k > 0 for k in counts):
        warn(report, "SNAPHU found no reliable connected components; unwrap.grd contains no accepted phase")
    if emit_los and (args.reference_pixel or args.reference_bbox):
        if args.reference_pixel:
            row, col = args.reference_pixel
            if not (0 <= row < len(y) and 0 <= col < len(x)):
                raise ValueError("--reference-pixel is outside the multilooked output grid")
            rr, cc = slice(row, row+1), slice(col, col+1)
        else:
            xmin, xmax, ymin, ymax = args.reference_bbox
            if not (xmin < xmax and ymin < ymax):
                raise ValueError("Invalid --reference-bbox")
            rr = slice(np.searchsorted(-y, -ymax), np.searchsorted(-y, -ymin, side="right"))
            cc = slice(np.searchsorted(x, xmin), np.searchsorted(x, xmax, side="right"))
        selected = (mask[rr, cc] != 0) & (labels[rr, cc] > 0) & np.isfinite(phase[rr, cc])
        components = np.unique(labels[rr, cc][selected])
        if len(components) != 1:
            raise ValueError("Reference area must contain reliable pixels in exactly one connected component")
        component = int(components[0])
        reference_phase = float(np.median(phase[rr, cc][selected]))
        with ExitStack() as stack:
            rw = GridWriter(outdir/"unwrap_ref.grd", x, y, crs, "Unwrapped phase relative to user reference area", "radian")
            lw = GridWriter(outdir/"los.grd", x, y, crs, "Relative LOS displacement, positive toward satellite", "m")
            stack.callback(rw.close); stack.callback(lw.close)
            for r in range(0, len(y), 256):
                good = (mask[r:r+256] != 0) & (labels[r:r+256] == component)
                relative = np.where(good, phase[r:r+256]-reference_phase, np.nan)
                rw[r:r+256] = relative
                lw[r:r+256] = -wavelength/(4*np.pi)*relative
        report["reference"] = {"component": component, "phase_rad": reference_phase,
                               "pixels": int(selected.sum()), "selection_is_user_supplied_not_verified_stable": True}
    elif emit_los:
        report['los'] = phase_to_los(outdir/'unwrap.grd', outdir/'los.grd', wavelength)
    del phase, labels, mask, wrapped
    if not args.keep_snaphu_files:
        for name in ("phase.bin", "corr.bin", "mask.bin", "unwrap.bin", "conncomp.bin"):
            (work/name).unlink()


def export_geotiffs(outdir, x, y, crs):
    import rasterio
    from rasterio.transform import Affine
    dx, dy = x[1]-x[0], y[1]-y[0]
    transform = Affine(dx, 0, x[0]-dx/2, 0, dy, y[0]-dy/2)
    for name in ("phase", "corr", "mask", "phasefilt", "unwrap", "conncomp", "unwrap_ref", "los", "display_amp"):
        if not (outdir/f"{name}.grd").exists():
            continue
        with netCDF4.Dataset(outdir/f"{name}.grd") as ds:
            ds.set_auto_mask(False)
            with rasterio.open(outdir/f"{name}.tif", "w", driver="GTiff", width=len(x), height=len(y), count=1,
                               dtype="float32", crs=crs.to_wkt(), transform=transform, nodata=np.nan,
                               compress="deflate", predictor=3, tiled=True, BIGTIFF="IF_SAFER") as dst:
                for r in range(0, len(y), 256):
                    h = min(256, len(y)-r)
                    dst.write(ds["z"][r:r+h], 1, window=rasterio.windows.Window(0, r, len(x), h))
                dst.update_tags(phase_convention="reference * conjugate(secondary)", quality_report="run.json")


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("reference", type=Path)
    p.add_argument("secondary", type=Path)
    p.add_argument("-o", "--output", type=Path, help="New or empty output directory")
    p.add_argument("--frequency", choices=("A", "B", "AB"), default="A", help="AB runs both subbands independently into separate subdirectories")
    p.add_argument("--polarization", choices=tuple(POL_BITS), default="HH")
    sampling = p.add_mutually_exclusive_group()
    sampling.add_argument("--looks", nargs=2, type=int, metavar=("Y", "X"), help="Explicit map-grid look factors; otherwise derived per band from output spacing")
    sampling.add_argument("--output-spacing", nargs=2, type=float, metavar=("X_METERS", "Y_METERS"),
                          help="Output east/north spacing in meters; each must be an integer multiple of the native spacing")
    p.add_argument("--posting", type=float, default=40.0, help="Target square grid spacing in meters when --looks is omitted; requires integer native-pixel ratios")
    roi = p.add_mutually_exclusive_group()
    roi.add_argument("--window", nargs=4, type=int, metavar=("ROW", "COL", "HEIGHT", "WIDTH"), help="Native north-up reference indices, zero-based")
    roi.add_argument("--bbox", nargs=4, type=float, metavar=("XMIN", "XMAX", "YMIN", "YMAX"), help="Pixel-center selection in product CRS meters")
    p.add_argument("--stage", choices=("intf", "filter", "unwrap"), default="intf", help="Run through this stage; each run starts from inputs")
    p.add_argument("--inspect", action="store_true", help="Inspect pair metadata and grid; no output or image data read")
    p.add_argument("--assume-flattened", action="store_true", help="Permit missing flattening flags if independently verified")
    p.add_argument("--skip-registration-check", action="store_true", help="Skip diagnostic amplitude cross-correlation")
    p.add_argument("--min-valid-fraction", type=float, default=0.8)
    p.add_argument("--block-size", type=int, default=128, help="Output tile size for HDF5 reading/multilooking")
    p.add_argument("--filter-wavelength", type=float, help="Gaussian wavelength in meters (GMTSAR w/5.3 sigma); overrides band defaults. 0 disables Gaussian processing")
    p.add_argument("--filter-wavelength-a", type=float, default=160., help="Frequency A Gaussian wavelength in meters")
    p.add_argument("--filter-wavelength-b", type=float, default=640., help="Frequency B Gaussian wavelength in meters")
    p.add_argument("--goldstein-mode", choices=("coherence", "fixed"), help="Default coherence, as local filter.csh; explicit --filter-alpha selects fixed mode")
    p.add_argument("--filter-alpha", type=float, help="Fixed Goldstein exponent, default 0.5 in fixed mode; incompatible with explicit coherence mode")
    p.add_argument("--filter-patch", type=int, default=32)
    p.add_argument("--filter-block-size", type=int, default=1024)
    p.add_argument("--phasefilt", help="Path to existing GMTSAR phasefilt")
    p.add_argument("--snaphu", help="Path to existing SNAPHU")
    p.add_argument("--correlation-threshold", "--threshold-snaphu", type=float, default=0.1,
                   help="SNAPHU coherence threshold; 0 skips unwrapping, following p2p_processing_nsr.csh")
    p.add_argument("--unwrap-mode", choices=("auto", "smooth", "deformation"), default="auto",
                   help="auto follows snaphu.csh: defomax=0 selects smooth, positive selects deformation")
    p.add_argument("--defomax", type=float, default=0.0, help="Maximum discontinuity, cycles; 0 selects smooth costs in auto mode")
    p.add_argument("--effective-looks", type=float, help="SNAPHU NCORRLOOKS; default is nominal Y*X, not measured ENL")
    reference = p.add_mutually_exclusive_group()
    reference.add_argument("--reference-pixel", nargs=2, type=int, metavar=("ROW", "COL"), help="Zero-based pixel in multilooked output; must be independently chosen as stable")
    reference.add_argument("--reference-bbox", nargs=4, type=float, metavar=("XMIN", "XMAX", "YMIN", "YMAX"), help="Stable reference area in product CRS meters")
    p.add_argument("--geotiff", action="store_true", help="Also export key products as georeferenced GeoTIFFs (requires rasterio)")
    p.add_argument("--pdf", action="store_true", help="Also draw GMTSAR p2p-style GMT PDFs and a multipage PDF per frequency")
    p.add_argument("--gmt", help="Path to GMT for --pdf")
    p.add_argument("--keep-snaphu-files", action="store_true")
    p.add_argument("--unwrap-tiles", nargs=2, type=int, default=[1, 1], metavar=("Y", "X"),
                   help="SNAPHU tile rows/columns; native overlap-aware assembly writes one full grid")
    p.add_argument("--unwrap-overlap", nargs=2, type=int, default=[512, 512], metavar=("Y", "X"),
                   help="SNAPHU overlap in output pixels; unused on axes with only one tile")
    p.add_argument("--unwrap-nproc", type=int, default=1, help="Concurrent SNAPHU tile workers")
    p.add_argument("--max-unwrap-pixels", type=int, default=16000000,
                   help="Pixel guard per SNAPHU solver tile; use --unwrap-tiles for full-scene grids")
    return p


def process_one(args, emit_inspection=True):
    args = argparse.Namespace(**vars(args))
    args.requested_stage = args.stage
    if args.stage == "unwrap" and args.correlation_threshold == 0:
        args.stage = "filter"
    if args.filter_wavelength is None:
        args.filter_wavelength = args.filter_wavelength_a if args.frequency == "A" else args.filter_wavelength_b
    if any(not np.isfinite(w) or (w != 0 and not 2 <= w <= 10000)
           for w in (args.filter_wavelength, args.filter_wavelength_a, args.filter_wavelength_b)):
        raise ValueError("Gaussian wavelengths must be 0 (disabled) or between 2 and 10000 meters")
    if args.goldstein_mode == "coherence" and args.filter_alpha is not None:
        raise ValueError("--filter-alpha sets a fixed exponent; do not combine it with --goldstein-mode coherence")
    args.goldstein_mode = args.goldstein_mode or ("fixed" if args.filter_alpha is not None else "coherence")
    args.filter_alpha = 0.5 if args.filter_alpha is None else args.filter_alpha
    if (args.looks is not None and min(args.looks) < 1) or args.block_size < 1:
        raise ValueError("Look factors and block size must be positive")
    if not np.isfinite(args.posting) or args.posting <= 0:
        raise ValueError("posting must be positive and finite")
    if args.output_spacing is not None and any(not np.isfinite(v) or v <= 0 for v in args.output_spacing):
        raise ValueError("output spacing must be positive and finite")
    if not 0 < args.min_valid_fraction <= 1 or not 0 <= args.correlation_threshold <= 1:
        raise ValueError("Invalid validity fraction or coherence threshold")
    if not 0 <= args.filter_alpha <= 1.5 or args.filter_patch < 8 or args.filter_patch & (args.filter_patch-1):
        raise ValueError("Filter alpha must be in [0, 1.5] and patch a power of two >=8")
    if args.filter_block_size < args.filter_patch or args.filter_block_size % args.filter_patch:
        raise ValueError("Filter block size must be a positive multiple of filter patch")
    if not np.isfinite(args.defomax) or args.defomax < 0:
        raise ValueError("defomax must be finite and nonnegative")
    if args.effective_looks is not None and (not np.isfinite(args.effective_looks) or args.effective_looks < 1):
        raise ValueError("effective looks must be finite and >=1")
    if args.max_unwrap_pixels < 16:
        raise ValueError("max unwrap pixels must be >=16")
    if not args.inspect and args.output is None:
        raise ValueError("--output is required unless --inspect is used")
    if (args.reference_pixel or args.reference_bbox) and args.stage != "unwrap" and not args.inspect:
        raise ValueError("A reference area is only used with --stage unwrap")
    report = {"version": VERSION, "status": "inspecting", "warnings": [], "phase_convention": "reference * conjugate(secondary)",
              "los_convention": "-wavelength/(4*pi) * phase; meters; positive toward satellite",
              "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              "created_utc": datetime.now(timezone.utc).isoformat()}
    start = time.monotonic()
    with h5py.File(args.reference, "r") as rh, h5py.File(args.secondary, "r") as sh:
        ref = Product(rh, args.reference, args.frequency, args.polarization)
        sec = Product(sh, args.secondary, args.frequency, args.polarization)
        if args.looks is None:
            dx, dy = args.output_spacing or (args.posting, args.posting)
            factors = np.array([dy/abs(ref.dy), dx/ref.dx])
            if np.any(factors < 1) or not np.allclose(factors, np.rint(factors), rtol=0, atol=1e-6):
                name = "--output-spacing" if args.output_spacing else "--posting"
                raise ValueError(f"{name} must be integer multiples of native X/Y spacings "
                                 f"({ref.dx:g}, {abs(ref.dy):g} m); finer grids are not resampled here")
            args.looks = [int(v) for v in np.rint(factors)]
        report["arguments"]["looks"] = args.looks
        region = pair_region(ref, sec, args, report)
        report["gaussian_filter"] = gaussian_plan(args, ref)
        report["processing_recipe"] = {
            "reference": "installed p2p_processing_nsr.csh, filter.csh, make_gaussian_filter.c, snaphu.csh",
            "coordinates": "GSLC map grid; native Gaussian weighting before output decimation",
            "radar_prefilter": "not copied: radar PRF/range-sampling preset kernels do not define a GSLC map-grid kernel",
            "snaphu_geometry": "legacy C-band geometry is not imported",
            "unimplemented_legacy_options": ["near_interp", "automatic landmask", "ionospheric correction", "phase gradient"]}
        ly, lx = args.looks
        x = ref.x[region.col:region.col+region.width].reshape(-1, lx).mean(axis=1)
        y = ref.y[region.row:region.row+region.height].reshape(-1, ly).mean(axis=1)
        report.update(reference_product=ref.meta, secondary_product=sec.meta, region=vars(region),
                      output_shape=[len(y), len(x)], output_spacing_xy_m=[ref.dx*lx, ref.dy*ly],
                      output_pixel_center_extent=[float(x[0]), float(x[-1]), float(y[-1]), float(y[0])],
                      wavelength_m=C/ref.frequency_hz, epsg=ref.epsg)
        if args.stage == "unwrap":
            report["unwrapping_plan"] = unwrap_plan(args, [len(y), len(x)])
        if report["gaussian_filter"]["enabled"] and max(ref.dx*lx, -ref.dy*ly) > args.filter_wavelength/2:
            warn(report, "Output spacing is coarse relative to the Gaussian wavelength; this sampling can alias fine fringes. "
                 "Use a finer --posting or a larger --filter-wavelength for quantitative processing.")
        if args.inspect:
            if emit_inspection:
                print(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False))
            return report
        if args.stage == "unwrap":
            enforce_unwrap_guard(args, report["unwrapping_plan"])
        # Check external dependencies before creating any processing outputs.
        binaries = {}
        if args.stage in ("filter", "unwrap"):
            binaries["phasefilt"] = executable("phasefilt", args.phasefilt)
        if args.stage == "unwrap":
            binaries["snaphu"] = executable("snaphu", args.snaphu)
        if args.geotiff:
            import rasterio  # noqa: F401
        if args.pdf:
            from plot_gslc_gmt import plot_dependencies
            binaries.update(plot_dependencies(args.gmt))
        report["executables"] = binaries
        outdir = args.output.resolve()
        if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
            raise ValueError(f"Output must be a new or empty directory: {outdir}")
        outdir.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(outdir/"processing.log", encoding="utf-8")
        file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        LOG.addHandler(file_handler)
        try:
            report["status"] = "running"
            (outdir/"run.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
            LOG.info("Output %d rows x %d columns, %.3f x %.3f m, EPSG:%d", len(y), len(x), ref.dx*lx, -ref.dy*ly, ref.epsg)
            if not args.skip_registration_check:
                report["registration_check"] = registration_check(ref, sec, region, report)
            else:
                warn(report, "Amplitude registration check skipped; grid alignment alone is not precision coregistration")
            LOG.info("Forming complex interferogram and sample coherence")
            report["stage_seconds"] = {}
            stage_start = time.monotonic()
            make_interferogram(ref, sec, region, args, outdir, x, y, report)
            report["stage_seconds"]["interferogram"] = time.monotonic()-stage_start
            if args.stage in ("filter", "unwrap"):
                stage_start = time.monotonic()
                filter_interferogram(args, outdir, x, y, ref.crs, report, binaries["phasefilt"])
                report["stage_seconds"]["goldstein"] = time.monotonic()-stage_start
            if args.stage == "unwrap":
                stage_start = time.monotonic()
                unwrap_interferogram(args, outdir, x, y, ref.crs, C/ref.frequency_hz, report, binaries["snaphu"])
                report["stage_seconds"]["unwrap_with_io"] = time.monotonic()-stage_start
            if args.pdf:
                from plot_gslc_gmt import plot_products
                stage_start = time.monotonic()
                report["plotting"] = plot_products(outdir, report, binaries)
                report["stage_seconds"]["pdf"] = time.monotonic()-stage_start
            if args.geotiff:
                stage_start = time.monotonic()
                export_geotiffs(outdir, x, y, ref.crs)
                report["stage_seconds"]["geotiff"] = time.monotonic()-stage_start
            report["status"] = "complete"
            report["elapsed_seconds"] = time.monotonic()-start
            report["outputs"] = sorted(p.name for p in outdir.iterdir() if p.is_file())
            LOG.info("Complete: %s (%.1f s)", outdir, report["elapsed_seconds"])
        except BaseException as exc:
            report["status"], report["error"] = "failed", str(exc)
            raise
        finally:
            (outdir/"run.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
            LOG.removeHandler(file_handler)
            file_handler.close()
    return report


def process(args):
    if args.frequency != "AB":
        return process_one(args)
    if args.window:
        raise ValueError("--frequency AB requires --bbox (map coordinates) instead of band-specific --window indices")
    if args.reference_pixel:
        raise ValueError("--frequency AB uses --reference-bbox so both bands share the same physical reference area")
    if not args.inspect and args.output is None:
        raise ValueError("--output is required unless --inspect is used")
    # Preflight both products before creating outputs, including the legacy B mask.
    planned = {}
    for band in ("A", "B"):
        options = argparse.Namespace(**vars(args))
        options.frequency, options.inspect = band, True
        planned[band] = process_one(options, emit_inspection=False)
        if planned[band]["arguments"]["stage"] == "unwrap" and not args.inspect:
            enforce_unwrap_guard(options, planned[band]["unwrapping_plan"])
    if args.inspect:
        report = {"version": VERSION, "status": "inspecting", "bands": planned}
        print(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False))
        return report
    outdir = args.output.resolve()
    if outdir.exists() and (not outdir.is_dir() or any(outdir.iterdir())):
        raise ValueError(f"Output must be a new or empty directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    report = {"version": VERSION, "status": "running", "bands": {},
              "note": "Independent A and B interferograms; not a split-spectrum ionospheric correction."}
    try:
        for band in ("A", "B"):
            options = argparse.Namespace(**vars(args))
            options.frequency = band
            options.output = outdir/f"frequency{band}_{args.polarization}"
            report["bands"][band] = process_one(options)
        report["status"] = "complete"
    except BaseException as exc:
        report["status"], report["error"] = "failed", str(exc)
        raise
    finally:
        (outdir/"batch.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    return report


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parser().parse_args(argv)
    try:
        process(args)
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        LOG.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        LOG.error("Interrupted; partial outputs are marked failed in run.json")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
