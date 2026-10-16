#!/usr/bin/env python3
"""Shared, side-effect-free contracts for the GOES01 physical validation.

The case deliberately exposes only the scientific claim supported by the archived
inputs: high-energy proton differential intensity at GOES-16.  In particular, no
density is derived from one broad SGPS channel and no directional response is
invented for the two telescope heads.
"""

from __future__ import annotations

import csv
import ctypes
import ctypes.util
import hashlib
import json
import math
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


CASE_ID = "GOES01"
SPACECRAFT = "GOES16"
EVENT_START = datetime(2017, 9, 11, 7, 30, tzinfo=timezone.utc)
EVENT_END_EXCLUSIVE = datetime(2017, 9, 14, 0, 0, tzinfo=timezone.utc)
P8A_INDEX = 8
# No cross-head covariance is documented.  Requiring the two delivered uncertainty
# intervals to overlap uses their stated magnitudes without inventing independence.
TWO_HEAD_ENVELOPE_MAX = 1.0
M2_TO_CM2 = 1.0e-4

# QUICK and COMPLETE use the same event, quality rules, physical model, and gates.
# Only the predeclared temporal and numerical resolution changes.  Epoch cadence is
# measured from EVENT_START using the SGPS five-minute record start, before AMPS is
# run and therefore before any model/reference residual can be known.
PROFILES: dict[str, dict[str, int]] = {
    "QUICK": {"cadence_minutes": 60, "ds_intervals": 24, "angular_samples": 96},
    "COMPLETE": {"cadence_minutes": 15, "ds_intervals": 48, "angular_samples": 384},
}

# These engineering gates are inherited unchanged from the SA01--SA04 observational
# workflows.  They are intentionally looser than an instrument-calibration claim:
# this is a global kinetic-transport regression using an independent upstream
# boundary, not an in-flight SGPS calibration exercise.
GATES = {
    "minimum_matched_epochs": 12,
    "median_abs_log10_ratio_max": 0.30,
    "factor_two_fraction_each_head_min": 0.60,
    "factor_two_fraction_pooled_min": 0.70,
    "maximum_unresolved_fraction": 0.01,
}

OBSERVATION_FILES = {
    "se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc":
        "f5da29fe3a7dc3092a311a191ee015d6cb97773c6425e64214db2d0dab5fafb2",
    "spacecraft/goes16_ephemeris_1min_20170910_20170913.csv":
        "f67f6e4a017cbedf71a7206670a85b1924d05916367fdeddfd6800642f7467e9",
}
MODEL_FILES = {
    "boundary_soho_erne_5min.dat":
        "96f2049eb05e3b444286dd3346e81aaf59ae0003457d931857447724e376f736",
    "ts05_driver_5min.txt":
        "97779580b058a3c15c7309f92ca5da6ee21b5494bc3cdbdafd67b034924445ac",
}

_QUOTED = re.compile(r'"([^"]+)"')
_ZONE_TIME = re.compile(r'T="[^"]*?((?:19|20)\d\d-\d\d-\d\dT\d\d:\d\d:\d\dZ?)')


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def iso_utc(value: datetime) -> str:
    """Render the half-minute model epochs without locale or host-timezone effects."""
    return value.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_utc(value: str) -> datetime:
    return datetime.fromisoformat(value.strip().replace("Z", "+00:00"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    """Publish complete strict JSON so interrupted runs never look successful."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty evidence table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path}: no rows")
    return rows


@dataclass(frozen=True)
class SpectrumZone:
    """One timestamped native Tecplot spectrum with its declared column order."""

    utc: str
    columns: tuple[str, ...]
    rows: np.ndarray


def _quoted_snapshot_epoch(line: str, path: Path) -> str:
    fields = _QUOTED.findall(line)
    if len(fields) != 1:
        raise ValueError(f"{path}: malformed SNAPSHOT_EPOCH_UTC")
    parse_utc(fields[0])
    return fields[0].rstrip("Z") + "Z"


def read_spectrum_zones(path: Path) -> list[SpectrumZone]:
    """Strictly read GRIDDED or GRIDLESS native AMPS spectrum zones.

    Declared widths and row counts are evidence, not hints. A snapshot may supply
    its UTC through AUXDATA, but conflicting zone and snapshot times are rejected.
    Invalid nodes are never dropped and interpolated across during comparison.
    """
    variables: tuple[str, ...] | None = None
    zones: list[SpectrumZone] = []
    current: str | None = None
    expected_rows = 0
    rows: list[list[float]] = []
    snapshot_epoch: str | None = None

    def finish() -> None:
        nonlocal current, expected_rows, rows
        if current is None:
            return
        if variables is None or len(rows) != expected_rows:
            raise ValueError(f"{path}: incomplete native spectrum zone {current}")
        values = np.asarray(rows, dtype=float)
        if np.any(~np.isfinite(values)):
            raise ValueError(f"{path}: nonfinite native spectrum value at {current}")
        zones.append(SpectrumZone(current, variables, values))
        current, expected_rows, rows = None, 0, []

    with path.open(encoding="utf-8") as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith("AUXDATA SNAPSHOT_EPOCH_UTC="):
                snapshot_epoch = _quoted_snapshot_epoch(line, path)
            elif line.startswith("VARIABLES"):
                variables = tuple(_QUOTED.findall(line))
                if not variables:
                    raise ValueError(f"{path}: cannot parse VARIABLES declaration")
            elif line.startswith("ZONE"):
                finish()
                match = _ZONE_TIME.search(line)
                count = re.search(r"\bI\s*=\s*(\d+)", line)
                if (match is None and snapshot_epoch is None) or count is None:
                    raise ValueError(f"{path}: malformed native spectrum ZONE")
                current = (match.group(1).rstrip("Z") + "Z"
                           if match is not None else snapshot_epoch)
                if snapshot_epoch is not None and current != snapshot_epoch:
                    raise ValueError(f"{path}: zone/snapshot epoch mismatch")
                expected_rows = int(count.group(1))
            elif current is not None and line and not line.startswith(("AUXDATA", "#", "!")):
                fields = [float(value) for value in line.split()]
                if variables is None or len(fields) != len(variables):
                    raise ValueError(f"{path}: native spectrum row width mismatch")
                rows.append(fields)
    finish()
    if not zones or len({zone.utc for zone in zones}) != len(zones):
        raise ValueError(f"{path}: no native spectrum zones or duplicate epochs")
    return zones


def _resolve_child(root: Path, child: str) -> Path:
    direct = root / child
    nested = root / "goes16_sgps_sep2017" / child
    if direct.is_file():
        return direct.resolve()
    if nested.is_file():
        return nested.resolve()
    raise FileNotFoundError(f"cannot locate GOES01 input {child!r} below {root}")


def resolve_observations(root: Path) -> tuple[Path, Path, Path]:
    """Resolve and authenticate the exact SGPS and ephemeris bytes under test."""
    resolved = root.expanduser().resolve()
    paths = {name: _resolve_child(resolved, name) for name in OBSERVATION_FILES}
    for name, expected in OBSERVATION_FILES.items():
        actual = sha256_file(paths[name])
        if actual != expected:
            raise ValueError(f"{paths[name]}: SHA-256 {actual} != frozen {expected}")
    event_root = paths[next(iter(OBSERVATION_FILES))].parent
    return event_root, paths[next(iter(OBSERVATION_FILES))], paths[
        "spacecraft/goes16_ephemeris_1min_20170910_20170913.csv"
    ]


def resolve_model_inputs(root: Path) -> tuple[Path, Path, Path]:
    """Accept the package root or model_inputs/SA04, then freeze processed inputs."""
    candidate = root.expanduser().resolve()
    probes = (candidate / "validation_data" / "processed" / "SA04",
              candidate.parents[1] / "validation_data" / "processed" / "SA04"
              if len(candidate.parents) > 1 else candidate)
    processed = next((path for path in probes if path.is_dir()), None)
    if processed is None:
        raise FileNotFoundError(
            f"{root}: expected package validation_data/processed/SA04 or model_inputs/SA04"
        )
    paths = {name: (processed / name).resolve() for name in MODEL_FILES}
    for name, expected in MODEL_FILES.items():
        if not paths[name].is_file():
            raise FileNotFoundError(paths[name])
        actual = sha256_file(paths[name])
        if actual != expected:
            raise ValueError(f"{paths[name]}: SHA-256 {actual} != frozen {expected}")
    return processed, paths["boundary_soho_erne_5min.dat"], paths["ts05_driver_5min.txt"]


class NetCDFReader:
    """Minimal read-only binding to libnetcdf for this fixed NOAA product schema.

    A full Python NetCDF package is unnecessary for eight primitive arrays.  Calling
    the installed C library removes an otherwise hidden virtual-environment
    dependency while retaining the library's own format validation and conversion.
    Every return code is checked and every expected variable shape is verified.
    """

    NC_NOWRITE = 0

    def __init__(self, path: Path):
        library = ctypes.util.find_library("netcdf")
        if not library:
            raise RuntimeError("GOES01 requires the system libnetcdf runtime")
        self.lib = ctypes.CDLL(library)
        self.lib.nc_strerror.restype = ctypes.c_char_p
        self.lib.nc_open.argtypes = (ctypes.c_char_p, ctypes.c_int,
                                     ctypes.POINTER(ctypes.c_int))
        self.lib.nc_close.argtypes = (ctypes.c_int,)
        self.lib.nc_inq_varid.argtypes = (ctypes.c_int, ctypes.c_char_p,
                                          ctypes.POINTER(ctypes.c_int))
        self.lib.nc_inq_varndims.argtypes = (ctypes.c_int, ctypes.c_int,
                                             ctypes.POINTER(ctypes.c_int))
        self.lib.nc_inq_vardimid.argtypes = (ctypes.c_int, ctypes.c_int,
                                             ctypes.POINTER(ctypes.c_int))
        self.lib.nc_inq_dimlen.argtypes = (ctypes.c_int, ctypes.c_int,
                                           ctypes.POINTER(ctypes.c_size_t))
        self.lib.nc_get_var_double.argtypes = (ctypes.c_int, ctypes.c_int,
                                               ctypes.POINTER(ctypes.c_double))
        self.lib.nc_get_var_uint.argtypes = (ctypes.c_int, ctypes.c_int,
                                             ctypes.POINTER(ctypes.c_uint))
        self.lib.nc_inq_attlen.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_char_p,
                                           ctypes.POINTER(ctypes.c_size_t))
        self.lib.nc_get_att_text.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_char_p,
                                             ctypes.c_char_p)
        handle = ctypes.c_int()
        self._check(self.lib.nc_open(os.fsencode(path), self.NC_NOWRITE,
                                     ctypes.byref(handle)), f"open {path}")
        self.handle = handle.value

    def _check(self, code: int, context: str) -> None:
        if code:
            message = self.lib.nc_strerror(code).decode("utf-8", "replace")
            raise ValueError(f"NetCDF {context}: {message}")

    def close(self) -> None:
        if getattr(self, "handle", -1) >= 0:
            self._check(self.lib.nc_close(self.handle), "close")
            self.handle = -1

    def __enter__(self) -> "NetCDFReader":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _variable(self, name: str) -> tuple[int, tuple[int, ...]]:
        variable = ctypes.c_int()
        self._check(self.lib.nc_inq_varid(self.handle, name.encode(),
                                          ctypes.byref(variable)), f"find {name}")
        count = ctypes.c_int()
        self._check(self.lib.nc_inq_varndims(self.handle, variable.value,
                                             ctypes.byref(count)), f"dimensions {name}")
        dimension_ids = (ctypes.c_int * count.value)()
        self._check(self.lib.nc_inq_vardimid(self.handle, variable.value, dimension_ids),
                    f"dimension ids {name}")
        shape: list[int] = []
        for dimension in dimension_ids:
            length = ctypes.c_size_t()
            self._check(self.lib.nc_inq_dimlen(self.handle, dimension,
                                               ctypes.byref(length)), f"dimension length {name}")
            shape.append(length.value)
        return variable.value, tuple(shape)

    def doubles(self, name: str, expected_tail: tuple[int, ...] = ()) -> np.ndarray:
        variable, shape = self._variable(name)
        if expected_tail and shape[-len(expected_tail):] != expected_tail:
            raise ValueError(f"{name}: shape {shape} does not end in {expected_tail}")
        values = np.empty(int(np.prod(shape)), dtype=np.float64)
        self._check(self.lib.nc_get_var_double(
            self.handle, variable, values.ctypes.data_as(ctypes.POINTER(ctypes.c_double))),
            f"read {name}")
        return values.reshape(shape)

    def uints(self, name: str, expected_tail: tuple[int, ...] = ()) -> np.ndarray:
        variable, shape = self._variable(name)
        if expected_tail and shape[-len(expected_tail):] != expected_tail:
            raise ValueError(f"{name}: shape {shape} does not end in {expected_tail}")
        values = np.empty(int(np.prod(shape)), dtype=np.uint32)
        self._check(self.lib.nc_get_var_uint(
            self.handle, variable, values.ctypes.data_as(ctypes.POINTER(ctypes.c_uint))),
            f"read {name}")
        return values.reshape(shape)

    def attribute(self, variable_name: str, attribute_name: str) -> str:
        variable, _ = self._variable(variable_name)
        length = ctypes.c_size_t()
        self._check(self.lib.nc_inq_attlen(self.handle, variable, attribute_name.encode(),
                                           ctypes.byref(length)),
                    f"attribute length {variable_name}:{attribute_name}")
        buffer = ctypes.create_string_buffer(length.value + 1)
        self._check(self.lib.nc_get_att_text(self.handle, variable, attribute_name.encode(), buffer),
                    f"attribute {variable_name}:{attribute_name}")
        return buffer.raw[:length.value].decode("utf-8").strip("\x00 ")


def load_sgps(path: Path) -> dict[str, np.ndarray]:
    """Load only documented arrays and reject unit/schema drift before selection."""
    with NetCDFReader(path) as reader:
        if reader.attribute("L2_SciData_TimeStamp", "units") != (
                "seconds since 2000-01-01 12:00:00"):
            raise ValueError("unexpected SGPS timestamp units")
        expected_flux_units = "protons/(cm^2 sr keV s)"
        if reader.attribute("AvgDiffProtonFlux", "units") != expected_flux_units:
            raise ValueError("unexpected SGPS differential-flux units")
        if reader.attribute("AvgDiffProtonFluxUncert", "units") != expected_flux_units:
            raise ValueError("unexpected SGPS differential-flux uncertainty units")
        if (reader.attribute("DiffProtonLowerEnergy", "units") != "keV" or
                reader.attribute("DiffProtonUpperEnergy", "units") != "keV"):
            raise ValueError("unexpected SGPS channel-bound units")
        channel_description = reader.attribute("DiffProtonLowerEnergy", "long_name")
        # P8A's index and the west/east sensor ordering come from this source
        # metadata. Binding both strings prevents a file with the same dimensions
        # but a reordered channel or sensor axis from being interpreted silently.
        if ("(P1, P2A, P2B, P3, P4, P5, P6, P7, P8A, P8B, P8C, P9, P10)"
                not in channel_description or
                "First unit (-X) looks Westward" not in channel_description or
                "Second unit (+X) looks Eastward" not in channel_description):
            raise ValueError("SGPS channel/sensor ordering metadata changed")
        result = {
            "time": reader.doubles("L2_SciData_TimeStamp"),
            "flux": reader.doubles("AvgDiffProtonFlux", (2, 13)),
            "uncertainty": reader.doubles("AvgDiffProtonFluxUncert", (2, 13)),
            "valid": reader.uints("DiffValidL1bSamplesInAvg", (2, 13)),
            "dtc": reader.uints("DiffDQFdtcSum", (2, 13)),
            "oob": reader.uints("DiffDQFoobSum", (2, 13)),
            "err": reader.uints("DiffDQFerrSum", (2, 13)),
            "records": reader.uints("L1bRecordsInAvg"),
            "yaw": reader.uints("YawFlipFlag"),
            "lower": reader.doubles("DiffProtonLowerEnergy", (2, 13)),
            "upper": reader.doubles("DiffProtonUpperEnergy", (2, 13)),
        }
    expected_shapes = {
        "time": (8640,), "flux": (8640, 2, 13), "uncertainty": (8640, 2, 13),
        "valid": (8640, 2, 13), "dtc": (8640, 2, 13), "oob": (8640, 2, 13),
        "err": (8640, 2, 13), "records": (8640,), "yaw": (8640,),
        "lower": (2, 13), "upper": (2, 13),
    }
    for name, expected in expected_shapes.items():
        if result[name].shape != expected:
            raise ValueError(f"SGPS {name} shape {result[name].shape} != frozen {expected}")
    return result


def ratio_metrics(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    ratios = np.asarray([float(row["model"]) / float(row["reference"]) for row in rows])
    if not len(ratios) or np.any(~np.isfinite(ratios)) or np.any(ratios < 0.0):
        raise ValueError("comparison ratios must be a nonempty finite nonnegative population")
    logs = np.asarray([math.log10(value) if value > 0.0 else -math.inf for value in ratios])
    median_abs = float(np.median(np.abs(logs)))
    median_bias = float(np.median(logs))
    return {
        "count": len(rows),
        "positive_model_count": int(np.count_nonzero(ratios > 0.0)),
        # Physical zeros stay in the ordered population as infinite residuals. They
        # are never dropped merely because logarithms cannot represent zero.
        "median_abs_log10_ratio": median_abs if math.isfinite(median_abs) else None,
        "median_log10_bias": median_bias if math.isfinite(median_bias) else None,
        "factor_two_fraction": float(np.mean(np.isfinite(logs) &
                                              (np.abs(logs) <= math.log10(2.0)))),
    }


def power_law_bin_average(energy: Iterable[float], intensity: Iterable[float],
                          lower: float, upper: float) -> float:
    """Integrate AMPS log-spaced spectral nodes with zero-preserving interpolation."""
    e = np.asarray(tuple(energy), dtype=float)
    j = np.asarray(tuple(intensity), dtype=float)
    if (len(e) < 2 or np.any(~np.isfinite(e)) or np.any(~np.isfinite(j)) or
            np.any(j < 0.0) or np.any(np.diff(e) <= 0.0) or
            lower < e[0] or upper > e[-1]):
        raise ValueError("model spectrum does not bracket the requested SGPS channel")
    edges = np.unique(np.concatenate(([lower], e[(e > lower) & (e < upper)], [upper])))
    total = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        index = max(0, min(int(np.searchsorted(e, math.sqrt(left * right)) - 1), len(e) - 2))
        e0, e1, j0, j1 = e[index], e[index + 1], j[index], j[index + 1]
        if j0 > 0.0 and j1 > 0.0:
            slope = math.log(j1 / j0) / math.log(e1 / e0)
            scale = j0 / (e0 ** slope)
            if abs(slope + 1.0) < 1.0e-12:
                total += scale * math.log(right / left)
            else:
                total += scale * (right ** (slope + 1.0) - left ** (slope + 1.0)) / (slope + 1.0)
        else:
            # A physical zero must stay a zero/linear boundary, never be replaced by
            # a logarithmic floor that would manufacture transmitted particles.
            left_j = j0 + (j1 - j0) * (left - e0) / (e1 - e0)
            right_j = j0 + (j1 - j0) * (right - e0) / (e1 - e0)
            total += 0.5 * (left_j + right_j) * (right - left)
    return float(total / (upper - lower))
