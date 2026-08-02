#!/usr/bin/env python3
"""Utilities for the C10 POES/MetOp SEM-2 reference-data pipeline.

This module contains the parts of C10 that are independent of AMPS:

* reading the historical NOAA/NCEI 16-second Level-2 SEM-2 files;
* converting sub-satellite positions to AACGM latitude and magnetic local time;
* separating the data into individual polar passes;
* finding the 50-percent-of-polar-cap-flux cutoff on each pass leg; and
* aggregating those individual measurements into the two-hour, one-hour-step
  windows used by the December-2006 POES cutoff-boundary literature.

The code intentionally keeps the measurement-level product.  A binned reference
CSV is convenient for a regression test, but the individual boundary crossings
are the actual observational reference and must not be discarded.

Historical-file support
-----------------------
The December 2006 Level-2 archive is available in both ASCII and CDF.  NOAA's
column names have varied slightly among products, so both readers use alias
sets and report the available columns when a required quantity cannot be found.
The ASCII reader is the preferred path because it has no binary-library
requirements and its fields are documented by NCEI's ``readme_16s_ascii.txt``.

Scientific boundary definition
------------------------------
Dmitriev et al. (2010) define a cutoff latitude as the invariant latitude where
an SEP channel falls to one half of its mean intensity in the polar cap.  C10
implements that definition separately on the inbound and outbound legs of each
polar pass.  The polar-cap plateau is estimated from samples above a configured
absolute AACGM latitude (75 degrees by default).  The crossing is interpolated
between the two adjacent 16-second records that bracket 0.5 times the plateau.

Important limitation
--------------------
MEPED P6-P9 are integral channels, not monochromatic rigidity measurements.  C10
assigns each channel the proton rigidity corresponding to its nominal lower
energy threshold.  This is a transparent and reproducible first-order mapping,
but it does not replace convolution with the full detector response and an
assumed SEP spectrum.  The channel, threshold, and mapping method are written
into every reference row and into the provenance manifest.
"""
from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
import io
import json
import math
import re
import statistics
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

PROTON_REST_MEV = 938.27208816

# The four omnidirectional channels used by Dmitriev et al. for the high-energy
# proton part of the cutoff analysis.  The rigidity is calculated from the
# nominal lower threshold and is therefore a channel label, not an assertion
# that the detector has a delta-function response at this rigidity.
CHANNEL_THRESHOLDS_MEV: Mapping[str, float] = {
    "P6": 16.0,
    "P7": 36.0,
    "P8": 70.0,
    "P9": 140.0,
}

CHANNEL_VARIABLE_ALIASES: Mapping[str, Tuple[str, ...]] = {
    "P6": ("mepomp6", "mep_omni_p6", "mep_omni_flux_p6", "omni_p6", "p6"),
    "P7": ("mepomp7", "mep_omni_p7", "mep_omni_flux_p7", "omni_p7", "p7"),
    "P8": ("mepomp8", "mep_omni_p8", "mep_omni_flux_p8", "omni_p8", "p8"),
    "P9": ("mepomp9", "mep_omni_p9", "mep_omni_flux_p9", "omni_p9", "p9"),
}

SATELLITE_ALIASES: Mapping[str, str] = {
    "n15": "NOAA-15", "noaa15": "NOAA-15", "noaa-15": "NOAA-15",
    "n16": "NOAA-16", "noaa16": "NOAA-16", "noaa-16": "NOAA-16",
    "n17": "NOAA-17", "noaa17": "NOAA-17", "noaa-17": "NOAA-17",
    "n18": "NOAA-18", "noaa18": "NOAA-18", "noaa-18": "NOAA-18",
    "m02": "MetOp-02", "metop02": "MetOp-02", "metop-02": "MetOp-02",
    "metopa": "MetOp-02", "metop-a": "MetOp-02",
}


@dataclass(frozen=True)
class Observation:
    """One valid 16-second SEM-2 Level-2 record."""

    time_utc: datetime
    satellite: str
    geographic_lat_deg: float
    geographic_lon_deg: float
    altitude_km: float
    aacgm_lat_deg: float
    mlt_hour: float
    l_value: Optional[float]
    flux_by_channel: Mapping[str, float]
    source_file: str
    source_sha256: str


@dataclass(frozen=True)
class BoundaryCrossing:
    """One measured half-polar-cap-flux crossing on one pass leg."""

    event_id: str
    satellite: str
    channel: str
    energy_threshold_mev: float
    assigned_rigidity_gv: float
    mapping_method: str
    hemisphere: str
    pass_id: str
    leg: str
    crossing_time_utc: datetime
    geographic_lat_deg: float
    geographic_lon_deg: float
    altitude_km: float
    aacgm_lat_deg: float
    mlt_hour: float
    polar_plateau_flux: float
    half_plateau_flux: float
    low_latitude_flux: Optional[float]
    boundary_uncertainty_deg: float
    n_polar_samples: int
    n_leg_samples: int
    quality_flags: str
    source_file: str
    source_sha256: str


@dataclass(frozen=True)
class ReferenceCell:
    """One window/channel/hemisphere/MLT-bin observational reference cell."""

    event_id: str
    interval_midpoint_utc: datetime
    interval_start_utc: datetime
    interval_end_utc: datetime
    rigidity_gv: float
    energy_threshold_mev: float
    channel: str
    hemisphere: str
    mlt_hour: float
    boundary_aacgm_lat_deg: Optional[float]
    sigma_deg: Optional[float]
    altitude_km: float
    n_crossings: int
    satellites: str
    missing: bool
    source: str
    source_ref: str
    notes: str


def normalize_name(name: str) -> str:
    """Return a conservative lowercase identifier for alias matching."""

    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())


def proton_rigidity_gv_from_kinetic_energy_mev(energy_mev: float) -> float:
    """Relativistic proton rigidity for charge number |Z|=1."""

    if energy_mev < 0.0:
        raise ValueError("kinetic energy must be non-negative")
    momentum_mev_c = math.sqrt(energy_mev * (energy_mev + 2.0 * PROTON_REST_MEV))
    return momentum_mev_c / 1000.0


def sha256_file(path: Path) -> str:
    """Calculate a streaming SHA-256 digest for provenance records."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def infer_satellite(path: Path) -> str:
    """Infer the spacecraft from a historical NCEI file name.

    Historical files normally contain ``n15``, ``n16``, ``n17``, ``n18``, or
    ``m02``.  A hard failure is preferable to silently mixing an unidentified
    spacecraft into the reference product.
    """

    token = normalize_name(path.name)
    for alias, canonical in SATELLITE_ALIASES.items():
        if normalize_name(alias) in token:
            return canonical
    raise ValueError(f"cannot infer satellite from file name: {path.name}")


def _open_text(path: Path) -> io.TextIOBase:
    """Open plain or gzip-compressed historical ASCII data as text."""

    if path.suffix.lower() == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _find_alias(columns: Sequence[str], aliases: Sequence[str], required: bool = True) -> Optional[str]:
    """Return the original column name matching one of ``aliases``."""

    normalized = {normalize_name(column): column for column in columns}
    for alias in aliases:
        key = normalize_name(alias)
        if key in normalized:
            return normalized[key]
    if required:
        raise KeyError(
            "none of aliases %s is present; available columns are: %s"
            % (", ".join(aliases), ", ".join(columns))
        )
    return None


def _float_or_none(value: object) -> Optional[float]:
    """Parse a numeric field while rejecting archive fill and impossible values."""

    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if not math.isfinite(result):
        return None
    # Historical SEM-2 products use negative fill values, and an old unpacker
    # could emit the unphysical count 19988488 for missing P8/P9 records.
    if result < 0.0 or result >= 1.0e7:
        return None
    return result


def _parse_ascii_header(lines: Sequence[str]) -> Tuple[int, List[str]]:
    """Locate the Level-2 column-name line in a small prefix of an ASCII file."""

    required = {"year", "sslat", "sslon"}
    for index, raw in enumerate(lines):
        stripped = raw.strip().lstrip("#").strip()
        if not stripped:
            continue
        columns = re.split(r"[\s,]+", stripped)
        normalized = {normalize_name(column) for column in columns}
        if required.issubset(normalized) and any(normalize_name(c) == "mepomp6" for c in columns):
            return index, columns
    raise ValueError(
        "could not locate the 16-second Level-2 header. Expected columns such "
        "as year, sslat, sslon, and mepomp6."
    )


def read_ascii_level2(path: Path, default_altitude_km: float = 850.0) -> List[Dict[str, object]]:
    """Read one historical NCEI Level-2 16-second ASCII file.

    The return value intentionally contains geographic coordinates and raw
    channel values only.  AACGM conversion is applied in a separate function so
    it can be tested and replaced independently.
    """

    with _open_text(path) as stream:
        raw_lines = stream.readlines()
    header_index, columns = _parse_ascii_header(raw_lines[:200])

    year_col = _find_alias(columns, ("year", "yr"))
    month_col = _find_alias(columns, ("mo", "month"))
    day_col = _find_alias(columns, ("dy", "day"))
    hour_col = _find_alias(columns, ("hr", "hour"))
    minute_col = _find_alias(columns, ("mi", "minute"))
    second_col = _find_alias(columns, ("second", "sec"))
    lat_col = _find_alias(columns, ("sslat", "satlat", "latitude"))
    lon_col = _find_alias(columns, ("sslon", "satlon", "longitude"))
    altitude_col = _find_alias(columns, ("alt", "altitude", "satalt", "height"), required=False)
    mlt_col = _find_alias(columns, ("mlt", "magneticlocaltime"), required=False)
    l_col = _find_alias(columns, ("lval", "lshell", "lvalue"), required=False)
    channel_columns = {
        channel: _find_alias(columns, aliases)
        for channel, aliases in CHANNEL_VARIABLE_ALIASES.items()
    }

    rows: List[Dict[str, object]] = []
    for line_number, raw in enumerate(raw_lines[header_index + 1 :], start=header_index + 2):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = re.split(r"[\s,]+", stripped)
        if len(fields) < len(columns):
            # Truncated daily files are documented in the archive.  Skipping a
            # malformed record and recording file-level provenance is safer than
            # shifting every later column by accident.
            continue
        record = dict(zip(columns, fields[: len(columns)]))
        try:
            second = float(record[second_col])
            whole_second = int(math.floor(second))
            microsecond = int(round((second - whole_second) * 1.0e6))
            when = datetime(
                int(float(record[year_col])), int(float(record[month_col])),
                int(float(record[day_col])), int(float(record[hour_col])),
                int(float(record[minute_col])), whole_second, microsecond,
                tzinfo=timezone.utc,
            )
            lat = float(record[lat_col])
            lon = float(record[lon_col])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number}: invalid time/position record: {exc}") from exc

        flux = {channel: _float_or_none(record[column]) for channel, column in channel_columns.items()}
        rows.append({
            "time_utc": when,
            "geographic_lat_deg": lat,
            "geographic_lon_deg": lon,
            "altitude_km": float(record[altitude_col]) if altitude_col and _float_or_none(record[altitude_col]) is not None else default_altitude_km,
            "archive_mlt_hour": _float_or_none(record[mlt_col]) if mlt_col else None,
            "l_value": _float_or_none(record[l_col]) if l_col else None,
            "flux_by_channel": flux,
        })
    return rows


def _read_cdf_variable(cdf: object, variable: str) -> object:
    """Small adapter around cdflib's variable getter for easier testing."""

    return cdf.varget(variable)  # type: ignore[attr-defined]


def read_cdf_level2(path: Path, default_altitude_km: float = 850.0) -> List[Dict[str, object]]:
    """Read one historical Level-2 CDF file using ``cdflib``.

    CDF naming differs among processing generations, so this reader discovers
    variables through the same alias mechanism as the ASCII reader.  For the
    December-2006 test the ASCII format remains the recommended and best-tested
    route; CDF support is provided for users who maintain a CDF mirror.
    """

    try:
        import cdflib  # type: ignore
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("CDF input requires cdflib and numpy") from exc

    cdf = cdflib.CDF(str(path))
    info = cdf.cdf_info()
    variables = list(getattr(info, "zVariables", []) or []) + list(getattr(info, "rVariables", []) or [])
    if not variables and isinstance(info, dict):
        variables = list(info.get("zVariables", [])) + list(info.get("rVariables", []))

    def alias(names: Sequence[str], required: bool = True) -> Optional[str]:
        return _find_alias(variables, names, required=required)

    epoch_name = alias(("epoch", "time", "datetime", "timestamp"), required=False)
    lat_name = alias(("sslat", "satlat", "latitude"))
    lon_name = alias(("sslon", "satlon", "longitude"))
    alt_name = alias(("alt", "altitude", "satalt", "height"), required=False)
    mlt_name = alias(("mlt", "magneticlocaltime"), required=False)
    l_name = alias(("lval", "lshell", "lvalue"), required=False)
    channel_names = {channel: alias(names) for channel, names in CHANNEL_VARIABLE_ALIASES.items()}

    if epoch_name:
        epochs = _read_cdf_variable(cdf, epoch_name)
        times = cdflib.cdfepoch.to_datetime(epochs)
        times = [
            value.astype("datetime64[us]").astype(datetime).replace(tzinfo=timezone.utc)
            if hasattr(value, "astype") else value.replace(tzinfo=timezone.utc)
            for value in times
        ]
    else:
        y = _read_cdf_variable(cdf, alias(("year", "yr")))
        mo = _read_cdf_variable(cdf, alias(("mo", "month")))
        d = _read_cdf_variable(cdf, alias(("dy", "day")))
        h = _read_cdf_variable(cdf, alias(("hr", "hour")))
        mi = _read_cdf_variable(cdf, alias(("mi", "minute")))
        sec = _read_cdf_variable(cdf, alias(("second", "sec")))
        times = [
            datetime(int(y[i]), int(mo[i]), int(d[i]), int(h[i]), int(mi[i]), int(float(sec[i])), tzinfo=timezone.utc)
            for i in range(len(y))
        ]

    latitudes = np.asarray(_read_cdf_variable(cdf, lat_name)).reshape(-1)
    longitudes = np.asarray(_read_cdf_variable(cdf, lon_name)).reshape(-1)
    altitudes = np.asarray(_read_cdf_variable(cdf, alt_name)).reshape(-1) if alt_name else None
    mlts = np.asarray(_read_cdf_variable(cdf, mlt_name)).reshape(-1) if mlt_name else None
    lvals = np.asarray(_read_cdf_variable(cdf, l_name)).reshape(-1) if l_name else None
    channels = {name: np.asarray(_read_cdf_variable(cdf, variable)).reshape(-1) for name, variable in channel_names.items()}

    rows: List[Dict[str, object]] = []
    for index, when in enumerate(times):
        rows.append({
            "time_utc": when,
            "geographic_lat_deg": float(latitudes[index]),
            "geographic_lon_deg": float(longitudes[index]),
            "altitude_km": float(altitudes[index]) if altitudes is not None and _float_or_none(altitudes[index]) is not None else default_altitude_km,
            "archive_mlt_hour": _float_or_none(mlts[index]) if mlts is not None else None,
            "l_value": _float_or_none(lvals[index]) if lvals is not None else None,
            "flux_by_channel": {channel: _float_or_none(values[index]) for channel, values in channels.items()},
        })
    return rows


def read_level2_file(path: Path, default_altitude_km: float = 850.0) -> List[Dict[str, object]]:
    """Dispatch to the ASCII or CDF reader based on the file suffix."""

    lower = path.name.lower()
    if lower.endswith((".txt", ".txt.gz", ".asc", ".asc.gz", ".dat", ".dat.gz")):
        return read_ascii_level2(path, default_altitude_km)
    if lower.endswith(".cdf"):
        return read_cdf_level2(path, default_altitude_km)
    raise ValueError(f"unsupported Level-2 file type: {path}")


def add_aacgm_coordinates(raw_rows: Sequence[Mapping[str, object]], path: Path) -> List[Observation]:
    """Convert geographic positions to AACGM latitude and MLT.

    ``aacgmv2.get_aacgm_coord`` returns AACGM latitude, AACGM longitude, and
    magnetic local time.  Records outside the library's supported domain are
    skipped instead of being assigned archive magnetic coordinates generated by
    a different field model; this keeps model and observation postprocessing
    internally consistent.
    """

    try:
        import aacgmv2  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("AACGM conversion requires aacgmv2") from exc

    satellite = infer_satellite(path)
    digest = sha256_file(path)
    observations: List[Observation] = []
    for row in raw_rows:
        when = row["time_utc"]
        lat = float(row["geographic_lat_deg"])
        lon = float(row["geographic_lon_deg"])
        altitude = float(row["altitude_km"])
        try:
            aacgm_lat, _aacgm_lon, mlt = aacgmv2.get_aacgm_coord(lat, lon, altitude, when)
        except Exception:
            continue
        if not all(math.isfinite(float(value)) for value in (aacgm_lat, mlt)):
            continue
        clean_flux = {
            channel: float(value)
            for channel, value in dict(row["flux_by_channel"]).items()
            if value is not None and math.isfinite(float(value)) and float(value) > 0.0
        }
        observations.append(Observation(
            time_utc=when,
            satellite=satellite,
            geographic_lat_deg=lat,
            geographic_lon_deg=lon,
            altitude_km=altitude,
            aacgm_lat_deg=float(aacgm_lat),
            mlt_hour=float(mlt) % 24.0,
            l_value=None if row.get("l_value") is None else float(row["l_value"]),
            flux_by_channel=clean_flux,
            source_file=path.name,
            source_sha256=digest,
        ))
    return observations


def load_observations(paths: Sequence[Path], default_altitude_km: float = 850.0) -> List[Observation]:
    """Read, convert, concatenate, and time-sort multiple daily files."""

    observations: List[Observation] = []
    for path in paths:
        raw = read_level2_file(path, default_altitude_km)
        observations.extend(add_aacgm_coordinates(raw, path))
    observations.sort(key=lambda row: (row.satellite, row.time_utc))
    return observations


def split_polar_passes(
    observations: Sequence[Observation],
    minimum_abs_lat_deg: float = 45.0,
    maximum_gap_seconds: float = 90.0,
    required_polar_abs_lat_deg: float = 75.0,
) -> List[List[Observation]]:
    """Split observations into contiguous single-hemisphere polar passes."""

    passes: List[List[Observation]] = []
    current: List[Observation] = []
    previous: Optional[Observation] = None

    def finish() -> None:
        nonlocal current
        if current and max(abs(row.aacgm_lat_deg) for row in current) >= required_polar_abs_lat_deg:
            passes.append(current)
        current = []

    for row in observations:
        if abs(row.aacgm_lat_deg) < minimum_abs_lat_deg:
            finish()
            previous = row
            continue
        hemisphere = 1 if row.aacgm_lat_deg >= 0.0 else -1
        discontinuity = False
        if previous is not None:
            discontinuity = (
                row.satellite != previous.satellite
                or (1 if previous.aacgm_lat_deg >= 0.0 else -1) != hemisphere
                or (row.time_utc - previous.time_utc).total_seconds() > maximum_gap_seconds
            )
        if discontinuity:
            finish()
        current.append(row)
        previous = row
    finish()
    return passes


def _interpolate_circular_hour(a: float, b: float, fraction: float) -> float:
    """Linearly interpolate MLT while respecting the 24/0-hour seam."""

    delta = ((b - a + 12.0) % 24.0) - 12.0
    return (a + fraction * delta) % 24.0


def _interpolate_crossing(
    a: Observation,
    b: Observation,
    fa: float,
    fb: float,
    target: float,
) -> Tuple[float, datetime, float, float, float, float, float]:
    """Interpolate location and time between two records bracketing ``target``."""

    if math.isclose(fa, fb):
        fraction = 0.5
    else:
        fraction = (target - fa) / (fb - fa)
    fraction = min(1.0, max(0.0, fraction))
    duration = b.time_utc - a.time_utc
    when = a.time_utc + timedelta(seconds=duration.total_seconds() * fraction)
    aacgm_lat = a.aacgm_lat_deg + fraction * (b.aacgm_lat_deg - a.aacgm_lat_deg)
    geo_lat = a.geographic_lat_deg + fraction * (b.geographic_lat_deg - a.geographic_lat_deg)
    lon_delta = ((b.geographic_lon_deg - a.geographic_lon_deg + 180.0) % 360.0) - 180.0
    geo_lon = ((a.geographic_lon_deg + fraction * lon_delta + 180.0) % 360.0) - 180.0
    altitude = a.altitude_km + fraction * (b.altitude_km - a.altitude_km)
    mlt = _interpolate_circular_hour(a.mlt_hour, b.mlt_hour, fraction)
    return fraction, when, geo_lat, geo_lon, altitude, aacgm_lat, mlt


def _find_leg_crossing(leg: Sequence[Observation], channel: str, target: float) -> Optional[Tuple[Observation, Observation, float, float]]:
    """Find the equatorward-most adjacent pair that brackets the half-flux level."""

    candidates: List[Tuple[float, Observation, Observation, float, float]] = []
    for first, second in zip(leg[:-1], leg[1:]):
        f1 = first.flux_by_channel.get(channel)
        f2 = second.flux_by_channel.get(channel)
        if f1 is None or f2 is None:
            continue
        if (f1 - target) * (f2 - target) <= 0.0 and not math.isclose(f1, f2):
            equatorward_lat = min(abs(first.aacgm_lat_deg), abs(second.aacgm_lat_deg))
            candidates.append((equatorward_lat, first, second, f1, f2))
    if not candidates:
        return None
    _, first, second, f1, f2 = min(candidates, key=lambda item: item[0])
    return first, second, f1, f2


def extract_boundary_crossings(
    observations: Sequence[Observation],
    event_id: str = "STORM_2006_12",
    minimum_abs_lat_deg: float = 45.0,
    polar_plateau_abs_lat_deg: float = 75.0,
    minimum_polar_samples: int = 4,
    minimum_plateau_to_low_ratio: float = 2.0,
    minimum_leg_samples: int = 4,
) -> List[BoundaryCrossing]:
    """Extract half-polar-cap-flux boundaries from all valid polar passes."""

    passes = split_polar_passes(
        observations,
        minimum_abs_lat_deg=minimum_abs_lat_deg,
        required_polar_abs_lat_deg=polar_plateau_abs_lat_deg,
    )
    crossings: List[BoundaryCrossing] = []

    for pass_index, polar_pass in enumerate(passes):
        peak_index = max(range(len(polar_pass)), key=lambda i: abs(polar_pass[i].aacgm_lat_deg))
        legs = {
            "INBOUND": polar_pass[: peak_index + 1],
            "OUTBOUND": polar_pass[peak_index:],
        }
        hemisphere = "N" if polar_pass[peak_index].aacgm_lat_deg >= 0.0 else "S"
        pass_id = "%s_%s_%04d" % (
            polar_pass[0].satellite.replace("-", ""),
            polar_pass[peak_index].time_utc.strftime("%Y%m%dT%H%M%S"),
            pass_index,
        )

        for channel, threshold_mev in CHANNEL_THRESHOLDS_MEV.items():
            polar_fluxes = [
                row.flux_by_channel[channel]
                for row in polar_pass
                if abs(row.aacgm_lat_deg) >= polar_plateau_abs_lat_deg
                and channel in row.flux_by_channel
            ]
            if len(polar_fluxes) < minimum_polar_samples:
                continue
            plateau = statistics.median(polar_fluxes)
            if not math.isfinite(plateau) or plateau <= 0.0:
                continue
            half_level = 0.5 * plateau
            low_fluxes = [
                row.flux_by_channel[channel]
                for row in polar_pass
                if abs(row.aacgm_lat_deg) <= minimum_abs_lat_deg + 5.0
                and channel in row.flux_by_channel
            ]
            low_flux = statistics.median(low_fluxes) if low_fluxes else None
            if low_flux is not None and plateau / max(low_flux, 1.0e-30) < minimum_plateau_to_low_ratio:
                continue

            for leg_name, leg in legs.items():
                if len(leg) < minimum_leg_samples:
                    continue
                bracket = _find_leg_crossing(leg, channel, half_level)
                if bracket is None:
                    continue
                first, second, f1, f2 = bracket
                (_fraction, when, geo_lat, geo_lon, altitude,
                 aacgm_lat, mlt) = _interpolate_crossing(first, second, f1, f2, half_level)
                spacing = abs(abs(second.aacgm_lat_deg) - abs(first.aacgm_lat_deg))
                uncertainty = max(0.20, 0.5 * spacing)
                flags: List[str] = []
                if low_flux is None:
                    flags.append("NO_LOW_LATITUDE_BACKGROUND")
                if channel in ("P8", "P9"):
                    flags.append("P8_P9_SUBCOMMUTATED_ARCHIVE_CHANNEL")

                crossings.append(BoundaryCrossing(
                    event_id=event_id,
                    satellite=first.satellite,
                    channel=channel,
                    energy_threshold_mev=threshold_mev,
                    assigned_rigidity_gv=proton_rigidity_gv_from_kinetic_energy_mev(threshold_mev),
                    mapping_method="NOMINAL_INTEGRAL_CHANNEL_LOWER_THRESHOLD",
                    hemisphere=hemisphere,
                    pass_id=pass_id,
                    leg=leg_name,
                    crossing_time_utc=when,
                    geographic_lat_deg=geo_lat,
                    geographic_lon_deg=geo_lon,
                    altitude_km=altitude,
                    aacgm_lat_deg=aacgm_lat,
                    mlt_hour=mlt,
                    polar_plateau_flux=plateau,
                    half_plateau_flux=half_level,
                    low_latitude_flux=low_flux,
                    boundary_uncertainty_deg=uncertainty,
                    n_polar_samples=len(polar_fluxes),
                    n_leg_samples=len(leg),
                    quality_flags=";".join(flags),
                    source_file=first.source_file,
                    source_sha256=first.source_sha256,
                ))
    crossings.sort(key=lambda row: (row.crossing_time_utc, row.channel, row.hemisphere, row.satellite, row.leg))
    return crossings


def nearest_mlt_bin(mlt_hour: float, bin_centers: Sequence[float]) -> float:
    """Return the circularly nearest configured MLT-bin center."""

    return min(bin_centers, key=lambda center: abs(((mlt_hour - center + 12.0) % 24.0) - 12.0))


def hourly_window_midpoints(start: datetime, end: datetime, step_hours: float = 1.0) -> List[datetime]:
    """Generate whole-step UTC window centers spanning an event interval."""

    first = start.replace(minute=0, second=0, microsecond=0)
    if first < start:
        first += timedelta(hours=step_hours)
    step = timedelta(hours=step_hours)
    result: List[datetime] = []
    current = first
    while current <= end:
        result.append(current)
        current += step
    return result


def aggregate_crossings(
    crossings: Sequence[BoundaryCrossing],
    event_start: datetime,
    event_end: datetime,
    window_hours: float = 2.0,
    step_hours: float = 1.0,
    mlt_bin_centers: Sequence[float] = tuple(float(v) for v in range(0, 24, 3)),
    minimum_crossings_per_cell: int = 1,
    source_ref: str = "Dmitriev et al. (2010), doi:10.1029/2010JA015380",
) -> List[ReferenceCell]:
    """Aggregate individual measured crossings into the C10 reference grid."""

    half_window = timedelta(hours=0.5 * window_hours)
    midpoints = hourly_window_midpoints(event_start, event_end, step_hours)
    cells: List[ReferenceCell] = []

    for midpoint in midpoints:
        start = midpoint - half_window
        end = midpoint + half_window
        for channel, threshold_mev in CHANNEL_THRESHOLDS_MEV.items():
            rigidity = proton_rigidity_gv_from_kinetic_energy_mev(threshold_mev)
            for hemisphere in ("N", "S"):
                for mlt_center in mlt_bin_centers:
                    selected = [
                        row for row in crossings
                        if row.channel == channel
                        and row.hemisphere == hemisphere
                        and start <= row.crossing_time_utc < end
                        and math.isclose(nearest_mlt_bin(row.mlt_hour, mlt_bin_centers), mlt_center)
                    ]
                    valid = len(selected) >= minimum_crossings_per_cell
                    latitudes = [abs(row.aacgm_lat_deg) for row in selected]
                    altitudes = [row.altitude_km for row in selected]
                    if valid:
                        boundary = statistics.median(latitudes)
                        if len(latitudes) > 1:
                            robust_sigma = 1.4826 * statistics.median(abs(value - boundary) for value in latitudes)
                        else:
                            robust_sigma = selected[0].boundary_uncertainty_deg
                        sigma = max(0.20, robust_sigma, statistics.fmean(row.boundary_uncertainty_deg for row in selected))
                        altitude = statistics.fmean(altitudes)
                    else:
                        boundary = sigma = None
                        altitude = statistics.fmean(altitudes) if altitudes else 850.0
                    satellites = ";".join(sorted({row.satellite for row in selected}))
                    cells.append(ReferenceCell(
                        event_id="STORM_2006_12",
                        interval_midpoint_utc=midpoint,
                        interval_start_utc=start,
                        interval_end_utc=end,
                        rigidity_gv=rigidity,
                        energy_threshold_mev=threshold_mev,
                        channel=channel,
                        hemisphere=hemisphere,
                        mlt_hour=float(mlt_center),
                        boundary_aacgm_lat_deg=boundary,
                        sigma_deg=sigma,
                        altitude_km=altitude,
                        n_crossings=len(selected),
                        satellites=satellites,
                        missing=not valid,
                        source="POES_NCEI_LEVEL2_16SEC",
                        source_ref=source_ref,
                        notes=(
                            "Median of independent inbound/outbound half-polar-cap-flux crossings; "
                            "rigidity is mapped from the channel nominal lower threshold."
                        ),
                    ))
    return cells


def _format_utc(value: datetime) -> str:
    return value.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_crossings_csv(crossings: Sequence[BoundaryCrossing], output: Path) -> None:
    """Write the primary, pass-level observational product."""

    output.parent.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, object]] = []
    for crossing in crossings:
        row = asdict(crossing)
        row["crossing_time_utc"] = _format_utc(crossing.crossing_time_utc)
        rows.append(row)
    fieldnames = list(rows[0].keys()) if rows else [field.name for field in BoundaryCrossing.__dataclass_fields__.values()]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


@contextlib.contextmanager
def _open_reference_text_output(output: Path):
    """Open a reference destination as deterministic UTF-8 text.

    Gzip normally embeds the current wall-clock time in its header.  Setting
    ``mtime=0`` and omitting the original filename makes repeated builds from
    identical measurements byte-for-byte reproducible, so the SHA-256 recorded
    by the reference summary is meaningful.
    """

    if output.name.lower().endswith(".gz"):
        with output.open("wb") as raw_stream:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_stream, mtime=0
            ) as compressed_stream:
                with io.TextIOWrapper(
                    compressed_stream, encoding="utf-8", newline=""
                ) as text_stream:
                    yield text_stream
    else:
        with output.open("w", encoding="utf-8", newline="") as text_stream:
            yield text_stream


def write_reference_csv(cells: Sequence[ReferenceCell], output: Path, manifest_sha256: str = "") -> None:
    """Write the windowed reference consumed by ``run_C10.py``.

    The production C10 reference is stored as
    ``reference_C10_poes_meped_boundary.csv.gz``.  A ``*.gz`` destination is
    written directly as deterministic UTF-8 gzip-compressed CSV; plain CSV
    remains supported for small developer fixtures.
    """

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "event_id", "interval_midpoint_utc", "interval_start_utc", "interval_end_utc",
        "rigidity_gv", "energy_threshold_mev", "channel", "hemisphere", "mlt_hour",
        "boundary_aacgm_lat_deg", "sigma_deg", "altitude_km", "n_crossings",
        "satellites", "missing", "source", "source_ref", "notes",
    ]
    with _open_reference_text_output(output) as stream:
        stream.write("# C10 real POES/MetOp SEM-2 MEPED cutoff-boundary reference\n")
        stream.write("# reference_kind=POES_NCEI_LEVEL2_16SEC\n")
        stream.write("# boundary_definition=50_percent_of_median_polar_cap_flux\n")
        stream.write("# rigidity_mapping=nominal_integral_channel_lower_threshold\n")
        if manifest_sha256:
            stream.write(f"# provenance_manifest_sha256={manifest_sha256}\n")
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for cell in cells:
            writer.writerow({
                "event_id": cell.event_id,
                "interval_midpoint_utc": _format_utc(cell.interval_midpoint_utc),
                "interval_start_utc": _format_utc(cell.interval_start_utc),
                "interval_end_utc": _format_utc(cell.interval_end_utc),
                "rigidity_gv": f"{cell.rigidity_gv:.9f}",
                "energy_threshold_mev": f"{cell.energy_threshold_mev:.3f}",
                "channel": cell.channel,
                "hemisphere": cell.hemisphere,
                "mlt_hour": f"{cell.mlt_hour:.1f}",
                "boundary_aacgm_lat_deg": "" if cell.boundary_aacgm_lat_deg is None else f"{cell.boundary_aacgm_lat_deg:.5f}",
                "sigma_deg": "" if cell.sigma_deg is None else f"{cell.sigma_deg:.5f}",
                "altitude_km": f"{cell.altitude_km:.3f}",
                "n_crossings": cell.n_crossings,
                "satellites": cell.satellites,
                "missing": "TRUE" if cell.missing else "FALSE",
                "source": cell.source,
                "source_ref": cell.source_ref,
                "notes": cell.notes,
            })


def write_manifest(
    source_files: Sequence[Path],
    crossings: Sequence[BoundaryCrossing],
    cells: Sequence[ReferenceCell],
    configuration: Mapping[str, object],
    output: Path,
) -> str:
    """Write a complete machine-readable provenance manifest and return its SHA-256."""

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "product": "C10_POES_NCEI_REFERENCE",
        "created_utc": _format_utc(datetime.now(timezone.utc)),
        "configuration": dict(configuration),
        "source_files": [
            {"path": str(path), "name": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
            for path in source_files
        ],
        "channel_mapping": {
            channel: {
                "nominal_lower_threshold_mev": energy,
                "assigned_rigidity_gv": proton_rigidity_gv_from_kinetic_energy_mev(energy),
            }
            for channel, energy in CHANNEL_THRESHOLDS_MEV.items()
        },
        "n_crossings": len(crossings),
        "n_reference_cells": len(cells),
        "n_nonmissing_reference_cells": sum(not cell.missing for cell in cells),
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")
    return sha256_file(output)
