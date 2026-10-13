#!/usr/bin/env python3
"""Shared, fail-closed utilities for the SA03 observational validation.

The module deliberately contains no AMPS launch side effects.  Preparation,
comparison, plotting, and tests import the same frozen event/channel definitions
and strict readers from here so that standalone modes cannot silently diverge.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shlex
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

import numpy as np


# These bounds come from the packaged SA03 case definition and provenance. They are
# deliberately literal and end-exclusive; deriving them from whichever rows happen
# to be present would allow a truncated data package to redefine the validation case.
START_UTC = "2017-09-06T00:00:00Z"
END_UTC_EXCLUSIVE = "2017-09-10T00:00:00Z"
# SA03 has exactly one eligible validation branch.  RBSP-A remains present in the
# immutable source package for audit, but its 2017 REPT proton data are overwhelmingly
# provider-flagged ``bias_st off``.  Treating it as a validation branch would violate
# the frozen Data_Quality==0 contract; this tuple drives every prepare/run/compare CLI
# default and therefore makes the exclusion structural rather than optional.
SPACECRAFT = ("RBSPB",)
PROTON_REST_MEV = 938.2720813
LIGHT_SPEED_CM_S = 2.99792458e10
M2_TO_CM2 = 1.0e-4

CHANNELS = (
    ("REPT_P1", 18.5, 24.0),
    ("REPT_P2", 24.0, 31.2),
    ("REPT_P3", 31.2, 40.6),
    ("REPT_P4", 40.6, 52.8),
    ("REPT_P5", 52.8, 68.6),
    ("REPT_P6", 68.6, 89.2),
    ("REPT_P7", 89.2, 116.0),
)

PROFILES: dict[str, dict[str, int]] = {
    "QUICK": {"cadence_minutes": 60, "ds_intervals": 24,
              "cutoff_nodes": 24, "angular_samples": 96},
    "STANDARD": {"cadence_minutes": 15, "ds_intervals": 48,
                 "cutoff_nodes": 50, "angular_samples": 384},
    "PUBLICATION": {"cadence_minutes": 5, "ds_intervals": 96,
                    "cutoff_nodes": 120, "angular_samples": 1536},
}

GATES = {
    "median_abs_log10_ratio_max": 0.30,
    "factor_two_fraction_each_min": 0.60,
    "factor_two_fraction_pooled_min": 0.70,
    "max_unresolved_fraction": 0.01,
    "min_matched_epochs": 12,
    "convergence_fractional_change_max": 0.02,
}

EXPECTED_OUTPUTS = (
    "gridless_points_spectrum.dat",
    "gridless_points_flux.dat",
    "gridless_points_density.dat",
    "gridless_termination_summary.dat",
    "standalone_run_manifest.json",
)

_QUOTED = re.compile(r'"([^"]+)"')
_ZONE_TIME = re.compile(r'T="[^"]*?((?:19|20)\d\d-\d\d-\d\dT\d\d:\d\d:\d\dZ?)')


def utc_now() -> str:
    """Use timezone-explicit UTC for provenance, independent of host timezone."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def parse_utc(value: str) -> datetime:
    """Parse ISO timestamps while retaining their explicit timezone information."""
    return datetime.fromisoformat(value.strip().replace("Z", "+00:00"))


def sha256_file(path: Path) -> str:
    """Hash bytes in bounded chunks, including large native numerical products."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    """Write strict JSON atomically; NaN/Infinity are never serialized."""

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
    """Refuse empty comparison artifacts, which otherwise resemble completed work."""
    if not rows:
        raise ValueError(f"refusing to write empty comparison table {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    """Load a nonempty table without silently coercing units or numeric strings."""
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def observation_epoch_inventory(processed: Path, spacecraft: str) -> tuple[set[str], dict[str, Any]]:
    """Return the immutable epochs that can support every SA03 observable.

    The packaged trajectory table is intentionally broader than the comparison
    population: it contains positions for five-minute bins in which only some REPT
    channels survived the frozen quality rules. Row counts alone therefore cannot
    establish comparison eligibility. An eligible epoch must contain each declared
    differential channel exactly once *and* the corresponding spectrum-derived
    density row. This routine validates that relationship before preparation can
    spend time modeling an epoch or comparison can select its reference population.

    Incomplete spectrum epochs remain valid audit evidence and are reported rather
    than treated as corrupt input. Conversely, a density row without exactly seven
    channels (or seven channels without its derived density row) is an inconsistent
    processed package and fails closed.
    """

    prefix = spacecraft.lower()
    spectrum_path = processed / f"{prefix}_rept_spectrum_5min.csv"
    density_path = processed / f"{prefix}_rept_density_5min.csv"
    channels_by_epoch: dict[str, set[int]] = {}
    expected_indices = set(range(len(CHANNELS)))

    for row in read_csv(spectrum_path):
        if row.get("spacecraft") != spacecraft:
            raise ValueError(f"{spectrum_path}: row has spacecraft {row.get('spacecraft')!r}")
        try:
            bounds = (float(row["energy_min_mev"]), float(row["energy_max_mev"]))
        except (KeyError, ValueError) as error:
            raise ValueError(f"{spectrum_path}: invalid REPT energy interval") from error
        matches = [
            index for index, (_, low, high) in enumerate(CHANNELS)
            if abs(bounds[0] - low) < 1.0e-8 and abs(bounds[1] - high) < 1.0e-8
        ]
        if len(matches) != 1:
            raise ValueError(f"{spectrum_path}: unknown REPT interval {bounds}")
        utc = row["utc"]
        present = channels_by_epoch.setdefault(utc, set())
        if matches[0] in present:
            raise ValueError(f"{spectrum_path}: duplicate epoch/channel {utc}/{CHANNELS[matches[0]][0]}")
        present.add(matches[0])

    complete = {utc for utc, present in channels_by_epoch.items() if present == expected_indices}
    incomplete = {
        utc: [CHANNELS[index][0] for index in sorted(expected_indices - present)]
        for utc, present in channels_by_epoch.items() if present != expected_indices
    }

    density_epochs: set[str] = set()
    for row in read_csv(density_path):
        if row.get("spacecraft") != spacecraft:
            raise ValueError(f"{density_path}: row has spacecraft {row.get('spacecraft')!r}")
        utc = row["utc"]
        if utc in density_epochs:
            raise ValueError(f"{density_path}: duplicate utc={utc}")
        density_epochs.add(utc)
        if int(row.get("complete_differential_channels", "-1")) != len(CHANNELS):
            raise ValueError(f"{density_path}: {utc} is not derived from all seven channels")

    if density_epochs != complete:
        missing_density = sorted(complete - density_epochs)
        unsupported_density = sorted(density_epochs - complete)
        raise ValueError(
            f"{spacecraft}: complete-spectrum/density epoch mismatch; "
            f"missing density={missing_density}, density without complete spectrum={unsupported_density}"
        )

    return complete, {
        "complete_epochs": len(complete),
        "incomplete_epochs": len(incomplete),
        "incomplete_epoch_missing_channels": dict(sorted(incomplete.items())),
        "density_epochs": len(density_epochs),
    }


def resolve_data_root(path: Path) -> tuple[Path, Path]:
    """Return the package and authenticated observation root for accepted layouts.

    The externally supplied SA03 location is the package's ``model_inputs``
    directory.  Its pre-rendered SA03 deck is an older, GRIDLESS STANDARD example,
    so the live QUICK/STANDARD workflows must not execute or copy it.  Accept that
    directory as an entry point, resolve its sibling ``validation_data``, and rebuild
    each requested GRIDDED profile from the hash-validated source tables instead.
    Package-root and direct ``validation_data`` forms remain useful for parity with
    SA01/SA02 and for relocatable extracted packages.
    """

    candidate = path.expanduser().resolve()
    if (candidate / "validation_data" / "processed" / "SA03").is_dir():
        return candidate, candidate / "validation_data"
    if (candidate.name == "model_inputs"
            and (candidate.parent / "validation_data" / "processed" / "SA03").is_dir()
            and (candidate / "SA03" / "RBSPB" / "prepared_run.json").is_file()):
        return candidate.parent, candidate.parent / "validation_data"
    if candidate.name == "validation_data" and (candidate / "processed" / "SA03").is_dir():
        return candidate.parent, candidate
    raise ValueError(
        "--data-root/--data-archive must be an extracted package root, its "
        f"model_inputs directory, or its validation_data directory; got {candidate}"
    )


def validate_spacecraft(values: Iterable[str]) -> tuple[str, ...]:
    result: list[str] = []
    for value in values:
        item = value.upper().replace("-", "")
        aliases = {"B": "RBSPB", "RBSPB": "RBSPB"}
        if item not in aliases:
            if item in ("A", "RBSPA"):
                raise ValueError(
                    "RBSPA is source-quality ineligible for SA03: the frozen "
                    "Data_Quality==0 rule leaves no complete seven-channel density epochs"
                )
            raise ValueError(f"unknown spacecraft {value!r}; SA03 permits RBSPB only")
        canonical = aliases[item]
        if canonical not in result:
            result.append(canonical)
    if not result:
        raise ValueError("at least one spacecraft is required")
    return tuple(result)


def validate_data_package(validation_root: Path) -> dict[str, Any]:
    """Audit immutable SA03 hashes, schemas, energy bins, times, and row counts."""

    processed = validation_root / "processed" / "SA03"
    raw = validation_root / "raw" / "SA03"
    provenance_path = processed / "provenance.json"
    manifest_path = raw / "source_manifest.json"
    if not provenance_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError("SA03 provenance/source manifest is missing")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if provenance.get("case_id") != "SA03" or manifest.get("case_id") != "SA03":
        raise ValueError("SA03 manifest identity mismatch")
    window = provenance.get("case_window", {})
    if window.get("start") != START_UTC or window.get("end_exclusive") != END_UTC_EXCLUSIVE:
        raise ValueError("SA03 event window differs from the frozen contract")

    # Preserve the scientific reason that this case has one branch.  Hash validation
    # protects the files themselves; these explicit provenance checks prevent a future
    # configuration edit from quietly promoting the archived, provider-flagged RBSP-A
    # rows into the comparison population.
    rbspa_audit = provenance.get("stages", {}).get("rbsp_rept", {}).get("rbspa", {})
    if (rbspa_audit.get("eligible_trajectory_samples") != 1
            or rbspa_audit.get("spectrum_rows") != 4
            or rbspa_audit.get("complete_density_rows") != 0
            or rbspa_audit.get("quality_flag_counts_at_l_ge_3p2", {}).get("2") != 23742):
        raise ValueError("SA03 RBSP-A source-quality exclusion differs from frozen provenance")

    checked: dict[str, str] = {}
    for group in ("raw_input_sha256", "processed_output_sha256"):
        for relative, expected in provenance[group].items():
            item = validation_root / relative
            if not item.is_file():
                raise FileNotFoundError(item)
            actual = sha256_file(item)
            if actual != expected:
                raise ValueError(f"SHA-256 mismatch: {item}")
            checked[relative] = actual

    raw_manifest_paths = {Path(item["path"]) for item in manifest.get("files", [])}
    actual_raw = {item.relative_to(validation_root) for item in raw.rglob("*") if item.is_file()}
    if actual_raw != raw_manifest_paths | {Path("raw/SA03/source_manifest.json")}:
        raise ValueError("raw/SA03 differs from its exact source manifest")

    bins_expected = [(name, low, high) for name, low, high in CHANNELS]
    spacecraft_summary: dict[str, Any] = {}
    for spacecraft in SPACECRAFT:
        prefix = spacecraft.lower()
        bins = read_csv(processed / f"{prefix}_rept_energy_bins.csv")
        actual_bins = [
            (f"REPT_P{int(row['channel'])}", float(row["energy_min_mev"]),
             float(row["energy_max_mev"])) for row in bins
        ]
        if len(actual_bins) != 7 or any(
            a[0] != b[0] or abs(a[1] - b[1]) > 1e-9 or abs(a[2] - b[2]) > 1e-9
            for a, b in zip(actual_bins, bins_expected)
        ):
            raise ValueError(f"{spacecraft}: REPT bin definition differs from contract")
        spectrum = read_csv(processed / f"{prefix}_rept_spectrum_5min.csv")
        trajectory = list(iter_data_lines(processed / f"{prefix}_trajectory_geo_5min.txt"))
        expected_samples = provenance["stages"]["rbsp_rept"][prefix]["eligible_trajectory_samples"]
        expected_rows = provenance["stages"]["rbsp_rept"][prefix]["spectrum_rows"]
        if len(trajectory) != expected_samples or len(spectrum) != expected_rows:
            raise ValueError(f"{spacecraft}: provenance row-count mismatch")
        for row in spectrum:
            epoch = parse_utc(row["utc"])
            if not (parse_utc(START_UTC) <= epoch < parse_utc(END_UTC_EXCLUSIVE)):
                raise ValueError(f"{spacecraft}: observation outside event window")
            if float(row["l_shell"]) < 3.2 or row.get("data_quality_rule") != "Data_Quality==0":
                raise ValueError(f"{spacecraft}: frozen L/quality rule violation")
        # A trajectory row proves that position and some qualified REPT data exist,
        # not that all seven comparison channels exist. Validate that relationship
        # explicitly so aggregate provenance counts cannot conceal a partial epoch.
        _, observation_inventory = observation_epoch_inventory(processed, spacecraft)
        spacecraft_summary[spacecraft] = {
            "trajectory_samples": len(trajectory), "spectrum_rows": len(spectrum),
            "observation_epoch_inventory": observation_inventory,
        }

    boundary = processed / "boundary_soho_erne_5min.dat"
    boundary_epochs: list[str] = []
    previous: datetime | None = None
    for line in iter_data_lines(boundary):
        fields = line.split()
        if len(fields) != 11:
            raise ValueError(f"{boundary}: expected UTC plus ten energies")
        epoch = parse_utc(fields[0])
        values = [float(value) for value in fields[1:]]
        if previous is not None and epoch <= previous:
            raise ValueError(f"{boundary}: duplicate or unordered epoch")
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise ValueError(f"{boundary}: nonpositive/nonfinite intensity")
        previous = epoch
        boundary_epochs.append(fields[0])
    expected_boundary = provenance["stages"]["soho_erne"]["complete_spectra"]
    if len(boundary_epochs) != expected_boundary:
        raise ValueError("boundary spectrum row count differs from provenance")

    return {
        "status": "PASS", "created_utc": utc_now(), "schema": provenance.get("schema"),
        "provenance_sha256": sha256_file(provenance_path),
        "source_manifest_sha256": sha256_file(manifest_path),
        "hashes_checked": len(checked), "boundary_rows": len(boundary_epochs),
        "spacecraft": spacecraft_summary,
    }


def iter_data_lines(path: Path) -> Iterator[str]:
    with path.open(encoding="utf-8") as stream:
        for raw in stream:
            line = raw.strip()
            if line and not line.startswith(("#", "!")):
                yield line


def relativistic_speed_cm_s(energy_mev: float) -> float:
    if not math.isfinite(energy_mev) or energy_mev <= 0.0:
        raise ValueError("proton kinetic energy must be finite and positive")
    gamma = 1.0 + energy_mev / PROTON_REST_MEV
    return LIGHT_SPEED_CM_S * math.sqrt(1.0 - gamma ** -2)


def isotropic_density_center_cm3(intensities: Sequence[float]) -> float:
    if len(intensities) != len(CHANNELS):
        raise ValueError("density requires all seven REPT channels")
    total = 0.0
    for intensity, (_, low, high) in zip(intensities, CHANNELS):
        if not math.isfinite(intensity) or intensity < 0.0:
            raise ValueError("density intensity must be finite and nonnegative")
        center = 0.5 * (low + high)
        total += intensity * (high - low) / relativistic_speed_cm_s(center)
    return 4.0 * math.pi * total


def _positive_power_integral(e0: float, j0: float, e1: float, j1: float,
                             lower: float, upper: float) -> float:
    slope = math.log(j1 / j0) / math.log(e1 / e0)
    coefficient = j0 / e0 ** slope
    if abs(slope + 1.0) < 1e-12:
        return coefficient * math.log(upper / lower)
    return coefficient * (upper ** (slope + 1.0) - lower ** (slope + 1.0)) / (slope + 1.0)


def _linear_integral(e0: float, j0: float, e1: float, j1: float,
                     lower: float, upper: float) -> float:
    """Integrate a linear transition used only when a segment touches zero."""

    slope = (j1 - j0) / (e1 - e0)
    return j0 * (upper - lower) + 0.5 * slope * ((upper - e0) ** 2 - (lower - e0) ** 2)


def spectrum_integral(energy: Iterable[float], intensity: Iterable[float],
                      lower: float, upper: float,
                      speed_weighted: bool = False) -> float:
    """Integrate a spectrum without deleting or bridging physical-zero nodes.

    Positive-positive segments use the exact power law implied by log-linear
    interpolation.  A segment touching zero uses a local linear transition, so a
    zero stays part of the support instead of being filtered and bridged.  For
    density reconstruction ``J/v`` is evaluated at a controlled fine log grid;
    the same zero-preserving interpolation rules remain in force.
    """

    e = np.asarray(tuple(energy), dtype=float)
    j = np.asarray(tuple(intensity), dtype=float)
    if e.ndim != 1 or j.ndim != 1 or len(e) != len(j) or len(e) < 2:
        raise ValueError("spectrum arrays must have equal one-dimensional length >= 2")
    if np.any(~np.isfinite(e)) or np.any(~np.isfinite(j)):
        raise ValueError("nonfinite modeled spectrum value")
    if np.any(e <= 0.0) or np.any(j < 0.0):
        raise ValueError("energies must be positive and intensities nonnegative")
    order = np.argsort(e)
    e, j = e[order], j[order]
    if np.any(np.diff(e) <= 0.0):
        raise ValueError("modeled spectral energy nodes must be unique")
    tolerance = 1e-10
    if lower < e[0] - tolerance or upper > e[-1] + tolerance or not upper > lower:
        raise ValueError(f"integration interval [{lower}, {upper}] outside model support [{e[0]}, {e[-1]}]")

    edges = [lower, *[float(value) for value in e if lower < value < upper], upper]
    total = 0.0
    for a, b in zip(edges, edges[1:]):
        middle = math.sqrt(a * b)
        index = max(0, min(int(np.searchsorted(e, middle) - 1), len(e) - 2))
        e0, e1, j0, j1 = map(float, (e[index], e[index + 1], j[index], j[index + 1]))
        if speed_weighted:
            # Dense deterministic quadrature isolates the center-speed observation
            # approximation from the native product/reconstruction check.
            grid = np.geomspace(a, b, 65)
            if j0 > 0.0 and j1 > 0.0:
                slope = math.log(j1 / j0) / math.log(e1 / e0)
                vals = j0 * (grid / e0) ** slope
            else:
                vals = j0 + (j1 - j0) * (grid - e0) / (e1 - e0)
            speeds = np.asarray([relativistic_speed_cm_s(float(value)) for value in grid])
            total += float(np.trapz(vals / speeds, grid))
        elif j0 == 0.0 or j1 == 0.0:
            total += _linear_integral(e0, j0, e1, j1, a, b)
        else:
            total += _positive_power_integral(e0, j0, e1, j1, a, b)
    return total


def bin_average(energy: Iterable[float], intensity: Iterable[float],
                lower: float, upper: float) -> float:
    return spectrum_integral(energy, intensity, lower, upper) / (upper - lower)


@dataclass(frozen=True)
class SpectrumZone:
    utc: str
    columns: tuple[str, ...]
    rows: np.ndarray


def _variables(line: str) -> tuple[str, ...]:
    """Keep the declared Tecplot column names/order instead of assuming a schema."""
    result = tuple(_QUOTED.findall(line))
    if not result:
        raise ValueError(f"cannot parse Tecplot VARIABLES line: {line}")
    return result


def read_spectrum_zones(path: Path) -> list[SpectrumZone]:
    """Validate declared widths/counts and bind each spectrum to an explicit epoch.

    GRIDLESS embeds UTC in the zone title; single-location Mode3D snapshots declare
    it in AUXDATA. Conflicting metadata, duplicate epochs and nonfinite nodes are
    errors. In particular, never drop invalid nodes then interpolate across them.
    """
    variables: tuple[str, ...] | None = None
    result: list[SpectrumZone] = []
    current: str | None = None
    expected = 0
    rows: list[list[float]] = []
    snapshot_epoch = None

    def finish() -> None:
        nonlocal current, expected, rows
        if current is None:
            return
        if variables is None or len(rows) != expected:
            raise ValueError(f"{path}: incomplete spectrum zone {current}")
        array = np.asarray(rows, dtype=float)
        if np.any(~np.isfinite(array)):
            raise ValueError(f"{path}: nonfinite spectrum value at {current}")
        result.append(SpectrumZone(current, variables, array))
        current, expected, rows = None, 0, []

    with path.open(encoding="utf-8") as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith("AUXDATA SNAPSHOT_EPOCH_UTC="):
                snapshot_epoch = _quoted_snapshot_epoch(line, path)
            elif line.startswith("VARIABLES"):
                variables = _variables(line)
            elif line.startswith("ZONE"):
                finish()
                match, count = _ZONE_TIME.search(line), re.search(r"\bI\s*=\s*(\d+)", line)
                if (match is None and snapshot_epoch is None) or count is None:
                    raise ValueError(f"{path}: malformed trajectory ZONE {line}")
                current = match.group(1).rstrip("Z") + "Z" if match else snapshot_epoch
                if snapshot_epoch and current != snapshot_epoch:
                    raise ValueError(f"{path}: zone/snapshot epoch mismatch")
                expected = int(count.group(1))
            elif current is not None and line and not line.startswith(("AUXDATA", "#", "!")):
                fields = [float(value) for value in line.split()]
                if variables is None or len(fields) != len(variables):
                    raise ValueError(f"{path}: spectrum row width mismatch")
                rows.append(fields)
    finish()
    if not result or len({zone.utc for zone in result}) != len(result):
        raise ValueError(f"{path}: no zones or duplicate epochs")
    return result


def _quoted_snapshot_epoch(line: str, path: Path) -> str:
    """Read a declared UTC epoch, never infer it from a filename or row order."""
    fields = _QUOTED.findall(line)
    if len(fields) != 1:
        raise ValueError(f"{path}: malformed snapshot epoch")
    parse_utc(fields[0])
    return fields[0].rstrip("Z") + "Z"


def read_point_table(path: Path) -> list[dict[str, str]]:
    """Read point products, accepting snapshot UTC only for one unambiguous row.

    Missing TimeUTC in a multi-location table cannot be repaired by assigning every
    row the snapshot epoch: that would conceal ambiguous location/epoch matching.
    Declared row counts remain mandatory completeness evidence.
    """
    variables: tuple[str, ...] | None = None
    rows: list[dict[str, str]] = []
    snapshot_epoch = None
    declared_rows = 0
    with path.open(encoding="utf-8") as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith("AUXDATA SNAPSHOT_EPOCH_UTC="):
                snapshot_epoch = _quoted_snapshot_epoch(line, path)
            elif line.startswith("VARIABLES"):
                variables = _variables(line)
            elif line.startswith("ZONE"):
                count = re.search(r"\bI\s*=\s*(\d+)", line)
                if count is None:
                    raise ValueError(f"{path}: missing point-zone row count")
                declared_rows += int(count.group(1))
            elif not line or line.startswith(("TITLE", "ZONE", "AUXDATA", "#", "!")):
                continue
            elif variables is not None:
                fields = shlex.split(line)
                if len(fields) != len(variables):
                    raise ValueError(f"{path}: row width {len(fields)} != {len(variables)}")
                row = dict(zip(variables, fields))
                for name, value in row.items():
                    if name != "TimeUTC" and not math.isfinite(float(value)):
                        raise ValueError(f"{path}: nonfinite {name}")
                rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no point rows")
    if declared_rows != len(rows):
        raise ValueError(f"{path}: incomplete point rows ({len(rows)} != {declared_rows})")
    # Mode3D's one-location trajectory snapshots store UTC as AUXDATA. A file
    # containing several untimed locations cannot be treated as a time series.
    if "TimeUTC" not in rows[0] and snapshot_epoch and len(rows) == 1:
        rows[0]["TimeUTC"] = snapshot_epoch
    if "TimeUTC" not in rows[0] or len({row["TimeUTC"] for row in rows}) != len(rows):
        raise ValueError(f"{path}: missing or duplicate TimeUTC")
    if snapshot_epoch and any(row["TimeUTC"].rstrip("Z") + "Z" != snapshot_epoch for row in rows):
        raise ValueError(f"{path}: row/snapshot epoch mismatch")
    return rows


def ratio_metrics(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Metrics whose denominator includes physical-zero model predictions."""

    if not rows:
        return {"count": 0, "positive_model_count": 0, "zero_model_count": 0,
                "median_abs_log10_ratio": None, "median_log10_bias": None,
                "rms_log10_ratio": None, "factor_two_fraction": 0.0}
    ratios: list[float] = []
    for row in rows:
        model, reference = float(row["model"]), float(row["reference"])
        if not math.isfinite(model) or not math.isfinite(reference) or model < 0.0 or reference <= 0.0:
            raise ValueError("comparison requires finite model >= 0 and positive reference")
        ratios.append(model / reference)
    logs = np.asarray([math.log10(value) if value > 0.0 else -math.inf for value in ratios])

    def finite_or_none(value: float) -> float | None:
        return float(value) if math.isfinite(float(value)) else None

    return {
        "count": len(ratios),
        "positive_model_count": sum(value > 0.0 for value in ratios),
        "zero_model_count": sum(value == 0.0 for value in ratios),
        "median_abs_log10_ratio": finite_or_none(float(np.median(np.abs(logs)))),
        "median_log10_bias": finite_or_none(float(np.median(logs))),
        "rms_log10_ratio": finite_or_none(float(np.sqrt(np.mean(logs * logs)))),
        "factor_two_fraction": float(np.mean(np.isfinite(logs) & (np.abs(logs) <= math.log10(2.0)))),
    }


def fractional_change(a: float, b: float, absolute_floor: float = 1e-30) -> float:
    """Symmetric relative change with an explicit absolute criterion near zero."""

    if not math.isfinite(a) or not math.isfinite(b):
        raise ValueError("convergence values must be finite")
    scale = max(abs(a), abs(b))
    if scale <= absolute_floor:
        return 0.0 if abs(a - b) <= absolute_floor else abs(a - b) / absolute_floor
    return abs(a - b) / scale


def command_text(command: Sequence[str]) -> str:
    return shlex.join([str(item) for item in command])
