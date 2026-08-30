#!/usr/bin/env python3
"""C15 — T05 driver interpolation and epoch reproducibility validation.

C15 performs independent Mode3D snapshots rather than relying on proposed
``TIME_SERIES`` input keywords.  At each selected epoch it runs T05 twice:

* once with the complete five-minute driver; and
* once with an independently materialized driver containing the analytically
  interpolated state as an exact row at the requested epoch.

The initialized fields and exact-rigidity access products must agree.  Additional
same-epoch repeats, scheduler variants, an epoch-invariant centered-dipole
control, a perturbed-driver sensitivity case, and field-continuity checks make
the test useful for both driver-reader and model-validation regressions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple


TEST_ID = "C15"
TEST_NAME = "T05 driver interpolation and epoch reproducibility"
RUNNER_SCHEMA_VERSION = 2
RUNNER_RELEASE = "2026-08-30"

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_TEMPLATE = SCRIPT_DIR / "AMPS_PARAM_C15_mode3d.in"
DEFAULT_DRIVER = SCRIPT_DIR / "data" / "ts05_driver_C15.txt"
DEFAULT_REFERENCE = SCRIPT_DIR / "reference_C15_driver_interpolation.csv"

EARTH_RADIUS_KM = 6371.2
PROTON_REST_MEV = 938.27208816
DRIVER_FIELDS = (
    "Bx", "By", "Bz", "Vx", "Vy", "Vz", "Np", "Temp", "DST",
    "IMFflag", "SWflag", "Tilt", "Pdyn", "W1", "W2", "W3", "W4",
    "W5", "W6",
)

# This spelling is part of the executable-facing contract, not just a display
# label.  Current AMPS driver-table parsing normalizes several legacy names
# (for example Np -> DEN_P) but requires the geomagnetic-index column to be
# named DST.  In particular, SYM-H/SYM_H is not accepted as an alias.  Keeping
# one canonical header here prevents generated reference drivers from drifting
# away from the checked-in full driver.
AMPS_DRIVER_HEADER = ("YYYY-MM-DDTHH:MM:SS",) + DRIVER_FIELDS

# Profiles control cost and temporal sampling only.  All profiles retain the
# exact interpolation, repeat, scheduler, dipole, and sensitivity contracts.
PROFILE_DEFAULTS = {
    "SMOKE": {
        "epochs": (
            "2012-05-17T05:55:00",
            "2012-05-17T05:57:30",
            "2012-05-17T06:00:00",
        ),
        "repeats": 2,
        "cross_schedulers": ("STATIC",),
        "shell_res_deg": 60.0,
        "mesh_res_earth_re": 0.75,
        "mesh_res_boundary_re": 2.0,
    },
    "ROUTINE": {
        "epochs": (
            "2012-05-17T05:55:00",
            "2012-05-17T05:57:30",
            "2012-05-17T06:00:00",
            "2012-05-17T06:02:30",
            "2012-05-17T06:05:00",
        ),
        "repeats": 2,
        "cross_schedulers": ("STATIC",),
        "shell_res_deg": 30.0,
        "mesh_res_earth_re": 0.50,
        "mesh_res_boundary_re": 1.50,
    },
    "THOROUGH": {
        "epochs": (
            "2012-05-17T05:55:00",
            "2012-05-17T05:56:15",
            "2012-05-17T05:57:30",
            "2012-05-17T05:58:45",
            "2012-05-17T06:00:00",
            "2012-05-17T06:01:15",
            "2012-05-17T06:02:30",
            "2012-05-17T06:03:45",
            "2012-05-17T06:05:00",
        ),
        "repeats": 3,
        "cross_schedulers": ("STATIC", "BLOCK_CYCLIC"),
        "shell_res_deg": 15.0,
        "mesh_res_earth_re": 0.35,
        "mesh_res_boundary_re": 1.0,
    },
}


@dataclass(frozen=True)
class DriverRow:
    """One complete AMPS T05 driver record."""

    epoch: datetime
    values: Tuple[float, ...]


@dataclass(frozen=True)
class InterpolatedDriver:
    """An independently calculated driver state and its interpolation metadata."""

    epoch: datetime
    values: Tuple[float, ...]
    left_epoch: datetime
    right_epoch: datetime
    weight: float

    @property
    def case_kind(self) -> str:
        return "exact" if self.left_epoch == self.right_epoch else "interpolated"


@dataclass
class RunCase:
    """One AMPS invocation in the C15 validation matrix."""

    case_id: str
    category: str
    model: str
    epoch: datetime
    driver_mode: str
    scheduler: str
    run_directory: Path
    paired_epoch: str = ""


@dataclass(frozen=True)
class Fingerprint:
    """Canonical binary64 fingerprint of all numeric rows in an output set."""

    path_names: Tuple[str, ...]
    numeric_rows: int
    numeric_values: int
    sha256: str


@dataclass
class RunRecord:
    """Execution and output-discovery result for one C15 case."""

    case_id: str
    category: str
    model: str
    epoch_utc: str
    driver_mode: str
    scheduler: str
    run_directory: str
    command: List[str]
    return_code: int
    field_fingerprint: Optional[Fingerprint]
    cutoff_fingerprint: Optional[Fingerprint]
    passed: bool
    error: str = ""


@dataclass(frozen=True)
class NumericComparison:
    """Whole-product numeric comparison used by reference equivalence gates."""

    values_compared: int
    mismatches: int
    max_abs_difference: float
    max_relative_difference: float
    rms_difference: float
    changed_value_relative_rms: float
    passed: bool
    error: str = ""


def parse_utc(text: str) -> datetime:
    """Parse ISO UTC and normalize to a timezone-naive UTC datetime."""
    value = datetime.fromisoformat(str(text).strip().replace("Z", "+00:00"))
    if value.tzinfo is not None:
        value = value.astimezone(timezone.utc).replace(tzinfo=None)
    return value


def format_utc(value: datetime) -> str:
    """Write the exact parser-compatible timestamp used in generated decks."""
    return value.replace(tzinfo=None).isoformat(timespec="seconds")


def epoch_label(value: datetime) -> str:
    """Return a sortable path label without punctuation."""
    return value.strftime("%Y%m%dT%H%M%S")


def sha256_file(path: Path) -> str:
    """Calculate a streaming checksum for driver/reference provenance."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_csv_values(text: str, cast, label: str) -> Tuple:
    """Parse a nonempty comma-separated option with contextual errors."""
    values = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        try:
            values.append(cast(token))
        except (TypeError, ValueError) as exc:
            raise argparse.ArgumentTypeError(
                "invalid %s value %r" % (label, token)) from exc
    if not values:
        raise argparse.ArgumentTypeError("%s list is empty" % label)
    return tuple(values)


def parse_epochs(text: str) -> Tuple[datetime, ...]:
    """Parse, sort, and de-duplicate requested snapshot epochs."""
    epochs = tuple(sorted(set(parse_csv_values(text, parse_utc, "epoch"))))
    if len(epochs) < 3:
        raise argparse.ArgumentTypeError(
            "C15 needs at least three epochs including a bracketed interpolation point")
    return epochs


def parse_schedulers(text: str) -> Tuple[str, ...]:
    """Parse cross-check scheduler names and reject unsupported spellings."""
    allowed = {"STATIC", "BLOCK_CYCLIC", "DYNAMIC"}
    schedulers = tuple(str(value).upper()
                       for value in parse_csv_values(text, str, "scheduler"))
    unknown = [value for value in schedulers if value not in allowed]
    if unknown:
        raise argparse.ArgumentTypeError("unsupported scheduler(s): %s" %
                                         ",".join(unknown))
    return tuple(dict.fromkeys(schedulers))


def parse_positive_floats(text: str) -> Tuple[float, ...]:
    """Parse a positive rigidity list."""
    values = tuple(float(value)
                   for value in parse_csv_values(text, float, "rigidity"))
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise argparse.ArgumentTypeError("rigidities must be finite and positive")
    return tuple(sorted(set(values)))


def kinetic_energy_mev_from_rigidity_gv(rigidity_gv: float) -> float:
    """Convert proton rigidity to kinetic energy for parser-required bounds."""
    momentum_mev = float(rigidity_gv) * 1000.0
    return math.hypot(momentum_mev, PROTON_REST_MEV) - PROTON_REST_MEV


def load_driver(path: Path, required_epochs: Sequence[datetime] = ()) -> Tuple[List[DriverRow], Dict[str, object]]:
    """Validate and load a standard timestamp-plus-19-value T05 driver.

    The selected interpolation window must have five-minute cadence, finite
    values, stable IMF/SW validity flags across interpolated intervals, and full
    coverage of every requested epoch.  Unknown/truncated input fails before any
    expensive field initialization begins.
    """
    rows: List[DriverRow] = []
    standard_header_seen = False
    with path.open(errors="replace") as stream:
        for line_number, raw in enumerate(stream, start=1):
            text = raw.strip()
            if not text:
                continue
            if text.startswith(("#", "!")):
                # AMPS discovers named columns from this comment.  Validate the
                # exact spelling here so a bad table fails immediately in the
                # lightweight runner instead of after a costly Mode3D launch.
                header_tokens = tuple(text.lstrip("#!").strip().split())
                if header_tokens and header_tokens[0] == AMPS_DRIVER_HEADER[0]:
                    if header_tokens != AMPS_DRIVER_HEADER:
                        raise ValueError(
                            "driver header must be exactly %r; got %r. "
                            "AMPS requires DST and does not accept SYM-H/SYM_H "
                            "for the T05 geomagnetic-index column" %
                            (" ".join(AMPS_DRIVER_HEADER),
                             " ".join(header_tokens)))
                    standard_header_seen = True
                continue
            fields = text.split()
            if len(fields) != 20:
                raise ValueError(
                    "driver line %d has %d fields; expected timestamp + 19 values" %
                    (line_number, len(fields)))
            epoch = parse_utc(fields[0])
            try:
                values = tuple(float(token) for token in fields[1:])
            except ValueError as exc:
                raise ValueError("driver line %d contains a nonnumeric value" %
                                 line_number) from exc
            if len(values) != len(DRIVER_FIELDS) or not all(
                    math.isfinite(value) for value in values):
                raise ValueError("driver line %d contains invalid values" % line_number)
            if rows and epoch <= rows[-1].epoch:
                raise ValueError("driver timestamps are not strictly increasing")
            rows.append(DriverRow(epoch, values))
    if len(rows) < 3:
        raise ValueError("driver must contain at least three records")
    if not standard_header_seen:
        raise ValueError(
            "driver is missing the canonical AMPS column header: %s" %
            " ".join(AMPS_DRIVER_HEADER))

    gaps = [(right.epoch - left.epoch).total_seconds()
            for left, right in zip(rows[:-1], rows[1:])]
    if not all(math.isclose(gap, 300.0, abs_tol=1.0e-9) for gap in gaps):
        raise ValueError("C15 driver must have uninterrupted five-minute cadence")
    if required_epochs and (min(required_epochs) < rows[0].epoch or
                            max(required_epochs) > rows[-1].epoch):
        raise ValueError("driver does not cover all selected C15 epochs")

    metadata = {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "record_count": len(rows),
        "first_epoch_utc": format_utc(rows[0].epoch),
        "last_epoch_utc": format_utc(rows[-1].epoch),
        "cadence_seconds": 300.0,
        "standard_header_seen": standard_header_seen,
    }
    return rows, metadata


def interpolate_driver(rows: Sequence[DriverRow], epoch: datetime) -> InterpolatedDriver:
    """Linearly interpolate all continuous driver values at one snapshot epoch.

    IMFflag and SWflag are categorical validity indicators.  Interpolation is
    allowed only when both bracket rows carry identical flag values; this avoids
    inventing a fractional flag in a custom driver.
    """
    for row in rows:
        if epoch == row.epoch:
            return InterpolatedDriver(epoch, row.values, row.epoch, row.epoch, 0.0)
    for left, right in zip(rows[:-1], rows[1:]):
        if left.epoch < epoch < right.epoch:
            width = (right.epoch - left.epoch).total_seconds()
            weight = (epoch - left.epoch).total_seconds() / width
            for flag_index in (9, 10):
                if left.values[flag_index] != right.values[flag_index]:
                    raise ValueError(
                        "cannot interpolate categorical %s across %s .. %s" %
                        (DRIVER_FIELDS[flag_index], format_utc(left.epoch),
                         format_utc(right.epoch)))
            values = tuple(
                a + weight * (b - a)
                for a, b in zip(left.values, right.values))
            return InterpolatedDriver(
                epoch, values, left.epoch, right.epoch, weight)
    raise ValueError("epoch %s is outside driver coverage" % format_utc(epoch))


def driver_reference_row(interpolated: InterpolatedDriver) -> Dict[str, object]:
    """Convert an interpolation result into the checked/output CSV schema."""
    row: Dict[str, object] = {
        "epoch_utc": format_utc(interpolated.epoch),
        "case_kind": ("exact" if interpolated.left_epoch == interpolated.right_epoch
                      else "midpoint" if math.isclose(interpolated.weight, 0.5)
                      else "interpolated"),
        "left_epoch_utc": format_utc(interpolated.left_epoch),
        "right_epoch_utc": format_utc(interpolated.right_epoch),
        "weight": interpolated.weight,
    }
    row.update(dict(zip(DRIVER_FIELDS, interpolated.values)))
    return row


def validate_checked_reference(path: Path, rows: Sequence[DriverRow]) -> None:
    """Ensure the committed default interpolation table matches the source driver."""
    with path.open(newline="") as stream:
        reference_rows = list(csv.DictReader(stream))
    if not reference_rows:
        raise ValueError("checked C15 interpolation reference is empty")
    for row in reference_rows:
        epoch = parse_utc(row["epoch_utc"])
        expected = driver_reference_row(interpolate_driver(rows, epoch))
        for key in ("case_kind", "left_epoch_utc", "right_epoch_utc"):
            if str(row[key]) != str(expected[key]):
                raise ValueError("reference mismatch at %s column %s" %
                                 (row["epoch_utc"], key))
        for key in ("weight",) + DRIVER_FIELDS:
            if not math.isclose(float(row[key]), float(expected[key]),
                                rel_tol=2.0e-14, abs_tol=2.0e-14):
                raise ValueError("reference mismatch at %s column %s" %
                                 (row["epoch_utc"], key))


def write_driver(path: Path, interpolated: InterpolatedDriver,
                 perturb: bool = False) -> Tuple[float, ...]:
    """Write a three-row exact-state driver used as an independent reference.

    The center row is exactly at the requested epoch.  Identical rows five
    minutes before and after satisfy the normal AMPS coverage/cadence contract
    while ensuring that the center record—not another interpolation—is selected.
    ``perturb=True`` changes physically relevant T05 inputs but leaves epoch,
    validity flags, and tilt intact for the same-epoch sensitivity guard.
    """
    values = list(interpolated.values)
    if perturb:
        values[1] += 3.0       # IMF By
        values[2] -= 5.0       # IMF Bz
        values[8] -= 50.0      # Dst geomagnetic disturbance index
        values[12] *= 1.5      # dynamic pressure
        for index in range(13, 19):
            values[index] = values[index] * 1.5 + 0.10
    header = (
        "# C15 independently materialized exact-state T05 driver.\n"
        "# " + " ".join(AMPS_DRIVER_HEADER) + "\n"
    )
    data_lines = []
    for delta in (-300, 0, 300):
        timestamp = interpolated.epoch + timedelta(seconds=delta)
        data_lines.append("%s %s" % (
            format_utc(timestamp),
            " ".join("%.17g" % value for value in values)))
    path.write_text(header + "\n".join(data_lines) + "\n")
    return tuple(values)


def replace_tokens(template_text: str, replacements: Mapping[str, str]) -> str:
    """Replace every named template token exactly once and reject leftovers."""
    output = template_text
    for name, value in replacements.items():
        token = "__%s__" % name
        if token not in output:
            raise ValueError("template token %s is missing" % token)
        output = output.replace(token, value)
    leftovers = sorted(set(re.findall(r"__[A-Z0-9_]+__", output)))
    if leftovers:
        raise ValueError("unresolved template token(s): %s" % ", ".join(leftovers))
    return output


def validate_template(path: Path) -> None:
    """Protect the parser-compatible C15 contract from speculative keywords."""
    text = path.read_text(errors="replace")
    required = {
        "vertical rigidity-list contract": r"^CUTOFF_SAMPLING\s+VERTICAL\s*$",
        "rigidity-list algorithm":
            r"^CUTOFF_SEARCH_ALGORITHM\s+RIGIDITY_LIST\s*$",
        "driver placeholder": r"^__DRIVER_DIRECTIVE__\s*$",
        "epoch placeholder": r"^EPOCH\s+__EPOCH__\s*$",
    }
    for label, pattern in required.items():
        if not re.search(pattern, text, re.MULTILINE):
            raise ValueError("C15 template lacks %s" % label)
    unsupported = (
        "TEMPORAL_MODE", "FIELD_UPDATE_DT", "TS_INPUT_MODE",
        "CUTOFF_UNRESOLVED_EXTENSION_PASSES", "CUTOFF_DEBUG_EXIT_TRACE",
    )
    for keyword in unsupported:
        if re.search(r"^\s*%s\b" % re.escape(keyword), text, re.MULTILINE):
            raise ValueError("C15 template contains unsupported keyword %s" % keyword)


def render_input(template: Path, destination: Path, case: RunCase,
                 driver_name: Optional[str], args: argparse.Namespace) -> None:
    """Render one complete input deck with fixed mesh/trajectory controls."""
    rigidities = tuple(sorted(args.rigidities_gv))
    replacements = {
        "RUN_ID": "C15_%s" % case.case_id,
        "CUTOFF_EMIN_MEV": "%.15g" % kinetic_energy_mev_from_rigidity_gv(
            0.8 * rigidities[0]),
        "CUTOFF_EMAX_MEV": "%.15g" % kinetic_energy_mev_from_rigidity_gv(
            1.2 * rigidities[-1]),
        "RIGIDITY_LIST_GV": ",".join("%.12g" % value for value in rigidities),
        "MAX_TRACE_TIME_S": "%.12g" % args.max_trace_time,
        "FIELD_MODEL": case.model,
        "EPOCH": format_utc(case.epoch),
        "DRIVER_DIRECTIVE": ("DRIVER_FILE                    %s" % driver_name
                             if driver_name else "! DIPOLE control: no DRIVER_FILE"),
        "MESH_RES_EARTH_RE": "%.12g" % args.mesh_res_earth_re,
        "MESH_RES_BOUNDARY_RE": "%.12g" % args.mesh_res_boundary_re,
        "DOMAIN_HALF_SIZE_KM": "%.12g" %
                               (args.domain_half_size_re * EARTH_RADIUS_KM),
        "SHELL_RES_DEG": "%.12g" % args.shell_res_deg,
        "DT_TRACE_S": "%.12g" % args.dt_trace,
        "MAX_STEPS": str(args.max_steps),
    }
    destination.write_text(replace_tokens(template.read_text(), replacements))


def build_command(case: RunCase, input_name: str, amps: Path,
                  args: argparse.Namespace) -> List[str]:
    """Build the exact one-line Mode3D invocation archived for each case."""
    return [
        args.mpirun, "-np", str(args.np), str(amps),
        "-mode", "3d", "-i", input_name,
        "-mode3d-field-eval", "MESH",
        "-mode3d-parallel", "THREADS",
        "-mode3d-threads", str(args.nt),
        "-mode3d-mpi-scheduler", case.scheduler,
        "-mode3d-mpi-dynamic-chunk", str(args.dynamic_chunk),
        "-mode3d-output-initialized",
        "-cutoff-search", "RIGIDITY_LIST",
        "-cutoff-rigidity-list-gv",
        ",".join("%.12g" % value for value in args.rigidities_gv),
        "-cutoff-access-abs-lat-min", "0",
        "-cutoff-access-abs-lat-max", "90",
        "-cutoff-trace-policy", "ACCURATE",
        "-mover", args.mover,
    ]


def run_process(command: Sequence[str], cwd: Path, log_path: Path) -> int:
    """Execute AMPS and flush every output line to both terminal and case log."""
    with log_path.open("w") as log:
        log.write("Command:\n  %s\n\n" % " ".join(command))
        log.flush()
        process = subprocess.Popen(
            list(command), cwd=str(cwd), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True)
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
            log.flush()
        return process.wait()


def numeric_records(path: Path) -> Iterator[Tuple[float, ...]]:
    """Yield finite all-numeric Tecplot records while ignoring headers."""
    with path.open(errors="replace") as stream:
        for raw in stream:
            text = raw.strip()
            if not text or text.startswith(("#", "!")):
                continue
            if text.upper().startswith((
                    "TITLE", "VARIABLES", "ZONE", "DATASETAUXDATA", "AUXDATA")):
                continue
            tokens = text.replace(",", " ").split()
            try:
                values = tuple(float(token) for token in tokens)
            except ValueError:
                continue
            if not values or not all(math.isfinite(value) for value in values):
                raise ValueError("non-finite numeric output in %s" % path)
            yield values


def numeric_values(paths: Sequence[Path]) -> Iterator[float]:
    """Flatten output files in stable filename/row/column order."""
    for path in sorted(paths, key=lambda item: item.name):
        for record in numeric_records(path):
            for value in record:
                yield value


def discover_field_files(run_dir: Path) -> List[Path]:
    """Find all Mode3D initialized-field data parts, excluding metadata files."""
    return sorted([
        path for path in run_dir.iterdir()
        if path.is_file() and path.name.startswith("amps_3d_initialized")
        and ".data.dat" in path.name
    ], key=lambda item: item.name)


def discover_cutoff_files(run_dir: Path) -> List[Path]:
    """Find the exact-rigidity shell product with compatible suffix handling."""
    primary = run_dir / "cutoff_3d_shells_access.dat"
    if primary.exists():
        return [primary]
    return sorted(run_dir.glob("cutoff_3d_shells_access*.dat*"),
                  key=lambda item: item.name)


def fingerprint_numeric_files(paths: Sequence[Path]) -> Fingerprint:
    """Canonicalize represented numeric values to exact binary64 hex strings."""
    if not paths:
        raise FileNotFoundError("no numeric output files found")
    digest = hashlib.sha256()
    rows = 0
    values_count = 0
    for path in sorted(paths, key=lambda item: item.name):
        digest.update(("FILE %s\n" % path.name).encode("utf-8"))
        for record in numeric_records(path):
            rows += 1
            values_count += len(record)
            digest.update(" ".join(value.hex() for value in record).encode("ascii"))
            digest.update(b"\n")
    if rows == 0:
        raise ValueError("no numeric records found")
    return Fingerprint(
        tuple(path.name for path in sorted(paths, key=lambda item: item.name)),
        rows, values_count, digest.hexdigest())


def compare_fingerprints(first: Optional[Fingerprint],
                         second: Optional[Fingerprint]) -> bool:
    """Require equal shapes and exact represented numeric values."""
    return bool(first is not None and second is not None and
                first.numeric_rows == second.numeric_rows and
                first.numeric_values == second.numeric_values and
                first.sha256 == second.sha256)


def compare_numeric_files(first_paths: Sequence[Path], second_paths: Sequence[Path],
                          rel_tol: float, abs_tol: float) -> NumericComparison:
    """Compare two complete output products with transparent residual metrics."""
    first = list(numeric_values(first_paths))
    second = list(numeric_values(second_paths))
    if len(first) != len(second):
        return NumericComparison(
            min(len(first), len(second)), abs(len(first) - len(second)),
            float("inf"), float("inf"), float("inf"), float("inf"), False,
            "numeric value counts differ: %d versus %d" %
            (len(first), len(second)))
    if not first:
        return NumericComparison(0, 1, float("inf"), float("inf"),
                                 float("inf"), float("inf"), False,
                                 "products contain no numeric values")
    mismatches = 0
    max_abs = 0.0
    max_rel = 0.0
    square_sum = 0.0
    changed_square = 0.0
    changed_scale_square = 0.0
    for left, right in zip(first, second):
        difference = abs(left - right)
        scale = max(abs(left), abs(right))
        relative = difference / scale if scale > 0.0 else difference
        max_abs = max(max_abs, difference)
        max_rel = max(max_rel, relative)
        square_sum += difference * difference
        if difference > 0.0:
            changed_square += difference * difference
            changed_scale_square += scale * scale
        if difference > abs_tol + rel_tol * scale:
            mismatches += 1
    changed_relative_rms = (
        math.sqrt(changed_square / changed_scale_square)
        if changed_scale_square > 0.0 else 0.0)
    return NumericComparison(
        len(first), mismatches, max_abs, max_rel,
        math.sqrt(square_sum / len(first)), changed_relative_rms,
        mismatches == 0)


def continuity_metrics(left_paths: Sequence[Path], middle_paths: Sequence[Path],
                       right_paths: Sequence[Path], weight: float) -> Dict[str, float]:
    """Measure field smoothness relative to linear interpolation in time.

    Coordinate and mesh-measure values are identical among cases and therefore
    cancel exactly.  Only time-varying field values contribute to the residuals.
    """
    left = list(numeric_values(left_paths))
    middle = list(numeric_values(middle_paths))
    right = list(numeric_values(right_paths))
    if not left or len(left) != len(middle) or len(left) != len(right):
        raise ValueError("continuity products have incompatible numeric shapes")
    endpoint_square = 0.0
    curvature_square = 0.0
    left_step_square = 0.0
    right_step_square = 0.0
    for a, b, c in zip(left, middle, right):
        endpoint_square += (c - a) ** 2
        expected = a + weight * (c - a)
        curvature_square += (b - expected) ** 2
        left_step_square += (b - a) ** 2
        right_step_square += (c - b) ** 2
    endpoint_rms = math.sqrt(endpoint_square / len(left))
    curvature_rms = math.sqrt(curvature_square / len(left))
    left_rate = math.sqrt(left_step_square / len(left)) / max(weight, 1.0e-30)
    right_rate = math.sqrt(right_step_square / len(left)) / max(1.0 - weight, 1.0e-30)
    smaller_rate = min(left_rate, right_rate)
    step_ratio = (max(left_rate, right_rate) / smaller_rate
                  if smaller_rate > 0.0 else
                  (1.0 if max(left_rate, right_rate) == 0.0 else float("inf")))
    return {
        "endpoint_rms_difference": endpoint_rms,
        "curvature_rms": curvature_rms,
        "normalized_curvature": (
            curvature_rms / endpoint_rms if endpoint_rms > 0.0 else float("inf")),
        "left_scaled_rate": left_rate,
        "right_scaled_rate": right_rate,
        "scaled_step_ratio": step_ratio,
    }


def case_paths(output_root: Path, epoch: datetime, kind: str) -> Path:
    """Centralize the stable, human-readable output layout."""
    if kind in ("full", "reference"):
        return output_root / "t05" / epoch_label(epoch) / kind
    return output_root / "checks" / kind


def build_cases(output_root: Path, epochs: Sequence[datetime], anchor: datetime,
                args: argparse.Namespace) -> List[RunCase]:
    """Build the complete matrix before touching any output directory."""
    cases: List[RunCase] = []
    for epoch in epochs:
        label = epoch_label(epoch)
        cases.append(RunCase(
            "t05_%s_full" % label, "driver_full", "T05", epoch, "full",
            args.scheduler, case_paths(output_root, epoch, "full"),
            paired_epoch=format_utc(epoch)))
        cases.append(RunCase(
            "t05_%s_reference" % label, "driver_reference", "T05", epoch,
            "reference", args.scheduler, case_paths(output_root, epoch, "reference"),
            paired_epoch=format_utc(epoch)))

    for repeat in range(2, args.repeats + 1):
        cases.append(RunCase(
            "repeat_%02d" % repeat, "repeat", "T05", anchor, "full",
            args.scheduler,
            case_paths(output_root, anchor, "repeat_%02d" % repeat),
            paired_epoch=format_utc(anchor)))
    for scheduler in args.cross_schedulers:
        if scheduler == args.scheduler:
            continue
        cases.append(RunCase(
            "scheduler_%s" % scheduler.lower(), "scheduler", "T05", anchor,
            "full", scheduler,
            case_paths(output_root, anchor, "scheduler_%s" % scheduler.lower()),
            paired_epoch=format_utc(anchor)))

    cases.append(RunCase(
        "dipole_first", "dipole_control", "DIPOLE", epochs[0], "none",
        args.scheduler, output_root / "dipole" / "first"))
    cases.append(RunCase(
        "dipole_last", "dipole_control", "DIPOLE", epochs[-1], "none",
        args.scheduler, output_root / "dipole" / "last"))
    cases.append(RunCase(
        "driver_perturbed", "driver_sensitivity", "T05", anchor, "perturbed",
        args.scheduler, case_paths(output_root, anchor, "driver_perturbed"),
        paired_epoch=format_utc(anchor)))
    return cases


def prepare_case(case: RunCase, template: Path, source_driver: Path,
                 driver_rows: Sequence[DriverRow], amps: Path,
                 args: argparse.Namespace) -> List[str]:
    """Create the case deck/driver and return its recorded command."""
    case.run_directory.mkdir(parents=True, exist_ok=True)
    driver_name: Optional[str] = None
    if case.model == "T05":
        driver_name = "ts05_driver_C15_case.txt"
        destination = case.run_directory / driver_name
        if case.driver_mode == "full":
            shutil.copy2(source_driver, destination)
        elif case.driver_mode in ("reference", "perturbed"):
            write_driver(
                destination, interpolate_driver(driver_rows, case.epoch),
                perturb=(case.driver_mode == "perturbed"))
        else:
            raise ValueError("unknown T05 driver mode %s" % case.driver_mode)
    input_path = case.run_directory / "AMPS_PARAM_C15.in"
    render_input(template, input_path, case, driver_name, args)
    command = build_command(case, input_path.name, amps, args)
    (case.run_directory / "C15_command.json").write_text(
        json.dumps(command, indent=2) + "\n")
    return command


def process_case(case: RunCase, command: List[str], args: argparse.Namespace) -> RunRecord:
    """Run or re-open one case and require both field and cutoff products."""
    return_code = 0
    if not args.skip_run and not args.dry_run:
        print("C15 %s command:\n  %s" % (case.case_id, " ".join(command)))
        return_code = run_process(
            command, case.run_directory, case.run_directory / "C15_amps.log")
    elif args.dry_run:
        print("C15 %s command:\n  %s" % (case.case_id, " ".join(command)))

    if args.dry_run:
        return RunRecord(
            case.case_id, case.category, case.model, format_utc(case.epoch),
            case.driver_mode, case.scheduler, str(case.run_directory), command,
            0, None, None, True)
    if return_code != 0:
        return RunRecord(
            case.case_id, case.category, case.model, format_utc(case.epoch),
            case.driver_mode, case.scheduler, str(case.run_directory), command,
            return_code, None, None, False,
            "AMPS exited with code %d" % return_code)

    errors = []
    try:
        field = fingerprint_numeric_files(discover_field_files(case.run_directory))
    except Exception as exc:
        field = None
        errors.append("initialized field: %s" % exc)
    try:
        cutoff = fingerprint_numeric_files(discover_cutoff_files(case.run_directory))
    except Exception as exc:
        cutoff = None
        errors.append("cutoff product: %s" % exc)
    return RunRecord(
        case.case_id, case.category, case.model, format_utc(case.epoch),
        case.driver_mode, case.scheduler, str(case.run_directory), command,
        return_code, field, cutoff, not errors, "; ".join(errors))


def record_to_dict(record: RunRecord) -> Dict[str, object]:
    """Convert nested fingerprint dataclasses into JSON-compatible dictionaries."""
    result = asdict(record)
    result["field_fingerprint"] = (
        asdict(record.field_fingerprint) if record.field_fingerprint else None)
    result["cutoff_fingerprint"] = (
        asdict(record.cutoff_fingerprint) if record.cutoff_fingerprint else None)
    return result


def comparison_row(label: str, first_case: RunCase, second_case: RunCase,
                   product: str, comparison: NumericComparison) -> Dict[str, object]:
    """Produce a stable diagnostic row for tolerant full-product comparisons."""
    row = {
        "comparison": label,
        "first_case": first_case.case_id,
        "second_case": second_case.case_id,
        "epoch_utc": format_utc(first_case.epoch),
        "product": product,
    }
    row.update(asdict(comparison))
    return row


def write_dict_csv(path: Path, rows: Sequence[Mapping[str, object]],
                   fieldnames: Optional[Sequence[str]] = None) -> None:
    """Write heterogeneous diagnostic rows with deterministic column ordering."""
    names = list(fieldnames or [])
    if not names:
        for row in rows:
            for name in row:
                if name not in names:
                    names.append(name)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def analyze_matrix(cases: Sequence[RunCase], records: Sequence[RunRecord],
                   driver_rows: Sequence[DriverRow], args: argparse.Namespace
                   ) -> Tuple[bool, Dict[str, List[Dict[str, object]]], List[str]]:
    """Evaluate every scientific gate and return complete diagnostics."""
    passed = all(record.passed for record in records) and bool(records)
    messages = ["%s: %s" % (record.case_id, record.error)
                for record in records if not record.passed]
    by_id = {case.case_id: case for case in cases}
    record_by_id = {record.case_id: record for record in records}
    full_by_epoch = {
        case.epoch: case for case in cases if case.category == "driver_full"}
    reference_by_epoch = {
        case.epoch: case for case in cases if case.category == "driver_reference"}

    equivalence_rows: List[Dict[str, object]] = []
    for epoch in sorted(full_by_epoch):
        first = full_by_epoch[epoch]
        second = reference_by_epoch[epoch]
        first_record = record_by_id[first.case_id]
        second_record = record_by_id[second.case_id]
        if not first_record.passed or not second_record.passed:
            passed = False
            messages.append("cannot compare driver equivalence at %s" % format_utc(epoch))
            continue
        for product, discover, rel_tol, abs_tol in (
                ("initialized_field", discover_field_files,
                 args.reference_rel_tol, args.reference_abs_tol),
                ("cutoff_access", discover_cutoff_files,
                 args.cutoff_reference_rel_tol, args.cutoff_reference_abs_tol)):
            comparison = compare_numeric_files(
                discover(first.run_directory), discover(second.run_directory),
                rel_tol, abs_tol)
            equivalence_rows.append(
                comparison_row("full_vs_materialized_driver", first, second,
                               product, comparison))
            if not comparison.passed:
                passed = False
                messages.append(
                    "%s driver equivalence failed for %s: %d mismatches, max_rel=%.3e" %
                    (format_utc(epoch), product, comparison.mismatches,
                     comparison.max_relative_difference))

    reproducibility_rows: List[Dict[str, object]] = []
    anchor = min(full_by_epoch, key=lambda value: abs(
        (value - args.anchor_epoch).total_seconds()))
    baseline_case = full_by_epoch[anchor]
    baseline_record = record_by_id[baseline_case.case_id]
    for case in cases:
        if case.category not in ("repeat", "scheduler"):
            continue
        record = record_by_id[case.case_id]
        field_equal = compare_fingerprints(
            baseline_record.field_fingerprint, record.field_fingerprint)
        cutoff_equal = compare_fingerprints(
            baseline_record.cutoff_fingerprint, record.cutoff_fingerprint)
        row = {
            "comparison": case.category,
            "baseline_case": baseline_case.case_id,
            "comparison_case": case.case_id,
            "epoch_utc": format_utc(anchor),
            "scheduler": case.scheduler,
            "field_exact": int(field_equal),
            "cutoff_exact": int(cutoff_equal),
            "passed": int(record.passed and field_equal and cutoff_equal),
        }
        reproducibility_rows.append(row)
        if not row["passed"]:
            passed = False
            messages.append("%s is not exactly reproducible" % case.case_id)

    dipole_cases = [case for case in cases if case.category == "dipole_control"]
    if len(dipole_cases) == 2:
        first_record = record_by_id[dipole_cases[0].case_id]
        second_record = record_by_id[dipole_cases[1].case_id]
        field_equal = compare_fingerprints(
            first_record.field_fingerprint, second_record.field_fingerprint)
        cutoff_equal = compare_fingerprints(
            first_record.cutoff_fingerprint, second_record.cutoff_fingerprint)
        row = {
            "comparison": "dipole_time_invariance",
            "baseline_case": dipole_cases[0].case_id,
            "comparison_case": dipole_cases[1].case_id,
            "epoch_utc": "%s..%s" %
                         (format_utc(dipole_cases[0].epoch),
                          format_utc(dipole_cases[1].epoch)),
            "scheduler": args.scheduler,
            "field_exact": int(field_equal),
            "cutoff_exact": int(cutoff_equal),
            "passed": int(first_record.passed and second_record.passed and
                          field_equal and cutoff_equal),
        }
        reproducibility_rows.append(row)
        if not row["passed"]:
            passed = False
            messages.append("centered DIPOLE changed with epoch")

    sensitivity_rows: List[Dict[str, object]] = []
    sensitivity_case = next(
        case for case in cases if case.category == "driver_sensitivity")
    sensitivity_record = record_by_id[sensitivity_case.case_id]
    if baseline_record.passed and sensitivity_record.passed:
        comparison = compare_numeric_files(
            discover_field_files(baseline_case.run_directory),
            discover_field_files(sensitivity_case.run_directory), 0.0, 0.0)
        sensitivity_passed = (
            comparison.mismatches > 0 and
            comparison.changed_value_relative_rms >= args.min_driver_sensitivity)
        sensitivity_rows.append({
            "baseline_case": baseline_case.case_id,
            "perturbed_case": sensitivity_case.case_id,
            "epoch_utc": format_utc(anchor),
            "changed_values": comparison.mismatches,
            "changed_value_relative_rms": comparison.changed_value_relative_rms,
            "minimum_required": args.min_driver_sensitivity,
            "passed": int(sensitivity_passed),
        })
        if not sensitivity_passed:
            passed = False
            messages.append(
                "perturbed DRIVER_FILE did not materially change the T05 field")
    else:
        passed = False
        messages.append("driver sensitivity comparison is unavailable")

    continuity_rows: List[Dict[str, object]] = []
    for epoch, middle_case in sorted(full_by_epoch.items()):
        interpolated = interpolate_driver(driver_rows, epoch)
        if interpolated.left_epoch == interpolated.right_epoch:
            continue
        left_case = full_by_epoch.get(interpolated.left_epoch)
        right_case = full_by_epoch.get(interpolated.right_epoch)
        if left_case is None or right_case is None:
            # Custom epoch lists may omit a bracket endpoint; equivalence still
            # validates interpolation, but continuity cannot be evaluated.
            passed = False
            messages.append(
                "continuity at %s lacks selected bracket endpoints" % format_utc(epoch))
            continue
        try:
            metrics = continuity_metrics(
                discover_field_files(left_case.run_directory),
                discover_field_files(middle_case.run_directory),
                discover_field_files(right_case.run_directory),
                interpolated.weight)
            row_passed = (
                metrics["endpoint_rms_difference"] > 0.0 and
                metrics["normalized_curvature"] <= args.max_normalized_curvature and
                metrics["scaled_step_ratio"] <= args.max_step_ratio)
        except Exception as exc:
            metrics = {
                "endpoint_rms_difference": float("nan"),
                "curvature_rms": float("nan"),
                "normalized_curvature": float("inf"),
                "left_scaled_rate": float("nan"),
                "right_scaled_rate": float("nan"),
                "scaled_step_ratio": float("inf"),
            }
            row_passed = False
            messages.append("continuity calculation failed at %s: %s" %
                            (format_utc(epoch), exc))
        row: Dict[str, object] = {
            "epoch_utc": format_utc(epoch),
            "left_epoch_utc": format_utc(interpolated.left_epoch),
            "right_epoch_utc": format_utc(interpolated.right_epoch),
            "weight": interpolated.weight,
            "max_normalized_curvature": args.max_normalized_curvature,
            "max_step_ratio": args.max_step_ratio,
            "passed": int(row_passed),
        }
        row.update(metrics)
        continuity_rows.append(row)
        if not row_passed:
            passed = False
            messages.append(
                "%s field continuity failed: curvature=%.3e step_ratio=%.3e" %
                (format_utc(epoch), metrics["normalized_curvature"],
                 metrics["scaled_step_ratio"]))

    return passed, {
        "equivalence": equivalence_rows,
        "reproducibility": reproducibility_rows,
        "sensitivity": sensitivity_rows,
        "continuity": continuity_rows,
    }, messages


def write_run_csv(path: Path, records: Sequence[RunRecord]) -> None:
    """Write compact fingerprints for every AMPS case."""
    fields = (
        "case_id", "category", "model", "epoch_utc", "driver_mode",
        "scheduler", "return_code", "field_numeric_rows", "field_sha256",
        "cutoff_numeric_rows", "cutoff_sha256", "passed", "error",
        "run_directory",
    )
    rows = []
    for record in records:
        rows.append({
            "case_id": record.case_id,
            "category": record.category,
            "model": record.model,
            "epoch_utc": record.epoch_utc,
            "driver_mode": record.driver_mode,
            "scheduler": record.scheduler,
            "return_code": record.return_code,
            "field_numeric_rows": (record.field_fingerprint.numeric_rows
                                   if record.field_fingerprint else ""),
            "field_sha256": (record.field_fingerprint.sha256
                             if record.field_fingerprint else ""),
            "cutoff_numeric_rows": (record.cutoff_fingerprint.numeric_rows
                                    if record.cutoff_fingerprint else ""),
            "cutoff_sha256": (record.cutoff_fingerprint.sha256
                              if record.cutoff_fingerprint else ""),
            "passed": record.passed,
            "error": record.error,
            "run_directory": record.run_directory,
        })
    write_dict_csv(path, rows, fields)


def write_summary(path: Path, passed: Optional[bool], args: argparse.Namespace,
                  records: Sequence[RunRecord], diagnostics: Mapping[str, Sequence[Mapping[str, object]]],
                  messages: Sequence[str]) -> None:
    """Write a concise human-readable result beside machine-readable products."""
    result_text = "DRY_RUN" if passed is None else ("PASS" if passed else "FAIL")
    lines = [
        "%s — %s" % (TEST_ID, TEST_NAME),
        "result: %s" % result_text,
        "profile: %s" % args.profile,
        "epochs: %s" % ",".join(format_utc(value) for value in args.epochs),
        "anchor epoch: %s" % format_utc(args.anchor_epoch),
        "MPI ranks / threads: %d / %d" % (args.np, args.nt),
        "primary scheduler: %s" % args.scheduler,
        "cross schedulers: %s" % ",".join(args.cross_schedulers),
        "AMPS cases: %d" % len(records),
        "driver equivalence rows: %d" % len(diagnostics.get("equivalence", [])),
        "reproducibility rows: %d" % len(diagnostics.get("reproducibility", [])),
        "continuity rows: %d" % len(diagnostics.get("continuity", [])),
        "",
    ]
    for record in records:
        lines.append("%-30s %-17s %s%s" %
                     (record.case_id, record.category,
                      "PASS" if record.passed else "FAIL",
                      (" — " + record.error) if record.error else ""))
    if messages:
        lines.extend(["", "Failure messages:"] + ["- " + item for item in messages])
    path.write_text("\n".join(lines) + "\n")


def validate_skip_configuration(existing: Mapping[str, object],
                                requested: Mapping[str, object]) -> None:
    """Prevent ``--skip-run`` from analyzing a different matrix by accident.

    Analysis tolerances may be changed deliberately during reprocessing, but the
    settings that determine which AMPS products exist must match the archived
    configuration exactly.
    """
    matrix_keys = (
        "epochs", "anchor_epoch", "driver_sha256", "np", "nt",
        "primary_scheduler", "cross_schedulers", "rigidities_GV",
        "repeats", "shell_res_deg", "mesh_res_earth_re",
        "mesh_res_boundary_re",
    )
    differences = []
    for key in matrix_keys:
        if existing.get(key) != requested.get(key):
            differences.append("%s: existing=%r requested=%r" %
                               (key, existing.get(key), requested.get(key)))
    if differences:
        raise ValueError(
            "--skip-run configuration does not match existing output:\n  " +
            "\n  ".join(differences))


def _write_synthetic_numeric(path: Path, values: Sequence[Sequence[float]]) -> None:
    """Write a small Tecplot-shaped product for runner/package self-tests."""
    lines = ['TITLE="C15 synthetic"', 'VARIABLES="x" "y" "z" "Bx"',
             'ZONE T="synthetic"']
    lines.extend(" ".join("%.17g" % value for value in row) for row in values)
    path.write_text("\n".join(lines) + "\n")


def self_test() -> int:
    """Exercise driver, reference, fingerprint, tolerance, and continuity guards."""
    try:
        validate_template(DEFAULT_TEMPLATE)
        driver_rows, _ = load_driver(DEFAULT_DRIVER)
        validate_checked_reference(DEFAULT_REFERENCE, driver_rows)
        midpoint = interpolate_driver(driver_rows, parse_utc("2012-05-17T05:57:30"))
        if (not math.isclose(midpoint.weight, 0.5) or
                not math.isclose(midpoint.values[0], -0.155, abs_tol=1.0e-15) or
                not math.isclose(midpoint.values[7], 10369.5, abs_tol=1.0e-12)):
            raise RuntimeError("midpoint interpolation is incorrect")
        exact = interpolate_driver(driver_rows, parse_utc("2012-05-17T06:00:00"))
        if exact.left_epoch != exact.right_epoch or exact.weight != 0.0:
            raise RuntimeError("exact-node selection is incorrect")

        with tempfile.TemporaryDirectory(prefix="c15_selftest_") as raw:
            root = Path(raw)
            materialized = root / "driver.txt"
            write_driver(materialized, midpoint)
            materialized_rows, _ = load_driver(materialized, [midpoint.epoch])
            selected = interpolate_driver(materialized_rows, midpoint.epoch)
            if any(not math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-15)
                   for a, b in zip(selected.values, midpoint.values)):
                raise RuntimeError("materialized exact driver changed values")

            # Reproduce the parser regression behind the original C15 failure:
            # SYM_H occupies the correct numeric position but is not a valid
            # named-column alias for the DST input required by T05.  The runner
            # must reject it before AMPS or MPI is launched.
            legacy = root / "driver_bad_symh.txt"
            legacy.write_text(materialized.read_text().replace(
                " Temp DST IMFflag ", " Temp SYM_H IMFflag "))
            try:
                load_driver(legacy, [midpoint.epoch])
            except ValueError as exc:
                if "requires DST" not in str(exc):
                    raise RuntimeError(
                        "legacy-header error did not explain the DST requirement")
            else:
                raise RuntimeError("unsupported SYM_H driver header was accepted")

            comparison_dirs = [root / name for name in ("a", "b", "c")]
            for directory in comparison_dirs:
                directory.mkdir()
            a, b, c = [directory / "field.dat" for directory in comparison_dirs]
            _write_synthetic_numeric(a, [(0, 0, 0, 1.0), (1, 0, 0, 2.0)])
            _write_synthetic_numeric(b, [(0, 0, 0, 1.0000000000000000),
                                         (1, 0, 0, 2.0)])
            _write_synthetic_numeric(c, [(0, 0, 0, 1.0001), (1, 0, 0, 2.0)])
            if not compare_fingerprints(
                    fingerprint_numeric_files([a]), fingerprint_numeric_files([b])):
                raise RuntimeError("equivalent numeric formatting changed fingerprint")
            if compare_fingerprints(
                    fingerprint_numeric_files([a]), fingerprint_numeric_files([c])):
                raise RuntimeError("changed numeric value was not detected")
            if not compare_numeric_files([a], [c], 2.0e-4, 0.0).passed:
                raise RuntimeError("tolerant numeric comparator rejected valid residual")
            if compare_numeric_files([a], [c], 1.0e-6, 0.0).passed:
                raise RuntimeError("tolerant numeric comparator accepted corruption")

            left = root / "left.dat"
            middle = root / "middle.dat"
            right = root / "right.dat"
            _write_synthetic_numeric(left, [(0, 0, 0, 1.0)])
            _write_synthetic_numeric(middle, [(0, 0, 0, 1.5)])
            _write_synthetic_numeric(right, [(0, 0, 0, 2.0)])
            metrics = continuity_metrics([left], [middle], [right], 0.5)
            if metrics["normalized_curvature"] != 0.0 or \
                    not math.isclose(metrics["scaled_step_ratio"], 1.0):
                raise RuntimeError("linear field failed continuity check")
            _write_synthetic_numeric(middle, [(0, 0, 0, 5.0)])
            bad = continuity_metrics([left], [middle], [right], 0.5)
            if bad["normalized_curvature"] <= 0.25:
                raise RuntimeError("discontinuous field passed continuity check")
    except Exception as exc:
        print("C15 self-test: FAIL: %s" % exc, file=sys.stderr)
        return 1
    print("C15 self-test: PASS")
    return 0


class ArgumentFormatter(argparse.ArgumentDefaultsHelpFormatter,
                        argparse.RawDescriptionHelpFormatter):
    """Preserve command examples while showing numeric defaults."""


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Define and validate the complete public runner interface."""
    parser = argparse.ArgumentParser(
        description=TEST_NAME,
        formatter_class=ArgumentFormatter,
        epilog=(
            "Recommended run from the AMPS repository root:\n"
            "  python3 srcEarth/test/C15/run_C15.py --profile ROUTINE "
            "--amps ./amps -np 4 -nt 16\n\n"
            "Package-only checks:\n"
            "  python3 srcEarth/test/C15/run_C15.py --self-test\n"
            "  python3 srcEarth/test/C15/run_C15.py --profile ROUTINE --dry-run"
        ))
    parser.add_argument("--profile", choices=tuple(PROFILE_DEFAULTS),
                        default="ROUTINE")
    parser.add_argument("--epochs", type=parse_epochs, default=None,
                        help="comma-separated snapshot epochs")
    parser.add_argument("--anchor-epoch", type=parse_utc,
                        default=parse_utc("2012-05-17T06:00:00"),
                        help="same-epoch repeat/scheduler/sensitivity epoch")
    parser.add_argument("--repeats", type=int, default=None,
                        help="total anchor executions including the primary run")
    parser.add_argument("--scheduler", choices=("STATIC", "BLOCK_CYCLIC", "DYNAMIC"),
                        default="DYNAMIC", help="primary cutoff scheduler")
    parser.add_argument("--cross-schedulers", type=parse_schedulers, default=None,
                        help="comma-separated scheduler equivalence cases")
    parser.add_argument("--dynamic-chunk", type=int, default=32)
    parser.add_argument("-np", "--np", type=int, default=4, help="MPI ranks")
    parser.add_argument("-nt", "--nt", type=int, default=16,
                        help="Mode3D cutoff worker threads")
    parser.add_argument("--amps", default="./amps")
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--mover", default="RK4")
    parser.add_argument("--template", default=str(DEFAULT_TEMPLATE))
    parser.add_argument("--driver", default=str(DEFAULT_DRIVER))
    parser.add_argument("--output-root", default="test_output/C15_t05_time")
    parser.add_argument("--rigidities-gv", type=parse_positive_floats,
                        default=(0.5, 1.0, 2.0, 5.0))
    parser.add_argument("--shell-res-deg", type=float, default=None)
    parser.add_argument("--domain-half-size-re", type=float, default=8.0)
    parser.add_argument("--mesh-res-earth-re", type=float, default=None)
    parser.add_argument("--mesh-res-boundary-re", type=float, default=None)
    parser.add_argument("--dt-trace", type=float, default=0.20)
    parser.add_argument("--max-steps", type=int, default=500000)
    parser.add_argument("--max-trace-time", type=float, default=60.0)
    parser.add_argument("--reference-rel-tol", type=float, default=5.0e-12)
    parser.add_argument("--reference-abs-tol", type=float, default=1.0e-15)
    parser.add_argument("--cutoff-reference-rel-tol", type=float, default=1.0e-8,
                        help="numeric tolerance for full/reference access products")
    parser.add_argument("--cutoff-reference-abs-tol", type=float, default=1.0e-10,
                        help="absolute tolerance for full/reference access products")
    parser.add_argument("--min-driver-sensitivity", type=float, default=1.0e-4)
    parser.add_argument("--max-normalized-curvature", type=float, default=0.25)
    parser.add_argument("--max-step-ratio", type=float, default=3.0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-run", action="store_true",
                        help="analyze the existing deterministic output matrix")
    parser.add_argument("--keep", action="store_true",
                        help="preserve existing output root before a new run")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--version", action="version",
                        version="C15 runner schema %d (%s)" %
                                (RUNNER_SCHEMA_VERSION, RUNNER_RELEASE))
    args = parser.parse_args(argv)

    defaults = PROFILE_DEFAULTS[args.profile]
    if args.epochs is None:
        args.epochs = tuple(parse_utc(value) for value in defaults["epochs"])
    if args.repeats is None:
        args.repeats = defaults["repeats"]
    if args.cross_schedulers is None:
        args.cross_schedulers = defaults["cross_schedulers"]
    if args.shell_res_deg is None:
        args.shell_res_deg = defaults["shell_res_deg"]
    if args.mesh_res_earth_re is None:
        args.mesh_res_earth_re = defaults["mesh_res_earth_re"]
    if args.mesh_res_boundary_re is None:
        args.mesh_res_boundary_re = defaults["mesh_res_boundary_re"]

    if args.repeats < 2:
        parser.error("--repeats must be >= 2")
    if args.np < 1 or args.nt < 1:
        parser.error("-np and -nt must be >= 1")
    if args.dynamic_chunk < 0:
        parser.error("--dynamic-chunk must be >= 0")
    if args.anchor_epoch not in args.epochs:
        parser.error("--anchor-epoch must be included in --epochs")
    positive = (
        args.shell_res_deg, args.domain_half_size_re, args.mesh_res_earth_re,
        args.mesh_res_boundary_re, args.dt_trace, args.max_trace_time,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in positive):
        parser.error("mesh, shell, domain, and time controls must be positive")
    if args.max_steps < 1:
        parser.error("--max-steps must be >= 1")
    for name in ("reference_rel_tol", "reference_abs_tol",
                 "cutoff_reference_rel_tol", "cutoff_reference_abs_tol",
                 "min_driver_sensitivity", "max_normalized_curvature",
                 "max_step_ratio"):
        value = getattr(args, name)
        if not math.isfinite(value) or value < 0.0:
            parser.error("--%s must be finite and nonnegative" %
                         name.replace("_", "-"))
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Prepare, execute, analyze, and archive the complete C15 matrix."""
    args = parse_args(argv)
    if args.self_test:
        return self_test()

    launch_dir = Path.cwd().resolve()
    template = Path(args.template).expanduser()
    if not template.is_absolute():
        template = (launch_dir / template).resolve() if template.exists() else DEFAULT_TEMPLATE
    driver = Path(args.driver).expanduser()
    if not driver.is_absolute():
        driver = (launch_dir / driver).resolve() if driver.exists() else DEFAULT_DRIVER
    output_root = Path(args.output_root).expanduser()
    if not output_root.is_absolute():
        output_root = (launch_dir / output_root).resolve()
    amps = Path(args.amps).expanduser()
    if not amps.is_absolute():
        amps = (launch_dir / amps).resolve()

    if not template.is_file() or not driver.is_file():
        print("C15 template or driver is missing", file=sys.stderr)
        return 2
    if not args.dry_run and not args.skip_run and not amps.is_file():
        print("AMPS executable not found: %s" % amps, file=sys.stderr)
        return 2
    try:
        validate_template(template)
        driver_rows, driver_metadata = load_driver(driver, args.epochs)
        validate_checked_reference(DEFAULT_REFERENCE, load_driver(DEFAULT_DRIVER)[0])
        for epoch in args.epochs:
            interpolate_driver(driver_rows, epoch)
    except Exception as exc:
        print("C15 input/reference validation failed: %s" % exc, file=sys.stderr)
        return 2

    if output_root.exists() and not args.keep and not args.skip_run:
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    configuration = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "runner_release": RUNNER_RELEASE,
        "profile": args.profile,
        "epochs": [format_utc(value) for value in args.epochs],
        "anchor_epoch": format_utc(args.anchor_epoch),
        "driver": driver_metadata,
        "driver_sha256": driver_metadata["sha256"],
        "source_reference_sha256": sha256_file(DEFAULT_REFERENCE),
        "np": args.np,
        "nt": args.nt,
        "primary_scheduler": args.scheduler,
        "cross_schedulers": list(args.cross_schedulers),
        "rigidities_GV": list(args.rigidities_gv),
        "repeats": args.repeats,
        "shell_res_deg": args.shell_res_deg,
        "mesh_res_earth_re": args.mesh_res_earth_re,
        "mesh_res_boundary_re": args.mesh_res_boundary_re,
        "reference_rel_tol": args.reference_rel_tol,
        "reference_abs_tol": args.reference_abs_tol,
        "cutoff_reference_rel_tol": args.cutoff_reference_rel_tol,
        "cutoff_reference_abs_tol": args.cutoff_reference_abs_tol,
        "min_driver_sensitivity": args.min_driver_sensitivity,
        "max_normalized_curvature": args.max_normalized_curvature,
        "max_step_ratio": args.max_step_ratio,
    }
    configuration_path = output_root / "C15_configuration.json"
    if args.skip_run:
        if not configuration_path.is_file():
            print("--skip-run requires existing %s" % configuration_path,
                  file=sys.stderr)
            return 2
        try:
            validate_skip_configuration(
                json.loads(configuration_path.read_text()), configuration)
        except Exception as exc:
            print("C15 reanalysis configuration error: %s" % exc,
                  file=sys.stderr)
            return 2
    else:
        configuration_path.write_text(json.dumps(configuration, indent=2) + "\n")

    cases = build_cases(output_root, args.epochs, args.anchor_epoch, args)
    commands: Dict[str, List[str]] = {}
    for case in cases:
        if not args.skip_run:
            command = prepare_case(case, template, driver, driver_rows, amps, args)
        else:
            command = build_command(case, "AMPS_PARAM_C15.in", amps, args)
        commands[case.case_id] = command

    interpolation_rows = [driver_reference_row(interpolate_driver(driver_rows, epoch))
                          for epoch in args.epochs]
    write_dict_csv(output_root / "C15_driver_reference_used.csv", interpolation_rows)
    shutil.copy2(SCRIPT_DIR / "reference_C15_acceptance_contract.csv",
                 output_root / "reference_C15_acceptance_contract.csv")

    records = [process_case(case, commands[case.case_id], args) for case in cases]
    (output_root / "C15_commands.json").write_text(json.dumps([
        {"case_id": case.case_id, "cwd": str(case.run_directory),
         "command": commands[case.case_id]}
        for case in cases
    ], indent=2) + "\n")
    write_run_csv(output_root / "C15_run_fingerprints.csv", records)

    if args.dry_run:
        diagnostics: Dict[str, List[Dict[str, object]]] = {
            "equivalence": [], "reproducibility": [],
            "sensitivity": [], "continuity": []}
        result_passed: Optional[bool] = None
        messages: List[str] = []
    else:
        result_passed, diagnostics, messages = analyze_matrix(
            cases, records, driver_rows, args)
        write_dict_csv(output_root / "C15_driver_equivalence.csv",
                       diagnostics["equivalence"])
        write_dict_csv(output_root / "C15_reproducibility.csv",
                       diagnostics["reproducibility"])
        write_dict_csv(output_root / "C15_driver_sensitivity.csv",
                       diagnostics["sensitivity"])
        write_dict_csv(output_root / "C15_continuity.csv",
                       diagnostics["continuity"])

    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "runner_release": RUNNER_RELEASE,
        "passed": result_passed,
        "dry_run": args.dry_run,
        "records": [record_to_dict(record) for record in records],
        "diagnostics": diagnostics,
        "messages": messages,
    }
    (output_root / "C15_result.json").write_text(json.dumps(result, indent=2) + "\n")
    write_summary(output_root / "C15_summary.txt", result_passed, args,
                  records, diagnostics, messages)
    print((output_root / "C15_summary.txt").read_text())

    if args.dry_run:
        return 0
    return 0 if result_passed else 1


if __name__ == "__main__":
    sys.exit(main())
