#!/usr/bin/env python3
"""C13 — GRIDLESS MPI scheduler correctness and load-distribution validation.

The test executes the same global-shell cutoff calculation through several
GRIDLESS MPI scheduling paths.  Parallel correctness is a hard requirement:
every case must conserve the AMPS task inventory and reproduce the one-rank
scientific product.  A centered dipole adds an independent analytical Størmer
reference; a driven T05 case exercises the production empirical-field path.

The runner intentionally separates three concepts that are often conflated:

* task closure (a hard correctness gate),
* task-count imbalance (diagnostic by default, optionally a local gate), and
* elapsed/per-rank time balance (reported diagnostically by default).

Task counts are not a proxy for trajectory cost, so C13 never claims that an
equal task count proves equal wall-clock load.  Use ``--max-rank-time-imbalance``
only when the tested AMPS executable emits explicit per-rank timing diagnostics.

Only Python's standard library is required.  ``--self-test`` validates the
package, parsers, comparison predicates, and deliberate negative cases without
an AMPS executable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import textwrap
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


TEST_ID = "C13"
TEST_NAME = "GRIDLESS MPI scheduler correctness and load distribution"
RUNNER_SCHEMA_VERSION = 3
RUNNER_RELEASE = "2026-08-31"

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_TEMPLATE = SCRIPT_DIR / "AMPS_PARAM_C13_gridless.in"
DEFAULT_DRIVER = SCRIPT_DIR / "data" / "ts05_driver_C13.txt"
DEFAULT_STORMER_REFERENCE = SCRIPT_DIR / "reference_C13_stormer_vertical.csv"
DEFAULT_ACCEPTANCE_REFERENCE = SCRIPT_DIR / "reference_C13_acceptance_contract.csv"

EARTH_RADIUS_KM = 6371.2
PROTON_REST_MEV = 938.27208816

# The dipole coefficient is calculated from the same physical constants used by
# C1/C12.  Keeping the expression here (rather than a rounded 14.9-GV literal)
# makes the executable reference and checked-in table independently auditable.
STORMER_R0_GV = (
    0.299792458 * 0.25 * 3.12e-5 * (EARTH_RADIUS_KM * 1000.0)
)
DEFAULT_STORMER_MIN_ABS_LAT_DEG = 30.0
DEFAULT_STORMER_MAX_ABS_LAT_DEG = 45.0
DEFAULT_STORMER_FLOOR_GV = 0.5

ALLOWED_MODELS = ("DIPOLE", "T05")
ALLOWED_SCHEDULERS = ("DYNAMIC", "BLOCK_CYCLIC", "STATIC")

# Profiles change computational cost, never the meaning of a gate.  ROUTINE is
# deliberately smaller than the full 10-degree campaign so it is suitable for
# regular validation; THOROUGH implements the complete planned rank envelope.
PROFILE_DEFAULTS = {
    "SMOKE": {
        "models": "DIPOLE",
        "ranks": "1,4",
        "shell_lon_res_deg": 60.0,
        "shell_lat_res_deg": 30.0,
        "cutoff_scan_n": 32,
        "max_trace_time": 60.0,
        "dipole_max_trace_time": 300.0,
        "chunks": "1,2",
    },
    "ROUTINE": {
        "models": "DIPOLE,T05",
        "ranks": "1,4,8",
        "shell_lon_res_deg": 30.0,
        "shell_lat_res_deg": 15.0,
        "cutoff_scan_n": 80,
        "max_trace_time": 300.0,
        "dipole_max_trace_time": 600.0,
        "chunks": "1,4",
    },
    "THOROUGH": {
        "models": "DIPOLE,T05",
        "ranks": "1,4,8,16",
        "shell_lon_res_deg": 10.0,
        "shell_lat_res_deg": 10.0,
        "cutoff_scan_n": 160,
        "max_trace_time": 600.0,
        "dipole_max_trace_time": 1200.0,
        "chunks": "1,8",
    },
}


@dataclass(frozen=True)
class CutoffRow:
    """One coordinate-keyed scientific record from PENUMBRA_SCAN output."""

    lon_deg: float
    lat_deg: float
    rc_lower_gv: float
    rc_effective_gv: float
    rc_upper_gv: float
    n_unresolved: int
    lower_bracket_unresolved: int
    upper_bracket_unresolved: int
    lower_below_range: int = 0
    lower_above_range: int = 0
    upper_below_range: int = 0
    upper_above_range: int = 0

    @property
    def key(self) -> Tuple[float, float]:
        """Return a stable grid key while treating 0 and 360 degrees equally."""
        lon = self.lon_deg % 360.0
        if math.isclose(lon, 360.0, abs_tol=1.0e-8) or math.isclose(
                lon, 0.0, abs_tol=1.0e-8):
            lon = 0.0
        lat = 0.0 if math.isclose(self.lat_deg, 0.0, abs_tol=1.0e-8) else self.lat_deg
        return round(lon, 7), round(lat, 7)


@dataclass(frozen=True)
class TaskDiagnostics:
    """AMPS task-distribution block parsed from one combined execution log."""

    expected_tasks: int
    reported_sum: int
    rank_tasks: Tuple[Tuple[int, int], ...]
    reported_min: int
    reported_average: float
    reported_max: int
    computed_sum: int
    computed_min: int
    computed_average: float
    computed_max: int
    task_imbalance: float
    closure_passed: bool
    error: str = ""


@dataclass(frozen=True)
class ProductComparison:
    """Result of comparing a parallel output with its one-rank baseline."""

    baseline_rows: int
    candidate_rows: int
    missing_keys: int
    extra_keys: int
    numeric_values_compared: int
    numeric_mismatches: int
    flag_mismatches: int
    max_absolute_difference: float
    max_relative_difference: float
    exact_fingerprint_match: bool
    passed: bool


@dataclass(frozen=True)
class StormerComparison:
    """Independent DIPOLE comparison against the vertical Størmer solution."""

    output_rows: int
    required_rows: int
    valid_rows: int
    passing_rows: int
    valid_fraction: float
    passing_fraction: float
    rmse_gv: float
    max_abs_error_gv: float
    tolerance_gv: float
    passed: bool


@dataclass
class RunCase:
    """Complete definition and later audit state for one AMPS invocation."""

    case_id: str
    model: str
    scheduler: str
    ranks: int
    dynamic_chunk: int
    workdir: Path
    input_file: Optional[Path] = None
    log_file: Optional[Path] = None
    command: Optional[List[str]] = None
    return_code: Optional[int] = None
    elapsed_seconds: Optional[float] = None
    output_file: Optional[Path] = None
    output_sha256: str = ""
    scientific_fingerprint: str = ""
    task_diagnostics: Optional[TaskDiagnostics] = None
    product_comparison: Optional[ProductComparison] = None
    stormer_comparison: Optional[StormerComparison] = None
    per_rank_times_seconds: Optional[Tuple[Tuple[int, float], ...]] = None
    rank_time_imbalance: Optional[float] = None
    speedup_vs_one_rank: Optional[float] = None
    parallel_efficiency: Optional[float] = None
    elapsed_ratio_to_best_same_rank: Optional[float] = None
    passed: bool = False
    error: str = ""


def parse_csv_names(text: str, label: str, allowed: Iterable[str]) -> List[str]:
    """Parse a case-insensitive comma list, preserving order and uniqueness."""
    allowed_set = {name.upper() for name in allowed}
    result: List[str] = []
    for token in str(text).replace(";", ",").split(","):
        name = token.strip().upper()
        if not name:
            continue
        if name not in allowed_set:
            raise SystemExit(
                "Unsupported %s %r; choose from %s" %
                (label, name, ",".join(sorted(allowed_set))))
        if name not in result:
            result.append(name)
    if not result:
        raise SystemExit("No values supplied for %s" % label)
    return result


def parse_csv_ints(text: str, label: str, minimum: int) -> List[int]:
    """Parse a positive/nonnegative integer list with deterministic ordering."""
    result: List[int] = []
    for token in str(text).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = int(token)
        except ValueError:
            raise SystemExit("Invalid %s value %r" % (label, token))
        if value < minimum:
            raise SystemExit("%s values must be >= %d" % (label, minimum))
        if value not in result:
            result.append(value)
    if not result:
        raise SystemExit("No values supplied for %s" % label)
    return sorted(result)


def fmt(value: float) -> str:
    """Write input numbers compactly without profile-dependent formatting."""
    return "%.15g" % float(value)


def rigidity_to_kinetic_mev(rigidity_gv: float) -> float:
    """Convert proton rigidity in GV to kinetic energy in MeV/nucleon."""
    momentum_mev_c = float(rigidity_gv) * 1000.0
    return math.sqrt(momentum_mev_c ** 2 + PROTON_REST_MEV ** 2) - PROTON_REST_MEV


def stormer_vertical_gv(lat_deg: float, altitude_km: float) -> float:
    """Analytical centered-dipole vertical cutoff at magnetic latitude."""
    radius_re = (EARTH_RADIUS_KM + float(altitude_km)) / EARTH_RADIUS_KM
    return STORMER_R0_GV * math.cos(math.radians(lat_deg)) ** 4 / radius_re ** 2


def expected_shell_task_count(args: argparse.Namespace) -> int:
    """Return the independent number of spatial shell tasks in the input deck.

    PENUMBRA_SCAN is scheduled as one GRIDLESS work item per shell coordinate;
    the rigidity scan is performed inside that work item.  Checking AMPS's
    reported ``totalTasks`` against this geometry-derived count prevents an
    internally self-consistent but incomplete task inventory from passing.
    """
    longitude_count = int(round(360.0 / args.shell_lon_res_deg))
    latitude_count = int(round(180.0 / args.shell_lat_res_deg)) + 1
    return longitude_count * latitude_count


def field_block(model: str, driver_name: str) -> str:
    """Return a complete parser-safe background-field block for one model."""
    if model == "DIPOLE":
        return "\n".join([
            "FIELD_MODEL                    DIPOLE",
            "DIPOLE_MOMENT                  1.0",
            "DIPOLE_TILT                    0.0",
        ])
    if model == "T05":
        # Static values are parser-required fallbacks.  The bundled driver row
        # at EPOCH is authoritative and exercises the production table reader.
        return "\n".join([
            "FIELD_MODEL                    T05",
            "EPOCH                          2012-05-17T06:00:00",
            "DRIVER_FILE                    %s" % driver_name,
            "DST                            -22.0",
            "PDYN                           1.03",
            "IMF_BX                         -0.34",
            "IMF_BY                         -8.71",
            "IMF_BZ                         1.52",
            "SW_VX                          -359.6",
            "SW_N                           4.09",
            "G1                             0.0",
            "G2                             0.0",
            "G3                             0.0",
            "T05_W1                         0.21",
            "T05_W2                         0.08",
            "T05_W3                         0.26",
            "T05_W4                         0.09",
            "T05_W5                         0.02",
            "T05_W6                         0.03",
        ])
    raise ValueError("Unsupported field model %s" % model)


def render_input(template: Path, destination: Path, case: RunCase,
                 driver_name: str, args: argparse.Namespace) -> None:
    """Render one immutable case input and reject every unresolved token."""
    text = template.read_text(errors="replace")
    replacements = {
        "__RUN_ID__": case.case_id,
        "__CUTOFF_EMIN_MEV__": fmt(rigidity_to_kinetic_mev(args.rigidity_min_gv)),
        "__CUTOFF_EMAX_MEV__": fmt(rigidity_to_kinetic_mev(args.rigidity_max_gv)),
        "__CUTOFF_SCAN_N__": str(args.cutoff_scan_n),
        "__MAX_TRACE_TIME__": fmt(
            args.dipole_max_trace_time if case.model == "DIPOLE"
            else args.max_trace_time),
        "__FIELD_BLOCK__": field_block(case.model, driver_name),
        "__SHELL_ALT_KM__": fmt(args.shell_alt_km),
        "__SHELL_LON_RES_DEG__": fmt(args.shell_lon_res_deg),
        "__SHELL_LAT_RES_DEG__": fmt(args.shell_lat_res_deg),
        "__DT_TRACE__": fmt(args.dt_trace),
        "__MAX_STEPS__": str(args.max_steps),
        "__MAX_TRACE_DISTANCE_RE__": fmt(args.max_trace_distance_re),
        "__GRIDLESS_MPI_SCHEDULER__": case.scheduler,
        "__GRIDLESS_MPI_DYNAMIC_CHUNK__": str(case.dynamic_chunk),
        "__GRIDLESS_THREADS__": str(args.nt),
    }
    for placeholder, replacement in replacements.items():
        text = text.replace(placeholder, replacement)
    remaining = sorted(set(re.findall(r"__[A-Z0-9_]+__", text)))
    if remaining:
        raise RuntimeError("Unresolved C13 template tokens: %s" % ", ".join(remaining))
    destination.write_text(text)


def validate_template(path: Path) -> None:
    """Protect C13 from reintroducing proposed/non-parser input keywords."""
    text = path.read_text(errors="replace")
    required = (
        "CUTOFF_SAMPLING", "CUTOFF_SEARCH_ALGORITHM", "PENUMBRA_SCAN",
        "GRIDLESS_MPI_SCHEDULER", "GRIDLESS_MPI_DYNAMIC_CHUNK",
        "GRIDLESS_PARALLEL", "GRIDLESS_THREADS", "__FIELD_BLOCK__",
    )
    missing = [token for token in required if token not in text]
    if missing:
        raise RuntimeError("C13 template lacks: %s" % ", ".join(missing))
    banned = (
        "CUTOFF_UNRESOLVED_EXTENSION_PASSES",
        "CUTOFF_UNRESOLVED_TIME_MULTIPLIER",
        "TS05_DST",
        "MESH_REFINEMENT_LEVELS",
    )
    active = "\n".join(
        line for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith(("!", "#")))
    present = [token for token in banned if re.search(r"\b%s\b" % token, active)]
    if present:
        raise RuntimeError("C13 template uses unsupported/proposed keywords: %s" %
                           ", ".join(present))


def validate_driver(path: Path) -> Dict[str, object]:
    """Validate the bundled T05 table, including the parser-critical DST name."""
    header: Optional[List[str]] = None
    rows: List[Tuple[str, List[float]]] = []
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            tokens = line[1:].strip().split()
            if tokens and tokens[0].upper().startswith("YYYY-MM-DD"):
                header = tokens
            continue
        parts = line.split()
        if len(parts) != 20:
            raise RuntimeError("C13 T05 driver row must contain epoch + 19 values: %s" % line)
        try:
            values = [float(value) for value in parts[1:]]
        except ValueError as exc:
            raise RuntimeError("Non-numeric C13 T05 driver row: %s" % exc)
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError("Non-finite C13 T05 driver value")
        rows.append((parts[0], values))
    if header is None:
        raise RuntimeError("C13 T05 driver lacks a named column header")
    if "DST" not in header or "SYM-H" in header or "SYM_H" in header:
        raise RuntimeError("C13 T05 driver must use the AMPS parser column name DST")
    if not rows or "2012-05-17T06:00:00" not in {row[0] for row in rows}:
        raise RuntimeError("C13 T05 driver does not cover its validation epoch")
    return {
        "path": str(path),
        "records": len(rows),
        "first_epoch": rows[0][0],
        "last_epoch": rows[-1][0],
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def load_stormer_reference(path: Path) -> List[Dict[str, str]]:
    """Load and independently verify every checked-in analytical reference row."""
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"altitude_km", "latitude_deg", "rc_stormer_gv", "required_for_pass"}
    if not rows or not required.issubset(rows[0]):
        raise RuntimeError("Malformed C13 Størmer reference table")
    seen = set()
    for row in rows:
        altitude = float(row["altitude_km"])
        latitude = float(row["latitude_deg"])
        stored = float(row["rc_stormer_gv"])
        computed = stormer_vertical_gv(latitude, altitude)
        if not math.isclose(stored, computed, rel_tol=5.0e-12, abs_tol=1.0e-14):
            raise RuntimeError(
                "C13 Størmer reference drift at lat=%g: stored=%.15g computed=%.15g" %
                (latitude, stored, computed))
        key = (altitude, latitude)
        if key in seen:
            raise RuntimeError("Duplicate C13 Størmer reference row %r" % (key,))
        expected_required = (
            DEFAULT_STORMER_MIN_ABS_LAT_DEG <= abs(latitude) <=
            DEFAULT_STORMER_MAX_ABS_LAT_DEG and
            computed >= DEFAULT_STORMER_FLOOR_GV)
        stored_required = row["required_for_pass"].strip().upper() in (
            "T", "TRUE", "1", "YES")
        if stored_required != expected_required:
            raise RuntimeError(
                "C13 Størmer required_for_pass drift at lat=%g" % latitude)
        seen.add(key)
    return rows


def validate_acceptance_reference(path: Path) -> List[Dict[str, str]]:
    """Ensure the human-readable gate inventory is complete and unambiguous."""
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    expected = {
        "execution_success", "task_closure", "rank_coverage",
        "task_inventory_geometry", "dynamic_task_count_imbalance",
        "dipole_stormer_reference",
        "dipole_scheduler_equivalence", "t05_scheduler_equivalence",
        "rank_timing_balance",
    }
    got = {row.get("gate_id", "") for row in rows}
    if got != expected:
        raise RuntimeError("C13 acceptance table IDs differ: missing=%s extra=%s" %
                           (sorted(expected - got), sorted(got - expected)))
    return rows


def make_cases(work_root: Path, models: Sequence[str], ranks: Sequence[int],
               schedulers: Sequence[str], chunks: Sequence[int]) -> List[RunCase]:
    """Build the minimal matrix covering rank, scheduler, and chunk dimensions.

    DYNAMIC/auto (chunk=0) runs at every requested rank.  Alternative
    schedulers and explicit DYNAMIC chunks run at the largest rank.  Thus each
    added case changes exactly one scheduling dimension relative to an already
    represented branch, while every model retains a one-rank baseline.
    """
    max_rank = max(ranks)
    identities = set()
    cases: List[RunCase] = []

    def add(model: str, scheduler: str, rank_count: int, chunk: int) -> None:
        identity = (model, scheduler, rank_count, chunk)
        if identity in identities:
            return
        identities.add(identity)
        case_id = "%s_%s_np%02d_chunk%03d" % (
            model.lower(), scheduler.lower(), rank_count, chunk)
        cases.append(RunCase(
            case_id=case_id,
            model=model,
            scheduler=scheduler,
            ranks=rank_count,
            dynamic_chunk=chunk,
            workdir=work_root / model.lower() / case_id,
        ))

    for model in models:
        for rank_count in ranks:
            add(model, "DYNAMIC", rank_count, 0)
        for scheduler in schedulers:
            if scheduler != "DYNAMIC":
                add(model, scheduler, max_rank, 0)
        for chunk in chunks:
            if chunk > 0:
                add(model, "DYNAMIC", max_rank, chunk)
    return cases


def compute_performance_metrics(cases: Sequence[RunCase],
                                models: Sequence[str]) -> None:
    """Attach diagnostic whole-job speedup and efficiency to completed cases.

    The one-rank DYNAMIC/auto job is the speedup baseline.  The best elapsed
    time at the same model/rank count provides a second, scheduler-local ratio.
    No threshold is applied here because portable performance gates require a
    controlled host, fixed process placement, and repeated timing samples.
    """
    for model in models:
        model_cases = [case for case in cases if case.model == model]
        baseline = next((case for case in model_cases if
                         case.scheduler == "DYNAMIC" and case.ranks == 1 and
                         case.dynamic_chunk == 0), None)
        baseline_time = baseline.elapsed_seconds if baseline else None
        best_by_rank = {
            rank_count: min(
                case.elapsed_seconds for case in model_cases
                if case.ranks == rank_count and case.elapsed_seconds is not None)
            for rank_count in sorted({case.ranks for case in model_cases})
            if any(case.elapsed_seconds is not None for case in model_cases
                   if case.ranks == rank_count)
        }
        for case in model_cases:
            if (baseline_time is not None and baseline_time > 0.0 and
                    case.elapsed_seconds is not None and case.elapsed_seconds > 0.0):
                case.speedup_vs_one_rank = baseline_time / case.elapsed_seconds
                case.parallel_efficiency = case.speedup_vs_one_rank / case.ranks
                case.elapsed_ratio_to_best_same_rank = (
                    case.elapsed_seconds / best_by_rank[case.ranks])


def command_for(case: RunCase, amps: Path, input_name: str,
                args: argparse.Namespace) -> List[str]:
    """Create the production one-line command recorded for a case."""
    command = [args.mpirun]
    for launcher_arg in args.mpirun_arg:
        command.append(launcher_arg)
    command.extend([
        "-np", str(case.ranks), str(amps),
        "-mode", "gridless", "-i", input_name,
        "-mover", args.mover,
        "-gridless-mpi-scheduler", case.scheduler,
        "-gridless-mpi-dynamic-chunk", str(case.dynamic_chunk),
        "-gridless-parallel", "THREADS",
        "-gridless-threads", str(args.nt),
        "-cutoff-search", "PENUMBRA_SCAN",
        "-cutoff-trace-policy", "ACCURATE",
    ])
    return command


def run_process(command: Sequence[str], workdir: Path,
                log_path: Path) -> Tuple[int, float]:
    """Run AMPS, mirror output for users, and measure whole-job elapsed time."""
    start = time.perf_counter()
    with log_path.open("w") as log:
        log.write("Command:\n  %s\n\n" % " ".join(command))
        log.flush()
        process = subprocess.Popen(
            list(command), cwd=str(workdir), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True)
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            log.write(line)
        return_code = process.wait()
    return return_code, time.perf_counter() - start


def parse_task_diagnostics(text: str, expected_ranks: int) -> TaskDiagnostics:
    """Parse and cross-check the last AMPS GRIDLESS task-distribution block."""
    starts = [match.start() for match in re.finditer(
        r"\[gridless\]\[MPI\]\s*Task distribution check:", text, re.IGNORECASE)]
    if not starts:
        raise RuntimeError("AMPS log lacks '[gridless][MPI] Task distribution check'")
    block = text[starts[-1]:]
    expected_match = re.search(r"totalTasks\s*\(expected\)\s*=\s*(\d+)", block)
    sum_match = re.search(r"sum\(all rank tasks\)\s*=\s*(\d+)", block)
    range_match = re.search(
        r"per-rank min/avg/max\s*=\s*(\d+)\s*/\s*([0-9.eE+\-]+)\s*/\s*(\d+)",
        block)
    if not (expected_match and sum_match and range_match):
        raise RuntimeError("Incomplete AMPS task-distribution summary")
    rank_pairs = [
        (int(rank), int(tasks))
        for rank, tasks in re.findall(r"^\s*rank\s+(\d+)\s*:\s*(\d+)\s+tasks\s*$",
                                      block, re.IGNORECASE | re.MULTILINE)
    ]
    ranks_seen = [rank for rank, _ in rank_pairs]
    if len(rank_pairs) != expected_ranks or sorted(ranks_seen) != list(range(expected_ranks)):
        raise RuntimeError(
            "Task block rank coverage is %s; expected ranks 0..%d" %
            (sorted(ranks_seen), expected_ranks - 1))
    if len(set(ranks_seen)) != len(ranks_seen):
        raise RuntimeError("Task block contains duplicate rank lines")

    expected_tasks = int(expected_match.group(1))
    reported_sum = int(sum_match.group(1))
    reported_min = int(range_match.group(1))
    reported_average = float(range_match.group(2))
    reported_max = int(range_match.group(3))
    values = [tasks for _, tasks in rank_pairs]
    computed_sum = sum(values)
    computed_min = min(values)
    computed_max = max(values)
    computed_average = float(computed_sum) / expected_ranks
    imbalance = (float(computed_max) / computed_average
                 if computed_average > 0.0 else math.inf)
    average_tolerance = max(0.051, 5.0e-4 * max(1.0, computed_average))
    closure = (
        expected_tasks > 0 and
        expected_tasks == reported_sum == computed_sum and
        reported_min == computed_min and reported_max == computed_max and
        math.isclose(reported_average, computed_average,
                     rel_tol=0.0, abs_tol=average_tolerance)
    )
    error = "" if closure else (
        "task closure mismatch: expected=%d reported_sum=%d computed_sum=%d "
        "reported_min/avg/max=%d/%g/%d computed=%d/%g/%d" %
        (expected_tasks, reported_sum, computed_sum, reported_min,
         reported_average, reported_max, computed_min, computed_average,
         computed_max))
    return TaskDiagnostics(
        expected_tasks=expected_tasks,
        reported_sum=reported_sum,
        rank_tasks=tuple(sorted(rank_pairs)),
        reported_min=reported_min,
        reported_average=reported_average,
        reported_max=reported_max,
        computed_sum=computed_sum,
        computed_min=computed_min,
        computed_average=computed_average,
        computed_max=computed_max,
        task_imbalance=imbalance,
        closure_passed=closure,
        error=error,
    )


def parse_per_rank_times(text: str, expected_ranks: int) -> Optional[Tuple[Tuple[int, float], ...]]:
    """Read optional explicit rank wall-time diagnostics without guessing.

    AMPS versions that do not emit a line containing both ``rank`` and one of
    ``wall time``, ``elapsed time``, or ``runtime`` simply return ``None``.  The
    parser deliberately does not reinterpret task counts as timing evidence.
    """
    patterns = (
        r"rank\s+(\d+).*?wall\s*time\s*[=:]\s*([0-9.eE+\-]+)\s*s",
        r"rank\s+(\d+).*?elapsed\s*time\s*[=:]\s*([0-9.eE+\-]+)\s*s",
        r"rank\s+(\d+).*?runtime\s*[=:]\s*([0-9.eE+\-]+)\s*s",
    )
    values: Dict[int, float] = {}
    for pattern in patterns:
        for rank_text, value_text in re.findall(pattern, text, re.IGNORECASE):
            value = float(value_text)
            if math.isfinite(value) and value >= 0.0:
                values[int(rank_text)] = value
    if not values:
        return None
    if sorted(values) != list(range(expected_ranks)):
        raise RuntimeError("Partial per-rank timing diagnostics: ranks=%s expected=0..%d" %
                           (sorted(values), expected_ranks - 1))
    return tuple(sorted(values.items()))


def normalize_variable_name(name: str) -> str:
    """Normalize Tecplot labels exactly once for strict named-column lookup."""
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def parse_penumbra_output(path: Path) -> List[CutoffRow]:
    """Parse PENUMBRA_SCAN output using column names, never fixed offsets."""
    variables: List[str] = []
    reading_variables = False
    numeric_rows: List[List[float]] = []
    for line_number, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        upper = line.upper()
        if upper.startswith("VARIABLES"):
            reading_variables = True
            variables.extend(normalize_variable_name(name)
                             for name in re.findall(r'"([^"]+)"', line))
            continue
        if reading_variables:
            # A quoted ZONE title is metadata, not another variable name.  Test
            # this boundary before collecting continuation-line quotations.
            if upper.startswith(("TITLE", "ZONE")):
                reading_variables = False
            else:
                quoted = re.findall(r'"([^"]+)"', line)
                if quoted:
                    variables.extend(normalize_variable_name(name) for name in quoted)
                    continue
                reading_variables = False
        if upper.startswith(("TITLE", "ZONE", "#", "!")):
            continue
        try:
            values = [float(token) for token in line.replace(",", " ").split()]
        except ValueError:
            continue
        if variables and len(values) == len(variables):
            numeric_rows.append(values)
        elif variables and values:
            raise RuntimeError(
                "%s numeric row %d has %d values for %d VARIABLES" %
                (path, line_number, len(values), len(variables)))

    required = {
        "lon_deg", "lat_deg", "rc_lower_gv", "rc_effective_gv", "rc_upper_gv",
        "n_unresolved", "lower_bracket_unresolved", "upper_bracket_unresolved",
    }
    missing = sorted(required - set(variables))
    if missing:
        raise RuntimeError("%s lacks required C13 columns: %s" %
                           (path, ", ".join(missing)))
    if not numeric_rows:
        raise RuntimeError("No complete numeric rows parsed from %s" % path)
    index = {name: variables.index(name) for name in variables}

    def optional_flag(values: Sequence[float], name: str) -> int:
        return int(round(values[index[name]])) if name in index else 0

    rows: List[CutoffRow] = []
    seen = set()
    for values in numeric_rows:
        row = CutoffRow(
            lon_deg=values[index["lon_deg"]],
            lat_deg=values[index["lat_deg"]],
            rc_lower_gv=values[index["rc_lower_gv"]],
            rc_effective_gv=values[index["rc_effective_gv"]],
            rc_upper_gv=values[index["rc_upper_gv"]],
            n_unresolved=int(round(values[index["n_unresolved"]])),
            lower_bracket_unresolved=int(round(values[index["lower_bracket_unresolved"]])),
            upper_bracket_unresolved=int(round(values[index["upper_bracket_unresolved"]])),
            lower_below_range=optional_flag(values, "lower_below_range"),
            lower_above_range=optional_flag(values, "lower_above_range"),
            upper_below_range=optional_flag(values, "upper_below_range"),
            upper_above_range=optional_flag(values, "upper_above_range"),
        )
        numeric = (row.lon_deg, row.lat_deg, row.rc_lower_gv,
                   row.rc_effective_gv, row.rc_upper_gv)
        if not all(math.isfinite(value) for value in numeric):
            raise RuntimeError("Non-finite C13 scientific value at key %r" % (row.key,))
        if row.key in seen:
            raise RuntimeError("Duplicate C13 output coordinate %r" % (row.key,))
        seen.add(row.key)
        rows.append(row)
    return sorted(rows, key=lambda row: row.key)


def find_penumbra_output(workdir: Path) -> Path:
    """Locate the documented output name, then compatible legacy spellings."""
    candidates = (
        workdir / "cutoff_gridless_shells_penumbra.dat",
        workdir / "cutoff_gridless_shells.dat",
        workdir / "cutoff_gridless_shells_dipole_compare.dat",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    matches = sorted(workdir.glob("cutoff*gridless*shell*penumbra*.dat"))
    if matches:
        return matches[0]
    raise RuntimeError("No GRIDLESS PENUMBRA_SCAN shell product in %s" % workdir)


def scientific_fingerprint(rows: Sequence[CutoffRow]) -> str:
    """Hash the canonical binary64/int scientific table, independent of row order."""
    digest = hashlib.sha256()
    for row in sorted(rows, key=lambda item: item.key):
        for value in (row.lon_deg % 360.0, row.lat_deg, row.rc_lower_gv,
                      row.rc_effective_gv, row.rc_upper_gv):
            digest.update(struct.pack("!d", float(value)))
        for value in (
                row.n_unresolved, row.lower_bracket_unresolved,
                row.upper_bracket_unresolved, row.lower_below_range,
                row.lower_above_range, row.upper_below_range,
                row.upper_above_range):
            digest.update(struct.pack("!q", int(value)))
    return digest.hexdigest()


def compare_products(baseline: Sequence[CutoffRow], candidate: Sequence[CutoffRow],
                     exact: bool, relative_tolerance: float,
                     absolute_tolerance: float) -> ProductComparison:
    """Compare coordinate-complete scientific products with strict flags."""
    base_map = {row.key: row for row in baseline}
    cand_map = {row.key: row for row in candidate}
    missing = sorted(set(base_map) - set(cand_map))
    extra = sorted(set(cand_map) - set(base_map))
    mismatches = 0
    flag_mismatches = 0
    compared = 0
    max_abs = 0.0
    max_rel = 0.0
    for key in sorted(set(base_map) & set(cand_map)):
        left = base_map[key]
        right = cand_map[key]
        for a, b in zip(
                (left.rc_lower_gv, left.rc_effective_gv, left.rc_upper_gv),
                (right.rc_lower_gv, right.rc_effective_gv, right.rc_upper_gv)):
            compared += 1
            difference = abs(a - b)
            scale = max(abs(a), abs(b), 1.0e-300)
            relative = difference / scale
            max_abs = max(max_abs, difference)
            max_rel = max(max_rel, relative)
            if exact:
                equal = struct.pack("!d", a) == struct.pack("!d", b)
            else:
                equal = difference <= absolute_tolerance + relative_tolerance * scale
            if not equal:
                mismatches += 1
        left_flags = (
            left.n_unresolved, left.lower_bracket_unresolved,
            left.upper_bracket_unresolved, left.lower_below_range,
            left.lower_above_range, left.upper_below_range, left.upper_above_range)
        right_flags = (
            right.n_unresolved, right.lower_bracket_unresolved,
            right.upper_bracket_unresolved, right.lower_below_range,
            right.lower_above_range, right.upper_below_range, right.upper_above_range)
        flag_mismatches += sum(a != b for a, b in zip(left_flags, right_flags))
    fingerprints_equal = scientific_fingerprint(baseline) == scientific_fingerprint(candidate)
    passed = (
        not missing and not extra and mismatches == 0 and flag_mismatches == 0 and
        (fingerprints_equal if exact else True)
    )
    return ProductComparison(
        baseline_rows=len(baseline), candidate_rows=len(candidate),
        missing_keys=len(missing), extra_keys=len(extra),
        numeric_values_compared=compared, numeric_mismatches=mismatches,
        flag_mismatches=flag_mismatches,
        max_absolute_difference=max_abs, max_relative_difference=max_rel,
        exact_fingerprint_match=fingerprints_equal, passed=passed)


def compare_stormer(rows: Sequence[CutoffRow], args: argparse.Namespace
                    ) -> Tuple[StormerComparison, List[Dict[str, object]]]:
    """Evaluate DIPOLE and retain a coordinate-level reference audit table."""
    # A finite rigidity scan cannot validate cutoffs below its first useful bin.
    # The default |latitude|<=60 region remains comfortably above that floor and
    # covers both low- and high-cutoff parts of the shell.
    required = [
        row for row in rows
        if args.stormer_min_abs_lat_deg <= abs(row.lat_deg) <=
        args.stormer_max_abs_lat_deg and
        stormer_vertical_gv(row.lat_deg, args.shell_alt_km) >= args.stormer_floor_gv
    ]
    scan_step = (
        (args.rigidity_max_gv - args.rigidity_min_gv) /
        max(1, args.cutoff_scan_n - 1))
    tolerance = args.stormer_abs_tolerance_gv + args.stormer_scan_step_factor * scan_step
    valid = []
    errors = []
    passing = 0
    details: List[Dict[str, object]] = []
    required_keys = {row.key for row in required}
    for row in rows:
        expected = stormer_vertical_gv(row.lat_deg, args.shell_alt_km)
        is_required = row.key in required_keys
        flags = (
            row.lower_bracket_unresolved, row.upper_bracket_unresolved,
            row.lower_below_range, row.lower_above_range,
            row.upper_below_range, row.upper_above_range)
        is_valid = is_required and row.n_unresolved == 0 and not any(flags)
        error = abs(row.rc_effective_gv - expected) if is_valid else math.nan
        allowed = tolerance + args.stormer_relative_tolerance * max(
            expected, args.stormer_floor_gv)
        point_passed = is_valid and error <= allowed
        if is_valid:
            valid.append(row)
            errors.append(error)
            passing += int(point_passed)
        if not is_required:
            if abs(row.lat_deg) < args.stormer_min_abs_lat_deg:
                exclusion = "below_moderate_latitude_band"
            elif abs(row.lat_deg) > args.stormer_max_abs_lat_deg:
                exclusion = "above_moderate_latitude_band"
            else:
                exclusion = "below_reference_floor"
        elif not is_valid:
            exclusion = "unresolved_or_out_of_range"
        else:
            exclusion = ""
        details.append({
            "lon_deg": row.lon_deg,
            "lat_deg": row.lat_deg,
            "altitude_km": args.shell_alt_km,
            "rc_stormer_gv": expected,
            "rc_effective_amps_gv": row.rc_effective_gv,
            "absolute_error_gv": error,
            "allowed_error_gv": allowed,
            "required_for_pass": int(is_required),
            "valid": int(is_valid),
            "passed": int(point_passed),
            "n_unresolved": row.n_unresolved,
            "lower_bracket_unresolved": row.lower_bracket_unresolved,
            "upper_bracket_unresolved": row.upper_bracket_unresolved,
            "lower_below_range": row.lower_below_range,
            "lower_above_range": row.lower_above_range,
            "upper_below_range": row.upper_below_range,
            "upper_above_range": row.upper_above_range,
            "exclusion_or_failure_reason": exclusion,
        })
    required_count = len(required)
    valid_count = len(valid)
    valid_fraction = float(valid_count) / required_count if required_count else 0.0
    passing_fraction = float(passing) / valid_count if valid_count else 0.0
    rmse = math.sqrt(sum(error * error for error in errors) / len(errors)) if errors else math.inf
    max_error = max(errors) if errors else math.inf
    passed = (
        required_count > 0 and
        valid_fraction >= args.stormer_required_valid_fraction and
        passing_fraction >= args.stormer_required_pass_fraction)
    metrics = StormerComparison(
        output_rows=len(rows), required_rows=required_count, valid_rows=valid_count,
        passing_rows=passing, valid_fraction=valid_fraction,
        passing_fraction=passing_fraction, rmse_gv=rmse,
        max_abs_error_gv=max_error, tolerance_gv=tolerance, passed=passed)
    return metrics, details


def case_record(case: RunCase, gridless_threads: int) -> Dict[str, object]:
    """Convert a case to JSON-safe provenance and audit content."""
    return {
        "case_id": case.case_id,
        "field_model": case.model,
        "scheduler": case.scheduler,
        "mpi_ranks": case.ranks,
        "gridless_threads": gridless_threads,
        "dynamic_chunk": case.dynamic_chunk,
        "workdir": str(case.workdir),
        "input_file": str(case.input_file or ""),
        "log_file": str(case.log_file or ""),
        "output_file": str(case.output_file or ""),
        "command": case.command or [],
        "return_code": case.return_code,
        "elapsed_seconds": case.elapsed_seconds,
        "output_sha256": case.output_sha256,
        "scientific_fingerprint": case.scientific_fingerprint,
        "task_diagnostics": (asdict(case.task_diagnostics)
                             if case.task_diagnostics else None),
        "product_comparison": (asdict(case.product_comparison)
                               if case.product_comparison else None),
        "stormer_comparison": (asdict(case.stormer_comparison)
                               if case.stormer_comparison else None),
        "per_rank_times_seconds": case.per_rank_times_seconds,
        "rank_time_imbalance": case.rank_time_imbalance,
        "speedup_vs_one_rank": case.speedup_vs_one_rank,
        "parallel_efficiency": case.parallel_efficiency,
        "elapsed_ratio_to_best_same_rank": case.elapsed_ratio_to_best_same_rank,
        "passed": case.passed,
        "error": case.error,
    }


def write_csv(path: Path, rows: Sequence[Mapping[str, object]],
              fieldnames: Sequence[str]) -> None:
    """Write deterministic summary tables, including a header for empty rows."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: object) -> object:
    """Recursively replace non-finite floats so result JSON stays standards-safe."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def write_outputs(work_root: Path, cases: Sequence[RunCase], result: Dict[str, object],
                  args: argparse.Namespace,
                  stormer_detail_rows: Sequence[Mapping[str, object]] = ()) -> None:
    """Write compact CSV views plus the complete machine-readable JSON result."""
    summary_rows = []
    rank_rows = []
    for case in cases:
        task = case.task_diagnostics
        comparison = case.product_comparison
        stormer = case.stormer_comparison
        summary_rows.append({
            "case_id": case.case_id,
            "field_model": case.model,
            "scheduler": case.scheduler,
            "mpi_ranks": case.ranks,
            "gridless_threads": args.nt,
            "dynamic_chunk": case.dynamic_chunk,
            "return_code": case.return_code,
            "elapsed_seconds": case.elapsed_seconds,
            "speedup_vs_one_rank": case.speedup_vs_one_rank if case.speedup_vs_one_rank is not None else "",
            "parallel_efficiency": case.parallel_efficiency if case.parallel_efficiency is not None else "",
            "elapsed_ratio_to_best_same_rank": (case.elapsed_ratio_to_best_same_rank
                                                  if case.elapsed_ratio_to_best_same_rank is not None else ""),
            "geometry_expected_tasks": expected_shell_task_count(args),
            "expected_tasks": task.expected_tasks if task else "",
            "computed_task_sum": task.computed_sum if task else "",
            "task_imbalance_max_over_avg": task.task_imbalance if task else "",
            "task_closure_passed": int(task.closure_passed) if task else 0,
            "numeric_mismatches": comparison.numeric_mismatches if comparison else "",
            "flag_mismatches": comparison.flag_mismatches if comparison else "",
            "max_relative_difference": comparison.max_relative_difference if comparison else "",
            "stormer_valid_fraction": stormer.valid_fraction if stormer else "",
            "stormer_passing_fraction": stormer.passing_fraction if stormer else "",
            "rank_time_imbalance": case.rank_time_imbalance if case.rank_time_imbalance is not None else "",
            "passed": int(case.passed),
            "error": case.error,
        })
        if task:
            for rank, tasks in task.rank_tasks:
                times = dict(case.per_rank_times_seconds or ())
                rank_rows.append({
                    "case_id": case.case_id,
                    "field_model": case.model,
                    "scheduler": case.scheduler,
                    "mpi_ranks": case.ranks,
                    "dynamic_chunk": case.dynamic_chunk,
                    "rank": rank,
                    "tasks": tasks,
                    "task_fraction": float(tasks) / task.computed_sum,
                    "rank_wall_seconds": times.get(rank, ""),
                })
    write_csv(work_root / "C13_summary.csv", summary_rows, (
        "case_id", "field_model", "scheduler", "mpi_ranks", "gridless_threads",
        "dynamic_chunk", "return_code", "elapsed_seconds", "speedup_vs_one_rank",
        "parallel_efficiency", "elapsed_ratio_to_best_same_rank",
        "geometry_expected_tasks", "expected_tasks",
        "computed_task_sum", "task_imbalance_max_over_avg", "task_closure_passed",
        "numeric_mismatches", "flag_mismatches", "max_relative_difference",
        "stormer_valid_fraction", "stormer_passing_fraction",
        "rank_time_imbalance", "passed", "error"))
    write_csv(work_root / "C13_rank_distribution.csv", rank_rows, (
        "case_id", "field_model", "scheduler", "mpi_ranks", "dynamic_chunk",
        "rank", "tasks", "task_fraction", "rank_wall_seconds"))
    write_csv(work_root / "C13_stormer_comparison.csv", stormer_detail_rows, (
        "lon_deg", "lat_deg", "altitude_km", "rc_stormer_gv",
        "rc_effective_amps_gv", "absolute_error_gv", "allowed_error_gv",
        "required_for_pass", "valid", "passed", "n_unresolved",
        "lower_bracket_unresolved", "upper_bracket_unresolved",
        "lower_below_range", "lower_above_range", "upper_below_range",
        "upper_above_range", "exclusion_or_failure_reason"))
    (work_root / "C13_result.json").write_text(
        json.dumps(json_safe(result), indent=2, allow_nan=False) + "\n")


def apply_profile(args: argparse.Namespace) -> None:
    """Fill only unspecified options from the selected profile."""
    defaults = PROFILE_DEFAULTS[args.profile]
    for name, value in defaults.items():
        if getattr(args, name) is None:
            setattr(args, name, value)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Define C13's complete CLI, including offline verification modes."""
    parser = argparse.ArgumentParser(
        description="C13 GRIDLESS MPI scheduler correctness/load-distribution test",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Recommended validation run from the AMPS repository root:
          python3 srcEarth/test/C13/run_C13.py --profile ROUTINE --amps ./amps -nt 1

        Full planned rank/shell campaign:
          python3 srcEarth/test/C13/run_C13.py --profile THOROUGH --amps ./amps -nt 1

        Package and command checks without AMPS:
          python3 srcEarth/test/C13/run_C13.py --self-test
          python3 srcEarth/test/C13/run_C13.py --profile SMOKE --dry-run

        Reanalyze an existing matrix:
          python3 srcEarth/test/C13/run_C13.py --profile ROUTINE --skip-run
        """))
    parser.add_argument("--profile", choices=sorted(PROFILE_DEFAULTS), default="ROUTINE")
    parser.add_argument("--models", default=None, help="comma-separated DIPOLE,T05")
    parser.add_argument("--ranks", default=None, help="MPI rank list; profile default")
    parser.add_argument("--schedulers", default="DYNAMIC,BLOCK_CYCLIC,STATIC")
    parser.add_argument("--chunks", default=None, help="explicit DYNAMIC chunks")
    parser.add_argument("-nt", type=int, default=1,
                        help="GRIDLESS threads/rank; 1 isolates MPI scheduling")
    parser.add_argument("--mover", default="BORIS",
                        choices=("BORIS", "RK2", "RK4", "RK6"))
    parser.add_argument("--amps", default="./amps")
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--mpirun-arg", action="append", default=[],
                        help="extra launcher argument; repeat as needed")
    parser.add_argument("--workdir", default="test_output/C13_gridless_scheduler")
    parser.add_argument("--shell-alt-km", type=float, default=500.0)
    parser.add_argument("--shell-lon-res-deg", type=float, default=None)
    parser.add_argument("--shell-lat-res-deg", type=float, default=None)
    parser.add_argument("--rigidity-min-gv", type=float, default=0.05)
    parser.add_argument("--rigidity-max-gv", type=float, default=15.0)
    parser.add_argument("--cutoff-scan-n", type=int, default=None)
    parser.add_argument("--dt-trace", type=float, default=0.20)
    parser.add_argument("--max-trace-time", type=float, default=None)
    parser.add_argument("--dipole-max-trace-time", type=float, default=None,
                        help="DIPOLE reference trace budget; profile default")
    parser.add_argument("--max-steps", type=int, default=2000000)
    parser.add_argument("--max-trace-distance-re", type=float, default=400.0)
    parser.add_argument("--t05-relative-tolerance", type=float, default=1.0e-7)
    parser.add_argument("--t05-absolute-tolerance-gv", type=float, default=1.0e-12)
    parser.add_argument("--max-dynamic-task-imbalance", type=float, default=0.0,
                        help="0 reports task-count imbalance only; >0 makes it a local gate")
    parser.add_argument("--max-rank-time-imbalance", type=float, default=0.0,
                        help="0 reports timing only; >0 requires explicit rank times")
    parser.add_argument("--stormer-min-abs-lat-deg", type=float,
                        default=DEFAULT_STORMER_MIN_ABS_LAT_DEG)
    parser.add_argument("--stormer-max-abs-lat-deg", type=float,
                        default=DEFAULT_STORMER_MAX_ABS_LAT_DEG)
    parser.add_argument("--stormer-floor-gv", type=float,
                        default=DEFAULT_STORMER_FLOOR_GV)
    parser.add_argument("--stormer-abs-tolerance-gv", type=float, default=0.10)
    parser.add_argument("--stormer-relative-tolerance", type=float, default=0.05)
    parser.add_argument("--stormer-scan-step-factor", type=float, default=1.50)
    parser.add_argument("--stormer-required-valid-fraction", type=float, default=0.95)
    parser.add_argument("--stormer-required-pass-fraction", type=float, default=0.95)
    parser.add_argument("--skip-run", action="store_true")
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--version", action="version",
                        version="C13 runner schema %d (%s)" %
                                (RUNNER_SCHEMA_VERSION, RUNNER_RELEASE))
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> Tuple[List[str], List[int], List[str], List[int]]:
    """Resolve list options and reject unsafe/degenerate validation requests."""
    models = parse_csv_names(args.models, "model", ALLOWED_MODELS)
    ranks = parse_csv_ints(args.ranks, "--ranks", 1)
    schedulers = parse_csv_names(args.schedulers, "scheduler", ALLOWED_SCHEDULERS)
    chunks = parse_csv_ints(args.chunks, "--chunks", 0)
    if "DYNAMIC" not in schedulers:
        raise SystemExit("C13 requires DYNAMIC in --schedulers as its rank-scaling baseline")
    if args.nt < 1:
        raise SystemExit("-nt must be >= 1")
    positive = (
        "shell_alt_km", "shell_lon_res_deg", "shell_lat_res_deg", "rigidity_min_gv",
        "rigidity_max_gv", "dt_trace", "max_trace_time", "dipole_max_trace_time",
        "max_trace_distance_re", "stormer_floor_gv")
    for name in positive:
        if getattr(args, name) <= 0.0:
            raise SystemExit("--%s must be positive" % name.replace("_", "-"))
    if args.rigidity_max_gv <= args.rigidity_min_gv:
        raise SystemExit("--rigidity-max-gv must exceed --rigidity-min-gv")
    if args.cutoff_scan_n < 3 or args.max_steps < 1:
        raise SystemExit("--cutoff-scan-n must be >=3 and --max-steps >=1")
    for name, total in (("shell_lon_res_deg", 360.0), ("shell_lat_res_deg", 180.0)):
        value = getattr(args, name)
        if not math.isclose(total / value, round(total / value),
                            rel_tol=0.0, abs_tol=1.0e-9):
            raise SystemExit("--%s must divide %g exactly" %
                             (name.replace("_", "-"), total))
    for name in ("t05_relative_tolerance", "t05_absolute_tolerance_gv",
                 "max_dynamic_task_imbalance", "max_rank_time_imbalance",
                 "stormer_abs_tolerance_gv",
                 "stormer_relative_tolerance", "stormer_scan_step_factor"):
        if getattr(args, name) < 0.0:
            raise SystemExit("--%s must be nonnegative" % name.replace("_", "-"))
    for name in ("stormer_required_valid_fraction", "stormer_required_pass_fraction"):
        if not 0.0 < getattr(args, name) <= 1.0:
            raise SystemExit("--%s must be in (0,1]" % name.replace("_", "-"))
    if not (0.0 <= args.stormer_min_abs_lat_deg <
            args.stormer_max_abs_lat_deg <= 90.0):
        raise SystemExit(
            "require 0 <= stormer-min-abs-lat-deg < "
            "stormer-max-abs-lat-deg <= 90")
    return models, ranks, schedulers, chunks


def self_test() -> int:
    """Exercise successful and deliberately broken package logic offline."""
    validate_template(DEFAULT_TEMPLATE)
    driver = validate_driver(DEFAULT_DRIVER)
    references = load_stormer_reference(DEFAULT_STORMER_REFERENCE)
    acceptance = validate_acceptance_reference(DEFAULT_ACCEPTANCE_REFERENCE)
    if driver["records"] != 5 or len(references) != 19 or len(acceptance) != 9:
        raise AssertionError("C13 checked-in reference counts changed unexpectedly")

    good_log = textwrap.dedent("""
        [gridless][MPI] Task distribution check:
          totalTasks (expected) = 100
          sum(all rank tasks)   = 100
          per-rank min/avg/max  = 24 / 25.0 / 26
            rank 0: 26 tasks
            rank 1: 24 tasks
            rank 2: 24 tasks
            rank 3: 26 tasks
        rank 0 wall time = 1.0 s
        rank 1 wall time = 1.1 s
        rank 2 wall time = 0.9 s
        rank 3 wall time = 1.0 s
    """)
    task = parse_task_diagnostics(good_log, 4)
    if not task.closure_passed or not math.isclose(task.task_imbalance, 1.04):
        raise AssertionError("C13 task parser rejected a valid closure block")
    times = parse_per_rank_times(good_log, 4)
    if times is None or len(times) != 4:
        raise AssertionError("C13 timing parser rejected complete explicit times")
    try:
        parse_task_diagnostics(good_log.replace("sum(all rank tasks)   = 100",
                                                "sum(all rank tasks)   = 99"), 4)
    except RuntimeError:
        raise AssertionError("Closure arithmetic should return a failed audit, not raise")
    broken = parse_task_diagnostics(
        good_log.replace("sum(all rank tasks)   = 100", "sum(all rank tasks)   = 99"), 4)
    if broken.closure_passed:
        raise AssertionError("C13 task parser accepted a corrupted reported sum")

    with tempfile.TemporaryDirectory(prefix="C13_selftest_") as temporary:
        root = Path(temporary)
        tecplot = root / "cutoff_gridless_shells_penumbra.dat"
        tecplot.write_text(textwrap.dedent('''
            TITLE = "C13 synthetic"
            VARIABLES = "lon_deg" "lat_deg" "Rc_lower_GV"
                        "Rc_effective_GV" "Rc_upper_GV" "n_unresolved"
                        "lower_bracket_unresolved" "upper_bracket_unresolved"
                        "lower_below_range" "lower_above_range"
                        "upper_below_range" "upper_above_range"
            ZONE T="alt_km=500"
            0 0 12.8 12.81 12.82 0 0 0 0 0 0 0
            30 0 12.8 12.81 12.82 0 0 0 0 0 0 0
        ''').strip() + "\n")
        rows = parse_penumbra_output(tecplot)
        if len(rows) != 2:
            raise AssertionError("C13 Tecplot parser row count mismatch")
        identical = compare_products(rows, list(reversed(rows)), True, 0.0, 0.0)
        if not identical.passed:
            raise AssertionError("C13 exact comparison depends on row order")
        changed = list(rows)
        changed[0] = CutoffRow(
            **{**asdict(changed[0]), "rc_effective_gv": changed[0].rc_effective_gv + 1.0e-3})
        if compare_products(rows, changed, True, 0.0, 0.0).passed:
            raise AssertionError("C13 exact comparison missed a perturbed cutoff")
        if not compare_products(rows, changed, False, 1.0e-3, 0.0).passed:
            raise AssertionError("C13 tolerant comparison rejected an in-tolerance value")

        args = parse_args(["--profile", "SMOKE"])
        apply_profile(args)
        models, ranks, schedulers, chunks = validate_args(args)
        exact_rc = stormer_vertical_gv(30.0, args.shell_alt_km)
        storm_rows = [
            CutoffRow(lon, 30.0, exact_rc, exact_rc, exact_rc, 0, 0, 0)
            for lon in (0.0, 60.0)
        ]
        storm_metrics, storm_details = compare_stormer(storm_rows, args)
        if not storm_metrics.passed or len(storm_details) != 2:
            raise AssertionError("C13 Størmer detail audit rejected exact controls")
        cases = make_cases(root / "runs", models, ranks, schedulers, chunks)
        if len(cases) != 6:
            raise AssertionError("C13 SMOKE matrix should contain six DIPOLE cases")
        if expected_shell_task_count(args) != 42:
            raise AssertionError("C13 SMOKE shell geometry should contain 42 tasks")
        for case in cases:
            case.elapsed_seconds = 100.0 / case.ranks
        compute_performance_metrics(cases, models)
        if not all(math.isclose(case.speedup_vs_one_rank or 0.0, case.ranks)
                   for case in cases):
            raise AssertionError("C13 speedup calculation failed ideal scaling")
        case = cases[0]
        case.workdir.mkdir(parents=True)
        render_input(DEFAULT_TEMPLATE, case.workdir / "AMPS_PARAM_C13.in",
                     case, DEFAULT_DRIVER.name, args)
        rendered = (case.workdir / "AMPS_PARAM_C13.in").read_text()
        if "__" in "\n".join(line for line in rendered.splitlines()
                               if not line.lstrip().startswith("!")):
            raise AssertionError("C13 rendered input retains a placeholder")

    print("C13 self-test PASS")
    print("  template: parser-safe scheduling/PENUMBRA_SCAN contract")
    print("  references: 19 analytical Størmer rows, 9 acceptance gates")
    print("  driver: %d rows with canonical DST column" % driver["records"])
    print("  parsers: task closure, optional rank timing, multiline Tecplot")
    print("  comparisons: exact DIPOLE and tolerant T05 negative paths")
    print("  metrics: shell task inventory and whole-job speedup/efficiency")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Prepare, execute, audit, summarize, and return a CI-friendly status."""
    args = parse_args(argv)
    if args.self_test:
        return self_test()
    apply_profile(args)
    models, ranks, schedulers, chunks = validate_args(args)
    validate_template(DEFAULT_TEMPLATE)
    stormer_reference = load_stormer_reference(DEFAULT_STORMER_REFERENCE)
    acceptance_reference = validate_acceptance_reference(DEFAULT_ACCEPTANCE_REFERENCE)
    driver_metadata = validate_driver(DEFAULT_DRIVER) if "T05" in models else None

    launch_dir = Path.cwd().resolve()
    work_root = Path(args.workdir).expanduser()
    if not work_root.is_absolute():
        work_root = (launch_dir / work_root).resolve()
    amps = Path(args.amps).expanduser()
    if not amps.is_absolute():
        amps = (launch_dir / amps).resolve()

    if args.skip_run:
        if not work_root.is_dir():
            raise SystemExit("--skip-run workdir does not exist: %s" % work_root)
    else:
        if work_root.exists() and not args.keep:
            shutil.rmtree(work_root)
        work_root.mkdir(parents=True, exist_ok=True)

    cases = make_cases(work_root, models, ranks, schedulers, chunks)
    print("C13 execution plan: %d cases (%s), ranks=%s, nt=%d" %
          (len(cases), ",".join(models), ",".join(map(str, ranks)), args.nt))
    if args.nt > 1 and "T05" in models:
        print(
            "C13 HYBRID warning: T05 will be evaluated concurrently by %d "
            "threads/rank. This is an intentional thread-safety stress test; "
            "large schedule-dependent cutoff or flag changes are hard failures. "
            "Use -nt 1 for the canonical MPI-only validation." % args.nt)

    # Preserve measured whole-job times during --skip-run reanalysis.  AMPS logs
    # contain task diagnostics but the runner itself measures elapsed time, so
    # the previous result JSON is the authoritative source for this diagnostic.
    previous_cases: Dict[str, Mapping[str, object]] = {}
    previous_result_path = work_root / "C13_result.json"
    if args.skip_run and previous_result_path.is_file():
        try:
            previous_result = json.loads(previous_result_path.read_text())
            previous_cases = {
                str(item.get("case_id")): item
                for item in previous_result.get("cases", [])
                if item.get("case_id")
            }
        except Exception as exc:
            print("C13 warning: cannot recover prior elapsed times: %s" % exc,
                  file=sys.stderr)

    # References are copied into results so an archived run remains intelligible
    # even when separated from the source tree.
    if not args.skip_run:
        shutil.copy2(DEFAULT_STORMER_REFERENCE,
                     work_root / DEFAULT_STORMER_REFERENCE.name)
        shutil.copy2(DEFAULT_ACCEPTANCE_REFERENCE,
                     work_root / DEFAULT_ACCEPTANCE_REFERENCE.name)

    for case in cases:
        case.workdir.mkdir(parents=True, exist_ok=True)
        case.input_file = case.workdir / "AMPS_PARAM_C13.in"
        case.log_file = case.workdir / "C13_amps.log"
        driver_name = DEFAULT_DRIVER.name
        if not args.skip_run:
            render_input(DEFAULT_TEMPLATE, case.input_file, case, driver_name, args)
            if case.model == "T05":
                shutil.copy2(DEFAULT_DRIVER, case.workdir / driver_name)
        case.command = command_for(case, amps, case.input_file.name, args)
        if case.case_id in previous_cases:
            prior = previous_cases[case.case_id]
            elapsed = prior.get("elapsed_seconds")
            if isinstance(elapsed, (int, float)) and math.isfinite(float(elapsed)):
                case.elapsed_seconds = float(elapsed)
            prior_rc = prior.get("return_code")
            if isinstance(prior_rc, int):
                case.return_code = prior_rc
        print("C13 %s command:\n  %s" % (case.case_id, " ".join(case.command)))
        if not args.skip_run and not args.dry_run:
            case.return_code, case.elapsed_seconds = run_process(
                case.command, case.workdir, case.log_file)

    if args.dry_run:
        result = {
            "test_id": TEST_ID,
            "test_name": TEST_NAME,
            "runner_schema_version": RUNNER_SCHEMA_VERSION,
            "runner_release": RUNNER_RELEASE,
            "profile": args.profile,
            "dry_run": True,
            "passed": None,
            "configuration": vars(args),
            "t05_driver": driver_metadata,
            "reference_rows": len(stormer_reference),
            "acceptance_gates": len(acceptance_reference),
            "cases": [case_record(case, args.nt) for case in cases],
        }
        # argparse's Path-free values are JSON-safe except no custom objects are
        # used here; retaining the exact CLI resolution aids reproduction.
        write_outputs(work_root, cases, result, args)
        print("C13 dry run complete: %d inputs under %s" % (len(cases), work_root))
        return 0

    output_rows: Dict[str, List[CutoffRow]] = {}
    stormer_detail_rows: List[Dict[str, object]] = []
    failures: List[str] = []
    for case in cases:
        if args.skip_run:
            case.return_code = 0 if case.return_code is None else case.return_code
        if case.return_code != 0:
            case.error = "AMPS exit code %s; see %s" % (case.return_code, case.log_file)
            failures.append("%s: %s" % (case.case_id, case.error))
            continue
        try:
            assert case.log_file is not None
            log_text = case.log_file.read_text(errors="replace")
            case.task_diagnostics = parse_task_diagnostics(log_text, case.ranks)
            case.per_rank_times_seconds = parse_per_rank_times(log_text, case.ranks)
            if case.per_rank_times_seconds:
                times = [value for _, value in case.per_rank_times_seconds]
                average = sum(times) / len(times)
                case.rank_time_imbalance = max(times) / average if average > 0.0 else math.inf
            case.output_file = find_penumbra_output(case.workdir)
            rows = parse_penumbra_output(case.output_file)
            output_rows[case.case_id] = rows
            case.output_sha256 = hashlib.sha256(case.output_file.read_bytes()).hexdigest()
            case.scientific_fingerprint = scientific_fingerprint(rows)
        except Exception as exc:
            case.error = str(exc)
            failures.append("%s: %s" % (case.case_id, case.error))

    # The one-rank DYNAMIC/auto case is the sole numerical baseline for each
    # model.  This prevents a majority of equally wrong parallel cases from
    # creating their own consensus reference.
    for model in models:
        baseline = next((case for case in cases if
                         case.model == model and case.scheduler == "DYNAMIC" and
                         case.ranks == 1 and case.dynamic_chunk == 0), None)
        if baseline is None or baseline.case_id not in output_rows:
            failures.append("%s: missing one-rank DYNAMIC baseline" % model)
            continue
        base_rows = output_rows[baseline.case_id]
        if model == "DIPOLE":
            baseline.stormer_comparison, stormer_detail_rows = compare_stormer(
                base_rows, args)
        for case in [item for item in cases if item.model == model and
                     item.case_id in output_rows]:
            case.product_comparison = compare_products(
                base_rows, output_rows[case.case_id], exact=(model == "DIPOLE"),
                relative_tolerance=args.t05_relative_tolerance,
                absolute_tolerance=args.t05_absolute_tolerance_gv)

    # Whole-job timing is informative even when explicit per-rank timers are not
    # available.  Report speedup and efficiency, but keep them diagnostic by
    # default because machine load and process placement are not portable gates.
    compute_performance_metrics(cases, models)

    for case in cases:
        if case.error:
            continue
        reasons = []
        if case.task_diagnostics is None or not case.task_diagnostics.closure_passed:
            reasons.append("task closure failed")
        elif case.task_diagnostics.expected_tasks != expected_shell_task_count(args):
            reasons.append("AMPS task inventory %d != shell geometry %d" %
                           (case.task_diagnostics.expected_tasks,
                            expected_shell_task_count(args)))
        if (case.scheduler == "DYNAMIC" and case.ranks > 1 and
                case.task_diagnostics is not None and
                args.max_dynamic_task_imbalance > 0.0 and
                case.task_diagnostics.task_imbalance > args.max_dynamic_task_imbalance):
            reasons.append("DYNAMIC task-count imbalance %.3f > %.3f" %
                           (case.task_diagnostics.task_imbalance,
                            args.max_dynamic_task_imbalance))
        if case.product_comparison is None or not case.product_comparison.passed:
            reasons.append("scientific product differs from one-rank baseline")
        if case.stormer_comparison is not None and not case.stormer_comparison.passed:
            reasons.append("analytical Størmer reference failed")
        if args.max_rank_time_imbalance > 0.0:
            if case.rank_time_imbalance is None:
                reasons.append("explicit per-rank timing diagnostics unavailable")
            elif case.rank_time_imbalance > args.max_rank_time_imbalance:
                reasons.append("rank-time imbalance %.3f > %.3f" %
                               (case.rank_time_imbalance,
                                args.max_rank_time_imbalance))
        if reasons:
            case.error = "; ".join(reasons)
            failures.append("%s: %s" % (case.case_id, case.error))
        else:
            case.passed = True

    passed = not failures and all(case.passed for case in cases)
    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "runner_release": RUNNER_RELEASE,
        "profile": args.profile,
        "passed": passed,
        "configuration": {
            "models": models,
            "ranks": ranks,
            "schedulers": schedulers,
            "chunks": chunks,
            "gridless_threads": args.nt,
            "mover": args.mover,
            "shell_alt_km": args.shell_alt_km,
            "shell_lon_res_deg": args.shell_lon_res_deg,
            "shell_lat_res_deg": args.shell_lat_res_deg,
            "rigidity_min_gv": args.rigidity_min_gv,
            "rigidity_max_gv": args.rigidity_max_gv,
            "cutoff_scan_n": args.cutoff_scan_n,
            "max_trace_time_s": args.max_trace_time,
            "dipole_max_trace_time_s": args.dipole_max_trace_time,
            "geometry_expected_tasks": expected_shell_task_count(args),
            "stormer_min_abs_lat_deg": args.stormer_min_abs_lat_deg,
            "stormer_max_abs_lat_deg": args.stormer_max_abs_lat_deg,
            "max_dynamic_task_imbalance": args.max_dynamic_task_imbalance,
            "max_rank_time_imbalance": args.max_rank_time_imbalance,
            "dipole_equivalence": "exact canonical binary64",
            "t05_relative_tolerance": args.t05_relative_tolerance,
        },
        "t05_driver": driver_metadata,
        "references": {
            "stormer_table": str(DEFAULT_STORMER_REFERENCE),
            "stormer_rows": len(stormer_reference),
            "acceptance_table": str(DEFAULT_ACCEPTANCE_REFERENCE),
            "acceptance_gates": len(acceptance_reference),
        },
        "failures": failures,
        "cases": [case_record(case, args.nt) for case in cases],
    }
    write_outputs(work_root, cases, result, args, stormer_detail_rows)

    print("\nC13 — %s" % TEST_NAME)
    print("result: %s" % ("PASS" if passed else "FAIL"))
    print("profile: %s; cases: %d; ranks: %s; nt: %d" %
          (args.profile, len(cases), ",".join(map(str, ranks)), args.nt))
    for case in cases:
        task_text = "n/a"
        if case.task_diagnostics:
            task_text = "%d tasks, max/avg=%.3f" % (
                case.task_diagnostics.computed_sum,
                case.task_diagnostics.task_imbalance)
        if case.speedup_vs_one_rank is not None:
            task_text += ", speedup=%.2fx" % case.speedup_vs_one_rank
        print("%-42s %-4s  %s" %
              (case.case_id, "PASS" if case.passed else "FAIL", task_text))
    if failures:
        print("\nFailure messages:")
        for failure in failures:
            print("- %s" % failure)
    print("\nC13 results: %s" % work_root)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
