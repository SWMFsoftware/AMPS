#!/usr/bin/env python3
"""Run the durable eight-case Stage-3 comet science-gate matrix.

This program is an *operational* runner, not the scientific evaluator.  It
freezes one already-configured executable plus its canonical input and
``picGlobal.dfn`` provenance, launches the approved mode/seed matrix one case
at a time, and records enough immutable state to recover after interruption.
Use ``evaluate_stage3_science_gate.py`` only after this runner reports all
eight cases complete.

Exit status contract
--------------------

``0``
    The requested operation succeeded.  For ``run``/``resume`` this means all
    eight operational runs and their required artifacts are complete; it does
    not mean that the scientific gate passed.
``2``
    Invalid command line, policy, or canonical input.
``3``
    Frozen provenance mismatch, archive-layout violation, active competing
    run, or dirty attempt directory.
``4``
    AMPS returned nonzero or a nominally successful run lacks required
    operational evidence.
``5``
    The runner received SIGINT/SIGTERM.  The attempt remains preserved and a
    later ``resume`` starts a new attempt for that case from iteration zero.

The science application has no validated production checkpoint/restart path.
Consequently this runner never reads particle restart data and never reuses a
partial ``PT`` tree.  Completed cases are skipped on resume, while an
interrupted or failed case receives a new, monotonically numbered attempt.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, TextIO, Tuple


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
sys.path.insert(0, str(THIS_DIR))

from stage3_gate_common import (  # noqa: E402
    GateError,
    MANIFEST_SCHEMA,
    atomic_write_json,
    expected_run_order,
    load_policy,
    parse_scalar_input_records,
    sha256_file,
    sha256_json,
    validate_canonical_input,
)


EXIT_SUCCESS = 0
EXIT_INVALID_POLICY = 2
EXIT_PREFLIGHT = 3
EXIT_OPERATIONAL = 4
EXIT_INTERRUPTED = 5

ATTEMPT_SCHEMA = "amps.comet.stage3-science-attempt.v1"
LOCK_FILE = ".runner.lock"
ATTEMPT_PATTERN = re.compile(r"^attempt-(\d{4})$")


class PreflightError(RuntimeError):
    """Frozen files, archive layout, or process ownership are unsafe."""


class OperationalError(RuntimeError):
    """A launched AMPS case failed or produced incomplete evidence."""


class RunInterrupted(RuntimeError):
    """The user or operating environment interrupted the active case."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_json(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PreflightError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PreflightError(f"JSON root is not an object: {path}")
    return value


def atomic_write_text(path: Path, text: str, mode: int = 0o644) -> None:
    """Publish a small text artifact with the same reboot-safe rename rule."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_copy(source: Path, destination: Path, mode: int) -> None:
    """Copy and fsync a frozen artifact before making its final name visible."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    try:
        with source.open("rb") as input_stream, os.fdopen(
            descriptor, "wb"
        ) as output_stream:
            shutil.copyfileobj(input_stream, output_stream, 1024 * 1024)
            output_stream.flush()
            os.fsync(output_stream.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, destination)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def command_version(command: str) -> str:
    """Capture launcher identity without making version formatting an API."""
    try:
        result = subprocess.run(
            [command, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise PreflightError(f"cannot execute MPI launcher {command}: {exc}") from exc
    first_lines = result.stdout.strip().splitlines()[:4]
    if result.returncode != 0 or not first_lines:
        raise PreflightError(f"MPI launcher version probe failed: {command}")
    return "\n".join(first_lines)


def validate_source_artifacts(
    executable: Path, input_path: Path, pic_global: Path, policy: Dict[str, Any]
) -> Dict[str, str]:
    """Validate and hash the three files that define the frozen runtime."""
    for label, path in (
        ("executable", executable),
        ("canonical input", input_path),
        ("picGlobal.dfn", pic_global),
    ):
        if not path.is_file():
            raise PreflightError(f"{label} is not a regular file: {path}")
    if not os.access(executable, os.X_OK):
        raise PreflightError(f"frozen executable is not executable: {executable}")

    try:
        validate_canonical_input(input_path, policy)
        records = parse_scalar_input_records(input_path)
    except OSError as exc:
        raise PreflightError(f"cannot read canonical input {input_path}: {exc}") from exc
    mode = records.get("ChemistryResidualSampling", "").upper()
    if mode not in {"BERNOULLI", "SYSTEMATIC"}:
        raise GateError("canonical input must explicitly select a residual mode")
    try:
        configured_seed = int(records["CometRandomSeed"])
    except (KeyError, ValueError) as exc:
        raise GateError("canonical input must contain an explicit valid CometRandomSeed") from exc
    if configured_seed < 1 or configured_seed > 2147483000:
        raise GateError("canonical CometRandomSeed is outside the accepted range")

    try:
        return {
            "executable": sha256_file(executable),
            "input": sha256_file(input_path),
            "pic_global": sha256_file(pic_global),
        }
    except OSError as exc:
        raise PreflightError(f"cannot hash a source artifact: {exc}") from exc


def frozen_paths(archive: Path) -> Dict[str, Path]:
    return {
        "executable": archive / "frozen" / "amps.stage3-science",
        "input": archive / "frozen" / "halley1d.input",
        "pic_global": archive / "frozen" / "picGlobal.dfn",
        "policy": archive / "policy.json",
    }


def verify_frozen(archive: Path, manifest: Dict[str, Any]) -> Dict[str, Path]:
    """Re-hash immutable evidence before every launch and after every run."""
    paths = frozen_paths(archive)
    expected = {
        "executable": manifest.get("frozen_executable_sha256"),
        "input": manifest.get("frozen_input_sha256"),
        "pic_global": manifest.get("frozen_pic_global_sha256"),
    }
    for name in ("executable", "input", "pic_global"):
        if not paths[name].is_file():
            raise PreflightError(f"missing frozen {name}: {paths[name]}")
        actual = sha256_file(paths[name])
        if not isinstance(expected[name], str) or actual != expected[name]:
            raise PreflightError(f"frozen {name} hash differs from RUN_MANIFEST.json")
    if not os.access(paths["executable"], os.X_OK):
        raise PreflightError("frozen executable lost its execute permission")

    try:
        policy = load_policy(paths["policy"])
    except GateError as exc:
        raise PreflightError(f"frozen policy is invalid: {exc}") from exc
    if sha256_json(policy) != manifest.get("policy_sha256"):
        raise PreflightError("frozen policy differs from RUN_MANIFEST.json")
    try:
        validate_canonical_input(paths["input"], policy)
    except (GateError, OSError) as exc:
        raise PreflightError(f"frozen canonical input is invalid: {exc}") from exc
    return paths


class ArchiveLock:
    """Hold one advisory lock for the complete scheduler invocation."""

    def __init__(self, archive: Path):
        self.path = archive / LOCK_FILE
        self.stream: TextIO | None = None

    def __enter__(self) -> "ArchiveLock":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.stream = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            self.stream.close()
            self.stream = None
            raise PreflightError(
                f"another Stage-3 runner owns {self.path}"
            ) from exc
        self.stream.seek(0)
        self.stream.truncate()
        self.stream.write(f"pid={os.getpid()} acquired={utc_now()}\n")
        self.stream.flush()
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if self.stream is not None:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_UN)
            self.stream.close()
            self.stream = None


def initialize_archive(
    archive: Path,
    policy_path: Path,
    executable: Path,
    input_path: Path,
    pic_global: Path,
    mpiexec: str,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Freeze provenance once; an existing manifest must use ``resume``."""
    policy = load_policy(policy_path)
    hashes = validate_source_artifacts(executable, input_path, pic_global, policy)

    archive.mkdir(parents=True, exist_ok=True)
    unexpected = [entry for entry in archive.iterdir() if entry.name != LOCK_FILE]
    if unexpected:
        raise PreflightError(
            f"run requires a new empty archive; found {unexpected[0]} (use resume)"
        )

    launcher = shutil.which(mpiexec)
    if launcher is None:
        raise PreflightError(f"MPI launcher is not on PATH: {mpiexec}")
    launcher_path = str(Path(launcher).resolve())
    paths = frozen_paths(archive)
    atomic_copy(executable, paths["executable"], 0o555)
    atomic_copy(input_path, paths["input"], 0o444)
    atomic_copy(pic_global, paths["pic_global"], 0o444)
    atomic_write_json(paths["policy"], policy)
    os.chmod(paths["policy"], 0o444)

    # Hash the copies, not merely their sources, before publishing the
    # manifest.  A failed or short copy can therefore never become evidence.
    copied_hashes = {
        "executable": sha256_file(paths["executable"]),
        "input": sha256_file(paths["input"]),
        "pic_global": sha256_file(paths["pic_global"]),
    }
    if copied_hashes != hashes:
        raise PreflightError("a frozen artifact changed while it was copied")

    manifest: Dict[str, Any] = {
        "schema": MANIFEST_SCHEMA,
        "created_utc": utc_now(),
        "updated_utc": utc_now(),
        "status": "initialized",
        "policy_sha256": sha256_json(policy),
        "policy_path": "policy.json",
        "frozen_executable": str(paths["executable"].relative_to(archive)),
        "frozen_executable_sha256": hashes["executable"],
        "frozen_input": str(paths["input"].relative_to(archive)),
        "frozen_input_sha256": hashes["input"],
        "frozen_pic_global": str(paths["pic_global"].relative_to(archive)),
        "frozen_pic_global_sha256": hashes["pic_global"],
        "source_artifacts": {
            "executable": str(executable.resolve()),
            "input": str(input_path.resolve()),
            "pic_global": str(pic_global.resolve()),
        },
        "mpi_launcher": launcher_path,
        "mpi_launcher_version": command_version(launcher_path),
        "host": platform.node(),
        "platform": platform.platform(),
        "python": sys.version.splitlines()[0],
        "environment": {"OMP_NUM_THREADS": "1", "LC_ALL": "C"},
        "execution_plan": expected_run_order(policy),
        "runs": [],
        "attempts": [],
    }
    atomic_write_json(archive / "RUN_MANIFEST.json", manifest)
    verify_frozen(archive, manifest)
    return policy, manifest


def load_archive(archive: Path) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    manifest_path = archive / "RUN_MANIFEST.json"
    if not manifest_path.is_file():
        raise PreflightError(f"archive has no RUN_MANIFEST.json: {archive}")
    manifest = load_json(manifest_path)
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise PreflightError("unknown Stage-3 science manifest schema")
    policy_path = frozen_paths(archive)["policy"]
    try:
        policy = load_policy(policy_path)
    except GateError as exc:
        raise PreflightError(f"frozen policy is invalid: {exc}") from exc
    if manifest.get("policy_sha256") != sha256_json(policy):
        raise PreflightError("manifest and frozen policy hashes differ")
    if manifest.get("execution_plan") != expected_run_order(policy):
        raise PreflightError("manifest execution plan differs from the approved policy")
    verify_frozen(archive, manifest)
    return policy, manifest


def process_start_ticks(pid: int) -> int | None:
    """Return Linux process start ticks, protecting against recycled PIDs."""
    try:
        fields = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8").split()
        return int(fields[21])
    except (OSError, ValueError, IndexError):
        return None


def process_matches(pid: object, ticks: object) -> bool:
    return (
        isinstance(pid, int)
        and isinstance(ticks, int)
        and process_start_ticks(pid) == ticks
    )


def expected_seed_directory(seed: int) -> str:
    return f"seed-{seed:010d}"


def scan_attempts(
    archive: Path, policy: Dict[str, Any], *, reject_live: bool = True
) -> List[Dict[str, Any]]:
    """Reject unowned files instead of deleting or guessing their meaning."""
    runs_root = archive / "runs"
    if not runs_root.exists():
        return []
    if not runs_root.is_dir():
        raise PreflightError(f"runs path is not a directory: {runs_root}")

    valid_seeds = {expected_seed_directory(seed): seed for seed in policy["seeds"]}
    attempts: List[Dict[str, Any]] = []
    for seed_path in sorted(runs_root.iterdir()):
        if not seed_path.is_dir() or seed_path.name not in valid_seeds:
            raise PreflightError(f"dirty or unexpected run entry: {seed_path}")
        seed = valid_seeds[seed_path.name]
        valid_modes = {mode.lower(): mode for mode in policy["modes"]}
        for mode_path in sorted(seed_path.iterdir()):
            if not mode_path.is_dir() or mode_path.name not in valid_modes:
                raise PreflightError(f"dirty or unexpected mode entry: {mode_path}")
            mode = valid_modes[mode_path.name]
            for attempt_path in sorted(mode_path.iterdir()):
                match = ATTEMPT_PATTERN.match(attempt_path.name)
                if not attempt_path.is_dir() or match is None:
                    raise PreflightError(f"dirty or unexpected attempt entry: {attempt_path}")
                status_path = attempt_path / "status.json"
                if not status_path.is_file():
                    raise PreflightError(f"attempt has no status.json: {attempt_path}")
                status = load_json(status_path)
                expected_relative = str(attempt_path.relative_to(archive))
                if (
                    status.get("schema") != ATTEMPT_SCHEMA
                    or status.get("seed") != seed
                    or status.get("mode") != mode
                    or status.get("attempt_path") != expected_relative
                    or status.get("attempt_number") != int(match.group(1))
                ):
                    raise PreflightError(f"attempt identity is inconsistent: {attempt_path}")
                if (
                    reject_live
                    and status.get("state") in {"launching", "running"}
                    and process_matches(
                        status.get("child_pid"), status.get("child_start_ticks")
                    )
                ):
                    raise PreflightError(f"attempt still has a live MPI launcher: {attempt_path}")
                attempts.append(status)
    return attempts


def validate_manifest_records(
    archive: Path,
    manifest: Dict[str, Any],
    policy: Dict[str, Any],
    attempts: List[Dict[str, Any]],
) -> None:
    plan = expected_run_order(policy)
    approved = {(item["mode"], item["seed"]) for item in plan}
    runs = manifest.get("runs")
    recorded_attempts = manifest.get("attempts")
    if not isinstance(runs, list) or not isinstance(recorded_attempts, list):
        raise PreflightError("manifest runs/attempts fields must be lists")
    seen = set()
    for run in runs:
        key = (run.get("mode"), run.get("seed"))
        if key not in approved or key in seen:
            raise PreflightError("manifest has a duplicate or unapproved completed run")
        seen.add(key)
        if run.get("status") != "completed" or run.get("return_code") != 0:
            raise PreflightError("manifest runs may contain only completed attempts")
        if (
            run.get("executable_sha256") != manifest.get("frozen_executable_sha256")
            or run.get("input_sha256") != manifest.get("frozen_input_sha256")
            or run.get("pic_global_sha256") != manifest.get("frozen_pic_global_sha256")
        ):
            raise PreflightError("completed run provenance differs from frozen artifacts")
        relative = run.get("attempt_path")
        if not isinstance(relative, str):
            raise PreflightError("completed run has no attempt_path")
        attempt_root = (archive / relative).resolve()
        try:
            attempt_root.relative_to(archive.resolve())
        except ValueError as exc:
            raise PreflightError("completed run attempt escapes the archive") from exc
        status_path = attempt_root / "status.json"
        if not status_path.is_file():
            raise PreflightError("completed run attempt status is missing")
        status = load_json(status_path)
        if status.get("state") != "completed" or status.get("return_code") != 0:
            raise PreflightError("completed run points to a non-completed attempt")
        checksum = attempt_root / str(status.get("artifact_manifest", ""))
        if (
            not checksum.is_file()
            or sha256_file(checksum) != run.get("artifact_manifest_sha256")
            or run.get("artifact_manifest_sha256")
            != status.get("artifact_manifest_sha256")
        ):
            raise PreflightError("completed run artifact manifest changed")
        verify_artifact_checksums(attempt_root, checksum)

    disk_paths = {item["attempt_path"] for item in attempts}
    manifest_paths = {item.get("attempt_path") for item in recorded_attempts}
    if None in manifest_paths or not manifest_paths.issubset(disk_paths):
        raise PreflightError("manifest refers to a missing attempt directory")
    # A crash can leave a valid atomic status on disk just before the top-level
    # manifest was refreshed.  Preserve it and let resume record it; do not
    # reject that recoverable ordering window.


def create_attempt(
    archive: Path, seed: int, mode: str, existing: Iterable[Dict[str, Any]]
) -> Tuple[Path, Dict[str, Any]]:
    numbers = [
        int(item["attempt_number"])
        for item in existing
        if item.get("seed") == seed and item.get("mode") == mode
    ]
    number = max(numbers, default=0) + 1
    mode_dir = archive / "runs" / expected_seed_directory(seed) / mode.lower()
    mode_dir.mkdir(parents=True, exist_ok=True)
    final = mode_dir / f"attempt-{number:04d}"
    if final.exists():
        raise PreflightError(f"refusing existing attempt directory: {final}")

    # Construct the identity file in a temporary sibling directory, then make
    # the complete attempt identity visible with one rename.  Resume therefore
    # never has to interpret a directory created only halfway before a reboot.
    preparing = mode_dir / f".attempt-{number:04d}.preparing-{os.getpid()}"
    if preparing.exists():
        raise PreflightError(f"dirty attempt-preparation directory: {preparing}")
    preparing.mkdir()
    status: Dict[str, Any] = {
        "schema": ATTEMPT_SCHEMA,
        "seed": seed,
        "mode": mode,
        "attempt_number": number,
        "attempt_path": str(final.relative_to(archive)),
        "state": "prepared",
        "created_utc": utc_now(),
        "restart_source": None,
        "start_iteration": 0,
        "child_pid": None,
        "child_start_ticks": None,
    }
    atomic_write_json(preparing / "status.json", status)
    os.replace(preparing, final)
    return final, status


def update_manifest_attempt(
    archive: Path, manifest: Dict[str, Any], status: Dict[str, Any]
) -> None:
    attempts = [
        item
        for item in manifest["attempts"]
        if item.get("attempt_path") != status["attempt_path"]
    ]
    attempts.append(
        {
            key: status.get(key)
            for key in (
                "attempt_path",
                "attempt_number",
                "seed",
                "mode",
                "state",
                "created_utc",
                "started_utc",
                "completed_utc",
                "return_code",
                "signal",
                "wall_seconds",
                "error",
            )
            if key in status
        }
    )
    manifest["attempts"] = sorted(
        attempts, key=lambda item: (item["seed"], item["mode"], item["attempt_number"])
    )
    manifest["updated_utc"] = utc_now()
    atomic_write_json(archive / "RUN_MANIFEST.json", manifest)


def add_completed_run(
    archive: Path,
    manifest: Dict[str, Any],
    policy: Dict[str, Any],
    run: Dict[str, Any],
) -> None:
    key = (run["mode"], run["seed"])
    retained = [
        item
        for item in manifest["runs"]
        if (item.get("mode"), item.get("seed")) != key
    ]
    retained.append(run)
    order = {
        (item["mode"], item["seed"]): index
        for index, item in enumerate(expected_run_order(policy))
    }
    manifest["runs"] = sorted(
        retained, key=lambda item: order[(item["mode"], item["seed"])]
    )
    manifest["status"] = (
        "completed" if len(manifest["runs"]) == 8 else "running"
    )
    manifest["updated_utc"] = utc_now()
    atomic_write_json(archive / "RUN_MANIFEST.json", manifest)


def parse_diagnostic_identity(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise OperationalError(f"cannot read diagnostic JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise OperationalError(f"diagnostic JSON root is not an object: {path}")
    if value.get("schema") != "amps.comet.stage0-diagnostics.v1":
        raise OperationalError(f"unknown diagnostic schema: {path}")
    return value


def validate_center_table(
    path: Path, seed: int, mode: str, ranks: int, sampling_length: int, output: int
) -> None:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            required = {
                "cell_index",
                "sampling_length",
                "output_index",
                "random_seed",
                "random_seed_configured",
                "residual_sampling",
                "mpi_ranks",
                "radius_m",
                "r_inner_m",
                "r_outer_m",
                "measure_m3",
                "number_density_m-3",
                "mean_model_particles",
                "block_index",
                "population_target_min_per_cell",
                "population_target_max_per_cell",
            }
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise OperationalError(f"incomplete center-cell schema: {path}")
            count = 0
            for count, row in enumerate(reader, start=1):
                if int(row["cell_index"]) != count - 1:
                    raise OperationalError(f"noncontiguous center-cell index: {path}")
                if (
                    int(row["sampling_length"]) != sampling_length
                    or int(row["output_index"]) != output
                    or int(row["random_seed"]) != seed
                    or int(row["random_seed_configured"]) != 1
                    or row["residual_sampling"].upper() != mode
                    or int(row["mpi_ranks"]) != ranks
                ):
                    raise OperationalError(f"center-cell provenance mismatch: {path}")
            if count == 0:
                raise OperationalError(f"empty center-cell table: {path}")
    except (OSError, ValueError, KeyError) as exc:
        raise OperationalError(f"cannot validate center-cell table {path}: {exc}") from exc


def validate_operational_artifacts(
    attempt: Path, seed: int, mode: str, policy: Dict[str, Any]
) -> Tuple[List[Path], Dict[str, Any]]:
    """Require complete, fresh, mode/seed-identifiable production evidence."""
    run_log = attempt / "run.log"
    error_log = attempt / "error.log"
    for path in (run_log, error_log):
        if not path.is_file():
            raise OperationalError(f"missing process log: {path}")
    combined = run_log.read_text(encoding="utf-8", errors="replace") + "\n" + error_log.read_text(
        encoding="utf-8", errors="replace"
    )
    if "Reading particle restart file" in combined:
        raise OperationalError("science attempt read a particle restart file")
    if f"The total number of the MPI processes={policy['mpi_ranks']}" not in combined:
        raise OperationalError("run log does not prove the approved MPI rank count")
    if f"comet random seed={seed}" not in combined:
        raise OperationalError("run log does not prove the requested random seed")

    diagnostics = attempt / "PT" / "diagnostics"
    if not diagnostics.is_dir():
        raise OperationalError(f"missing diagnostic directory: {diagnostics}")
    expected_outputs = list(range(policy["required_converged_windows"] + 1))
    found_reports: Dict[int, Path] = {}
    for path in diagnostics.glob("comet_stage0.out=*.json"):
        match = re.fullmatch(r"comet_stage0\.out=(\d+)\.json", path.name)
        if match is not None:
            found_reports[int(match.group(1))] = path
    if sorted(found_reports) != expected_outputs:
        raise OperationalError(
            f"complete diagnostic windows are {sorted(found_reports)}, expected {expected_outputs}"
        )

    required: List[Path] = [run_log, error_log, diagnostics / "cell_measures.csv"]
    if not required[-1].is_file():
        raise OperationalError("missing exact cell-measure table")
    for output in expected_outputs:
        path = found_reports[output]
        report = parse_diagnostic_identity(path)
        expected_iteration = (output + 1) * policy["sampling_length"]
        if (
            report.get("iteration") != expected_iteration
            or report.get("output_index") != output
            or report.get("sampling_length") != policy["sampling_length"]
            or report.get("random_seed") != seed
            or report.get("random_seed_configured") is not True
            or report.get("mpi_ranks") != policy["mpi_ranks"]
            or report.get("residual_sampling") != mode
        ):
            raise OperationalError(f"diagnostic provenance mismatch: {path}")
        required.append(path)
        csv_report = diagnostics / f"comet_stage0.out={output}.csv"
        if not csv_report.is_file():
            raise OperationalError(f"missing diagnostic CSV: {csv_report}")
        required.append(csv_report)
        for species in policy["important_species"]:
            table = diagnostics / f"stage3_cells.{species}.out={output}.csv"
            if not table.is_file():
                raise OperationalError(f"missing center-cell table: {table}")
            validate_center_table(
                table,
                seed,
                mode,
                policy["mpi_ranks"],
                policy["sampling_length"],
                output,
            )
            required.append(table)

    # Fresh iteration-zero execution is established jointly by the empty
    # attempt directory at launch, absence of a restart read, and the exact
    # out=0/iteration=sampling_length record validated above.
    metadata = {
        "complete_output_indices": expected_outputs,
        "final_iteration": policy["maximum_iterations"],
        "fresh_start_verified": True,
        "required_artifact_count": len(required),
    }
    return required, metadata


def write_artifact_checksums(attempt: Path, artifacts: Iterable[Path]) -> Tuple[Path, str]:
    rows = []
    for path in sorted(set(artifacts)):
        rows.append(f"{sha256_file(path)}  {path.relative_to(attempt)}")
    destination = attempt / "artifacts.sha256"
    atomic_write_text(destination, "\n".join(rows) + "\n", 0o444)
    return destination, sha256_file(destination)


def verify_artifact_checksums(attempt: Path, checksum_path: Path) -> None:
    """Verify every file named by an immutable attempt checksum manifest.

    Hashing only ``artifacts.sha256`` would detect edits to the manifest but
    not later edits to the output files it describes.  Resume and recovery use
    this stricter check so a completed attempt can never silently become a
    different scientific data set between invocations.
    """
    try:
        rows = checksum_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise PreflightError(
            f"cannot read artifact checksum manifest {checksum_path}: {exc}"
        ) from exc
    if not rows:
        raise PreflightError(f"empty artifact checksum manifest: {checksum_path}")

    root = attempt.resolve()
    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", row)
        if match is None:
            raise PreflightError(
                f"malformed artifact checksum row {line_number}: {checksum_path}"
            )
        expected, relative = match.groups()
        if relative in seen:
            raise PreflightError(
                f"duplicate artifact checksum path {relative}: {checksum_path}"
            )
        seen.add(relative)
        candidate = (attempt / relative).resolve()
        try:
            candidate.relative_to(root)
        except ValueError as exc:
            raise PreflightError(
                f"artifact checksum path escapes its attempt: {relative}"
            ) from exc
        if not candidate.is_file():
            raise PreflightError(f"completed artifact is missing: {candidate}")
        if sha256_file(candidate) != expected:
            raise PreflightError(f"completed artifact changed: {candidate}")


class SignalForwarder:
    """Forward termination to the complete MPI launcher process group."""

    def __init__(self, process: subprocess.Popen[bytes]):
        self.process = process
        self.received: int | None = None
        self.previous: Dict[int, Any] = {}

    def __enter__(self) -> "SignalForwarder":
        for signum in (signal.SIGINT, signal.SIGTERM):
            self.previous[signum] = signal.getsignal(signum)
            signal.signal(signum, self._handle)
        return self

    def _handle(self, signum: int, frame: object) -> None:
        self.received = signum
        if self.process.poll() is None:
            try:
                os.killpg(self.process.pid, signum)
            except ProcessLookupError:
                pass

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        for signum, previous in self.previous.items():
            signal.signal(signum, previous)


def run_attempt(
    archive: Path,
    attempt: Path,
    status: Dict[str, Any],
    manifest: Dict[str, Any],
    policy: Dict[str, Any],
) -> Dict[str, Any]:
    paths = verify_frozen(archive, manifest)
    seed = int(status["seed"])
    mode = str(status["mode"])
    command = [
        str(manifest["mpi_launcher"]),
        "-n",
        str(policy["mpi_ranks"]),
        str(paths["executable"]),
        "--random-seed",
        str(seed),
        "--chemistry-residual-sampling",
        mode.lower(),
    ]
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(policy["omp_threads"])
    environment["LC_ALL"] = "C"
    environment["LANG"] = "C"

    # No AMPS-created file may predate launch.  status.json is the sole file
    # intentionally present in a prepared attempt.
    unexpected = [entry for entry in attempt.iterdir() if entry.name != "status.json"]
    if unexpected:
        raise PreflightError(f"prepared attempt directory is dirty: {unexpected[0]}")

    run_log = attempt / "run.log"
    error_log = attempt / "error.log"
    started_monotonic = time.monotonic()
    status.update(
        {
            "state": "launching",
            "started_utc": utc_now(),
            "command": command,
            "cwd": str(attempt.resolve()),
            "environment": {"OMP_NUM_THREADS": "1", "LC_ALL": "C", "LANG": "C"},
            "executable_sha256": manifest["frozen_executable_sha256"],
            "input_sha256": manifest["frozen_input_sha256"],
            "pic_global_sha256": manifest["frozen_pic_global_sha256"],
        }
    )
    atomic_write_json(attempt / "status.json", status)
    update_manifest_attempt(archive, manifest, status)

    with run_log.open("wb") as stdout, error_log.open("wb") as stderr:
        try:
            process = subprocess.Popen(
                command,
                cwd=attempt,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
        except OSError as exc:
            status.update(
                {
                    "state": "failed",
                    "completed_utc": utc_now(),
                    "return_code": None,
                    "error": f"cannot launch MPI process: {exc}",
                    "wall_seconds": time.monotonic() - started_monotonic,
                }
            )
            atomic_write_json(attempt / "status.json", status)
            update_manifest_attempt(archive, manifest, status)
            raise OperationalError(status["error"]) from exc

        status.update(
            {
                "state": "running",
                "child_pid": process.pid,
                "child_start_ticks": process_start_ticks(process.pid),
            }
        )
        atomic_write_json(attempt / "status.json", status)
        update_manifest_attempt(archive, manifest, status)

        with SignalForwarder(process) as forwarder:
            return_code = process.wait()

    elapsed = time.monotonic() - started_monotonic
    if forwarder.received is not None:
        status.update(
            {
                "state": "interrupted",
                "completed_utc": utc_now(),
                "return_code": return_code,
                "signal": signal.Signals(forwarder.received).name,
                "wall_seconds": elapsed,
            }
        )
        atomic_write_json(attempt / "status.json", status)
        update_manifest_attempt(archive, manifest, status)
        raise RunInterrupted(
            f"attempt interrupted by {signal.Signals(forwarder.received).name}"
        )
    if return_code != 0:
        status.update(
            {
                "state": "failed",
                "completed_utc": utc_now(),
                "return_code": return_code,
                "wall_seconds": elapsed,
                "error": f"AMPS/MPI returned {return_code}",
            }
        )
        atomic_write_json(attempt / "status.json", status)
        update_manifest_attempt(archive, manifest, status)
        raise OperationalError(status["error"])

    try:
        artifacts, evidence = validate_operational_artifacts(
            attempt, seed, mode, policy
        )
        checksum_path, checksum_hash = write_artifact_checksums(attempt, artifacts)
        verify_frozen(archive, manifest)
    except (GateError, PreflightError, OperationalError, OSError, ValueError) as exc:
        status.update(
            {
                "state": "invalid_artifacts",
                "completed_utc": utc_now(),
                "return_code": return_code,
                "wall_seconds": elapsed,
                "error": str(exc),
            }
        )
        atomic_write_json(attempt / "status.json", status)
        update_manifest_attempt(archive, manifest, status)
        if isinstance(exc, PreflightError):
            raise
        raise OperationalError(str(exc)) from exc

    status.update(
        {
            "state": "completed",
            "completed_utc": utc_now(),
            "return_code": 0,
            "wall_seconds": elapsed,
            "artifact_manifest": str(checksum_path.relative_to(attempt)),
            "artifact_manifest_sha256": checksum_hash,
            "operational_evidence": evidence,
        }
    )
    atomic_write_json(attempt / "status.json", status)
    update_manifest_attempt(archive, manifest, status)
    return {
        "mode": mode,
        "seed": seed,
        "status": "completed",
        "return_code": 0,
        "mpi_ranks": policy["mpi_ranks"],
        "omp_threads": policy["omp_threads"],
        "attempt_path": str(attempt.relative_to(archive)),
        "command": command,
        "started_utc": status["started_utc"],
        "completed_utc": status["completed_utc"],
        "wall_seconds": elapsed,
        "executable_sha256": manifest["frozen_executable_sha256"],
        "input_sha256": manifest["frozen_input_sha256"],
        "pic_global_sha256": manifest["frozen_pic_global_sha256"],
        "artifact_manifest_sha256": checksum_hash,
        **evidence,
    }


def recover_attempt_records(
    archive: Path,
    manifest: Dict[str, Any],
    attempts: List[Dict[str, Any]],
) -> None:
    """Reconcile status files written just before a crash with the manifest."""
    known = {item.get("attempt_path") for item in manifest["attempts"]}
    for status in attempts:
        if status["attempt_path"] not in known:
            update_manifest_attempt(archive, manifest, status)
            known.add(status["attempt_path"])


def recover_completed_runs(
    archive: Path,
    manifest: Dict[str, Any],
    policy: Dict[str, Any],
    attempts: Sequence[Dict[str, Any]],
) -> None:
    """Promote a fully validated case if a crash preceded manifest promotion.

    ``run_attempt`` deliberately publishes the completed attempt status before
    adding the compact run record to ``RUN_MANIFEST.json``.  A power loss in
    that narrow interval must not waste a many-hour science run.  Promotion is
    allowed only when the immutable checksum manifest and every frozen hash
    still match; ambiguous duplicate completions are rejected.
    """
    recorded = {(item["mode"], item["seed"]) for item in manifest["runs"]}
    candidates: Dict[Tuple[str, int], List[Dict[str, Any]]] = {}
    for status in attempts:
        key = (status["mode"], status["seed"])
        if status.get("state") == "completed" and key not in recorded:
            candidates.setdefault(key, []).append(status)
    for key, values in candidates.items():
        if len(values) != 1:
            raise PreflightError(
                f"multiple unrecorded completed attempts exist for {key}"
            )
        status = values[0]
        attempt = archive / status["attempt_path"]
        checksum = attempt / str(status.get("artifact_manifest", ""))
        if (
            status.get("return_code") != 0
            or status.get("executable_sha256")
            != manifest["frozen_executable_sha256"]
            or status.get("input_sha256") != manifest["frozen_input_sha256"]
            or status.get("pic_global_sha256")
            != manifest["frozen_pic_global_sha256"]
            or not checksum.is_file()
            or sha256_file(checksum) != status.get("artifact_manifest_sha256")
            or not isinstance(status.get("operational_evidence"), dict)
        ):
            raise PreflightError(
                f"unrecorded completed attempt failed recovery validation: {attempt}"
            )
        verify_artifact_checksums(attempt, checksum)
        run = {
            "mode": status["mode"],
            "seed": status["seed"],
            "status": "completed",
            "return_code": 0,
            "mpi_ranks": policy["mpi_ranks"],
            "omp_threads": policy["omp_threads"],
            "attempt_path": status["attempt_path"],
            "command": status["command"],
            "started_utc": status["started_utc"],
            "completed_utc": status["completed_utc"],
            "wall_seconds": status["wall_seconds"],
            "executable_sha256": manifest["frozen_executable_sha256"],
            "input_sha256": manifest["frozen_input_sha256"],
            "pic_global_sha256": manifest["frozen_pic_global_sha256"],
            "artifact_manifest_sha256": status["artifact_manifest_sha256"],
            **status["operational_evidence"],
        }
        add_completed_run(archive, manifest, policy, run)


def execute_remaining(
    archive: Path, policy: Dict[str, Any], manifest: Dict[str, Any]
) -> None:
    attempts = scan_attempts(archive, policy)
    validate_manifest_records(archive, manifest, policy, attempts)
    recover_attempt_records(archive, manifest, attempts)
    recover_completed_runs(archive, manifest, policy, attempts)
    complete = {(item["mode"], item["seed"]) for item in manifest["runs"]}

    for item in expected_run_order(policy):
        key = (item["mode"], item["seed"])
        if key in complete:
            continue
        # Re-scan immediately before allocation so an unexpected file created
        # while the archive lock is held is treated as corruption, not reused.
        attempts = scan_attempts(archive, policy)
        attempt_path, status = create_attempt(
            archive, item["seed"], item["mode"], attempts
        )
        update_manifest_attempt(archive, manifest, status)
        print(
            f"START Stage-3 science case mode={item['mode']} seed={item['seed']} "
            f"attempt={attempt_path.relative_to(archive)}",
            flush=True,
        )
        completed = run_attempt(
            archive, attempt_path, status, manifest, policy
        )
        add_completed_run(archive, manifest, policy, completed)
        complete.add(key)
        print(
            f"PASS Stage-3 operational case mode={item['mode']} seed={item['seed']} "
            f"wall={completed['wall_seconds']:.1f}s",
            flush=True,
        )

    if len(manifest["runs"]) != 8:
        raise OperationalError("scheduler ended without eight completed cases")
    print(f"PASS Stage-3 operational matrix archive={archive}")


def plan_payload(
    archive: Path,
    policy_path: Path,
    executable: Path,
    input_path: Path,
    pic_global: Path,
    mpiexec: str,
) -> Dict[str, Any]:
    policy = load_policy(policy_path)
    hashes = validate_source_artifacts(executable, input_path, pic_global, policy)
    launcher = shutil.which(mpiexec)
    if launcher is None:
        raise PreflightError(f"MPI launcher is not on PATH: {mpiexec}")
    return {
        "archive": str(archive.resolve()),
        "policy": str(policy_path.resolve()),
        "policy_sha256": sha256_json(policy),
        "frozen_source_paths": {
            "executable": str(executable.resolve()),
            "input": str(input_path.resolve()),
            "pic_global": str(pic_global.resolve()),
        },
        "hashes": hashes,
        "mpi_launcher": str(Path(launcher).resolve()),
        "mpi_ranks": policy["mpi_ranks"],
        "omp_threads": policy["omp_threads"],
        "environment": {"OMP_NUM_THREADS": "1", "LC_ALL": "C", "LANG": "C"},
        "run_order": expected_run_order(policy),
        "restart_policy": "new-attempt-from-iteration-zero",
    }


def print_status(archive: Path, as_json: bool) -> None:
    policy, manifest = load_archive(archive)
    attempts = scan_attempts(archive, policy, reject_live=False)
    validate_manifest_records(archive, manifest, policy, attempts)
    completed = {(item["mode"], item["seed"]) for item in manifest["runs"]}
    payload = {
        "archive": str(archive.resolve()),
        "manifest_status": manifest.get("status"),
        "completed_cases": len(completed),
        "total_cases": 8,
        "attempt_count": len(attempts),
        "remaining": [
            item
            for item in expected_run_order(policy)
            if (item["mode"], item["seed"]) not in completed
        ],
        "attempt_states": [
            {
                "mode": item["mode"],
                "seed": item["seed"],
                "attempt_number": item["attempt_number"],
                "state": item["state"],
                "attempt_path": item["attempt_path"],
            }
            for item in attempts
        ],
    }
    if as_json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return
    print(
        f"Stage-3 science runner: {payload['completed_cases']}/8 complete, "
        f"attempts={payload['attempt_count']} status={payload['manifest_status']}"
    )
    for item in payload["remaining"]:
        print(f"  pending mode={item['mode']} seed={item['seed']}")
    for item in payload["attempt_states"]:
        print(
            f"  attempt mode={item['mode']} seed={item['seed']} "
            f"number={item['attempt_number']} state={item['state']}"
        )


def add_source_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--executable", type=Path, default=REPO_ROOT / "amps",
        help="configured executable to freeze (default: repository root amps)",
    )
    parser.add_argument(
        "--input", type=Path, default=REPO_ROOT / "input" / "halley1d.input",
        help="canonical input snapshot",
    )
    parser.add_argument(
        "--pic-global", type=Path, default=REPO_ROOT / "build" / "pic" / "picGlobal.dfn",
        help="generated compile-configuration snapshot",
    )
    parser.add_argument(
        "--mpiexec", default="mpiexec", help="MPI launcher name or path"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    default_policy = THIS_DIR / "stage3_science_policy.json"

    plan = subparsers.add_parser("plan", help="validate and print the immutable run plan")
    plan.add_argument("archive", type=Path)
    plan.add_argument("--policy", type=Path, default=default_policy)
    add_source_arguments(plan)

    run = subparsers.add_parser("run", help="create a new archive and execute all cases")
    run.add_argument("archive", type=Path)
    run.add_argument("--policy", type=Path, default=default_policy)
    add_source_arguments(run)

    resume = subparsers.add_parser(
        "resume", help="continue missing cases; interrupted cases restart at iteration zero"
    )
    resume.add_argument("archive", type=Path)

    status = subparsers.add_parser("status", help="report archive progress without launching")
    status.add_argument("archive", type=Path)
    status.add_argument("--json", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        archive = arguments.archive.resolve()
        if arguments.command == "plan":
            payload = plan_payload(
                archive,
                arguments.policy.resolve(),
                arguments.executable.resolve(),
                arguments.input.resolve(),
                arguments.pic_global.resolve(),
                arguments.mpiexec,
            )
            print(json.dumps(payload, indent=2, sort_keys=True))
            return EXIT_SUCCESS
        if arguments.command == "status":
            print_status(archive, arguments.json)
            return EXIT_SUCCESS
        if arguments.command == "run":
            # Validate policy and sources before creating the archive.  A typo
            # must not leave a directory that a later invocation could mistake
            # for an initialized evidence root.
            policy_path = arguments.policy.resolve()
            executable = arguments.executable.resolve()
            input_path = arguments.input.resolve()
            pic_global = arguments.pic_global.resolve()
            policy = load_policy(policy_path)
            validate_source_artifacts(executable, input_path, pic_global, policy)
            with ArchiveLock(archive):
                policy, manifest = initialize_archive(
                    archive,
                    policy_path,
                    executable,
                    input_path,
                    pic_global,
                    arguments.mpiexec,
                )
                execute_remaining(archive, policy, manifest)
            return EXIT_SUCCESS
        if arguments.command == "resume":
            with ArchiveLock(archive):
                policy, manifest = load_archive(archive)
                execute_remaining(archive, policy, manifest)
            return EXIT_SUCCESS
        raise AssertionError("unreachable command")
    except GateError as exc:
        print(f"run_stage3_science_gate.py: invalid policy/input: {exc}", file=sys.stderr)
        return EXIT_INVALID_POLICY
    except PreflightError as exc:
        print(f"run_stage3_science_gate.py: preflight error: {exc}", file=sys.stderr)
        return EXIT_PREFLIGHT
    except OperationalError as exc:
        print(f"run_stage3_science_gate.py: operational failure: {exc}", file=sys.stderr)
        return EXIT_OPERATIONAL
    except RunInterrupted as exc:
        print(f"run_stage3_science_gate.py: interrupted: {exc}", file=sys.stderr)
        return EXIT_INTERRUPTED
    except KeyboardInterrupt:
        # The normal active-run handler converts SIGINT into RunInterrupted.
        # This fallback covers interruption during preflight or between cases.
        print("run_stage3_science_gate.py: interrupted", file=sys.stderr)
        return EXIT_INTERRUPTED


if __name__ == "__main__":
    raise SystemExit(main())
