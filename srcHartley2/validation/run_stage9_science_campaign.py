#!/usr/bin/env python3
"""Plan, run, inspect, and resume the Stage-9 Halley science campaign.

The campaign consists of one ANALOG and one ENHANCED run for each of four
documented seeds.  All eight cases use one frozen Stage-9-enabled executable,
input, core configuration, and generated chemistry policy set.  Cases execute
sequentially at exactly 12 MPI ranks and always start from iteration zero.

``plan`` creates the archive and freezes provenance but never launches AMPS.
``run`` executes a pristine planned archive. ``resume`` preserves completed
cases and starts a new attempt for each interrupted/failed case; partial PT
trees are never reused. ``status`` is read-only.  Before every ENHANCED long
case the runner executes the frozen production-NTC transported-daughter
fixture in an isolated subdirectory.  The summary builder authenticates that
log, so no observer is inserted into the general collision hot path.

The template policy deliberately contains ``DECISION_REQUIRED`` science
thresholds.  Plan refuses it until those values are explicitly resolved.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
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
from typing import Any, Callable, Dict, Iterable, List, Mapping, Sequence, TextIO, Tuple


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
sys.path.insert(0, str(THIS_DIR))

from stage9_run_summary import (  # noqa: E402
    SummaryError,
    atomic_write_json,
    build_summary,
    sha256_file,
)


STATE_SCHEMA = "amps.comet.stage9-campaign-state.v1"
POLICY_SCHEMA = "amps.comet.stage9-campaign-policy.v1"
MANIFEST_SCHEMA = "amps.comet.stage9-science-manifest.v1"
ATTEMPT_SCHEMA = "amps.comet.stage9-attempt.v1"
HEADER_NAMES = (
    "CometDaughterPacketPolicy.generated.h",
    "CometChannelImportancePolicy.generated.h",
    "CometTotalReactionForcingPolicy.generated.h",
    "CometChemistry.generated.h",
)
SEEDS = (104729, 1299709, 15485863, 32452843)
ARMS = ("analog", "enhanced")


class PolicyError(RuntimeError):
    pass


class PreflightError(RuntimeError):
    pass


class OperationalError(RuntimeError):
    pass


class Interrupted(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_object(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PreflightError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PreflightError(f"JSON root is not an object: {path}")
    return value


def sha256_json(value: object) -> str:
    encoded = (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()
    return hashlib.sha256(encoded).hexdigest()


def atomic_copy(source: Path, destination: Path, mode: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    try:
        with source.open("rb") as input_stream, os.fdopen(descriptor, "wb") as output_stream:
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


def atomic_write_text(path: Path, text: str, mode: int = 0o644) -> None:
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


def exact_keys(value: Mapping[str, Any], expected: Iterable[str], label: str) -> None:
    actual, wanted = set(value), set(expected)
    if actual != wanted:
        raise PolicyError(
            f"{label} keys differ: missing={sorted(wanted-actual)} extra={sorted(actual-wanted)}"
        )


def finite(value: object, label: str, *, minimum: float = 0.0) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        if value == "DECISION_REQUIRED":
            raise PolicyError(f"{label} remains DECISION_REQUIRED")
        raise PolicyError(f"{label} is not numeric")
    result = float(value)
    if not math.isfinite(result) or result < minimum:
        raise PolicyError(f"{label} is outside its valid range")
    return result


def validate_policy(value: Dict[str, Any]) -> Dict[str, Any]:
    required = {
        "schema",
        "seeds",
        "mpi_ranks",
        "omp_threads",
        "sampling_length",
        "burn_in_iterations",
        "maximum_iterations",
        "minimum_steady_windows",
        "chemistry_validation_events",
        "gate_policy",
    }
    exact_keys(value, required, "campaign policy")
    if value["schema"] != POLICY_SCHEMA:
        raise PolicyError("unknown Stage-9 campaign policy schema")
    if value["seeds"] != list(SEEDS):
        raise PolicyError("Stage-9 campaign requires the documented four seeds")
    expected_int = {
        "mpi_ranks": 12,
        "omp_threads": 1,
        "sampling_length": 10000,
        "burn_in_iterations": 10000,
        "maximum_iterations": 40000,
        "minimum_steady_windows": 3,
    }
    for key, expected in expected_int.items():
        if value[key] != expected:
            raise PolicyError(f"campaign {key} must equal {expected}")
    if not isinstance(value["chemistry_validation_events"], int) or value[
        "chemistry_validation_events"
    ] < 1000:
        raise PolicyError("chemistry_validation_events must be an integer >=1000")

    gate = value["gate_policy"]
    if not isinstance(gate, dict):
        raise PolicyError("gate_policy is not an object")
    gate_keys = {
        "minimum_seeds",
        "mpi_ranks",
        "omp_threads",
        "burn_in_iterations",
        "minimum_steady_windows",
        "maximum_buffer_utilization",
        "confidence_multiplier",
        "minimum_forced_daughter_collisions",
        "flux_closure",
        "primary_species",
        "daughter_species",
    }
    exact_keys(gate, gate_keys, "gate_policy")
    mirrored = {
        "minimum_seeds": 4,
        "mpi_ranks": 12,
        "omp_threads": 1,
        "burn_in_iterations": 10000,
        "minimum_steady_windows": 3,
    }
    for key, expected in mirrored.items():
        if gate[key] != expected:
            raise PolicyError(f"gate_policy {key} must equal {expected}")
    maximum_buffer = finite(gate["maximum_buffer_utilization"], "maximum_buffer_utilization")
    if not 0.0 < maximum_buffer < 1.0:
        raise PolicyError("maximum_buffer_utilization must be in (0,1)")
    if finite(gate["confidence_multiplier"], "confidence_multiplier") <= 0.0:
        raise PolicyError("confidence_multiplier must be positive")
    if not isinstance(gate["minimum_forced_daughter_collisions"], int) or gate[
        "minimum_forced_daughter_collisions"
    ] < 1:
        raise PolicyError("minimum_forced_daughter_collisions must be positive")
    flux = gate["flux_closure"]
    if not isinstance(flux, dict):
        raise PolicyError("flux_closure policy is not an object")
    exact_keys(flux, {"method", "relative_tolerance"}, "flux_closure policy")
    method = flux["method"]
    if method == "DECISION_REQUIRED":
        raise PolicyError("flux_closure method remains DECISION_REQUIRED")
    if method == "independent_seed_ci_contains_zero":
        if flux["relative_tolerance"] is not None:
            raise PolicyError(
                "flux_closure relative_tolerance must be null for independent_seed_ci_contains_zero"
            )
    elif method == "relative_tolerance":
        tolerance = finite(
            flux["relative_tolerance"], "flux_closure relative_tolerance"
        )
        if tolerance <= 0.0:
            raise PolicyError("flux_closure relative_tolerance must be positive")
    else:
        raise PolicyError(f"unknown flux_closure method: {method}")
    primary, daughters = gate["primary_species"], gate["daughter_species"]
    if not isinstance(primary, dict) or not primary or not isinstance(daughters, dict) or not daughters:
        raise PolicyError("gate_policy requires primary and daughter species")
    for species, settings in primary.items():
        if not isinstance(settings, dict):
            raise PolicyError(f"primary {species} policy is not an object")
        exact_keys(settings, {"relative_tolerance", "combined_se_multiplier"}, f"primary {species}")
        finite(settings["relative_tolerance"], f"{species} relative_tolerance")
        finite(settings["combined_se_multiplier"], f"{species} combined_se_multiplier")
    daughter_keys = {
        "relative_tolerance",
        "combined_se_multiplier",
        "minimum_ess_fraction",
        "maximum_zero_cell_fraction",
        "minimum_variance_cost_gain",
        "minimum_gain_bin_fraction",
    }
    for species, settings in daughters.items():
        if not isinstance(settings, dict):
            raise PolicyError(f"daughter {species} policy is not an object")
        exact_keys(settings, daughter_keys, f"daughter {species}")
        for key in daughter_keys:
            finite(settings[key], f"{species} {key}")
        for key in ("minimum_ess_fraction", "maximum_zero_cell_fraction", "minimum_gain_bin_fraction"):
            if float(settings[key]) > 1.0:
                raise PolicyError(f"{species} {key} exceeds one")
        if float(settings["minimum_variance_cost_gain"]) <= 1.0:
            raise PolicyError(f"{species} minimum_variance_cost_gain must exceed one")
    return value


def load_policy(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PolicyError(f"cannot read campaign policy {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PolicyError("campaign policy root is not an object")
    return validate_policy(value)


def parse_input(path: Path) -> Dict[str, str]:
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*(?:\([^)]*\))?)\s*=\s*(.*?)\s*$")
    result: Dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.split("!", 1)[0])
        if match:
            key, raw = match.groups()
            result[key] = raw.strip()
    return result


def validate_input(path: Path, policy: Dict[str, Any]) -> None:
    records = parse_input(path)
    required = {
        "Stage9ScienceDiagnostics": "on",
        "RadialSamplingLength": str(policy["sampling_length"]),
        "RadialMaxIterations": str(policy["maximum_iterations"]),
        "RadialConvergenceWindows": str(policy["minimum_steady_windows"]),
        "CometDiagnostics": "on",
        "ChemistryMode": "on",
        "ChemistryResidualSampling": "bernoulli",
        "ParticleWeightCorrectionWindowMode": "on",
        "ChemistryDaughterPacketMode": "on",
        "ChemistryChannelImportanceMode": "on",
        "ChemistryTotalReactionForcingMode": "on",
    }
    for key, expected in required.items():
        if records.get(key, "").lower() != expected.lower():
            raise PolicyError(f"Stage-9 input {key} must equal {expected}")


def execution_plan(policy: Dict[str, Any]) -> List[Dict[str, Any]]:
    result: List[Dict[str, Any]] = []
    for index, seed in enumerate(policy["seeds"]):
        order = ARMS if index % 2 == 0 else tuple(reversed(ARMS))
        result.extend({"seed": seed, "estimator_mode": arm} for arm in order)
    return result


def frozen_paths(archive: Path) -> Dict[str, Any]:
    root = archive / "frozen"
    return {
        "executable": root / "amps.stage9-science",
        "input": root / "halley1d.stage9.input",
        "pic_global": root / "picGlobal.stage9.dfn",
        "policy": root / "stage9_policy.json",
        "headers": {name: root / name for name in HEADER_NAMES},
    }


def frozen_record(path: Path, archive: Path) -> Dict[str, str]:
    return {"path": str(path.relative_to(archive)), "sha256": sha256_file(path)}


def verify_frozen(archive: Path, state: Dict[str, Any]) -> Dict[str, Any]:
    paths = frozen_paths(archive)
    records = state.get("frozen")
    if not isinstance(records, dict):
        raise PreflightError("campaign state lacks frozen provenance")
    simple = {"executable": paths["executable"], "input": paths["input"], "pic_global": paths["pic_global"]}
    for key, path in simple.items():
        record = records.get(key)
        if not path.is_file() or not isinstance(record, dict) or sha256_file(path) != record.get("sha256"):
            raise PreflightError(f"frozen {key} changed or is missing")
    for name, path in paths["headers"].items():
        record = records.get("headers", {}).get(name)
        if not path.is_file() or not isinstance(record, dict) or sha256_file(path) != record.get("sha256"):
            raise PreflightError(f"frozen generated header changed or is missing: {name}")
    policy = load_policy(paths["policy"])
    if sha256_json(policy) != state.get("policy_sha256"):
        raise PreflightError("frozen campaign policy differs from state")
    if not os.access(paths["executable"], os.X_OK):
        raise PreflightError("frozen executable is not executable")
    return paths


class ArchiveLock:
    def __init__(self, archive: Path):
        self.path = archive / ".stage9-runner.lock"
        self.stream: TextIO | None = None

    def __enter__(self) -> "ArchiveLock":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.stream = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise PreflightError("another Stage-9 campaign process owns the archive") from exc
        self.stream.seek(0)
        self.stream.truncate()
        self.stream.write(f"pid={os.getpid()} acquired={utc_now()}\n")
        self.stream.flush()
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if self.stream is not None:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_UN)
            self.stream.close()


def initialize_archive(
    archive: Path,
    policy_path: Path,
    executable: Path,
    input_path: Path,
    pic_global: Path,
    header_dir: Path,
    mpiexec: str,
) -> Dict[str, Any]:
    policy = load_policy(policy_path)
    validate_input(input_path, policy)
    sources = {"executable": executable, "input": input_path, "pic_global": pic_global}
    sources["headers"] = {name: header_dir / name for name in HEADER_NAMES}  # type: ignore[assignment]
    for key, path in (("executable", executable), ("input", input_path), ("pic_global", pic_global)):
        if not path.is_file():
            raise PreflightError(f"source {key} is missing: {path}")
    if not os.access(executable, os.X_OK):
        raise PreflightError(f"source executable is not executable: {executable}")
    for name, path in sources["headers"].items():  # type: ignore[union-attr]
        if not path.is_file():
            raise PreflightError(f"generated header is missing: {path}")
    launcher = shutil.which(mpiexec)
    if launcher is None:
        raise PreflightError(f"MPI launcher is not on PATH: {mpiexec}")
    archive.mkdir(parents=True, exist_ok=True)
    unexpected = [entry for entry in archive.iterdir() if entry.name != ".stage9-runner.lock"]
    if unexpected:
        raise PreflightError(f"plan requires a new empty archive: {unexpected[0]}")

    paths = frozen_paths(archive)
    atomic_copy(executable, paths["executable"], 0o555)
    atomic_copy(input_path, paths["input"], 0o444)
    atomic_copy(pic_global, paths["pic_global"], 0o444)
    atomic_write_json(paths["policy"], policy)
    os.chmod(paths["policy"], 0o444)
    for name, source in sources["headers"].items():  # type: ignore[union-attr]
        atomic_copy(source, paths["headers"][name], 0o444)

    frozen = {
        "executable": frozen_record(paths["executable"], archive),
        "input": frozen_record(paths["input"], archive),
        "pic_global": frozen_record(paths["pic_global"], archive),
        "headers": {
            name: frozen_record(path, archive) for name, path in paths["headers"].items()
        },
    }
    state = {
        "schema": STATE_SCHEMA,
        "created_utc": utc_now(),
        "updated_utc": utc_now(),
        "status": "planned",
        "policy_sha256": sha256_json(policy),
        "frozen": frozen,
        "mpi_launcher": str(Path(launcher).resolve()),
        "host": platform.node(),
        "platform": platform.platform(),
        "execution_plan": execution_plan(policy),
        "attempts": [],
        "completed": [],
    }
    atomic_write_json(archive / "CAMPAIGN_STATE.json", state)
    verify_frozen(archive, state)
    return state


def load_archive(archive: Path) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    state = load_object(archive / "CAMPAIGN_STATE.json")
    if state.get("schema") != STATE_SCHEMA:
        raise PreflightError("unknown Stage-9 campaign-state schema")
    paths = verify_frozen(archive, state)
    policy = load_policy(paths["policy"])
    if state.get("execution_plan") != execution_plan(policy):
        raise PreflightError("campaign execution plan differs from policy")
    return policy, state


def process_start_ticks(pid: int) -> int | None:
    try:
        return int(Path(f"/proc/{pid}/stat").read_text(encoding="utf-8").split()[21])
    except (OSError, ValueError, IndexError):
        return None


def process_alive(pid: object, ticks: object) -> bool:
    return isinstance(pid, int) and isinstance(ticks, int) and process_start_ticks(pid) == ticks


def update_state(archive: Path, state: Dict[str, Any]) -> None:
    state["updated_utc"] = utc_now()
    atomic_write_json(archive / "CAMPAIGN_STATE.json", state)


def attempt_root(archive: Path, seed: int, arm: str, number: int) -> Path:
    return archive / "runs" / f"seed-{seed:010d}" / arm / f"attempt-{number:04d}"


def reconcile(archive: Path, policy: Dict[str, Any], state: Dict[str, Any]) -> None:
    """Recover atomic attempt statuses written just before a reboot."""
    known = {item["attempt_path"] for item in state["attempts"]}
    for status_path in sorted((archive / "runs").glob("seed-*/[a-z]*/attempt-*/status.json")) if (archive / "runs").exists() else ():
        status = load_object(status_path)
        relative = str(status_path.parent.relative_to(archive))
        if status.get("schema") != ATTEMPT_SCHEMA or status.get("attempt_path") != relative:
            raise PreflightError(f"inconsistent attempt identity: {status_path}")
        if relative not in known:
            state["attempts"].append(status)
            known.add(relative)
    for index, recorded in enumerate(list(state["attempts"])):
        status_path = archive / recorded["attempt_path"] / "status.json"
        status = load_object(status_path)
        if status.get("state") in {"launching", "running", "collision_preflight"}:
            if process_alive(status.get("child_pid"), status.get("child_start_ticks")):
                raise PreflightError(f"attempt still has a live launcher: {status_path.parent}")
            status.update(state="interrupted", completed_utc=utc_now(), error="launcher not live during resume")
            atomic_write_json(status_path, status)
        state["attempts"][index] = status
    completed_keys = {(row["method"], row["seed"]) for row in state["completed"]}
    for status in state["attempts"]:
        if status.get("state") != "completed":
            continue
        method = "ANALOG" if status["estimator_mode"] == "analog" else "FORCED"
        key = (method, status["seed"])
        summary = archive / status["attempt_path"] / "stage9_summary.json"
        if not summary.is_file() or sha256_file(summary) != status.get("summary_sha256"):
            raise PreflightError(f"completed summary changed: {summary}")
        checksum = archive / status["attempt_path"] / str(status.get("artifact_manifest", ""))
        if not checksum.is_file() or sha256_file(checksum) != status.get(
            "artifact_manifest_sha256"
        ):
            raise PreflightError(f"completed artifact manifest changed: {checksum}")
        verify_checksums(archive / status["attempt_path"], checksum)
        if key not in completed_keys:
            state["completed"].append(
                {
                    "method": method,
                    "seed": status["seed"],
                    "summary": str(summary.relative_to(archive)),
                    "sha256": status["summary_sha256"],
                }
            )
            completed_keys.add(key)
    update_state(archive, state)


class SignalForwarder:
    def __init__(self, process: subprocess.Popen[bytes]):
        self.process = process
        self.received: int | None = None
        self.previous: Dict[int, Any] = {}

    def __enter__(self) -> "SignalForwarder":
        for signum in (signal.SIGINT, signal.SIGTERM):
            self.previous[signum] = signal.getsignal(signum)
            signal.signal(signum, self.handle)
        return self

    def handle(self, signum: int, frame: object) -> None:
        self.received = signum
        if self.process.poll() is None:
            try:
                os.killpg(self.process.pid, signum)
            except ProcessLookupError:
                pass

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        for signum, handler in self.previous.items():
            signal.signal(signum, handler)


def launch(
    command: List[str],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
    environment: Dict[str, str],
    on_start: Callable[[int, int | None], None] | None = None,
) -> Tuple[int, float, int | None]:
    started = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            env=environment,
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        if on_start is not None:
            on_start(process.pid, process_start_ticks(process.pid))
        with SignalForwarder(process) as forwarder:
            return_code = process.wait()
    return return_code, time.monotonic() - started, forwarder.received


def write_checksums(attempt: Path, artifacts: Iterable[Path]) -> Tuple[Path, str]:
    rows = [f"{sha256_file(path)}  {path.relative_to(attempt)}" for path in sorted(set(artifacts))]
    destination = attempt / "artifacts.sha256"
    atomic_write_text(destination, "\n".join(rows) + "\n", 0o444)
    return destination, sha256_file(destination)


def verify_checksums(attempt: Path, manifest: Path) -> None:
    root = attempt.resolve()
    seen = set()
    for line_number, line in enumerate(
        manifest.read_text(encoding="utf-8").splitlines(), start=1
    ):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            raise PreflightError(f"malformed artifact checksum row {line_number}: {manifest}")
        expected, relative = match.groups()
        if relative in seen:
            raise PreflightError(f"duplicate artifact checksum path: {relative}")
        seen.add(relative)
        path = (attempt / relative).resolve()
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise PreflightError(f"artifact checksum escapes attempt: {relative}") from exc
        if not path.is_file() or sha256_file(path) != expected:
            raise PreflightError(f"completed artifact changed: {path}")
    if not seen:
        raise PreflightError(f"empty artifact checksum manifest: {manifest}")


def publish_manifest(archive: Path, policy: Dict[str, Any], state: Dict[str, Any]) -> None:
    if len(state["completed"]) != 8:
        return
    expected = {(arm, seed) for arm in ("ANALOG", "FORCED") for seed in policy["seeds"]}
    actual = {(row["method"], row["seed"]) for row in state["completed"]}
    if actual != expected:
        raise PreflightError("completed matrix is not balanced")
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "policy": policy["gate_policy"],
        "runs": sorted(state["completed"], key=lambda row: (row["seed"], row["method"])),
    }
    atomic_write_json(archive / "RUN_MANIFEST.json", manifest)
    state["status"] = "completed"
    update_state(archive, state)


def create_attempt(archive: Path, state: Dict[str, Any], seed: int, arm: str) -> Tuple[Path, Dict[str, Any]]:
    prior = [
        item
        for item in state["attempts"]
        if item.get("seed") == seed and item.get("estimator_mode") == arm
    ]
    root = attempt_root(archive, seed, arm, len(prior) + 1)
    if root.exists():
        raise PreflightError(f"new attempt path already exists: {root}")
    root.mkdir(parents=True)
    status = {
        "schema": ATTEMPT_SCHEMA,
        "attempt_path": str(root.relative_to(archive)),
        "attempt_number": len(prior) + 1,
        "seed": seed,
        "estimator_mode": arm,
        "state": "prepared",
        "created_utc": utc_now(),
        "restart_used": False,
        "fresh_start_verified": False,
    }
    atomic_write_json(root / "status.json", status)
    state["attempts"].append(status)
    update_state(archive, state)
    return root, status


def execute_attempt(
    archive: Path,
    policy: Dict[str, Any],
    state: Dict[str, Any],
    attempt: Path,
    status: Dict[str, Any],
) -> None:
    paths = verify_frozen(archive, state)
    launcher = state["mpi_launcher"]
    ranks = policy["mpi_ranks"]
    seed = status["seed"]
    arm = status["estimator_mode"]
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", LC_ALL="C", LANG="C")
    base = [launcher, "-n", str(ranks), str(paths["executable"]), "--random-seed", str(seed), "--chemistry-estimator-mode", arm]

    if arm == "enhanced":
        preflight = attempt / "collision-preflight"
        preflight.mkdir()
        command = base + [
            "--validate-transported-daughter-collision",
            "--chemistry-validation-events",
            str(policy["chemistry_validation_events"]),
        ]
        status.update(state="collision_preflight", collision_command=command)
        atomic_write_json(attempt / "status.json", status)
        def collision_started(pid: int, ticks: int | None) -> None:
            status.update(child_pid=pid, child_start_ticks=ticks)
            atomic_write_json(attempt / "status.json", status)
        return_code, elapsed, received = launch(
            command,
            preflight,
            preflight / "run.log",
            preflight / "error.log",
            environment,
            collision_started,
        )
        if received is not None:
            status.update(state="interrupted", return_code=return_code, collision_wall_seconds=elapsed)
            atomic_write_json(attempt / "status.json", status)
            raise Interrupted("collision preflight was interrupted")
        if return_code != 0:
            status.update(state="failed", return_code=return_code, collision_wall_seconds=elapsed)
            atomic_write_json(attempt / "status.json", status)
            raise OperationalError("transported-daughter collision preflight failed")
        status["collision_wall_seconds"] = elapsed

    if (attempt / "PT").exists():
        raise PreflightError("long-run attempt already contains PT before launch")
    command = base
    status.update(
        state="running",
        command=command,
        started_utc=utc_now(),
        fresh_start_verified=True,
        restart_used=False,
    )
    atomic_write_json(attempt / "status.json", status)
    def long_run_started(pid: int, ticks: int | None) -> None:
        status.update(child_pid=pid, child_start_ticks=ticks)
        atomic_write_json(attempt / "status.json", status)
    return_code, elapsed, received = launch(
        command,
        attempt,
        attempt / "run.log",
        attempt / "error.log",
        environment,
        long_run_started,
    )
    if received is not None:
        status.update(state="interrupted", return_code=return_code, wall_seconds=elapsed)
        atomic_write_json(attempt / "status.json", status)
        raise Interrupted("long science run was interrupted")
    if return_code != 0:
        status.update(state="failed", return_code=return_code, wall_seconds=elapsed)
        atomic_write_json(attempt / "status.json", status)
        raise OperationalError(f"AMPS/MPI returned {return_code}")
    combined = (attempt / "run.log").read_text(encoding="utf-8", errors="replace") + (attempt / "error.log").read_text(encoding="utf-8", errors="replace")
    if "Reading particle restart file" in combined:
        raise OperationalError("science attempt read a particle restart")
    status.update(state="summarizing", return_code=0, wall_seconds=elapsed, completed_utc=utc_now())
    atomic_write_json(attempt / "status.json", status)
    try:
        summary, artifacts = build_summary(archive, attempt)
    except SummaryError as exc:
        status.update(state="invalid_artifacts", error=str(exc))
        atomic_write_json(attempt / "status.json", status)
        raise OperationalError(str(exc)) from exc
    summary_path = attempt / "stage9_summary.json"
    atomic_write_json(summary_path, summary)
    artifacts.append(summary_path)
    checksum_path, checksum_hash = write_checksums(attempt, artifacts)
    status.update(
        state="completed",
        summary_sha256=sha256_file(summary_path),
        artifact_manifest=str(checksum_path.relative_to(attempt)),
        artifact_manifest_sha256=checksum_hash,
    )
    atomic_write_json(attempt / "status.json", status)
    for index, item in enumerate(state["attempts"]):
        if item["attempt_path"] == status["attempt_path"]:
            state["attempts"][index] = status
            break
    method = "ANALOG" if arm == "analog" else "FORCED"
    state["completed"].append(
        {
            "method": method,
            "seed": seed,
            "summary": str(summary_path.relative_to(archive)),
            "sha256": status["summary_sha256"],
        }
    )
    state["status"] = "running"
    update_state(archive, state)


def execute_remaining(archive: Path, policy: Dict[str, Any], state: Dict[str, Any]) -> None:
    reconcile(archive, policy, state)
    completed = {(item["method"], item["seed"]) for item in state["completed"]}
    for item in state["execution_plan"]:
        method = "ANALOG" if item["estimator_mode"] == "analog" else "FORCED"
        key = (method, item["seed"])
        if key in completed:
            continue
        attempt, status = create_attempt(
            archive, state, item["seed"], item["estimator_mode"]
        )
        print(
            f"START Stage-9 method={method} seed={item['seed']} attempt={attempt.relative_to(archive)}",
            flush=True,
        )
        execute_attempt(archive, policy, state, attempt, status)
        completed.add(key)
        print(f"PASS Stage-9 operational method={method} seed={item['seed']}", flush=True)
    publish_manifest(archive, policy, state)
    if len(state["completed"]) != 8:
        raise OperationalError("campaign ended without eight completed cases")
    print(f"PASS Stage-9 operational campaign archive={archive}")


def status_payload(archive: Path) -> Dict[str, Any]:
    policy, state = load_archive(archive)
    # Status is read-only: overlay the newest atomically published per-attempt
    # records without reconciling or changing campaign state.
    attempts = []
    for recorded in state["attempts"]:
        path = archive / recorded["attempt_path"] / "status.json"
        attempts.append(load_object(path) if path.is_file() else recorded)
    completed = {(item["method"], item["seed"]) for item in state["completed"]}
    return {
        "archive": str(archive),
        "status": state["status"],
        "completed_cases": len(completed),
        "total_cases": 8,
        "attempt_count": len(attempts),
        "remaining": [
            item
            for item in state["execution_plan"]
            if (("ANALOG" if item["estimator_mode"] == "analog" else "FORCED"), item["seed"])
            not in completed
        ],
        "attempts": [
            {
                key: item.get(key)
                for key in ("seed", "estimator_mode", "attempt_number", "state", "attempt_path")
            }
            for item in attempts
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan = subparsers.add_parser("plan", help="freeze a new campaign without launching AMPS")
    plan.add_argument("archive", type=Path)
    plan.add_argument("--policy", type=Path, default=THIS_DIR / "stage9_science_policy.template.json")
    plan.add_argument("--executable", type=Path, default=REPO_ROOT / "amps")
    plan.add_argument("--input", type=Path, default=REPO_ROOT / "input" / "halley1d.input")
    plan.add_argument("--pic-global", type=Path, default=REPO_ROOT / "build" / "pic" / "picGlobal.dfn")
    plan.add_argument("--header-dir", type=Path, default=REPO_ROOT / "build" / "main")
    plan.add_argument("--mpiexec", default="mpiexec")
    for name in ("run", "resume"):
        command = subparsers.add_parser(name, help=f"{name} the frozen sequential campaign")
        command.add_argument("archive", type=Path)
    status = subparsers.add_parser("status", help="report progress without launching")
    status.add_argument("archive", type=Path)
    status.add_argument("--json", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    archive = arguments.archive.resolve()
    try:
        if arguments.command == "plan":
            with ArchiveLock(archive):
                state = initialize_archive(
                    archive,
                    arguments.policy.resolve(),
                    arguments.executable.resolve(),
                    arguments.input.resolve(),
                    arguments.pic_global.resolve(),
                    arguments.header_dir.resolve(),
                    arguments.mpiexec,
                )
            print(json.dumps(status_payload(archive), indent=2, sort_keys=True))
            return 0
        if arguments.command == "status":
            payload = status_payload(archive)
            if arguments.json:
                print(json.dumps(payload, indent=2, sort_keys=True))
            else:
                print(
                    f"Stage-9 campaign: {payload['completed_cases']}/8 complete "
                    f"attempts={payload['attempt_count']} status={payload['status']}"
                )
                for item in payload["remaining"]:
                    print(f"  pending mode={item['estimator_mode']} seed={item['seed']}")
            return 0
        with ArchiveLock(archive):
            policy, state = load_archive(archive)
            if arguments.command == "run" and state["attempts"]:
                raise PreflightError("run requires a pristine plan; use resume")
            execute_remaining(archive, policy, state)
        return 0
    except PolicyError as exc:
        print(f"run_stage9_science_campaign.py: policy/input error: {exc}", file=sys.stderr)
        return 2
    except PreflightError as exc:
        print(f"run_stage9_science_campaign.py: preflight error: {exc}", file=sys.stderr)
        return 3
    except (OperationalError, OSError, subprocess.SubprocessError) as exc:
        print(f"run_stage9_science_campaign.py: operational failure: {exc}", file=sys.stderr)
        return 4
    except (Interrupted, KeyboardInterrupt) as exc:
        print(f"run_stage9_science_campaign.py: interrupted: {exc}", file=sys.stderr)
        return 5


if __name__ == "__main__":
    raise SystemExit(main())
