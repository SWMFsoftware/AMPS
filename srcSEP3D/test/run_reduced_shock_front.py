#!/usr/bin/env python3
"""Run the complete reduced shock-front plus ambient verification campaign.

This is an evidence orchestrator, not a physics implementation and not a
replacement for any assertion in the underlying executables.  It runs:

* the dependency-light shared RSH assertion harness;
* the shared-library architecture boundary;
* the portable srcSEP3D RSHAPP factory/adapter tests;
* explicit RSH24--RSH27 native smoke cases on one and four MPI ranks; and
* uninterrupted, same-rank restart, and one-to-four-rank restart comparisons;
* the separate four-rank RSH24--RSH28 low-corona-to-1-AU campaign.

Every child has a fresh log.  Native JSON is parsed as the result authority;
stdout alone can never create a PASS.  The final text and JSON summaries list
every FAIL/ERROR and its complete execution log.  A reduced PASS does not
qualify BG3D-4 or claim a physical downstream CME volume.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import selectors
import shlex
import signal
import subprocess
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


AMPS_ROOT = Path(__file__).resolve().parents[2]
EXPECTED_ROOT = Path("/home/vtenishe/Mars2/AMPS")
STATUSES = ("PASS", "FAIL", "SKIP", "ERROR")
SMOKE_IDS = ("RSH24", "RSH25", "RSH26", "RSH27")
LONG_IDS = SMOKE_IDS + ("RSH28",)
RESTART_IDS = SMOKE_IDS + ("RSH29",)


class RunnerError(RuntimeError):
    """A setup/report failure that cannot be called a scientific FAIL."""


@dataclass
class TestResult:
    scope: str
    test_id: str
    status: str
    message: str
    log: str
    report: str = ""


@dataclass
class Operation:
    phase: str
    command: List[str]
    cwd: str
    log: str
    returncode: int
    elapsed_seconds: float
    error: str = ""


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Run all reduced shock-front shared, portable and native tests.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    result.add_argument("--list", action="store_true",
                        help="list phases and stable test IDs, then exit")
    result.add_argument("--dry-run", action="store_true",
                        help="print resolved commands without executing them")
    result.add_argument("--output-dir", type=Path,
                        help="fresh evidence directory; existing paths are rejected")
    result.add_argument("--amps", type=Path, default=AMPS_ROOT / "amps",
                        help="configured srcSEP3D AMPS executable")
    result.add_argument(
        "--mpi-launch-prefix", default="mpiexec -n {ranks}",
        help="shell-like argv prefix; must contain the literal {ranks}")
    result.add_argument("--skip-native", action="store_true",
                        help="run shared/portable gates and report native cases SKIP")
    result.add_argument("--rebuild-native", action="store_true",
                        help="perform the mandated clean srcSEP3D native rebuild first")
    result.add_argument("--timeout", type=float, default=900.0,
                        help="shared/portable/smoke timeout in seconds")
    result.add_argument("--long-timeout", type=float, default=7200.0,
                        help="four-rank 1-AU campaign timeout in seconds")
    result.add_argument("--progress-interval", type=float, default=30.0,
                        help="heartbeat interval while a child is quiet")
    result.add_argument("--verbose", action="store_true",
                        help="stream complete child output as well as saving logs")
    return result


def list_tests() -> None:
    print("shared-rsh: RSH00--RSH24, RSH28--RSH29, RSH31--RSH32, RSH36--RSH40 (54 assertions)")
    print("architecture: ARCHCSWC01")
    print("portable-srcsep3d: RSHAPP01 RSHAPP02 RSHAPP03 RSHAPP04")
    print("native-smoke-1: " + " ".join(SMOKE_IDS))
    print("native-smoke-4: " + " ".join(SMOKE_IDS))
    print("native-restart-reference-4: " + " ".join(SMOKE_IDS))
    print("native-restart-checkpoint-4: " + " ".join(SMOKE_IDS))
    print("native-restart-resume-4: " + " ".join(RESTART_IDS))
    print("native-restart-checkpoint-1: " + " ".join(SMOKE_IDS))
    print("native-restart-repartition-4: " + " ".join(RESTART_IDS))
    print("native-1au-4: " + " ".join(LONG_IDS))
    print("Expected complete selected-profile total: 93 PASS, 0 FAIL, 0 SKIP, 0 ERROR")


def command_text(command: Sequence[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def stop_process(process: Optional[subprocess.Popen]) -> None:
    if process is None:
        return
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        elif process.poll() is None:
            process.kill()
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass


def execute(phase: str, command: Sequence[str], cwd: Path, log: Path,
            timeout: float, progress_interval: float, verbose: bool) -> Tuple[Operation, str]:
    """Execute one phase while preserving raw combined stdout/stderr.

    The runner prints result-looking lines and periodic heartbeats by default;
    --verbose streams everything.  A new POSIX process group lets a timeout
    terminate an MPI launcher and its local children rather than leaking ranks
    that could contaminate a later phase.
    """
    log.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    print(f"[reduced-runner] START {phase} log={log}", flush=True)
    print("[reduced-runner] RUN " + command_text(command), flush=True)
    process: Optional[subprocess.Popen] = None
    error = ""
    output_chunks: List[str] = []
    selector = selectors.DefaultSelector()
    next_heartbeat = started + progress_interval
    try:
        environment = dict(os.environ, PYTHONUNBUFFERED="1")
        process = subprocess.Popen(
            list(command), cwd=str(cwd), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=False, bufsize=0, env=environment,
            start_new_session=(os.name == "posix"))
        assert process.stdout is not None
        selector.register(process.stdout, selectors.EVENT_READ)
        with log.open("wb") as stream:
            header = ("phase=" + phase + "\ncwd=" + str(cwd) +
                      "\ncommand=" + command_text(command) + "\n\n").encode()
            stream.write(header); stream.flush()
            pending = ""
            while process.poll() is None or selector.get_map():
                now = time.monotonic()
                if now - started >= timeout:
                    raise subprocess.TimeoutExpired(list(command), timeout)
                wait = max(0.0, min(0.25, next_heartbeat-now,
                                    timeout-(now-started)))
                for key, _ in selector.select(wait):
                    chunk = os.read(key.fileobj.fileno(), 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        continue
                    stream.write(chunk); stream.flush()
                    text = chunk.decode("utf-8", errors="replace")
                    output_chunks.append(text)
                    if verbose:
                        print(text, end="", flush=True)
                    else:
                        pending += text
                        lines = pending.split("\n"); pending = lines.pop()
                        for line in lines:
                            if (re.match(r"^(PASS|FAIL|SKIP|ERROR|SUMMARY|\[[A-Z])", line)
                                    or " Summary:" in line):
                                print(line, flush=True)
                now = time.monotonic()
                if now >= next_heartbeat:
                    print(f"[reduced-runner] RUNNING {phase} elapsed={now-started:.1f}s",
                          flush=True)
                    next_heartbeat = now + progress_interval
            returncode = process.wait()
    except (OSError, subprocess.TimeoutExpired) as failure:
        error = str(failure)
        returncode = 2
        stop_process(process)
        if not log.exists():
            log.write_text("launch failure: " + error + "\n", encoding="utf-8")
    finally:
        selector.close()
        if process is not None and process.stdout is not None:
            process.stdout.close()
    elapsed = time.monotonic() - started
    print(f"[reduced-runner] END {phase} exit={returncode} elapsed={elapsed:.1f}s" +
          (" error=" + error if error else ""), flush=True)
    operation = Operation(phase, list(command), str(cwd), str(log), returncode,
                          elapsed, error)
    text = "".join(output_chunks)
    if not text and log.exists():
        text = log.read_text(encoding="utf-8", errors="replace")
    return operation, text


def infrastructure_error(scope: str, message: str, log: Path) -> TestResult:
    return TestResult(scope, "SETUP", "ERROR", message, str(log))


def parse_shared(text: str, operation: Operation) -> List[TestResult]:
    rows: List[TestResult] = []
    occurrences: Dict[str, int] = defaultdict(int)
    for line in text.splitlines():
        match = re.match(r"^(PASS|FAIL|SKIP|ERROR)\s+(RSH[A-Z0-9]+)\s*(.*)$", line)
        if not match:
            continue
        status, base_id, message = match.groups()
        occurrences[base_id] += 1
        test_id = base_id + "#" + str(occurrences[base_id])
        rows.append(TestResult("shared-rsh", test_id, status, message,
                               operation.log))
    summary = re.findall(
        r"^SUMMARY PASS=(\d+) FAIL=(\d+) SKIP=(\d+) ERROR=(\d+)\s*$",
        text, flags=re.MULTILINE)
    if len(summary) != 1:
        raise RunnerError("shared RSH output has no unique SUMMARY line")
    expected = tuple(int(value) for value in summary[0])
    observed = Counter(row.status for row in rows)
    actual = tuple(observed[name] for name in STATUSES)
    if expected != actual:
        raise RunnerError(f"shared RSH summary {expected} disagrees with parsed rows {actual}")
    expected_exit = 2 if observed["ERROR"] else (1 if observed["FAIL"] else 0)
    if operation.returncode != expected_exit:
        raise RunnerError("shared RSH process exit disagrees with reported status")
    return rows


def parse_architecture(text: str, operation: Operation) -> List[TestResult]:
    matches = re.findall(r"^\[(ARCHCSWC01)\]\s+(PASS|FAIL|SKIP|ERROR)\s*(.*)$",
                         text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise RunnerError("architecture output has no unique ARCHCSWC01 result")
    test_id, status, message = matches[0]
    expected_exit = 2 if status == "ERROR" else (1 if status == "FAIL" else 0)
    if operation.returncode != expected_exit:
        raise RunnerError("architecture process exit disagrees with result")
    return [TestResult("architecture", test_id, status, message, operation.log)]


def parse_portable(text: str, operation: Operation) -> List[TestResult]:
    matches = re.findall(r"^\[(RSHAPP\d+)\]\s+(PASS|FAIL|SKIP|ERROR).*?message=\"(.*)\"$",
                         text, flags=re.MULTILINE)
    ids = [row[0] for row in matches]
    expected_ids = ["RSHAPP01", "RSHAPP02", "RSHAPP03", "RSHAPP04"]
    if ids != expected_ids:
        raise RunnerError(f"portable RSHAPP results {ids!r} do not equal {expected_ids!r}")
    rows = [TestResult("portable-srcsep3d", test_id, status, message,
                       operation.log) for test_id, status, message in matches]
    observed = Counter(row.status for row in rows)
    expected_exit = 2 if observed["ERROR"] else (1 if observed["FAIL"] else 0)
    if operation.returncode != expected_exit:
        raise RunnerError("portable RSHAPP process exit disagrees with results")
    return rows


def read_json(path: Path) -> dict:
    def pairs(items: Iterable[Tuple[str, object]]) -> dict:
        result: dict = {}
        for key, value in items:
            if key in result:
                raise RunnerError("duplicate JSON key in native report: " + key)
            result[key] = value
        return result
    def nonfinite(value: str) -> None:
        raise RunnerError("nonfinite native JSON value: " + value)
    def floating(value: str) -> float:
        number = float(value)
        if not math.isfinite(number):
            raise RunnerError("out-of-range native JSON number: " + value)
        return number
    try:
        return json.loads(path.read_text(encoding="utf-8"),
                          object_pairs_hook=pairs, parse_constant=nonfinite,
                          parse_float=floating)
    except (OSError, ValueError) as error:
        raise RunnerError("missing or malformed native JSON: " + str(error)) from error


def parse_native(scope: str, report: Path, operation: Operation,
                 ids: Sequence[str], ranks: int,
                 completed_steps: int) -> List[TestResult]:
    payload = read_json(report)
    if payload.get("schema") != "srcsep-component-tests-v1":
        raise RunnerError("native report has an unsupported schema")
    if (payload.get("mpi_ranks") != ranks or
            payload.get("completed_steps") != completed_steps):
        raise RunnerError("native rank/step receipt differs from the requested campaign")
    providers = payload.get("providers")
    if not isinstance(providers, dict) or providers.get("background") != "runtime-model" or providers.get("shock") != "none":
        raise RunnerError("native report did not use the reduced runtime background with shock.authority=none")
    particle = payload.get("particles")
    if not isinstance(particle, dict) or particle.get("global_count") != 0 or particle.get("global_injected") != 0:
        raise RunnerError("native report is not a zero-particle execution")
    reduced = payload.get("reduced_front")
    if (not isinstance(reduced, dict) or reduced.get("selected") is not True or
            not isinstance(reduced.get("event_identity"), str) or
            not reduced["event_identity"] or
            not isinstance(reduced.get("generation"), int) or
            reduced["generation"] <= 0 or
            reduced.get("generation") != reduced.get("ambient_generation")):
        raise RunnerError("native report has no coherent selected reduced-front epoch")
    source_rows = payload.get("results")
    if not isinstance(source_rows, list):
        raise RunnerError("native report has no results array")
    by_id: Dict[str, dict] = {}
    for row in source_rows:
        if not isinstance(row, dict) or not isinstance(row.get("id"), str):
            raise RunnerError("native report contains a malformed result")
        if row["id"] in by_id:
            raise RunnerError("native report contains duplicate result " + row["id"])
        by_id[row["id"]] = row
    if set(by_id) != set(ids):
        raise RunnerError(f"native result IDs {sorted(by_id)!r} do not equal {sorted(ids)!r}")
    rows: List[TestResult] = []
    for test_id in ids:
        row = by_id[test_id]
        status = row.get("status")
        if status not in STATUSES:
            raise RunnerError("native result has invalid status: " + test_id)
        message = str(row.get("message", ""))
        rows.append(TestResult(scope, test_id, status, message, operation.log,
                               str(report)))
    observed = Counter(row.status for row in rows)
    expected_exit = 2 if observed["ERROR"] else (1 if observed["FAIL"] else 0)
    if operation.returncode != expected_exit:
        raise RunnerError("native process exit disagrees with JSON results")
    return rows


def launcher(prefix: str, ranks: int) -> List[str]:
    if "{ranks}" not in prefix:
        raise RunnerError("--mpi-launch-prefix must contain the literal {ranks}")
    return [part.replace("{ranks}", str(ranks)) for part in shlex.split(prefix)]


def native_command(args: argparse.Namespace, ids: Sequence[str], ranks: int,
                   deck: Path, steps: int, directory: Path,
                   restart: Optional[Path] = None) -> Tuple[List[str], Path]:
    report = directory / "native.json"
    command = launcher(args.mpi_launch_prefix, ranks) + [str(args.amps.resolve())]
    for test_id in ids:
        command.extend(("--test", test_id))
    command.extend((
        "--test-input", str(deck), "--test-steps", str(steps),
        "--expect-mpi-ranks", str(ranks),
        "--output-dir", str(directory / "products"),
        "--test-json", str(report),
        "--artifact-directory", str(directory / "artifacts")))
    if restart is not None:
        command.extend(("--restart", str(restart)))
    return command, report


def write_restart_deck(output: Path) -> Path:
    """Create one evidence-local deck with a checkpoint at ticks 5 and 10.

    Checkpoint cadence participates in the physics fingerprint, so every leg
    uses these exact bytes.  Output paths are relocation-only.  The model
    asset becomes absolute because the generated deck lives under the evidence
    directory rather than beside ``handoff_smoke.event``.
    """
    source = AMPS_ROOT / "srcSEP3D/examples/shock-front/handoff_smoke.in"
    text = source.read_text(encoding="utf-8")
    replacements = (
        ("model_asset = handoff_smoke.event",
         "model_asset = " + str((source.parent / "handoff_smoke.event").resolve())),
        ("checkpoint_cadence_steps = 0", "checkpoint_cadence_steps = 5"),
        ("output_path = test_output/reduced-front/handoff-smoke/restart.chk",
         "output_path = restart.chk"),
    )
    for old, new in replacements:
        if text.count(old) != 1:
            raise RunnerError("restart-deck source has no unique assignment: " + old)
        text = text.replace(old, new)
    marker = ("# Generated by run_reduced_shock_front.py.  This uses the same "
              "front/ambient physics as handoff_smoke.in; only checkpoint "
              "cadence and relocation-only paths differ.\n")
    target = output / "restart-smoke.in"
    target.write_text(marker + text, encoding="utf-8")
    return target


def compare_restart_receipts(reference_path: Path, resumed_path: Path,
                             source_ranks: int, current_ranks: int) -> str:
    """Compare independently executed native states at the same final time.

    Owner hashes are rank-independent XOR+sum reductions over positions and
    complete ambient samples, after the native current/previous slots have
    passed byte readback.  The front fingerprint covers every surface label,
    geometry and shock classification.  Thus equality is stronger than
    comparing two apex radii or two provider-generated expectations.
    """
    reference = read_json(reference_path)
    resumed = read_json(resumed_path)
    failures: List[str] = []

    def same(path: str, left: object, right: object) -> None:
        if left != right:
            failures.append(f"{path}: uninterrupted={left!r} resumed={right!r}")

    same("configuration_fingerprint", reference.get("configuration_fingerprint"),
         resumed.get("configuration_fingerprint"))
    same("completed_steps", reference.get("completed_steps"),
         resumed.get("completed_steps"))
    mesh_a = reference.get("runtime_mesh_background", {})
    mesh_b = resumed.get("runtime_mesh_background", {})
    for key in ("owner_fingerprint_xor", "owner_fingerprint_sum"):
        same("runtime_mesh_background." + key, mesh_a.get(key), mesh_b.get(key))
    for key in ("owner_fields_match", "ghost_fields_match", "provider_matches"):
        if mesh_a.get(key) is not True or mesh_b.get(key) is not True:
            failures.append("runtime_mesh_background." + key + " is not true in both runs")
    if current_ranks > 1:
        if not isinstance(mesh_b.get("ghost_cells_checked"), int) or mesh_b["ghost_cells_checked"] <= 0:
            failures.append("resumed run checked no received MPI ghost cells")
        same("runtime_mesh_background.ghost_cells_checked",
             mesh_a.get("ghost_cells_checked"), mesh_b.get("ghost_cells_checked"))

    front_a = reference.get("reduced_front", {})
    front_b = resumed.get("reduced_front", {})
    for key in ("selected", "event_identity", "generation", "ambient_generation",
                "state_fingerprint", "epoch_s", "phase", "apex_radius_m",
                "apex_speed_m_s", "accepted_area_m2",
                "numerical_failure_area_m2", "geometric_endpoint_reached",
                "apex_shock_accepted", "endpoint_time_s",
                "endpoint_observer_status", "endpoint_observer_shock_accepted"):
        same("reduced_front." + key, front_a.get(key), front_b.get(key))

    same("particles", reference.get("particles"), resumed.get("particles"))
    particles = resumed.get("particles", {})
    if (particles.get("allocation_requested_zero") is not True or
            particles.get("global_count") != 0 or
            particles.get("global_injected") != 0):
        failures.append("resumed native particle/source counts are not exactly zero")

    restart = resumed.get("restart", {})
    if restart.get("configured") is not True:
        failures.append("resumed receipt is not marked as a native restart")
    if restart.get("source_ranks") != source_ranks:
        failures.append("restart source rank count differs from checkpoint creator")
    if resumed.get("mpi_ranks") != current_ranks:
        failures.append("resumed MPI rank count differs from requested rank count")
    if restart.get("input_tick") != 5 or restart.get("input_background_generation") != 6:
        failures.append("restart did not begin from the frozen tick-5/generation-6 boundary")
    same("restart.checkpoint_sequence", reference.get("restart", {}).get("checkpoint_sequence"),
         restart.get("checkpoint_sequence"))

    if failures:
        raise RunnerError("native restart equivalence failed: " + "; ".join(failures))
    mode = "same-rank" if source_ranks == current_ranks else "changed-rank"
    return (mode + " native checkpoint/resume exactly matches the uninterrupted "
            "tick-10 front, shock classification, owner fields, ghosts, epochs, and zero-particle state")


def require_fresh_output(path: Path) -> None:
    if path.exists():
        raise RunnerError("output directory already exists; refusing overwrite: " + str(path))
    path.mkdir(parents=True)


def rebuild_native(args: argparse.Namespace, output: Path,
                   operations: List[Operation]) -> Optional[TestResult]:
    """Apply the mandatory clean native rebuild sequence exactly once.

    The fixed root check and symlink rejection ensure ``rm -rf -- build`` can
    only name the generated Mars2 root build directory.  Any observed build,
    compiler, MPI or AMPS process aborts before removal.
    """
    log = output / "native-rebuild.log"
    if AMPS_ROOT.resolve() != EXPECTED_ROOT.resolve() or Path.cwd().resolve() != AMPS_ROOT.resolve():
        message = ("native rebuild requires cwd and resolved root exactly " +
                   str(EXPECTED_ROOT))
        log.write_text(message + "\n", encoding="utf-8")
        return infrastructure_error("native-rebuild", message, log)
    build = AMPS_ROOT / "build"
    if build.is_symlink():
        message = "refusing to remove symlinked root build directory"
        log.write_text(message + "\n", encoding="utf-8")
        return infrastructure_error("native-rebuild", message, log)
    process_command = ["ps", "-C", "make", "-C", "gmake", "-C", "g++",
                       "-C", "gcc", "-C", "cc1plus", "-C", "mpiexec",
                       "-C", "mpirun", "-C", "amps", "-o",
                       "pid=,ppid=,stat=,etime=,args="]
    operation, text = execute("native-preflight", process_command, AMPS_ROOT,
                              output / "native-preflight.log", args.timeout,
                              args.progress_interval, args.verbose)
    operations.append(operation)
    if operation.returncode not in (0, 1):
        return infrastructure_error("native-rebuild", "process audit failed",
                                    Path(operation.log))
    active = [line for line in text.splitlines() if line.strip() and not line.startswith("phase=") and not line.startswith("cwd=") and not line.startswith("command=")]
    if active:
        return infrastructure_error(
            "native-rebuild", "active build/test process blocks root build removal: " +
            " | ".join(active), Path(operation.log))
    commands = (
        ["rm", "-rf", "--", "build"],
        ["./Config.pl", "-application=sep3d"],
        ["./ampsConfig.pl", "-input", "sep3d.input", "-no-compile"],
        ["make", "-C", "srcSEP3D", "prepare-production"],
        ["make", "-j16", "amps"],
    )
    with log.open("w", encoding="utf-8") as manifest:
        manifest.write("confirmed_root=" + str(AMPS_ROOT) + "\n")
        manifest.write("selected_application=srcSEP3D\n")
        manifest.write("process_audit=" + operation.log + "\n")
    for index, command in enumerate(commands, start=1):
        phase = f"native-rebuild-{index}"
        operation, _ = execute(phase, command, AMPS_ROOT,
                               output / (phase + ".log"), args.long_timeout,
                               args.progress_interval, args.verbose)
        operations.append(operation)
        with log.open("a", encoding="utf-8") as manifest:
            manifest.write(phase + "=" + operation.log +
                           " exit=" + str(operation.returncode) + "\n")
        if operation.returncode != 0:
            return infrastructure_error("native-rebuild",
                                        "native rebuild failed during " + phase,
                                        Path(operation.log))
    return None


def write_summary(output: Path, results: List[TestResult],
                  operations: List[Operation], started: str) -> int:
    counts = Counter(row.status for row in results)
    totals = {name.lower(): counts[name] for name in STATUSES}
    failed = [row for row in results if row.status in ("FAIL", "ERROR")]
    skipped = [row for row in results if row.status == "SKIP"]
    completed = datetime.now(timezone.utc).isoformat()
    payload = {
        "schema": "srcsep3d-reduced-shock-runner-v1",
        "started_utc": started,
        "completed_utc": completed,
        "amps_root": str(AMPS_ROOT),
        "qualification": "reduced-front-plus-ambient-only; BG3D-4 remains unqualified",
        "totals": totals,
        "results": [asdict(row) for row in results],
        "failed_tests": [asdict(row) for row in failed],
        "skipped_tests": [asdict(row) for row in skipped],
        "operations": [asdict(operation) for operation in operations],
    }
    (output / "summary.json").write_text(json.dumps(payload, indent=2) + "\n",
                                          encoding="utf-8")
    lines = [
        "Reduced shock-front plus ambient test summary",
        "evidence_directory=" + str(output),
        "qualification=reduced-front-plus-ambient-only; BG3D-4 remains unqualified",
        ("SUMMARY PASS={pass} FAIL={fail} SKIP={skip} ERROR={error}".format(**totals)),
    ]
    if failed:
        lines.append("FAILED TESTS:")
        for row in failed:
            lines.append(f"  {row.status} {row.scope}/{row.test_id}: {row.message}")
            lines.append("    log=" + row.log)
            if row.report:
                lines.append("    report=" + row.report)
    else:
        lines.append("FAILED TESTS: none")
    if skipped:
        lines.append("SKIPPED TESTS:")
        for row in skipped:
            lines.append(f"  SKIP {row.scope}/{row.test_id}: {row.message}")
            lines.append("    log=" + row.log)
    lines.append("PHASE LOGS:")
    for operation in operations:
        lines.append(f"  {operation.phase}: {operation.log}")
    summary = "\n".join(lines) + "\n"
    (output / "summary.txt").write_text(summary, encoding="utf-8")
    print("\n" + summary, end="", flush=True)
    return 1 if counts["FAIL"] or counts["ERROR"] else 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv)
    if args.list:
        list_tests(); return 0
    if args.timeout <= 0 or args.long_timeout <= 0 or args.progress_interval <= 0:
        parser().error("timeouts and progress interval must be positive")
    if args.rebuild_native and args.skip_native:
        parser().error("--rebuild-native and --skip-native are mutually exclusive")
    output = (args.output_dir or
              AMPS_ROOT / "test_output" / "reduced-front" / "runner" /
              (utc_stamp() + "-" + str(os.getpid()))).expanduser().resolve()
    args.amps = args.amps.expanduser().resolve()
    if args.dry_run:
        print("AMPS_ROOT=" + str(AMPS_ROOT))
        print("output_dir=" + str(output))
        print("shared: make -C src/models/sep_corona_swcme -j16 shock-front-test")
        print("architecture: make -C src/models/sep_corona_swcme check-architecture")
        print("portable-build: make -C srcSEP3D -j16 test/stage1")
        print("portable: cwd=srcSEP3D ./test/stage1 --test-group RSHAPP")
        for name, ids, ranks, deck, steps in (
                ("native-smoke-1", SMOKE_IDS, 1, "handoff_smoke.in", 10),
                ("native-smoke-4", SMOKE_IDS, 4, "handoff_smoke.in", 10),
                ("native-1au-4", LONG_IDS, 4, "corona_to_1au.in", 340)):
            directory = output / name
            command, _ = native_command(
                args, ids, ranks,
                AMPS_ROOT / "srcSEP3D/examples/shock-front" / deck,
                steps, directory)
            print(name + ": " + command_text(command))
        restart_deck = output / "restart-smoke.in"
        restart_cases = (
            ("native-restart-reference-4", SMOKE_IDS, 4, 10, None),
            ("native-restart-checkpoint-4", SMOKE_IDS, 4, 5, None),
            ("native-restart-resume-4", RESTART_IDS, 4, 5,
             output / "native-restart-checkpoint-4/restart.chk"),
            ("native-restart-checkpoint-1", SMOKE_IDS, 1, 5, None),
            ("native-restart-repartition-4", RESTART_IDS, 4, 5,
             output / "native-restart-checkpoint-1/restart.chk"),
        )
        for name, ids, ranks, steps, restart in restart_cases:
            directory = output / name
            command, _ = native_command(args, ids, ranks, restart_deck, steps,
                                         directory, restart)
            print(name + ": cwd=" + str(directory) + " " + command_text(command))
        return 0

    require_fresh_output(output)
    started = datetime.now(timezone.utc).isoformat()
    operations: List[Operation] = []
    results: List[TestResult] = []
    runner_log = output / "runner.log"
    runner_log.write_text("AMPS_ROOT=" + str(AMPS_ROOT) + "\n", encoding="utf-8")

    phases = (
        ("shared-rsh", ["make", "-C", "src/models/sep_corona_swcme", "-j16",
                        "shock-front-test"], AMPS_ROOT, parse_shared),
        ("architecture", ["make", "-C", "src/models/sep_corona_swcme",
                          "check-architecture"], AMPS_ROOT, parse_architecture),
    )
    for phase, command, cwd, reader in phases:
        operation, text = execute(phase, command, cwd, output/(phase + ".log"),
                                  args.timeout, args.progress_interval, args.verbose)
        operations.append(operation)
        try:
            results.extend(reader(text, operation))
        except RunnerError as error:
            results.append(infrastructure_error(phase, str(error), Path(operation.log)))

    build_operation, _ = execute(
        "portable-build", ["make", "-C", "srcSEP3D", "-j16", "test/stage1"],
        AMPS_ROOT, output/"portable-build.log", args.timeout,
        args.progress_interval, args.verbose)
    operations.append(build_operation)
    if build_operation.returncode != 0:
        results.append(infrastructure_error(
            "portable-srcsep3d", "portable stage1 build failed",
            Path(build_operation.log)))
    else:
        operation, text = execute(
            "portable-rshapp", ["./test/stage1", "--test-group", "RSHAPP"],
            AMPS_ROOT/"srcSEP3D", output/"portable-rshapp.log", args.timeout,
            args.progress_interval, args.verbose)
        operations.append(operation)
        try:
            results.extend(parse_portable(text, operation))
        except RunnerError as error:
            results.append(infrastructure_error(
                "portable-srcsep3d", str(error), Path(operation.log)))

    native_cases = (
        ("native-smoke-1", SMOKE_IDS, 1, "handoff_smoke.in", 10, args.timeout),
        ("native-smoke-4", SMOKE_IDS, 4, "handoff_smoke.in", 10, args.timeout),
        ("native-1au-4", LONG_IDS, 4, "corona_to_1au.in", 340, args.long_timeout),
    )
    restart_case_ids = (
        ("native-restart-reference-4", SMOKE_IDS),
        ("native-restart-checkpoint-4", SMOKE_IDS),
        ("native-restart-resume-4", RESTART_IDS),
        ("native-restart-checkpoint-1", SMOKE_IDS),
        ("native-restart-repartition-4", RESTART_IDS),
    )
    if args.skip_native:
        for scope, ids, _, _, _, _ in native_cases:
            for test_id in ids:
                results.append(TestResult(scope, test_id, "SKIP",
                    "native execution explicitly disabled by --skip-native",
                    str(runner_log)))
        for scope, ids in restart_case_ids:
            for test_id in ids:
                results.append(TestResult(scope, test_id, "SKIP",
                    "native restart execution explicitly disabled by --skip-native",
                    str(runner_log)))
    else:
        rebuild_error = rebuild_native(args, output, operations) if args.rebuild_native else None
        if rebuild_error is not None:
            results.append(rebuild_error)
        elif not args.amps.is_file():
            results.append(infrastructure_error(
                "native-amps", "configured AMPS executable is absent: " + str(args.amps),
                runner_log))
        else:
            def run_case(scope: str, ids: Sequence[str], ranks: int,
                         deck: Path, requested_steps: int,
                         completed_steps: int, timeout: float,
                         cwd: Path = AMPS_ROOT,
                         restart: Optional[Path] = None
                         ) -> Tuple[List[TestResult], Path]:
                directory = output / scope
                directory.mkdir(parents=True)
                (directory / "artifacts").mkdir()
                command, report = native_command(
                    args, ids, ranks, deck, requested_steps, directory, restart)
                operation, _ = execute(scope, command, cwd,
                                       directory/"execution.log", timeout,
                                       args.progress_interval, args.verbose)
                operations.append(operation)
                try:
                    rows = parse_native(scope, report, operation, ids, ranks,
                                        completed_steps)
                    results.extend(rows)
                    return rows, report
                except RunnerError as error:
                    results.append(infrastructure_error(scope, str(error),
                                                        Path(operation.log)))
                    return [], report

            example_root = AMPS_ROOT / "srcSEP3D/examples/shock-front"
            for scope, ids, ranks, deck_name, steps, timeout in native_cases[:2]:
                run_case(scope, ids, ranks, example_root/deck_name, steps,
                         steps, timeout)

            # Native restart qualification uses one generated deck for every
            # leg, because checkpoint cadence is physics-fingerprinted.  Each
            # process has its own cwd so the relative restart.chk is retained
            # beside that leg's raw log and cannot overwrite another receipt.
            restart_deck = write_restart_deck(output)
            reference_rows, reference_report = run_case(
                "native-restart-reference-4", SMOKE_IDS, 4, restart_deck,
                10, 10, args.timeout, output/"native-restart-reference-4")
            checkpoint4_rows, checkpoint4_report = run_case(
                "native-restart-checkpoint-4", SMOKE_IDS, 4, restart_deck,
                5, 5, args.timeout, output/"native-restart-checkpoint-4")
            checkpoint4 = output/"native-restart-checkpoint-4/restart.chk"
            if checkpoint4_rows and not checkpoint4.is_file():
                results.append(infrastructure_error(
                    "native-restart-checkpoint-4",
                    "tick-5 native checkpoint was not written", checkpoint4_report))
            resume4_rows: List[TestResult] = []
            resume4_report = output/"native-restart-resume-4/native.json"
            if reference_rows and checkpoint4.is_file():
                resume4_rows, resume4_report = run_case(
                    "native-restart-resume-4", RESTART_IDS, 4, restart_deck,
                    5, 10, args.timeout, output/"native-restart-resume-4",
                    checkpoint4)
            if resume4_rows:
                try:
                    message = compare_restart_receipts(
                        reference_report, resume4_report, 4, 4)
                    for row in resume4_rows:
                        if row.test_id == "RSH29": row.message = message
                except RunnerError as error:
                    for row in resume4_rows:
                        if row.test_id == "RSH29":
                            row.status = "FAIL"; row.message = str(error)

            checkpoint1_rows, checkpoint1_report = run_case(
                "native-restart-checkpoint-1", SMOKE_IDS, 1, restart_deck,
                5, 5, args.timeout, output/"native-restart-checkpoint-1")
            checkpoint1 = output/"native-restart-checkpoint-1/restart.chk"
            if checkpoint1_rows and not checkpoint1.is_file():
                results.append(infrastructure_error(
                    "native-restart-checkpoint-1",
                    "one-rank tick-5 native checkpoint was not written",
                    checkpoint1_report))
            repartition_rows: List[TestResult] = []
            repartition_report = output/"native-restart-repartition-4/native.json"
            if reference_rows and checkpoint1.is_file():
                repartition_rows, repartition_report = run_case(
                    "native-restart-repartition-4", RESTART_IDS, 4,
                    restart_deck, 5, 10, args.timeout,
                    output/"native-restart-repartition-4", checkpoint1)
            if repartition_rows:
                try:
                    message = compare_restart_receipts(
                        reference_report, repartition_report, 1, 4)
                    for row in repartition_rows:
                        if row.test_id == "RSH29": row.message = message
                except RunnerError as error:
                    for row in repartition_rows:
                        if row.test_id == "RSH29":
                            row.status = "FAIL"; row.message = str(error)

            # Preserve the separate long campaign and its negative
            # non-forward 1-AU outcome after restart qualification.
            scope, ids, ranks, deck_name, steps, timeout = native_cases[2]
            run_case(scope, ids, ranks, example_root/deck_name, steps, steps,
                     timeout)

    return write_summary(output, results, operations, started)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RunnerError as error:
        print("ERROR: " + str(error), file=sys.stderr)
        raise SystemExit(2)
