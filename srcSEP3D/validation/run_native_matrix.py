#!/usr/bin/env python3
"""Run a fail-closed srcSEP3D native/MPI qualification matrix.

This command is an evidence collector, not an MPI simulator. Every matrix
entry is executed by the user-supplied launcher against a configured, linked
AMPS executable. Before starting an expensive matrix, the runner asks that
executable for ``--list-tests`` and refuses to continue unless every required
native case is advertised. This prevents an ordinary production binary, an
old binary, or the AMPS-independent ``test/stage1`` executable from being
mistaken for a linked qualification host.

The exact executable, input deck, profile, commands, thread counts, reports,
and stdout logs are hashed or recorded. A missing/malformed native report is
an ERROR; it can never inherit PASS from a process return code or a stale file.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import datetime as _datetime
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
PROFILE_FILE = ROOT / "validation" / "native_profiles.json"
NATIVE_REPORT_SCHEMA = "srcsep-component-tests-v1"
MATRIX_REPORT_SCHEMA = "srcsep3d-native-matrix-v2"


class MatrixError(RuntimeError):
    """A structural qualification error, distinct from a failed test case."""


@dataclass(frozen=True)
class CommandResult:
    """Lossless result of one launcher invocation.

    Timeout and operating-system launch errors are converted to structured
    results so a partially completed matrix can still publish useful evidence.
    The discovery preflight promotes either condition to ``MatrixError``.
    """

    returncode: int
    stdout: str
    elapsed_seconds: float
    execution_error: Optional[str] = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise MatrixError(f"cannot read {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise MatrixError(f"invalid JSON in {path}: {error}") from error


def _write_text_atomic(path: Path, text: str) -> None:
    """Publish one text artifact without exposing a partially written file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    _write_text_atomic(path, json.dumps(payload, indent=2) + "\n")


def _positive_integer_list(value: Any, field: str) -> List[int]:
    if not isinstance(value, list) or not value:
        raise MatrixError(f"native profile {field} must be a nonempty array")
    result: List[int] = []
    for item in value:
        # ``bool`` is an ``int`` subclass in Python and must be rejected here.
        if isinstance(item, bool) or not isinstance(item, int) or item < 1:
            raise MatrixError(
                f"native profile {field} contains a non-positive integer")
        if item in result:
            raise MatrixError(f"native profile {field} contains duplicate {item}")
        result.append(item)
    return result


def _load_profile(name: str) -> Dict[str, Any]:
    payload = _read_json(PROFILE_FILE)
    if not isinstance(payload, dict) or payload.get("schema") != \
            "srcsep3d-native-profiles-v1":
        raise MatrixError("native profile registry has an unsupported schema")
    profiles = payload.get("profiles")
    if not isinstance(profiles, list):
        raise MatrixError("native profile registry has no profiles array")
    matches = [item for item in profiles
               if isinstance(item, dict) and item.get("name") == name]
    if len(matches) != 1:
        raise MatrixError(f"native profile '{name}' is absent or duplicated")

    source = matches[0]
    required = source.get("required_cases")
    if not isinstance(required, list) or not required:
        raise MatrixError("native profile required_cases must be nonempty")
    cases: List[str] = []
    for item in required:
        case_id = str(item).strip().upper()
        if not case_id or case_id in cases:
            raise MatrixError(
                "native profile required_cases contains an empty or duplicate ID")
        cases.append(case_id)
    return {
        "name": name,
        "ranks": _positive_integer_list(source.get("ranks"), "ranks"),
        "threads": _positive_integer_list(source.get("threads"), "threads"),
        "required_cases": cases,
    }


def _launcher_arguments(template: str, ranks: int) -> List[str]:
    """Expand a launcher as argv, never through a shell.

    Requiring ``{ranks}`` prevents a nominal multi-rank profile from silently
    running every row with the same rank count. Literal braces can be escaped
    as ``{{`` and ``}}`` using normal ``str.format`` syntax.
    """

    if "{ranks}" not in template:
        raise MatrixError("--launcher must contain the literal {ranks} field")
    try:
        expanded = template.format(ranks=ranks)
    except (KeyError, IndexError, ValueError) as error:
        raise MatrixError(f"invalid --launcher format: {error}") from error
    try:
        arguments = shlex.split(expanded)
    except ValueError as error:
        raise MatrixError(f"cannot parse --launcher as argv: {error}") from error
    if not arguments:
        raise MatrixError("--launcher expands to an empty command")
    return arguments


def _execute(command: Sequence[str], environment: Mapping[str, str],
             timeout: float) -> CommandResult:
    started = time.monotonic()
    try:
        completed = subprocess.run(
            list(command), cwd=str(ROOT), env=dict(environment), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            timeout=timeout)
        return CommandResult(completed.returncode, completed.stdout,
                             time.monotonic() - started)
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        message = f"command timed out after {timeout:g} seconds"
        return CommandResult(124, output, time.monotonic() - started, message)
    except OSError as error:
        return CommandResult(127, "", time.monotonic() - started,
                             f"cannot execute command: {error}")


def _advertised_cases(executable: Path, launcher: str,
                      timeout: float) -> tuple[List[str], Dict[str, Any]]:
    """Discover linked callbacks with one rank before running the matrix."""

    command = [*_launcher_arguments(launcher, 1), str(executable),
               "--list-tests"]
    result = _execute(command, os.environ, timeout)
    if result.execution_error is not None or result.returncode != 0:
        detail = result.execution_error or result.stdout[-2000:]
        raise MatrixError(
            "linked executable rejected --list-tests; this usually means it "
            f"is not a native qualification host: {detail}")
    advertised = sorted({
        line.split("|", 1)[0].strip().upper()
        for line in result.stdout.splitlines()
        if "|" in line and line.split("|", 1)[0].strip()
    })
    return advertised, {
        "command": command,
        "returncode": result.returncode,
        "elapsed_seconds": result.elapsed_seconds,
        "stdout": result.stdout,
    }


def _native_result(report: Path, case_id: str,
                   process: CommandResult) -> tuple[str, str, List[Any], List[Any]]:
    """Validate one native report and reconcile it with the process status."""

    if process.execution_error is not None:
        return "ERROR", process.execution_error, [], []
    if not report.is_file():
        return ("ERROR",
                f"native JSON report is absent (process exit {process.returncode})",
                [], [])
    try:
        payload = _read_json(report)
    except MatrixError as error:
        return "ERROR", str(error), [], []
    if not isinstance(payload, dict) or payload.get("schema") != \
            NATIVE_REPORT_SCHEMA:
        return "ERROR", "native JSON report has an unsupported schema", [], []
    rows = payload.get("results")
    if not isinstance(rows, list):
        return "ERROR", "native JSON report contains no results array", [], []
    matches = [row for row in rows if isinstance(row, dict) and
               str(row.get("id", "")).upper() == case_id]
    if len(matches) != 1:
        return ("ERROR",
                f"native JSON report contains {len(matches)} records for {case_id}",
                [], [])
    row = matches[0]
    status = str(row.get("status", "ERROR")).upper()
    message = str(row.get("message", ""))
    metrics = row.get("metrics", [])
    artifacts = row.get("artifacts", [])
    if status not in ("PASS", "FAIL", "SKIP", "ERROR"):
        return "ERROR", f"native report has invalid status '{status}'", [], []
    if not isinstance(metrics, list) or not isinstance(artifacts, list):
        return "ERROR", "native report metrics/artifacts must be arrays", [], []

    # The shared test harness has a documented 0/1/2 exit contract. Checking
    # it here prevents a report copied from a different invocation from being
    # accepted merely because its case ID happens to match.
    expected_exit = 2 if status == "ERROR" else (1 if status == "FAIL" else 0)
    if process.returncode != expected_exit:
        return ("ERROR",
                f"native status {status} requires process exit {expected_exit}, "
                f"observed {process.returncode}", metrics, artifacts)
    return status, message, metrics, artifacts


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a real linked srcSEP3D AMPS rank/thread matrix")
    parser.add_argument("--amps", type=Path, required=True,
                        help="configured linked native-test executable")
    parser.add_argument("--test-input", type=Path, required=True,
                        help="complete immutable srcSEP3D input deck")
    parser.add_argument("--profile", choices=("small", "medium", "production"),
                        default="small")
    parser.add_argument(
        "--launcher", default="mpiexec -n {ranks}",
        help="direct argv launcher template containing {ranks}")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=3600.0,
                        help="timeout for each discovery or case command")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.timeout <= 0.0:
        raise MatrixError("--timeout must be positive")

    executable = args.amps.expanduser().resolve()
    test_input = args.test_input.expanduser().resolve()
    output_root = args.output_dir.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise MatrixError(
            f"linked AMPS executable is absent or not executable: {executable}")
    if not test_input.is_file():
        raise MatrixError(f"native test input is absent or not a file: {test_input}")

    profile = _load_profile(args.profile)
    # Parse every declared rank now so a malformed launcher cannot fail halfway
    # through an otherwise expensive qualification campaign.
    for ranks in profile["ranks"]:
        _launcher_arguments(args.launcher, ranks)

    advertised, discovery = _advertised_cases(
        executable, args.launcher, args.timeout)
    missing = [case for case in profile["required_cases"]
               if case not in advertised]
    if missing:
        raise MatrixError(
            "linked executable does not advertise required native case(s): " +
            ", ".join(missing))

    output_root.mkdir(parents=True, exist_ok=True)
    _write_text_atomic(output_root / "discovery.stdout.log",
                       str(discovery.pop("stdout")))
    records: List[Dict[str, Any]] = []
    for ranks in profile["ranks"]:
        for threads in profile["threads"]:
            for case_id in profile["required_cases"]:
                artifact_dir = output_root / f"r{ranks}-t{threads}" / case_id
                artifact_dir.mkdir(parents=True, exist_ok=True)
                report = artifact_dir / "native.json"
                # Never consume a report left by an earlier interrupted run.
                # Only this exact per-case evidence file is replaced; other
                # native artifacts remain available to the callback.
                report.unlink(missing_ok=True)
                command = [
                    *_launcher_arguments(args.launcher, ranks),
                    str(executable),
                    "--test", case_id,
                    "--test-input", str(test_input),
                    "--test-json", str(report),
                    "--artifact-directory", str(artifact_dir),
                ]
                environment = dict(os.environ)
                environment["OMP_NUM_THREADS"] = str(threads)
                process = _execute(command, environment, args.timeout)
                stdout_log = artifact_dir / "native.stdout.log"
                _write_text_atomic(stdout_log, process.stdout)
                status, message, metrics, artifacts = _native_result(
                    report, case_id, process)
                records.append({
                    "case": case_id,
                    "ranks": ranks,
                    "threads": threads,
                    "status": status,
                    "message": message,
                    "returncode": process.returncode,
                    "elapsed_seconds": process.elapsed_seconds,
                    "command": command,
                    "native_report": str(report),
                    "native_report_sha256": (
                        _sha256(report) if report.is_file() else None),
                    "stdout_log": str(stdout_log),
                    "stdout_sha256": _sha256(stdout_log),
                    "stdout_tail": process.stdout[-2000:],
                    "metrics": metrics,
                    "artifacts": artifacts,
                })
                print(f"[{case_id} r={ranks} t={threads}] {status} {message}",
                      flush=True)

    totals = {status: sum(row["status"] == status for row in records)
              for status in ("PASS", "FAIL", "SKIP", "ERROR")}
    summary: Dict[str, Any] = {
        "schema": MATRIX_REPORT_SCHEMA,
        "generated_utc": _datetime.datetime.now(
            _datetime.timezone.utc).isoformat(),
        "profile": profile,
        "profile_registry": str(PROFILE_FILE),
        "profile_registry_sha256": _sha256(PROFILE_FILE),
        "executable": str(executable),
        "executable_sha256": _sha256(executable),
        "test_input": str(test_input),
        "test_input_sha256": _sha256(test_input),
        "launcher": args.launcher,
        "discovery": discovery,
        "advertised_cases": advertised,
        "totals": totals,
        "records": records,
    }
    _write_json_atomic(output_root / "native-matrix.json", summary)

    # A SKIP is permitted by the component harness but cannot close a required
    # qualification cell. Structural ERROR retains exit 2; any other non-PASS
    # matrix outcome is incomplete and returns exit 1.
    if totals["ERROR"]:
        return 2
    if totals["FAIL"] or totals["SKIP"]:
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MatrixError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
