"""Utilities shared by individual AMPS core-test runners."""

from __future__ import annotations

import argparse
import json
import os
import queue
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence


STATUSES = {"PASS", "FAIL", "SKIP", "ERROR"}
EXIT_STATUS = {"PASS": 0, "SKIP": 0, "FAIL": 1, "ERROR": 2}


class RunnerError(RuntimeError):
    """Harness/configuration failure, distinct from a tested assertion failure."""

    pass


def load_context(description: str) -> Dict[str, Any]:
    """Read the global runner's versioned, resolved execution receipt."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--context", required=True, type=Path,
                        help="resolved context JSON produced by the global runner")
    arguments = parser.parse_args()
    try:
        with arguments.context.open(encoding="utf-8") as stream:
            context = json.load(stream)
    except (OSError, ValueError) as error:
        raise RunnerError("cannot read runner context {}: {}".format(arguments.context, error))
    if not isinstance(context, dict) or context.get("schema_version") != 1:
        raise RunnerError("runner context must be a schema_version 1 object")
    for key in ("test", "repo_root", "work_dir", "artifact_dir", "result_file",
                "paths", "tools", "execution"):
        if key not in context:
            raise RunnerError("runner context is missing '{}'".format(key))
    return context


def command(context: Mapping[str, Any], name: str) -> List[str]:
    """Return a defensive copy of a configured tool argv vector."""

    try:
        value = context["tools"][name]
    except KeyError:
        raise RunnerError("tool '{}' is not configured".format(name))
    if not isinstance(value, list) or not value or not all(isinstance(item, str) for item in value):
        raise RunnerError("tool '{}' must resolve to a nonempty argv list".format(name))
    return list(value)


def run_command(
    argv: Sequence[str],
    cwd: Path,
    timeout: Optional[int] = None,
    env: Optional[Mapping[str, str]] = None,
) -> subprocess.CompletedProcess:
    """Run without a shell and relay output while retaining a complete copy.

    A reader thread prevents ``readline`` from defeating the timeout when a
    child is quiet.  Lines are flushed immediately so the global runner can
    announce a native case before that case completes, while the returned
    ``CompletedProcess.stdout`` remains the exact merged transcript used by
    registry parsers and artifacts.  The argument vector avoids quoting
    differences and command injection when paths, launcher wrappers, or
    compiler options vary between machines.
    """

    merged = dict(os.environ)
    if env:
        merged.update({str(key): str(value) for key, value in env.items()})
    print("$ {}".format(" ".join(_quote(item) for item in argv)), flush=True)
    started = time.monotonic()
    # When invoked by the global runner, this adapter already owns a dedicated
    # POSIX process group. Nested make/mpiexec commands must inherit it so the
    # global timeout can terminate the entire tree. A directly invoked adapter
    # has no such parent contract and therefore creates/owns a group itself.
    globally_managed_group = (
        os.name == "posix" and bool(os.environ.get("AMPS_TEST_CONTEXT")))
    process = subprocess.Popen(
        list(argv), cwd=str(cwd), env=merged, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1,
        start_new_session=(os.name == "posix" and not globally_managed_group))
    lines: List[str] = []
    events: "queue.Queue[object]" = queue.Queue()
    completed_marker = object()

    def read_output() -> None:
        try:
            assert process.stdout is not None
            for line in process.stdout:
                events.put(line)
        finally:
            events.put(completed_marker)

    reader = threading.Thread(target=read_output, name="amps-test-output-reader")
    reader.daemon = True
    reader.start()
    deadline = None if timeout is None else started + timeout
    timed_out = False
    try:
        while True:
            wait = 0.1
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0.0:
                    timed_out = True
                    break
                wait = min(wait, remaining)
            try:
                item = events.get(timeout=wait)
            except queue.Empty:
                continue
            if item is completed_marker:
                break
            line = str(item)
            lines.append(line)
            print(line, end="", flush=True)
    except BaseException:
        # A direct adapter owns the complete child process group. Under the
        # global runner the adapter itself is the group owner, so terminate the
        # direct make/mpiexec child and let the global exception boundary reap
        # any remaining descendants with the adapter group. Never leave a
        # native lock/process alive after Ctrl-C or a Python control failure.
        if os.name == "posix" and not globally_managed_group:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
        reader.join(timeout=1.0)
        raise

    if timed_out:
        # A directly invoked adapter owns its group and can terminate the full
        # make/mpiexec tree. Under the global runner, kill the direct child and
        # let the global group owner perform final descendant cleanup when the
        # adapter returns its timeout error.
        if os.name == "posix" and not globally_managed_group:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
        reader.join(timeout=1.0)
        while True:
            try:
                item = events.get_nowait()
            except queue.Empty:
                break
            if item is completed_marker:
                continue
            line = str(item)
            lines.append(line)
            print(line, end="", flush=True)
        raise RunnerError("command timed out after {} seconds".format(timeout))

    try:
        if deadline is None:
            return_code = process.wait()
        else:
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                raise subprocess.TimeoutExpired(list(argv), timeout)
            return_code = process.wait(timeout=remaining)
    except subprocess.TimeoutExpired:
        if os.name == "posix" and not globally_managed_group:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
        reader.join(timeout=1.0)
        raise RunnerError("command timed out after {} seconds".format(timeout))
    except BaseException:
        if os.name == "posix" and not globally_managed_group:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        process.wait()
        reader.join(timeout=1.0)
        raise
    reader.join(timeout=1.0)
    result = subprocess.CompletedProcess(
        list(argv), return_code, stdout="".join(lines), stderr=None)
    print("[command] exit={} seconds={:.3f}".format(
        result.returncode, time.monotonic() - started), flush=True)
    return result


def write_result(
    context: Mapping[str, Any],
    status: str,
    message: str,
    cases: Optional[List[Dict[str, Any]]] = None,
    evidence: Optional[List[str]] = None,
    metrics: Optional[Dict[str, Any]] = None,
) -> int:
    """Atomically publish a schema-1 result and return its mandated exit code."""

    status = status.upper()
    if status not in STATUSES:
        raise RunnerError("invalid result status '{}'".format(status))
    result = {
        "schema_version": 1,
        "id": context["test"]["id"],
        "status": status,
        "message": " ".join(str(message).split()),
        "cases": list(cases or []),
        "evidence": list(evidence or []),
        "metrics": dict(metrics or {}),
    }
    target = Path(str(context["result_file"]))
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(str(temporary), str(target))
    print("[{}] {} {}".format(result["id"], status, result["message"]), flush=True)
    return EXIT_STATUS[status]


def error_result(context: Mapping[str, Any], error: BaseException) -> int:
    """Publish an infrastructure failure without losing requested case IDs.

    The global protocol validates the exact selected case registry even when a
    suite cannot start its compiler, parser, or native subprocess.  Emitting
    only a suite-level ERROR would therefore replace the useful root cause
    with a secondary ``registry mismatch`` and make the requested cases vanish
    from the final inventory.  Attribute the same control-plane failure to
    every requested identity instead.  This does *not* claim that each tested
    assertion ran: the explicit ERROR state records precisely the opposite.

    ``selected_cases`` is normally resolved by the global scheduler.  The
    expected registry is the fail-safe fallback for a directly constructed
    context.  Malformed selections fall back once more to the expected list so
    error publication itself remains best-effort and cannot conceal the first
    exception.
    """

    message = "{}: {}".format(type(error).__name__, error)
    test = context.get("test", {})
    expected = test.get("expected_cases", []) if isinstance(test, Mapping) else []
    expected_is_valid = (
        isinstance(expected, list) and bool(expected) and
        all(isinstance(identifier, str) and identifier for identifier in expected) and
        len(set(expected)) == len(expected))
    selected = test.get("selected_cases", expected) if isinstance(test, Mapping) else []
    # A malformed selection is itself a control-plane error, but publishing
    # duplicate, foreign, or zero rows would make the global validator replace
    # the original exception with another registry diagnostic.  Fall back to
    # the complete known-good owner registry so the first exception remains
    # the visible reason for every case.  Reverse/shuffled subsets are valid;
    # only membership and uniqueness, not canonical order, are required here.
    selected_is_valid = (
        expected_is_valid and isinstance(selected, list) and bool(selected) and
        all(isinstance(identifier, str) and identifier for identifier in selected) and
        len(set(selected)) == len(selected) and
        all(identifier in expected for identifier in selected))
    if not selected_is_valid:
        selected = expected if expected_is_valid else []
    rows = [case(identifier, "ERROR", message) for identifier in selected]
    return write_result(context, "ERROR", message, cases=rows)


def case(case_id: str, status: str, message: str = "") -> Dict[str, Any]:
    status = status.upper()
    if status not in STATUSES:
        raise RunnerError("invalid case status '{}'".format(status))
    return {"id": case_id, "status": status, "message": " ".join(message.split())}


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _quote(value: str) -> str:
    import shlex
    return shlex.quote(str(value))


def main_guard(function) -> None:
    """Execute a runner main function with a useful last-resort diagnostic.

    Individual runners should translate tested assertion failures themselves.
    Reaching this guard means their control plane failed, so exit 2 classifies
    the outcome as ERROR even if no result receipt could be written.
    """

    try:
        sys.exit(function())
    except RunnerError as error:
        print("runner error: {}".format(error), file=sys.stderr)
        sys.exit(2)
