#!/usr/bin/env python3
"""Shared fail-closed runtime helpers for executable Earth validation tests.

The scientific runners under ``srcEarth/test`` accept an ``--amps`` path so a
reviewed executable can be tested without rebuilding it.  A filesystem check is
not sufficient: one AMPS checkout can be rebuilt for another application while
leaving the same executable pathname in place.  C9/C10 experienced exactly that
failure mode when a historical absolute path began naming the SEP application.

This module deliberately checks *application identity*, not numerical results.
It does not alter a solver command, input deck, reference product, tolerance, or
pass/fail threshold.  Its other responsibility is atomic status/provenance I/O,
so an interrupted or rejected launch cannot leave an older result looking like
the result of the new invocation.

``report_amps_failure`` turns a failed AMPS launch into a self-contained
diagnosis: the decoded exit status (for example ``139 = 128 + 11``, SIGSEGV in
an MPI rank), the crash messages MPI/UCX/the C++ runtime wrote into the AMPS
log, the log tail, and the host limits that commonly matter for a crash.  It
only reports; the caller still returns the unchanged nonzero status, so a
failed run stays a failed test.
"""

from __future__ import annotations

import collections
import hashlib
import json
import os
import re
import signal
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, TextIO, Tuple

try:  # POSIX only; the Earth validation runs on Linux
    import resource
except ImportError:  # pragma: no cover
    resource = None  # type: ignore[assignment]


EARTH_HELP_BANNER = "AMPS Earth: Gridless Energetic Particle Solver"
EARTH_HELP_REQUIRED_TOKENS = (
    EARTH_HELP_BANNER,
    "-mode <3d|gridless>",
    "No solver is run.",
)


class EarthExecutableError(RuntimeError):
    """Raised when an ``--amps`` target is absent or is not the Earth driver."""


def utc_now() -> str:
    """Return a stable UTC timestamp for JSON provenance."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_file(path: Path) -> str:
    """Hash a file without loading a multi-megabyte executable into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def probe_earth_executable(
    executable: Path,
    *,
    cwd: Optional[Path] = None,
    timeout_seconds: float = 30.0,
) -> Dict[str, Any]:
    """Resolve and positively identify an AMPS Earth executable.

    The probe invokes only ``-h``.  The Earth driver documents that this exits
    before model initialization, so the check does not initialize MPI, Geopack,
    SPICE, a mesh, or a physical calculation.  Requiring several Earth-specific
    help tokens prevents an executable for SEP or another AMPS application from
    reaching a costly multi-rank launch and failing after it has created a
    misleading partial output tree.
    """

    path = executable.expanduser()
    if not path.is_absolute():
        path = ((cwd or Path.cwd()) / path).resolve()
    else:
        path = path.resolve()

    if not path.is_file():
        raise EarthExecutableError("AMPS executable not found: %s" % path)
    if not os.access(str(path), os.X_OK):
        raise EarthExecutableError("AMPS executable is not executable: %s" % path)

    try:
        completed = subprocess.run(
            [str(path), "-h"],
            cwd=str((cwd or Path.cwd()).resolve()),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        raise EarthExecutableError(
            "AMPS executable identity probe timed out after %.1f s: %s"
            % (timeout_seconds, path)
        ) from exc
    except OSError as exc:
        raise EarthExecutableError(
            "AMPS executable identity probe could not start %s: %s" % (path, exc)
        ) from exc

    output = completed.stdout or ""
    missing = [token for token in EARTH_HELP_REQUIRED_TOKENS if token not in output]
    if completed.returncode != 0 or missing:
        details = []
        if completed.returncode != 0:
            details.append("help exit=%d" % completed.returncode)
        if missing:
            details.append("missing Earth help token(s): %s" % ", ".join(missing))
        raise EarthExecutableError(
            "AMPS executable is not the required Earth driver: %s (%s)"
            % (path, "; ".join(details))
        )

    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "identity_probe": "-h",
        "identity_banner": EARTH_HELP_BANNER,
        "identity_probe_return_code": completed.returncode,
    }


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Atomically replace a JSON status/result file in its destination directory."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".%s.tmp.%d" % (path.name, os.getpid()))
    try:
        temporary.write_text(json.dumps(dict(payload), indent=2) + "\n", encoding="utf-8")
        os.replace(str(temporary), str(path))
    finally:
        # ``os.replace`` removes the temporary path.  The cleanup matters only
        # when serialization or replacement raises midway through the write.
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def archive_previous_result(path: Path, *, run_started_utc: str) -> Optional[Path]:
    """Move a prior top-level result aside before a new physical calculation.

    The archive is recoverable and never overwrites an earlier archive.  Branch
    directories are managed by C9/C10 themselves, but their top-level result can
    otherwise survive a failed first launch and be mistaken for current output.
    """

    if not path.exists():
        return None
    token = run_started_utc.replace("-", "").replace(":", "")
    candidate = path.with_name("%s.previous-%s%s" % (path.stem, token, path.suffix))
    suffix = 1
    while candidate.exists():
        candidate = path.with_name(
            "%s.previous-%s-%d%s" % (path.stem, token, suffix, path.suffix)
        )
        suffix += 1
    os.replace(str(path), str(candidate))
    return candidate


def initial_run_status(
    *,
    schema: str,
    test_id: str,
    output_root: Path,
    argv: Sequence[str],
    execution_mode: str,
) -> Dict[str, Any]:
    """Create the fail-closed status object written before solver execution."""

    started = utc_now()
    return {
        "schema": schema,
        "test_id": test_id,
        "status": "STARTING",
        "started_utc": started,
        "finished_utc": None,
        "execution_mode": execution_mode,
        "output_root": str(output_root.resolve()),
        "argv": list(argv),
        "executable": None,
        "archived_previous_result": None,
        "result_current": False,
        "return_code": None,
        "message": "validation invocation has not completed",
    }


def finish_run_status(
    status: Dict[str, Any],
    *,
    state: str,
    return_code: int,
    message: str,
    result_current: bool,
) -> Dict[str, Any]:
    """Finalize a status object without interpreting any scientific threshold."""

    status.update(
        {
            "status": state,
            "finished_utc": utc_now(),
            "result_current": bool(result_current),
            "return_code": int(return_code),
            "message": message,
        }
    )
    return status


# --------------------------------------------------------------------------
# Failure diagnosis of an AMPS launch
# --------------------------------------------------------------------------

# Prefix of the one-line diagnosis.  test_runner.py copies the last such line
# of a test log into its reports, and the nightly reader picks it as the last
# line containing "FAIL", so it is printed as the final line of the report.
FAILURE_CAUSE_PREFIX = "FAILURE CAUSE:"

# Messages written when a process dies abnormally, most specific first.  The
# first group names the rank and signal (Open MPI mpirun, MPICH/Hydra and Intel
# MPI); the second is written by the dying process itself (UCX and Open MPI
# signal handlers, the C++ runtime, assertions, MPI_Abort); the third is a
# generic kernel/shell wording.
CRASH_MESSAGE_PATTERNS: Tuple[Tuple[int, "re.Pattern[str]"], ...] = (
    (0, re.compile(r"exited on signal \d+|KILLED BY SIGNAL|BAD TERMINATION OF ONE OF YOUR APPLICATION")),
    (1, re.compile(
        r"Caught signal \d+|Process received signal|Signal: .*\(\d+\)|Failing at address|"
        r"terminate called|what\(\):|std::bad_alloc|Assertion .* failed|"
        r"MPI_ABORT was invoked|MPI_Abort\(|[Oo]ut of memory|oom-kill"
    )),
    (2, re.compile(
        r"Segmentation fault|Bus error|Floating point exception|Illegal instruction|"
        r"\bAborted\b|\bKilled\b|core dumped"
    )),
)

# Largest number 128 + N decoded as a signal (SIGRTMAX is 64 on Linux).
_MAX_SIGNAL = 64


def describe_signal(number: int) -> str:
    """Return ``signal N (NAME: description)`` for a signal number."""

    try:
        name = signal.Signals(number).name
    except ValueError:
        name = "SIG%d" % number
    strsignal = getattr(signal, "strsignal", None)
    description = (strsignal(number) if strsignal else None) or "unknown signal"
    return "signal %d (%s: %s)" % (number, name, description)


def describe_exit_status(returncode: Optional[int], *, compact: bool = False) -> str:
    """Describe a process exit status exactly, without guessing its cause.

    ``subprocess`` reports ``-N`` when the process it started was itself killed
    by signal N.  A shell, and an MPI launcher such as Open MPI ``mpirun`` (seen
    here with Open MPI 4.1.7: a rank's SIGSEGV gives ``mpirun`` status 139),
    instead exit with ``128 + N`` when a process they ran was killed by signal N.
    A program may also exit with such a value deliberately, so ``128 + N`` is
    reported as the convention it is.  test_runner.py keeps an identical copy
    (``describe_exit_code``) so that it stays self-contained; UTR-F09 checks
    that the two agree.  ``compact`` gives the short form used in one-line
    causes, e.g. ``exit 139 = 128 + 11 (SIGSEGV: Segmentation fault)``.
    """

    if returncode is None:
        if compact:
            return "no exit status"
        return "no exit status (the process did not start or was not waited for)"
    if returncode == 0:
        return "exit 0" if compact else "exit status 0"
    if returncode < 0:
        return "killed by %s" % describe_signal(-returncode)
    if 128 < returncode <= 128 + _MAX_SIGNAL:
        number = returncode - 128
        if compact:
            return "exit %d = 128 + %d (%s)" % (
                returncode, number, describe_signal(number).split(" (", 1)[1].rstrip(")"))
        return (
            "exit status %d = 128 + %d: by the shell/MPI-launcher convention, "
            "a process was killed by %s"
            % (returncode, number, describe_signal(number))
        )
    return ("exit %d" if compact else "exit status %d") % returncode


def scan_crash_messages(
    path: Path, *, tail_lines: int = 30, max_messages: int = 20,
) -> Tuple[List[Tuple[int, int, str]], List[str], Optional[Tuple[int, int, str]]]:
    """Return ``(messages, tail, primary)`` for a log file.

    ``messages`` are the last ``max_messages`` crash messages as
    ``(priority, line_no, text)``, ``tail`` the last ``tail_lines`` lines, and
    ``primary`` the first message of the most specific kind: in an MPI job the
    first process to fail usually takes the others down, so later messages are
    mostly consequences.  The file is streamed once, so a multi-gigabyte AMPS
    log costs no memory.  A missing or unreadable log gives ``([], [], None)``.
    """

    messages: "collections.deque[Tuple[int, int, str]]" = collections.deque(maxlen=max_messages)
    tail: "collections.deque[str]" = collections.deque(maxlen=max(tail_lines, 0))
    primary: Optional[Tuple[int, int, str]] = None
    try:
        with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
            for line_no, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if tail.maxlen:
                    tail.append(line)
                for priority, pattern in CRASH_MESSAGE_PATTERNS:
                    if pattern.search(line):
                        message = (priority, line_no, line.strip())
                        messages.append(message)
                        if primary is None or priority < primary[0]:
                            primary = message
                        break
    except OSError:
        return [], [], None
    return list(messages), list(tail), primary


def _format_limit(value: int) -> str:
    if resource is not None and value == resource.RLIM_INFINITY:
        return "unlimited"
    if value >= 1024 and value % 1024 == 0:
        return "%d KiB" % (value // 1024)
    return "%d bytes" % value


def failure_environment() -> List[str]:
    """Facts about this process that matter for a crash, as ``name: value`` lines.

    They describe the wrapper process; mpirun and the ranks it starts on this
    host inherit the resource limits and, unless the launcher changes them, the
    environment variables.
    """

    lines = ["host: %s" % socket.gethostname()]
    if resource is not None:
        for label, limit in (("stack size limit (ulimit -s)", resource.RLIMIT_STACK),
                             ("core file limit (ulimit -c)", resource.RLIMIT_CORE)):
            soft, _hard = resource.getrlimit(limit)
            lines.append("%s: %s" % (label, _format_limit(soft)))
    for name in ("OMP_NUM_THREADS", "OMP_STACKSIZE", "OMP_PROC_BIND", "OMP_PLACES"):
        lines.append("%s=%s" % (name, os.environ.get(name, "<unset>")))
    return lines


def report_amps_failure(
    returncode: Optional[int],
    command: Sequence[str],
    workdir: Path,
    log_path: Path,
    *,
    label: str,
    tail_lines: int = 30,
    stream: Optional[TextIO] = None,
) -> str:
    """Print and log a diagnosis of a failed AMPS launch; return its cause line.

    The report goes to ``stream`` (default: stderr, which the test runner
    captures in the test log) and is appended to ``log_path``.  Its last line is
    ``FAILURE CAUSE: <label>: <exit status>; <most specific crash message>``.
    Nothing is inferred beyond what the exit status and the log state: when the
    log holds no crash message, the cause says so.
    """

    log_path = Path(log_path)
    status = describe_exit_status(returncode)
    messages, tail, primary = scan_crash_messages(log_path, tail_lines=tail_lines)
    if primary is not None:
        detail = primary[2]
    elif not log_path.exists():
        detail = "the AMPS log %s does not exist" % log_path
    else:
        detail = "no crash message (signal, rank, backtrace, abort) in the AMPS log"
    cause = "%s %s: AMPS %s; %s" % (
        FAILURE_CAUSE_PREFIX, label, describe_exit_status(returncode, compact=True), detail)

    rule = "=" * 80
    lines = [
        rule,
        "AMPS FAILURE REPORT: %s" % label,
        "  status:   %s" % status,
        "  command:  %s" % " ".join(str(part) for part in command),
        "  workdir:  %s" % workdir,
        "  AMPS log: %s" % log_path,
    ]
    lines += ["  %s" % item for item in failure_environment()]
    if messages:
        lines.append("  crash messages in the AMPS log (last %d shown, with line numbers):" % len(messages))
        lines += ["    L%d: %s" % (line_no, text) for _priority, line_no, text in messages]
    else:
        lines.append("  crash messages in the AMPS log: none found")
    if tail:
        lines.append("  last %d lines of the AMPS log:" % len(tail))
        lines += ["    | %s" % line for line in tail]
    lines += [rule, cause]
    text = "\n".join(lines) + "\n"

    # AMPS output was teed to stdout; flush it so the report follows it in the
    # test log instead of overtaking buffered lines.
    sys.stdout.flush()
    out = stream if stream is not None else sys.stderr
    out.write(text)
    out.flush()
    try:
        with log_path.open("a", encoding="utf-8") as log:
            log.write("\n" + text)
    except OSError as exc:
        out.write("note: could not append the failure report to %s: %s\n" % (log_path, exc))
    return cause

