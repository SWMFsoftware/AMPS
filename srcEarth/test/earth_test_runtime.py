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
import math
import os
import re
import signal
import socket
import subprocess
import sys
from dataclasses import dataclass, field
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


EXTREMUM_TIME_PASS = "PASS"
EXTREMUM_TIME_FAIL = "FAIL"
EXTREMUM_TIME_NOT_EVALUATED = "NOT_EVALUATED"


@dataclass
class ExtremumTimeEvaluation:
    """Outcome of a storm-extremum timing criterion (evaluate_extremum_time())."""

    status: str
    reason: str
    observed_extremum_time: Optional[datetime] = None
    observed_extremum_separation_sigma: Optional[float] = None
    modeled_extremum_time: Optional[datetime] = None
    modeled_extremum_candidates: List[datetime] = field(default_factory=list)
    modeled_extremum_margin: Optional[float] = None
    time_error_minutes: Optional[float] = None
    argextremum_time_error_minutes: Optional[float] = None


def evaluate_extremum_time(
    observed: Sequence[Tuple[datetime, float, Optional[float], Optional[float]]],
    modeled: Sequence[Tuple[datetime, float]],
    *,
    extremum: str,
    tie_tolerance: float,
    min_separation_sigma: float,
    max_error_minutes: float,
) -> ExtremumTimeEvaluation:
    """Compare the time of a storm extremum, model against observation.

    Used by C9 (minimum cutoff latitude) and C10 (maximum boundary expansion).
    ``extremum`` is "min" or "max".  ``observed`` holds (interval time UTC,
    value, sigma_up, sigma_down), where sigma_up/sigma_down are the one-sided
    uncertainties toward larger/smaller values; ``modeled`` holds (interval
    time UTC, value).  Values and sigmas share one unit (deg).

    Taking the argmin/argmax of a few sampled intervals and comparing the two
    times decides the result by noise when two intervals are nearly equal: the
    time error can only take the gaps between sampled intervals, so a rounding
    change swaps it between 0 and several hours (C9, 2026-11-19: 0.004 deg
    between two minima, 0 vs 188 min).  The criterion is therefore stated in
    terms of what each side resolves:

    * Model: every interval within ``tie_tolerance`` of the extreme modeled
      value is a candidate; the time error is measured to the nearest
      candidate.  A real timing error still fails, because the observed time is
      not a candidate when the model's extremum exceeds its value there by more
      than the tolerance.
    * Observation: the time is scored only when the observed extremum is beyond
      the next interval by at least ``min_separation_sigma`` combined sigma (the
      extremum's sigma toward the next value and the next value's sigma toward
      the extremum, in quadrature).  Otherwise the status is NOT_EVALUATED,
      which callers report as such and never as PASS.

    Failure behavior: no modeled values gives FAIL; no or one observed value,
    or a missing/zero uncertainty, gives NOT_EVALUATED with the reason.
    """
    if extremum not in ("min", "max"):
        raise ValueError("extremum must be 'min' or 'max', got %r" % (extremum,))
    sign = 1.0 if extremum == "min" else -1.0
    word = "minimum" if extremum == "min" else "maximum"

    if not modeled:
        return ExtremumTimeEvaluation(EXTREMUM_TIME_FAIL, "no modeled values for the %s" % word)

    # Work with sign*value so that the extremum is always the smallest value.
    ordered_model = sorted(modeled, key=lambda item: sign * item[1])
    extreme = sign * ordered_model[0][1]
    margin = (sign * ordered_model[1][1] - extreme) if len(ordered_model) > 1 else None
    candidates = sorted(t for t, v in modeled if sign * v - extreme <= tie_tolerance)
    out = ExtremumTimeEvaluation(
        EXTREMUM_TIME_NOT_EVALUATED, "",
        modeled_extremum_time=ordered_model[0][0],
        modeled_extremum_candidates=candidates,
        modeled_extremum_margin=margin,
    )

    if not observed:
        out.reason = "no observed values for the %s" % word
        return out
    ordered_obs = sorted(observed, key=lambda item: sign * item[1])
    obs_time, obs_value, obs_up, obs_down = ordered_obs[0]
    out.observed_extremum_time = obs_time

    def minutes(a: datetime, b: datetime) -> float:
        return abs((a - b).total_seconds()) / 60.0

    out.time_error_minutes = min(minutes(t, obs_time) for t in candidates)
    out.argextremum_time_error_minutes = minutes(ordered_model[0][0], obs_time)

    if len(ordered_obs) < 2:
        out.reason = "fewer than two observed intervals; the %s is not resolved" % word
        return out
    _, next_value, next_up, next_down = ordered_obs[1]
    # Toward the next value: a minimum's upward and the next value's downward
    # sigma; for a maximum the opposite.
    sigma_a = obs_up if extremum == "min" else obs_down
    sigma_b = next_down if extremum == "min" else next_up
    if sigma_a is None or sigma_b is None:
        out.reason = "observed uncertainty missing; the separation of the %s is unknown" % word
        return out
    combined = math.hypot(float(sigma_a), float(sigma_b))
    if combined <= 0.0:
        out.reason = "observed uncertainty is zero; the separation of the %s is unknown" % word
        return out
    separation = sign * (next_value - obs_value) / combined
    out.observed_extremum_separation_sigma = separation
    if separation < min_separation_sigma:
        out.reason = ("observed %s is %.2f sigma beyond the next interval (< %.2f sigma); "
                      "the observation does not resolve which interval holds it"
                      % (word, separation, min_separation_sigma))
        return out

    error = out.time_error_minutes
    if error <= max_error_minutes:
        out.status = EXTREMUM_TIME_PASS
        out.reason = ("nearest modeled %s candidate is %.1f min from the observed %s "
                      "(<= %.1f min)" % (word, error, word, max_error_minutes))
    else:
        out.status = EXTREMUM_TIME_FAIL
        out.reason = ("nearest modeled %s candidate is %.1f min from the observed %s "
                      "(> %.1f min)" % (word, error, word, max_error_minutes))
    return out


# Options that only say how to treat an existing work directory.  They are left out
# of the configuration key, so "--skip-run" with the options of an earlier run
# finds that run's directory.
WORKDIR_KEY_IGNORED_OPTIONS = ("--skip-run", "--keep")


def configuration_workdir(base: str, argv: Optional[Sequence[str]] = None) -> str:
    """Return the default work directory of one test configuration.

    Several entries of srcEarth/test/list run the same wrapper with different
    options.  The wrappers delete their default work directory when they start,
    so with a fixed default two entries that the runner starts in parallel
    (``test_runner.py -j N``) delete and overwrite each other's files.  The C4,
    C5 and C14 entries failed that way, and AMPS then crashed with SIGSEGV
    writing into the deleted directory.

    The returned directory is ``base`` when no options are given, so a bare
    invocation and the documented ``--skip-run`` examples keep their paths.
    Otherwise it is ``base`` + ``_`` + the first 10 hex digits of the SHA-1 of
    the options (in the order given, NUL-separated), leaving out
    WORKDIR_KEY_IGNORED_OPTIONS.  The same options give the same directory;
    different options give different directories.  An explicit --workdir is
    used as given and bypasses this function.
    """
    tokens = list(sys.argv[1:] if argv is None else argv)
    key = [t for t in tokens if t not in WORKDIR_KEY_IGNORED_OPTIONS]
    if not key:
        return base
    digest = hashlib.sha1("\0".join(key).encode("utf-8")).hexdigest()[:10]
    return "%s_%s" % (base.rstrip("/"), digest)


def record_workdir_configuration(workdir: Path, argv: Optional[Sequence[str]] = None) -> None:
    """Write the options that selected a configuration work directory into it.

    The directory name carries only a hash, so the options are kept beside the
    results as ``workdir_configuration.json`` for a reader of the directory.
    """
    tokens = list(sys.argv[1:] if argv is None else argv)
    payload = {
        "argv": tokens,
        "key_options": [t for t in tokens if t not in WORKDIR_KEY_IGNORED_OPTIONS],
        "ignored_options": list(WORKDIR_KEY_IGNORED_OPTIONS),
        "written_utc": utc_now(),
    }
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "workdir_configuration.json").write_text(json.dumps(payload, indent=2) + "\n")


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

