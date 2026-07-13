#!/usr/bin/env python3
"""
Concurrent pass/fail test runner for AMPS validation commands.

Purpose
-------
The validation list contains shell commands marked with an expected result:

  ! this is a comment

  P command expected to return exit code 0
  last pass: optional git commit id for this individual test

  F command expected to return nonzero exit code
  last pass: optional git commit id for this individual test

Blank lines are ignored.  Lines beginning with '!' after leading whitespace are
ignored.  The first non-whitespace token on a test line must be P or F; the rest
of the line is executed as the command.  A test command may optionally be
followed by a metadata line beginning with 'last pass:'.  Metadata lines are
associated with the immediately preceding test command and are not executed.

Use cases
---------
1. Run the full validation list and report mismatches between expected and
   actual P/F states.
2. Run many short serial tests concurrently without changing the individual
   test wrappers.
3. Run MPI/OpenMP AMPS tests concurrently while avoiding physical memory
   oversubscription, by watching real available memory instead of a static
   per-command CPU estimate.  This matters for commands such as

       srcEarth/test/C3/run_C3.py -np 4 -nt 16

   where several such commands running at once can exhaust physical memory
   even though nothing here changes how many ranks or threads any one of
   them uses.
4. Keep per-test provenance by updating each test's ``last pass:`` commit id
   after successful command execution.

Scheduling algorithm
--------------------
Concurrency is governed by two independent things:

* ``-j``/``--jobs`` is a hard ceiling on how many commands may be running at
  once.

* A physical-memory ramp-up gate.  Commands are never rewritten: every test
  always runs with exactly the ``-np``/``-nt`` (or any other options) it was
  written with.  Instead of estimating CPU cost ahead of time, the runner
  watches real memory headroom and paces new launches against it:

    1. Start a test.
    2. Wait ``--memory-check-interval`` seconds (default 60).
    3. Check available physical memory. If more than
       ``--min-free-memory-fraction`` (default 0.5, i.e. half) of total
       physical memory is still available, start another test.
    4. Repeat until ``-j`` tests are running concurrently.

  The only launches exempt from this gate are ones where nothing of this
  runner's is currently running: the very first test of the run, and any
  later launch that happens to start right after the pool has fully drained
  back to zero.  In both cases there is no in-flight test to pace against, so
  waiting out a stale timer would just be dead time; the runner launches
  immediately (after a best-effort memory reading for display) and the
  pacing clock restarts from that launch.  Every other launch -- ramping up
  toward ``-j`` while other tests are still running, or backfilling a slot
  next to tests that are still running -- waits out the interval-and-headroom
  check described above.  If memory is too tight, the runner keeps waiting
  and rechecking on that cadence rather than giving up: unlike a fixed
  CPU-slot pool, memory headroom can recover on its own as other tests
  finish or the kernel reclaims cache, so there is nothing to deadlock on.

  One consequence worth knowing: concurrency can only build up between tests
  that are actually running at the same time the gate re-checks.  A test
  that finishes faster than ``--memory-check-interval`` never gets a sibling
  launched next to it -- by the time the next launch is allowed, it has
  already exited, so that pair effectively runs serially regardless of
  ``-j``.  For AMPS's MPI/OpenMP tests (tens of seconds to tens of minutes)
  this is rarely an issue at the 60s default; a validation list dominated by
  short/cheap tests will see less concurrency than ``-j`` suggests unless
  ``--memory-check-interval`` is lowered to something shorter than those
  tests' typical runtime.

An earlier version of this runner estimated CPU cost from ``-np``/``-nt`` and
used it both to cap concurrency and to silently lower a test's ``-nt`` so more
tests would fit an assumed CPU budget. That changed what was actually being
validated -- a test's reference P/F result reflects the process/thread count
it was authored and last validated with, not whatever a scheduling formula
computed for a given ``-j`` -- and the CPU estimate never modeled memory,
which is what concurrently running MPI ranks actually contend for regardless
of how many OpenMP threads each rank uses. This version tracks only real
memory headroom and never touches a test's command line. Pass
``--no-memory-gate`` to fall back to plain ``-j``-limited concurrency with no
pacing or memory checks at all (every test launches as soon as a ``-j`` slot
is free) -- also needed on systems without ``/proc/meminfo``.

Memory detection
----------------
Available memory is read from ``/proc/meminfo``'s ``MemAvailable`` field,
which already accounts for reclaimable page cache/buffers and is a much
better estimate of "how much can a new allocation actually get" than raw free
memory. Kernels too old to report ``MemAvailable`` fall back to ``MemFree +
Buffers + Cached``. This is Linux-specific; on a system without
``/proc/meminfo`` the runner exits immediately with an explanatory error
before running anything, unless ``--no-memory-gate`` is given.

Thread environment policy
-------------------------
By default the runner does not alter ``OMP_NUM_THREADS`` or BLAS thread-count
environment variables. Passing ``--set-thread-env`` sets those variables for
each command to its own parsed ``-nt`` value, which can reduce hidden nested
threading in libraries. Existing values are preserved unless
``--overwrite-thread-env`` is also supplied. This never changes what ``-nt``
a command is launched with; it only sets environment variables to match the
value the command already declares.

Progress and result handling
----------------------------
The runner reports progress in three places. First, when a command has
passed the ``-j`` job-count limit and the memory gate and is actually being
launched, it prints a ``[START]`` line together with the memory headroom
observed at that moment. Second, if the memory gate is currently blocking the
next launch, the runner prints a ``[WAITING]`` line on the same cadence as
``--memory-check-interval`` showing the most recently observed free-memory
fraction, so a terminal that is simply waiting out the gate does not look
like a hung runner. Third, as soon as each command finishes, the runner
prints the historical completion line beginning with ``[OK]`` or
``[MISMATCH]``; that line is intentionally kept compatible with the original
runner output.

The process exit status is converted to actual P/F. Timeouts and launch
errors are treated as actual F. The full JSON/CSV reports include each
test's parsed ``-np``/``-nt`` and the physical-memory reading observed at the
moment it was launched. Two triage reports are also written:

  <report-prefix>_to_address.txt/.csv
      tests whose actual P/F result does not match the reference P/F marker;
      these are the tests that require investigation.
  <report-prefix>_actual_failed.txt/.csv
      all commands that returned actual F, including expected-failure tests.

Use ``--update-last-pass`` to rewrite the ``last pass:`` line of each test whose
command actually exited 0 on this run. Existing metadata lines are updated in
place; missing metadata lines are inserted below the corresponding test command.
This is independent of whether the result matched the reference marker: an
F-marked test that unexpectedly exits 0 still gets its ``last pass:`` commit id
updated and is also reported as an unexpected pass.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import json
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, List, Optional


@dataclass
class TestCase:
    index: int
    line_no: int
    expected: str  # "P" or "F"
    command: str
    last_pass: Optional[str] = None
    last_pass_line_no: Optional[int] = None


@dataclass
class TestResult:
    index: int
    line_no: int
    expected: str
    actual: str
    matched_reference: bool
    exit_code: Optional[int]
    timed_out: bool
    elapsed_s: float
    command: str
    log_file: str
    np_value: int = 1
    nt_value: int = 1
    np_nt_details: str = ""
    mem_total_bytes: Optional[int] = None
    mem_available_bytes_at_launch: Optional[int] = None
    mem_available_fraction_at_launch: Optional[float] = None
    last_pass: Optional[str] = None
    last_pass_line_no: Optional[int] = None


@dataclass
class TestPlan:
    """Per-test info the dispatcher and ``run_one_test`` need.

    There is no launch-time command rewriting: ``np_value``/``nt_value`` are
    parsed only for reporting and for ``--set-thread-env``. The command that
    actually runs is always ``test.command``, unmodified.
    """

    test: TestCase
    np_value: int
    nt_value: int
    np_nt_details: str


LAST_PASS_RE = re.compile(r"^(?P<prefix>\s*last\s+pass\s*:\s*)(?P<value>.*)$", re.IGNORECASE)


def _parse_positive_int(text: str) -> Optional[int]:
    """Parse an integer-like value and reject zero/negative values.

    Some wrapper scripts pass values that look like floats, e.g. ``"4.0"``.
    Accepting ``int(float(x))`` makes the parser tolerant without allowing
    invalid nonpositive values.
    """
    try:
        value = int(float(str(text)))
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


def _option_values(argv: List[str], names: Iterable[str]) -> List[str]:
    """Return values for options that can appear as ``--x v`` or ``--x=v``."""
    names = set(names)
    values: List[str] = []
    i = 0
    while i < len(argv):
        token = argv[i]
        if token in names and i + 1 < len(argv):
            values.append(argv[i + 1])
            i += 2
            continue
        for name in names:
            prefix = name + "="
            if token.startswith(prefix):
                values.append(token[len(prefix):])
                break
        i += 1
    return values


def _first_positive_option(argv: List[str], names: Iterable[str]) -> Optional[int]:
    """Return the first positive integer value found for any option name."""
    for value in _option_values(argv, names):
        parsed = _parse_positive_int(value)
        if parsed is not None:
            return parsed
    return None


def parse_np_nt(command: str, default_np: int, default_nt: int) -> tuple[int, int, str]:
    """Parse MPI rank count and thread count out of a command, for reporting only.

    This used to feed a CPU-slot scheduling estimate; it no longer does.
    Nothing in this module uses the returned values to change scheduling or
    to rewrite a command. They are kept so reports stay informative and so
    ``--set-thread-env`` still has an ``-nt`` value to propagate into
    thread-count environment variables.
    """
    try:
        argv = shlex.split(command, posix=True)
    except ValueError:
        argv = command.split()

    np_value = _first_positive_option(argv, ["-np", "-n", "--np", "--n", "--mpi-ranks"])
    nt_value = _first_positive_option(
        argv,
        [
            "-nt",
            "--nt",
            "--threads",
            "--mode3d-threads",
            "--density-threads",
            "--omp-threads",
        ],
    )

    if np_value is None:
        np_value = max(1, default_np)
    if nt_value is None:
        nt_value = max(1, default_nt)

    details = f"np={np_value}, nt={nt_value}"
    return np_value, nt_value, details


def set_thread_env(base_env: dict, nt: int, only_if_missing: bool) -> dict:
    """Return environment for one command with optional thread-count controls.

    AMPS wrappers normally receive ``-nt`` explicitly, but imported numerical
    libraries can still start their own thread pools. Setting these variables
    is optional because some HPC environments already configure them
    carefully. ``only_if_missing`` preserves existing settings unless the
    user requested ``--overwrite-thread-env``.
    """
    env = base_env.copy()
    if nt < 1:
        return env
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        if only_if_missing and env.get(name):
            continue
        env[name] = str(nt)
    return env


def parse_test_file(path: Path) -> List[TestCase]:
    tests: List[TestCase] = []
    pending_metadata_test: Optional[TestCase] = None

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()

            if not line:
                # Allow a blank line between a command and its metadata, but do
                # not otherwise give blank lines semantic meaning.
                continue

            if line.startswith("!"):
                # Section comments are not metadata. Reset the pending command
                # so a later stray 'last pass:' line is reported clearly.
                pending_metadata_test = None
                continue

            last_pass_match = LAST_PASS_RE.match(raw.rstrip("\n"))
            if last_pass_match:
                if pending_metadata_test is None:
                    raise ValueError(
                        f"Invalid test-list line {line_no}: 'last pass:' must "
                        "immediately follow a P/F command"
                    )
                if pending_metadata_test.last_pass_line_no is not None:
                    raise ValueError(
                        f"Invalid test-list line {line_no}: duplicate 'last pass:' "
                        f"metadata for test from line {pending_metadata_test.line_no}"
                    )
                pending_metadata_test.last_pass = last_pass_match.group("value").strip()
                pending_metadata_test.last_pass_line_no = line_no
                pending_metadata_test = None
                continue

            parts = line.split(maxsplit=1)
            if len(parts) != 2 or parts[0].upper() not in {"P", "F"}:
                raise ValueError(
                    f"Invalid test-list line {line_no}: expected 'P <command>' "
                    f"or 'F <command>', got: {raw.rstrip()}"
                )

            test = TestCase(
                index=len(tests) + 1,
                line_no=line_no,
                expected=parts[0].upper(),
                command=parts[1].strip(),
            )
            tests.append(test)
            pending_metadata_test = test

    return tests


def safe_log_name(test: TestCase) -> str:
    # Keep names deterministic and readable while avoiding shell metacharacters.
    words = shlex.split(test.command, posix=True) if test.command.strip() else []
    base = Path(words[0]).name if words else "command"
    base = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in base)
    return f"test_{test.index:03d}_line_{test.line_no}_{base}.log"


def read_proc_meminfo(path: str = "/proc/meminfo") -> dict[str, int]:
    """Parse /proc/meminfo into a dict of field name -> value in bytes."""
    info: dict[str, int] = {}
    unit_multiplier = {"b": 1, "kb": 1024, "mb": 1024**2, "gb": 1024**3}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if ":" not in line:
                continue
            key, _, rest = line.partition(":")
            tokens = rest.split()
            if not tokens:
                continue
            try:
                value = int(tokens[0])
            except ValueError:
                continue
            unit = tokens[1].lower() if len(tokens) > 1 else "kb"
            info[key.strip()] = value * unit_multiplier.get(unit, 1024)
    return info


def read_memory_totals() -> tuple[int, int]:
    """Return (total_bytes, available_bytes) for physical memory.

    Uses ``/proc/meminfo``'s ``MemAvailable`` when present (Linux 3.14+),
    which already accounts for reclaimable page cache/buffers and is a far
    better estimate of real headroom than ``MemFree`` alone. Falls back to
    ``MemFree + Buffers + Cached`` on older kernels. Raises ``RuntimeError``
    with an actionable message if memory cannot be read at all (e.g.
    non-Linux), so callers can fail fast instead of silently mis-scheduling.
    """
    try:
        info = read_proc_meminfo()
    except OSError as exc:
        raise RuntimeError(
            "could not read /proc/meminfo to track available physical memory "
            "(memory tracking currently requires Linux); pass --no-memory-gate "
            "to run with plain -j concurrency and no memory checks"
        ) from exc

    total = info.get("MemTotal")
    if not total:
        raise RuntimeError("/proc/meminfo did not report a usable MemTotal")

    available = info.get("MemAvailable")
    if available is None:
        available = info.get("MemFree", 0) + info.get("Buffers", 0) + info.get("Cached", 0)

    return total, available


def format_bytes(n: int) -> str:
    """Human-readable binary size, e.g. 134744072192 -> '125.5 GiB'."""
    value = float(n)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024.0:
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{value:.1f} TiB"


def build_test_plan(test: TestCase, *, default_np: int, default_nt: int) -> TestPlan:
    np_value, nt_value, details = parse_np_nt(test.command, default_np, default_nt)
    return TestPlan(test=test, np_value=np_value, nt_value=nt_value, np_nt_details=details)


async def run_one_test(
    plan: TestPlan,
    *,
    workdir: Path,
    log_dir: Path,
    timeout_s: Optional[float],
    use_shell: bool,
    base_env: dict,
    set_thread_env_vars: bool,
    preserve_thread_env: bool,
    print_start: bool,
    mem_total_bytes: Optional[int],
    mem_available_bytes: Optional[int],
) -> TestResult:
    """Run one already-scheduled command, exactly as written in the test list.

    This coroutine does not wait on any throttle itself; ``run_all_tests``
    decides when a ``-j`` slot and the memory gate both allow a launch, then
    creates this task. Nothing here -- or anywhere else in this module --
    ever rewrites ``plan.test.command``: it always runs with the same
    ``-np``/``-nt``/other options it was authored with.
    """
    test = plan.test
    log_path = log_dir / safe_log_name(test)
    start = time.monotonic()
    exit_code: Optional[int] = None
    timed_out = False

    env = set_thread_env(base_env, plan.nt_value, preserve_thread_env) if set_thread_env_vars else base_env

    mem_fraction = (
        mem_available_bytes / mem_total_bytes
        if (mem_total_bytes and mem_available_bytes is not None)
        else None
    )

    if print_start:
        if mem_fraction is not None:
            mem_note = (
                f"mem_avail={mem_fraction:.0%} "
                f"({format_bytes(mem_available_bytes)}/{format_bytes(mem_total_bytes)}), "
            )
        else:
            mem_note = ""
        print(
            f"[START] #{test.index:03d} line {test.line_no}: "
            f"{mem_note}command={test.command}",
            flush=True,
        )

    if mem_fraction is not None:
        mem_line = (
            f"# Memory at launch: {mem_fraction:.1%} available "
            f"({format_bytes(mem_available_bytes)} / {format_bytes(mem_total_bytes)})\n"
        )
    else:
        mem_line = "# Memory at launch: not tracked (--no-memory-gate)\n"

    with log_path.open("wb") as log:
        header = (
            f"# Test {test.index} from line {test.line_no}\n"
            f"# Expected: {test.expected}\n"
            f"# Command: {test.command}\n"
            f"# Workdir: {workdir}\n"
            f"# {plan.np_nt_details}\n"
            f"{mem_line}"
            f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"{'=' * 80}\n"
        )
        log.write(header.encode("utf-8", errors="replace"))
        log.flush()

        try:
            if use_shell:
                # start_new_session=True makes timeout cleanup kill the whole
                # process group, including mpirun-launched children when they
                # remain in the same session.
                proc = await asyncio.create_subprocess_shell(
                    test.command,
                    cwd=str(workdir),
                    stdout=log,
                    stderr=asyncio.subprocess.STDOUT,
                    env=env,
                    start_new_session=True,
                )
            else:
                argv = shlex.split(test.command)
                if not argv:
                    raise ValueError("empty command")
                proc = await asyncio.create_subprocess_exec(
                    *argv,
                    cwd=str(workdir),
                    stdout=log,
                    stderr=asyncio.subprocess.STDOUT,
                    env=env,
                    start_new_session=True,
                )

            try:
                exit_code = await asyncio.wait_for(proc.wait(), timeout=timeout_s)
            except asyncio.TimeoutError:
                timed_out = True
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    await asyncio.wait_for(proc.wait(), timeout=10.0)
                except asyncio.TimeoutError:
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    await proc.wait()
                exit_code = proc.returncode

        except Exception as exc:  # launch failure is also a failed test
            log.write(f"\nERROR launching command: {exc}\n".encode("utf-8", errors="replace"))
            exit_code = None

        elapsed = time.monotonic() - start
        footer = (
            f"\n{'=' * 80}\n"
            f"# Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"# Elapsed seconds: {elapsed:.3f}\n"
            f"# Timed out: {timed_out}\n"
            f"# Exit code: {exit_code}\n"
        )
        log.write(footer.encode("utf-8", errors="replace"))

    actual = "P" if (exit_code == 0 and not timed_out) else "F"
    matched = actual == test.expected

    return TestResult(
        index=test.index,
        line_no=test.line_no,
        expected=test.expected,
        actual=actual,
        matched_reference=matched,
        exit_code=exit_code,
        timed_out=timed_out,
        elapsed_s=elapsed,
        command=test.command,
        log_file=str(log_path),
        np_value=plan.np_value,
        nt_value=plan.nt_value,
        np_nt_details=plan.np_nt_details,
        mem_total_bytes=mem_total_bytes,
        mem_available_bytes_at_launch=mem_available_bytes,
        mem_available_fraction_at_launch=mem_fraction,
        last_pass=test.last_pass,
        last_pass_line_no=test.last_pass_line_no,
    )


async def run_all_tests(
    tests: Iterable[TestCase],
    *,
    jobs: int,
    workdir: Path,
    log_dir: Path,
    timeout_s: Optional[float],
    use_shell: bool,
    default_np: int,
    default_nt: int,
    set_thread_env_vars: bool,
    preserve_thread_env: bool,
    print_start: bool,
    memory_gate_enabled: bool,
    min_free_memory_fraction: float,
    memory_check_interval: float,
) -> List[TestResult]:
    """Schedule tests behind a -j job-count cap and a physical-memory ramp-up gate.

    A launch into an empty pool -- the first test of the run, or any later
    test that starts once nothing else is running -- always launches
    immediately. Every other launch attempt -- ramping up toward ``jobs``
    concurrent tests while others are still running, or backfilling a slot
    next to tests that are still running -- must wait at least
    ``memory_check_interval`` seconds since the previous launch attempt and
    then observe more than ``min_free_memory_fraction`` of total physical
    memory still available. Tests are launched in list order; there is no
    per-test resource estimate to reorder around, since every test is
    subject to exactly the same gate regardless of its own ``-np``/``-nt``.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    base_env = os.environ.copy()

    pending: List[TestPlan] = [
        build_test_plan(test, default_np=default_np, default_nt=default_nt) for test in tests
    ]
    running: dict[asyncio.Task[TestResult], TestPlan] = {}
    results: List[TestResult] = []

    last_check_time: Optional[float] = None
    last_known_total: Optional[int] = None
    last_known_available: Optional[int] = None

    def try_gate() -> bool:
        """Return True if a new test may launch right now.

        A launch into an empty pool -- nothing of this runner's currently
        running -- always launches immediately and unconditionally: this
        covers the first test of the run, and any later test that starts
        once everything before it has finished. There is no in-flight
        launch to pace against in that case. Every other attempt must wait
        at least ``memory_check_interval`` seconds since the previous
        attempt (successful or not), then observe more than
        ``min_free_memory_fraction`` of total physical memory free. A
        transient failure to read /proc/meminfo during one of those later
        checks is treated as "not enough information to launch yet" rather
        than aborting the run, and is retried on the same cadence.

        Whenever this performs an actual memory read, it records the
        reading in the enclosing ``last_known_*`` variables so callers (the
        [START] line and TestResult) can report it even on calls that
        return False. It also prints a [WAITING] line itself whenever a
        real check comes back insufficient -- this is the only place that
        happens, since a blocked launch can occur equally whether other
        tests are currently running (the common case once the pool is warm)
        or not, and only this function knows when a check was actually
        performed rather than skipped for still being inside the pacing
        window.
        """
        nonlocal last_check_time, last_known_total, last_known_available
        if not memory_gate_enabled:
            return True

        now = time.monotonic()
        if not running:
            # Nothing of ours is currently running, so there is no ongoing
            # ramp-up to pace against: launch immediately. This always covers
            # the very first test of the run, and also re-covers any later
            # point where the pool has fully drained before the next test
            # starts -- there is no in-flight memory footprint left over from
            # a stale pacing timer to wait out. Still take a best-effort
            # reading for display; it does not gate this launch.
            last_check_time = now
            try:
                last_known_total, last_known_available = read_memory_totals()
            except RuntimeError:
                pass
            return True

        if last_check_time is not None and (now - last_check_time) < memory_check_interval:
            return False

        last_check_time = now
        try:
            total, available = read_memory_totals()
        except RuntimeError as exc:
            if print_start:
                print(
                    f"[WAITING] memory check failed ({exc}); will retry in "
                    f"{memory_check_interval:.0f}s ({len(running)} test(s) currently running)",
                    flush=True,
                )
            return False
        last_known_total, last_known_available = total, available
        if available > min_free_memory_fraction * total:
            return True
        if print_start:
            print(
                f"[WAITING] {available / total:.0%} of physical memory free "
                f"(need > {min_free_memory_fraction:.0%}); holding off starting "
                f"another test, rechecking in {memory_check_interval:.0f}s "
                f"({len(running)} test(s) currently running)",
                flush=True,
            )
        return False

    while pending or running:
        while pending and len(running) < jobs:
            if not try_gate():
                break
            plan = pending.pop(0)
            task = asyncio.create_task(
                run_one_test(
                    plan,
                    workdir=workdir,
                    log_dir=log_dir,
                    timeout_s=timeout_s,
                    use_shell=use_shell,
                    base_env=base_env,
                    set_thread_env_vars=set_thread_env_vars,
                    preserve_thread_env=preserve_thread_env,
                    print_start=print_start,
                    mem_total_bytes=last_known_total,
                    mem_available_bytes=last_known_available,
                )
            )
            running[task] = plan

        if not running:
            # Defensive fallback, not expected in normal operation: try_gate()
            # always admits a launch into an empty pool unconditionally, so
            # reaching here with pending tests left would mean that invariant
            # broke. Don't crash the run over a scheduling bug -- wait a
            # moment and retry -- but say so loudly, since it would indicate
            # something worth investigating rather than ordinary memory
            # pressure (which is already reported from inside try_gate).
            print(
                "[WAITING] scheduler made no progress this cycle although the "
                "pool is empty (unexpected); retrying",
                flush=True,
            )
            await asyncio.sleep(memory_check_interval if memory_gate_enabled else 1.0)
            continue

        if pending and len(running) < jobs and memory_gate_enabled and last_check_time is not None:
            wait_timeout = max(0.1, memory_check_interval - (time.monotonic() - last_check_time))
        else:
            wait_timeout = None

        done, _pending_tasks = await asyncio.wait(
            running.keys(), timeout=wait_timeout, return_when=asyncio.FIRST_COMPLETED
        )

        for task in done:
            plan = running.pop(task)
            result = await task
            results.append(result)
            status = "OK" if result.matched_reference else "MISMATCH"
            print(
                f"[{status}] #{result.index:03d} line {result.line_no}: "
                f"expected {result.expected}, actual {result.actual}, "
                f"exit={result.exit_code}, {result.elapsed_s:.1f}s",
                flush=True,
            )

    results.sort(key=lambda r: r.index)
    return results


def write_reports(results: List[TestResult], report_prefix: Path) -> tuple[Path, Path, Path, Path, Path, Path]:
    json_path = report_prefix.with_suffix(".json")
    csv_path = report_prefix.with_suffix(".csv")
    to_address_txt_path = report_prefix.with_name(report_prefix.name + "_to_address.txt")
    to_address_csv_path = report_prefix.with_name(report_prefix.name + "_to_address.csv")
    actual_failed_txt_path = report_prefix.with_name(report_prefix.name + "_actual_failed.txt")
    actual_failed_csv_path = report_prefix.with_name(report_prefix.name + "_actual_failed.csv")

    with json_path.open("w", encoding="utf-8") as f:
        json.dump([asdict(r) for r in results], f, indent=2)
        f.write("\n")

    fieldnames = [
        "index",
        "line_no",
        "expected",
        "actual",
        "matched_reference",
        "exit_code",
        "timed_out",
        "elapsed_s",
        "command",
        "log_file",
        "np_value",
        "nt_value",
        "np_nt_details",
        "mem_total_bytes",
        "mem_available_bytes_at_launch",
        "mem_available_fraction_at_launch",
        "last_pass",
        "last_pass_line_no",
    ]

    def write_csv(path: Path, rows: List[TestResult]) -> None:
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in rows:
                row = asdict(r)
                row["elapsed_s"] = f"{r.elapsed_s:.6f}"
                if r.mem_available_fraction_at_launch is not None:
                    row["mem_available_fraction_at_launch"] = f"{r.mem_available_fraction_at_launch:.6f}"
                writer.writerow(row)

    def write_txt(path: Path, title: str, rows: List[TestResult]) -> None:
        with path.open("w", encoding="utf-8") as f:
            f.write(title + "\n")
            f.write("=" * len(title) + "\n")
            f.write(f"Count: {len(rows)}\n\n")
            if not rows:
                f.write("None.\n")
                return
            for r in rows:
                if r.mem_total_bytes and r.mem_available_bytes_at_launch is not None:
                    mem_line = (
                        f"  memory:   {r.mem_available_fraction_at_launch:.0%} available at launch "
                        f"({format_bytes(r.mem_available_bytes_at_launch)} / {format_bytes(r.mem_total_bytes)})\n"
                    )
                else:
                    mem_line = "  memory:   not tracked (--no-memory-gate)\n"
                f.write(
                    f"#{r.index:03d} line {r.line_no}\n"
                    f"  expected: {r.expected}\n"
                    f"  actual:   {r.actual}\n"
                    f"  exit:     {r.exit_code}\n"
                    f"  timeout:  {r.timed_out}\n"
                    f"  elapsed:  {r.elapsed_s:.3f} s\n"
                    f"  {r.np_nt_details}\n"
                    f"{mem_line}"
                    f"  log:      {r.log_file}\n"
                    f"  command:  {r.command}\n\n"
                )

    to_address = [r for r in results if not r.matched_reference]
    actual_failed = [r for r in results if r.actual == "F"]

    write_csv(csv_path, results)
    write_csv(to_address_csv_path, to_address)
    write_csv(actual_failed_csv_path, actual_failed)
    write_txt(to_address_txt_path, "Tests to address: actual result differs from reference P/F", to_address)
    write_txt(actual_failed_txt_path, "Actual failed commands: exit code was nonzero or timed out", actual_failed)

    return (
        json_path,
        csv_path,
        to_address_txt_path,
        to_address_csv_path,
        actual_failed_txt_path,
        actual_failed_csv_path,
    )


def get_git_commit_id(workdir: Path) -> str:
    """Return the current git commit hash for the requested working tree."""
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(workdir),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            "could not determine the current git commit id; run from inside a "
            "git work tree or pass --commit-id explicitly"
        ) from exc

    commit_id = completed.stdout.strip()
    if not commit_id:
        raise RuntimeError("git rev-parse HEAD returned an empty commit id")
    return commit_id


def update_last_pass_entries(test_file: Path, tests: List[TestCase], results: List[TestResult], commit_id: str) -> int:
    """
    Update the test-list 'last pass:' metadata for tests whose command actually
    exited 0 on this run.

    This tracks the literal exit status of each command (actual == "P"),
    independent of the test's P/F marker in the list. An F-marked test that
    unexpectedly exits 0 still gets its 'last pass:' commit id updated, since
    the command did in fact pass; it will also appear in the to-address report
    as an unexpected pass so the marker itself can be revisited separately.
    """
    passed_by_line = {r.line_no for r in results if r.actual == "P"}
    tests_by_line = {t.line_no: t for t in tests}

    update_existing_line: dict[int, str] = {}
    insert_after_line: dict[int, str] = {}

    for line_no in passed_by_line:
        test = tests_by_line[line_no]
        if test.last_pass_line_no is None:
            insert_after_line[test.line_no] = commit_id
        else:
            update_existing_line[test.last_pass_line_no] = commit_id

    if not update_existing_line and not insert_after_line:
        return 0

    original_lines = test_file.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    new_lines: List[str] = []

    for line_no, raw in enumerate(original_lines, start=1):
        if line_no in update_existing_line:
            newline = "\n" if raw.endswith("\n") else ""
            stripped = raw.rstrip("\n")
            match = LAST_PASS_RE.match(stripped)
            if match:
                indent_match = re.match(r"^(\s*)", match.group("prefix"))
                indent = indent_match.group(1) if indent_match else ""
                prefix = f"{indent}last pass: "
            else:
                prefix = "last pass: "
            new_lines.append(f"{prefix}{update_existing_line[line_no]}{newline}")
        else:
            new_lines.append(raw)

        if line_no in insert_after_line:
            new_lines.append(f"last pass: {insert_after_line[line_no]}\n")

    tmp_path = test_file.with_name(test_file.name + ".tmp")
    tmp_path.write_text("".join(new_lines), encoding="utf-8")
    tmp_path.replace(test_file)
    return len(update_existing_line) + len(insert_after_line)


def print_final_summary(
    results: List[TestResult],
    *,
    json_path: Path,
    csv_path: Path,
    to_address_txt_path: Path,
    to_address_csv_path: Path,
    actual_failed_txt_path: Path,
    actual_failed_csv_path: Path,
) -> None:
    n = len(results)
    n_expected_pass = sum(r.expected == "P" for r in results)
    n_expected_fail = sum(r.expected == "F" for r in results)
    n_actual_pass = sum(r.actual == "P" for r in results)
    n_actual_fail = sum(r.actual == "F" for r in results)
    n_match = sum(r.matched_reference for r in results)
    n_mismatch = n - n_match
    n_timeout = sum(r.timed_out for r in results)
    to_address = [r for r in results if not r.matched_reference]
    actual_failed = [r for r in results if r.actual == "F"]

    print("\nSummary")
    print("-------")
    print(f"Total tests:        {n}")
    print(f"Expected pass/fail: {n_expected_pass}/{n_expected_fail}")
    print(f"Actual pass/fail:   {n_actual_pass}/{n_actual_fail}")
    print(f"Matched reference:  {n_match}")
    print(f"Mismatches:         {n_mismatch}")
    print(f"Timeouts:           {n_timeout}")
    print(f"JSON report:        {json_path}")
    print(f"CSV report:         {csv_path}")
    print(f"To-address text:    {to_address_txt_path}")
    print(f"To-address CSV:     {to_address_csv_path}")
    print(f"Actual-failed text: {actual_failed_txt_path}")
    print(f"Actual-failed CSV:  {actual_failed_csv_path}")

    if to_address:
        print("\nTests to address: actual result differs from reference P/F")
        for r in to_address:
            reason = "unexpected failure" if r.expected == "P" and r.actual == "F" else "unexpected pass"
            print(
                f"  #{r.index:03d} line {r.line_no}: {reason}; "
                f"expected {r.expected}, actual {r.actual}, exit={r.exit_code}, log={r.log_file}\n"
                f"      {r.command}"
            )
    else:
        print("\nTests to address: none; all actual results match the reference P/F markers.")

    # This second list is useful when the reference file intentionally contains
    # expected-failure tests. It shows every command that physically failed,
    # even when that failure matched an expected F marker.
    if actual_failed:
        print("\nAll actual failed commands, including expected F tests:")
        for r in actual_failed:
            marker = "expected" if r.expected == "F" else "unexpected"
            print(
                f"  #{r.index:03d} line {r.line_no}: {marker} actual F; "
                f"exit={r.exit_code}, log={r.log_file}\n"
                f"      {r.command}"
            )


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run P/F commands from a test list concurrently and compare with reference pass/fail markers."
    )
    parser.add_argument("test_file", help="Path to the test list file")
    parser.add_argument(
        "jobs_positional",
        nargs="?",
        type=int,
        help="Number of tests to run concurrently; can also be set with -j/--jobs",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        help="Number of tests to run concurrently; overrides the positional jobs value",
    )
    # Memory-aware scheduling controls. These are runner-level options rather
    # than test-list syntax so existing validation lists keep working
    # unchanged, and no test's -np/-nt is ever rewritten because of them.
    parser.add_argument(
        "--min-free-memory-fraction",
        type=float,
        default=0.5,
        help=(
            "Before starting another concurrent test, require more than this "
            "fraction of total physical memory (from /proc/meminfo's "
            "MemAvailable) to still be free. Default: 0.5 (more than half)."
        ),
    )
    parser.add_argument(
        "--memory-check-interval",
        type=float,
        default=60.0,
        help=(
            "Seconds to wait after a test launch attempt before considering "
            "another (giving a just-started test time to reach its "
            "steady-state memory footprint), and how often to recheck while "
            "waiting for memory to free up. A launch into an empty pool "
            "(the first test of the run, or any later test that starts once "
            "nothing else is running) is exempt and always starts "
            "immediately, since there is no in-flight test to pace against. "
            "Default: 60."
        ),
    )
    parser.add_argument(
        "--no-memory-gate",
        action="store_true",
        help=(
            "Disable the physical-memory ramp-up gate; launch tests up to -j "
            "immediately with no pacing or memory checks, like plain "
            "job-count-limited concurrency. Also needed on systems without "
            "/proc/meminfo."
        ),
    )
    parser.add_argument(
        "--default-np",
        type=int,
        default=1,
        help=(
            "MPI rank count assumed for commands without -np/--np; used only "
            "for reporting and --set-thread-env, not for scheduling. Default: 1"
        ),
    )
    parser.add_argument(
        "--default-nt",
        type=int,
        default=1,
        help=(
            "Thread count assumed for commands without -nt/--nt; used only "
            "for reporting and --set-thread-env, not for scheduling. Default: 1"
        ),
    )
    parser.add_argument(
        "--set-thread-env",
        action="store_true",
        help=(
            "Set OMP_NUM_THREADS/OPENBLAS_NUM_THREADS/MKL_NUM_THREADS/NUMEXPR_NUM_THREADS "
            "for each launched command to its own parsed -nt value"
        ),
    )
    parser.add_argument(
        "--overwrite-thread-env",
        action="store_true",
        help="With --set-thread-env, overwrite existing thread-count environment variables",
    )
    parser.add_argument(
        "--workdir",
        default=".",
        help="Working directory used to launch all commands; default: current directory",
    )
    parser.add_argument(
        "--log-dir",
        default="test_runner_logs",
        help="Directory for per-test stdout/stderr logs; default: test_runner_logs",
    )
    parser.add_argument(
        "--report-prefix",
        default="test_runner_report",
        help="Prefix for JSON/CSV reports; default: test_runner_report",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="Optional timeout per test in seconds; timeout is treated as actual F",
    )
    parser.add_argument(
        "--no-shell",
        action="store_true",
        help="Execute commands without a shell using shlex.split(); default uses the shell",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse and list tests without executing them",
    )
    parser.add_argument(
        "--no-start-progress",
        action="store_true",
        help=(
            "Do not print [START]/[WAITING] lines when commands launch or the "
            "memory gate blocks the next launch; completion [OK]/[MISMATCH] "
            "lines are still printed"
        ),
    )
    parser.add_argument(
        "--update-last-pass",
        action="store_true",
        help=(
            "After running tests, update or insert the 'last pass:' git commit id "
            "for each test whose command actually exited 0 this run, regardless "
            "of its P/F marker"
        ),
    )
    parser.add_argument(
        "--commit-id",
        default=None,
        help=(
            "Commit id to write with --update-last-pass; default: current "
            "'git rev-parse HEAD' in --workdir"
        ),
    )

    args = parser.parse_args(argv)

    jobs = args.jobs if args.jobs is not None else args.jobs_positional
    if jobs is None:
        jobs = 1
    if jobs < 1:
        parser.error("jobs must be >= 1")
    if args.default_np < 1:
        parser.error("--default-np must be >= 1")
    if args.default_nt < 1:
        parser.error("--default-nt must be >= 1")
    if not (0.0 <= args.min_free_memory_fraction <= 1.0):
        parser.error("--min-free-memory-fraction must be between 0 and 1")
    if args.memory_check_interval < 0:
        parser.error("--memory-check-interval must be >= 0")

    memory_gate_enabled = not args.no_memory_gate

    test_file = Path(args.test_file).expanduser().resolve()
    workdir = Path(args.workdir).expanduser().resolve()
    log_dir = Path(args.log_dir).expanduser().resolve()
    report_prefix = Path(args.report_prefix).expanduser().resolve()

    try:
        tests = parse_test_file(test_file)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if not tests:
        print(f"No tests found in {test_file}")
        return 0

    # Fail fast, before running anything, if memory can't be tracked at all.
    mem_total_at_start: Optional[int] = None
    mem_available_at_start: Optional[int] = None
    if memory_gate_enabled:
        try:
            mem_total_at_start, mem_available_at_start = read_memory_totals()
        except RuntimeError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2

    print(f"Parsed {len(tests)} tests from {test_file}")
    print(f"Concurrent jobs: {jobs}")
    if memory_gate_enabled:
        print(
            f"Memory gate: start next test only when > {args.min_free_memory_fraction:.0%} "
            f"of physical memory is free, rechecked every {args.memory_check_interval:.0f}s "
            "after each launch attempt (launches into an empty pool are exempt)"
        )
        print(
            f"Physical memory: {format_bytes(mem_total_at_start)} total, "
            f"{format_bytes(mem_available_at_start)} available now "
            f"({mem_available_at_start / mem_total_at_start:.0%})"
        )
    else:
        print("Memory gate: disabled (--no-memory-gate); launching up to -j immediately")
    print(f"Workdir: {workdir}")
    print(f"Execution mode: {'no shell' if args.no_shell else 'shell'}")

    if args.dry_run:
        # Dry-run lists the tests and their parsed np/nt, but actual
        # concurrency now depends on runtime memory use, so it cannot be
        # previewed the way the old CPU-slot plan could be.
        print("\nDry run test list:")
        for t in tests:
            plan = build_test_plan(t, default_np=args.default_np, default_nt=args.default_nt)
            suffix = f"; last pass: {t.last_pass}" if t.last_pass else "; last pass: <none>"
            print(
                f"  #{t.index:03d} line {t.line_no}: expected {t.expected}: "
                f"{plan.np_nt_details}: {t.command}{suffix}"
            )
        if memory_gate_enabled:
            print(
                f"\nActual concurrency depends on runtime memory use and cannot be "
                f"previewed here: once a test is already running, the next one only "
                f"starts once > {args.min_free_memory_fraction:.0%} of physical memory "
                f"is free, checked every {args.memory_check_interval:.0f}s (a launch "
                f"into an idle pool always starts immediately)."
            )
        return 0

    try:
        results = asyncio.run(
            run_all_tests(
                tests,
                jobs=jobs,
                workdir=workdir,
                log_dir=log_dir,
                timeout_s=args.timeout,
                use_shell=not args.no_shell,
                default_np=args.default_np,
                default_nt=args.default_nt,
                set_thread_env_vars=args.set_thread_env,
                preserve_thread_env=not args.overwrite_thread_env,
                print_start=not args.no_start_progress,
                memory_gate_enabled=memory_gate_enabled,
                min_free_memory_fraction=args.min_free_memory_fraction,
                memory_check_interval=args.memory_check_interval,
            )
        )
    except KeyboardInterrupt:
        print("Interrupted by user", file=sys.stderr)
        return 130

    (
        json_path,
        csv_path,
        to_address_txt_path,
        to_address_csv_path,
        actual_failed_txt_path,
        actual_failed_csv_path,
    ) = write_reports(results, report_prefix)
    print_final_summary(
        results,
        json_path=json_path,
        csv_path=csv_path,
        to_address_txt_path=to_address_txt_path,
        to_address_csv_path=to_address_csv_path,
        actual_failed_txt_path=actual_failed_txt_path,
        actual_failed_csv_path=actual_failed_csv_path,
    )

    if args.update_last_pass:
        try:
            commit_id = args.commit_id or get_git_commit_id(workdir)
            n_updated = update_last_pass_entries(test_file, tests, results, commit_id)
        except Exception as exc:
            print(f"ERROR updating last-pass metadata: {exc}", file=sys.stderr)
            return 2
        print(
            f"Updated {n_updated} last-pass entr{'y' if n_updated == 1 else 'ies'} "
            f"in {test_file} to commit {commit_id}"
        )

    return 0 if all(r.matched_reference for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
