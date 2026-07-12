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
3. Run MPI/OpenMP AMPS tests concurrently while avoiding CPU oversubscription.
   This is important for commands such as

       srcEarth/test/C3/run_C3.py -np 4 -nt 16

   where one command can already consume 64 CPU execution slots.
4. Keep per-test provenance by updating each test's ``last pass:`` commit id
   after successful command execution.

Scheduling algorithm
--------------------
The runner uses two independent throttles:

* Command throttle: ``-j/--jobs`` is a traditional asyncio semaphore limiting
  how many test commands can be active at once.  This preserves the old behavior
  and prevents launching an unbounded number of subprocesses.

* Resource throttle: in ``--resource-mode auto`` each command is assigned an
  estimated slot weight.  The estimate is parsed from common AMPS/MPI wrapper
  options, especially ``-np``/``--np`` and ``-nt``/``--nt``.  A command with
  ``-np 4 -nt 16`` requests ``4*16 = 64`` slots.  The weighted slot limiter only
  launches a command when the current free-slot count is large enough.

The two throttles are intentionally combined.  For example, ``-j 8 --max-slots
128`` allows up to eight active commands, but only if their estimated total CPU
slot use stays below 128.  This permits safe concurrency among small tests while
preventing several large MPI/OpenMP jobs from starting together and failing due
to thread/process oversubscription.

Slot detection and fallbacks
----------------------------
If ``--max-slots`` is not supplied, the runner tries to infer the allocation from
common scheduler variables such as ``SLURM_CPUS_ON_NODE``, ``PBS_NP``, ``NSLOTS``,
or from the number of entries in ``PBS_NODEFILE``.  If no scheduler allocation is
visible, it falls back to ``os.cpu_count()``.  The environment variable
``TEST_RUNNER_MAX_SLOTS`` can be used to override this detection without changing
command lines.

The slot estimator is deliberately conservative and non-invasive: it never
rewrites a test command.  If a command does not contain recognizable MPI/thread
options, ``--default-np`` and ``--default-nt`` are used, both defaulting to 1.
If an individual command requests more slots than the configured pool, it is
allowed to acquire the whole pool and run alone rather than deadlocking forever.

Thread environment policy
-------------------------
By default the runner does not alter ``OMP_NUM_THREADS`` or BLAS thread-count
environment variables.  Passing ``--set-thread-env`` sets those variables for
each command to the estimated ``-nt`` value, which can reduce hidden nested
threading in libraries.  Existing values are preserved unless
``--overwrite-thread-env`` is also supplied.

Result handling
---------------
The process exit status is converted to actual P/F.  Timeouts and launch errors
are treated as actual F.  The full JSON/CSV reports include the slot estimate and
acquired slot count for each test.  Two triage reports are also written:

  <report-prefix>_to_address.txt/.csv
      tests whose actual P/F result does not match the reference P/F marker;
      these are the tests that require investigation.
  <report-prefix>_actual_failed.txt/.csv
      all commands that returned actual F, including expected-failure tests.

Use ``--update-last-pass`` to rewrite the ``last pass:`` line of each test whose
command actually exited 0 on this run.  Existing metadata lines are updated in
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
    # requested_slots is the raw estimate from the command line, while
    # effective_slots is the number actually acquired from the slot limiter.
    # They differ only when a single command asks for more slots than the
    # configured pool; in that case the command receives the full pool and runs
    # alone.
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
    requested_slots: int = 1
    effective_slots: int = 1
    slot_details: str = ""
    last_pass: Optional[str] = None
    last_pass_line_no: Optional[int] = None


LAST_PASS_RE = re.compile(r"^(?P<prefix>\s*last\s+pass\s*:\s*)(?P<value>.*)$", re.IGNORECASE)


class WeightedSlotLimiter:
    """Async weighted semaphore used to cap total estimated CPU slots.

    ``asyncio.Semaphore`` only models one identical token per job.  That is not
    enough for the AMPS validation suite because two commands can be very
    different in cost: one may be a quick Python-only check, while another may
    launch ``mpirun -np 4`` with 16 OpenMP threads per rank.  This class models
    a pool of weighted slots instead.  Each test requests N slots; the coroutine
    blocks until at least N slots are available, then subtracts N from the pool.

    The class is intentionally minimal: it lives in this single script, uses one
    ``asyncio.Condition`` for all waiters, and relies on the caller to release
    exactly the number of slots returned by ``acquire``.
    """

    def __init__(self, capacity: int):
        if capacity < 1:
            raise ValueError("slot capacity must be >= 1")
        self.capacity = capacity
        self._available = capacity
        self._cond = asyncio.Condition()

    async def acquire(self, slots: int) -> int:
        """Reserve up to ``slots`` CPU slots and return the acquired weight.

        A command can request more slots than ``capacity`` if, for example, the
        test list says ``-np 8 -nt 32`` but the current developer laptop exposes
        only 16 CPUs.  Refusing to run such a command would be inconvenient, and
        waiting for 256 free slots would deadlock.  Therefore a too-large command
        is capped to ``capacity`` and runs exclusively.
        """
        effective_slots = min(max(1, int(slots)), self.capacity)
        async with self._cond:
            # Condition.wait() is used in a loop because multiple tests may wake
            # up when one job finishes; only the first ones that fit should start.
            while self._available < effective_slots:
                await self._cond.wait()
            self._available -= effective_slots
        return effective_slots

    async def release(self, slots: int) -> None:
        """Return previously acquired slots and wake blocked launchers."""
        effective_slots = min(max(1, int(slots)), self.capacity)
        async with self._cond:
            self._available += effective_slots
            # Clamp defensively so a future bug in caller release accounting does
            # not leave the limiter with more slots than its configured capacity.
            if self._available > self.capacity:
                self._available = self.capacity
            self._cond.notify_all()


def _parse_positive_int(text: str) -> Optional[int]:
    """Parse an integer-like value and reject zero/negative values.

    Scheduler variables sometimes arrive as strings and some wrapper scripts pass
    values that look like floats, e.g. ``"4.0"``.  Accepting ``int(float(x))``
    makes the parser tolerant without allowing invalid nonpositive slot counts.
    """
    try:
        value = int(float(str(text)))
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


def _option_values(argv: List[str], names: Iterable[str]) -> List[str]:
    """Return values for options that can appear as ``--x v`` or ``--x=v``.

    This helper intentionally does not use argparse because the command belongs
    to an arbitrary test runner, not to this script.  It performs a lightweight
    scan that is sufficient for the AMPS wrappers and direct ``mpirun`` commands
    in the validation list.
    """
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


def estimate_test_slots(command: str, default_np: int, default_nt: int) -> tuple[int, int, int, str]:
    """Estimate MPI ranks, threads per rank, and total CPU slots for one test.

    The validation list can contain heterogeneous commands: Python wrappers,
    direct ``mpirun`` invocations, shell fragments, or serial tools.  The runner
    needs a cheap estimate before launching the command.  The convention used by
    AMPS tests is that MPI rank count is carried by ``-np``/``--np`` or related
    aliases, and OpenMP/thread count is carried by ``-nt``/``--nt`` or specific
    AMPS options such as ``--mode3d-threads``.

    The estimate is

        slots = max(1, np * nt)

    where missing values fall back to ``--default-np`` and ``--default-nt``.  The
    function returns both the numeric estimate and a short string that is copied
    into logs/reports so failed tests can be triaged together with their resource
    assumptions.
    """
    try:
        # shlex gives shell-like tokenization without executing the command.  It
        # correctly handles quoted paths/arguments in most test-list entries.
        argv = shlex.split(command, posix=True)
    except ValueError:
        # If a command has unmatched quotes, let the test itself fail later, but
        # still produce a conservative slot estimate instead of aborting parsing.
        argv = command.split()

    # Recognize both direct mpirun-style rank options and wrapper aliases.
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

    slots = max(1, np_value * nt_value)
    details = f"np={np_value}, nt={nt_value}, slots={slots}"
    return np_value, nt_value, slots, details


def detect_default_max_slots() -> int:
    """Pick a conservative total slot count for the current batch run.

    Priority order:
    1. Explicit ``TEST_RUNNER_MAX_SLOTS`` override.
    2. Common scheduler environment variables.
    3. ``PBS_NODEFILE`` line count, which is common on Pleiades-style systems.
    4. Local CPU count as a last resort for workstation/developer use.
    """
    for name in (
        "TEST_RUNNER_MAX_SLOTS",
        "SLURM_CPUS_ON_NODE",
        "SLURM_CPUS_PER_TASK",
        "PBS_NP",
        "LSB_DJOB_NUMPROC",
        "NSLOTS",
    ):
        value = os.environ.get(name)
        parsed = _parse_positive_int(value) if value else None
        if parsed is not None:
            return parsed

    # PBS_NODEFILE often contains one line per allocated CPU slot.
    pbs_nodefile = os.environ.get("PBS_NODEFILE")
    if pbs_nodefile:
        try:
            with open(pbs_nodefile, "r", encoding="utf-8", errors="replace") as f:
                n = sum(1 for line in f if line.strip())
            if n > 0:
                return n
        except OSError:
            pass

    return os.cpu_count() or 1


def set_thread_env(base_env: dict, nt: int, only_if_missing: bool) -> dict:
    """Return environment for one command with optional thread-count controls.

    AMPS wrappers normally receive ``-nt`` explicitly, but imported numerical
    libraries can still start their own thread pools.  Setting these variables is
    optional because some HPC environments already configure them carefully.
    ``only_if_missing`` preserves existing settings unless the user requested
    ``--overwrite-thread-env``.
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
                # Section comments are not metadata.  Reset the pending command
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


async def run_one_test(
    test: TestCase,
    *,
    job_semaphore: asyncio.Semaphore,
    slot_limiter: Optional[WeightedSlotLimiter],
    workdir: Path,
    log_dir: Path,
    timeout_s: Optional[float],
    use_shell: bool,
    base_env: dict,
    default_np: int,
    default_nt: int,
    set_thread_env_vars: bool,
    preserve_thread_env: bool,
) -> TestResult:
    """Run one command after passing command-count and slot-count gates.

    The order is intentional: first acquire the simple job semaphore, then the
    weighted slot limiter.  This keeps the number of coroutines waiting on CPU
    slots bounded by ``-j`` while still preventing overcommit among the active
    candidates.  Slot acquisition happens before the subprocess is created, so a
    test never starts unless its estimated resources are available.
    """
    _np, nt, requested_slots, slot_details = estimate_test_slots(test.command, default_np, default_nt)
    effective_slots = requested_slots

    async with job_semaphore:
        if slot_limiter is not None:
            effective_slots = await slot_limiter.acquire(requested_slots)
        try:
            log_path = log_dir / safe_log_name(test)
            start = time.monotonic()
            exit_code: Optional[int] = None
            timed_out = False
            # Build a per-command environment after the slot estimate is known.
            # The base environment is never mutated, so concurrently running tests
            # cannot accidentally change each other's thread settings.
            env = set_thread_env(base_env, nt, preserve_thread_env) if set_thread_env_vars else base_env

            with log_path.open("wb") as log:
                header = (
                    f"# Test {test.index} from line {test.line_no}\n"
                    f"# Expected: {test.expected}\n"
                    f"# Command: {test.command}\n"
                    f"# Workdir: {workdir}\n"
                    f"# Resource estimate: {slot_details}; acquired={effective_slots}\n"
                    f"# Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                    f"{'=' * 80}\n"
                )
                log.write(header.encode("utf-8", errors="replace"))
                log.flush()

                try:
                    if use_shell:
                        # start_new_session=True makes timeout cleanup kill the
                        # whole process group, including mpirun-launched children
                        # when they remain in the same session.
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
                requested_slots=requested_slots,
                effective_slots=effective_slots,
                slot_details=slot_details,
                last_pass=test.last_pass,
                last_pass_line_no=test.last_pass_line_no,
            )
        finally:
            # Release slots even if command launch fails, the test times out, or
            # the Python wrapper raises.  Without this finally block, one failed
            # launch could starve the remaining queued tests.
            if slot_limiter is not None:
                await slot_limiter.release(effective_slots)


async def run_all_tests(
    tests: Iterable[TestCase],
    *,
    jobs: int,
    workdir: Path,
    log_dir: Path,
    timeout_s: Optional[float],
    use_shell: bool,
    resource_mode: str,
    max_slots: int,
    default_np: int,
    default_nt: int,
    set_thread_env_vars: bool,
    preserve_thread_env: bool,
) -> List[TestResult]:
    """Schedule all tests concurrently and return results in list order.

    All tasks are created immediately so asyncio can keep the pipeline full.
    Actual subprocess launch is gated inside ``run_one_test`` by both the job
    semaphore and, when enabled, the weighted slot limiter.  Results are printed
    as soon as each test finishes, then sorted by test index before being written
    to reports.
    """
    job_semaphore = asyncio.Semaphore(jobs)
    slot_limiter = WeightedSlotLimiter(max_slots) if resource_mode == "auto" else None
    log_dir.mkdir(parents=True, exist_ok=True)

    # Capture the parent environment once.  Each test gets either this dict or a
    # copied/thread-limited variant derived from it.
    base_env = os.environ.copy()
    tasks = [
        asyncio.create_task(
            run_one_test(
                test,
                job_semaphore=job_semaphore,
                slot_limiter=slot_limiter,
                workdir=workdir,
                log_dir=log_dir,
                timeout_s=timeout_s,
                use_shell=use_shell,
                base_env=base_env,
                default_np=default_np,
                default_nt=default_nt,
                set_thread_env_vars=set_thread_env_vars,
                preserve_thread_env=preserve_thread_env,
            )
        )
        for test in tests
    ]

    results: List[TestResult] = []
    for task in asyncio.as_completed(tasks):
        result = await task
        results.append(result)
        status = "OK" if result.matched_reference else "MISMATCH"
        print(
            f"[{status}] #{result.index:03d} line {result.line_no}: "
            f"expected {result.expected}, actual {result.actual}, "
            f"exit={result.exit_code}, {result.elapsed_s:.1f}s, "
            f"slots={result.effective_slots}/{result.requested_slots}"
        )
        sys.stdout.flush()

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
        "log_file",
        "requested_slots",
        "effective_slots",
        "slot_details",
        "last_pass",
        "last_pass_line_no",
        "command",
    ]

    def write_csv(path: Path, rows: List[TestResult]) -> None:
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in rows:
                row = asdict(r)
                row["elapsed_s"] = f"{r.elapsed_s:.6f}"
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
                f.write(
                    f"#{r.index:03d} line {r.line_no}\n"
                    f"  expected: {r.expected}\n"
                    f"  actual:   {r.actual}\n"
                    f"  exit:     {r.exit_code}\n"
                    f"  timeout:  {r.timed_out}\n"
                    f"  elapsed:  {r.elapsed_s:.3f} s\n"
                    f"  slots:    {r.effective_slots}/{r.requested_slots} ({r.slot_details})\n"
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
    independent of the test's P/F marker in the list.  An F-marked test that
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
    # expected-failure tests.  It shows every command that physically failed,
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
    # Resource-aware scheduling controls.  These are intentionally runner-level
    # options rather than test-list syntax so existing validation lists keep
    # working unchanged.
    parser.add_argument(
        "--resource-mode",
        choices=["auto", "off"],
        default="auto",
        help=(
            "Resource-aware scheduling mode. auto reserves estimated CPU slots "
            "from -np/-nt before launching each command; off keeps the old pure "
            "command-count behavior. Default: auto"
        ),
    )
    parser.add_argument(
        "--max-slots",
        type=int,
        default=None,
        help=(
            "Maximum total CPU slots to run concurrently in --resource-mode auto; "
            "default: TEST_RUNNER_MAX_SLOTS, scheduler allocation variables, "
            "PBS_NODEFILE count, or os.cpu_count()"
        ),
    )
    parser.add_argument(
        "--default-np",
        type=int,
        default=1,
        help="MPI rank estimate for commands that do not contain -np/--np; default: 1",
    )
    parser.add_argument(
        "--default-nt",
        type=int,
        default=1,
        help="thread estimate for commands that do not contain -nt/--nt; default: 1",
    )
    parser.add_argument(
        "--set-thread-env",
        action="store_true",
        help=(
            "Set OMP_NUM_THREADS/OPENBLAS_NUM_THREADS/MKL_NUM_THREADS/NUMEXPR_NUM_THREADS "
            "for each launched command to the estimated -nt value"
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
    # Resolve the total slot budget once at startup.  All test commands in this
    # invocation share the same budget.
    max_slots = args.max_slots if args.max_slots is not None else detect_default_max_slots()
    if max_slots < 1:
        parser.error("--max-slots must be >= 1")

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

    print(f"Parsed {len(tests)} tests from {test_file}")
    print(f"Concurrent jobs: {jobs}")
    print(f"Resource mode: {args.resource_mode}")
    if args.resource_mode == "auto":
        print(f"Max CPU slots: {max_slots}")
    print(f"Workdir: {workdir}")
    print(f"Execution mode: {'no shell' if args.no_shell else 'shell'}")

    if args.dry_run:
        # Dry-run is also useful for checking the scheduler's interpretation of
        # resource requests before launching expensive MPI jobs.
        print("\nDry run test list:")
        for t in tests:
            _np, _nt, slots, slot_details = estimate_test_slots(t.command, args.default_np, args.default_nt)
            suffix = f"; last pass: {t.last_pass}" if t.last_pass else "; last pass: <none>"
            print(
                f"  #{t.index:03d} line {t.line_no}: expected {t.expected}: "
                f"slots={slots} ({slot_details}): {t.command}{suffix}"
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
                resource_mode=args.resource_mode,
                max_slots=max_slots,
                default_np=args.default_np,
                default_nt=args.default_nt,
                set_thread_env_vars=args.set_thread_env,
                preserve_thread_env=not args.overwrite_thread_env,
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
