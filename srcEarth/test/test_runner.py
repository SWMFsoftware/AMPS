#!/usr/bin/env python3
"""
Concurrent pass/fail test runner.

Input file format:
  ! this is a comment

  P command expected to return exit code 0
  F command expected to return nonzero exit code

Blank lines are ignored.  Lines beginning with '!' after leading whitespace are
ignored.  The first non-whitespace token on a test line must be P or F; the rest
of the line is executed as the command.

By default commands are executed through the user's shell so that ordinary test
commands, environment variables, and shell wrappers work as written.  Use
--no-shell if the command lines should be parsed with shlex and executed without
an intermediate shell.

In addition to the full JSON/CSV reports, the runner writes:
  <report-prefix>_to_address.txt/.csv
      tests whose actual P/F result does not match the reference P/F marker;
      these are the tests that require investigation.
  <report-prefix>_actual_failed.txt/.csv
      all commands that returned actual F, including expected-failure tests.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import json
import os
import shlex
import signal
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, List, Optional


@dataclass(frozen=True)
class TestCase:
    index: int
    line_no: int
    expected: str  # "P" or "F"
    command: str


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


def parse_test_file(path: Path) -> List[TestCase]:
    tests: List[TestCase] = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()

            if not line:
                continue
            if line.startswith("!"):
                continue

            parts = line.split(maxsplit=1)
            if len(parts) != 2 or parts[0].upper() not in {"P", "F"}:
                raise ValueError(
                    f"Invalid test-list line {line_no}: expected 'P <command>' "
                    f"or 'F <command>', got: {raw.rstrip()}"
                )

            tests.append(
                TestCase(
                    index=len(tests) + 1,
                    line_no=line_no,
                    expected=parts[0].upper(),
                    command=parts[1].strip(),
                )
            )

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
    semaphore: asyncio.Semaphore,
    workdir: Path,
    log_dir: Path,
    timeout_s: Optional[float],
    use_shell: bool,
    env: dict,
) -> TestResult:
    async with semaphore:
        log_path = log_dir / safe_log_name(test)
        start = time.monotonic()
        exit_code: Optional[int] = None
        timed_out = False

        with log_path.open("wb") as log:
            header = (
                f"# Test {test.index} from line {test.line_no}\n"
                f"# Expected: {test.expected}\n"
                f"# Command: {test.command}\n"
                f"# Workdir: {workdir}\n"
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
        )


async def run_all_tests(
    tests: Iterable[TestCase],
    *,
    jobs: int,
    workdir: Path,
    log_dir: Path,
    timeout_s: Optional[float],
    use_shell: bool,
) -> List[TestResult]:
    semaphore = asyncio.Semaphore(jobs)
    log_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    tasks = [
        asyncio.create_task(
            run_one_test(
                test,
                semaphore=semaphore,
                workdir=workdir,
                log_dir=log_dir,
                timeout_s=timeout_s,
                use_shell=use_shell,
                env=env,
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
            f"exit={result.exit_code}, {result.elapsed_s:.1f}s"
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

    args = parser.parse_args(argv)

    jobs = args.jobs if args.jobs is not None else args.jobs_positional
    if jobs is None:
        jobs = 1
    if jobs < 1:
        parser.error("jobs must be >= 1")

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
    print(f"Workdir: {workdir}")
    print(f"Execution mode: {'no shell' if args.no_shell else 'shell'}")

    if args.dry_run:
        print("\nDry run test list:")
        for t in tests:
            print(f"  #{t.index:03d} line {t.line_no}: expected {t.expected}: {t.command}")
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

    return 0 if all(r.matched_reference for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
