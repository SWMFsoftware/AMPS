#!/usr/bin/env python3
"""Run and summarize the srcSEP3D parallel-diffusion / particle-mover tests.

This runner executes the focused test suites added or updated for the shared
parallel-diffusion library's input-section parser and mover-facing pointer,
the srcSEP3D schema-5 [parallel_diffusion] binding, and the deck-selected
runtime mover dispatch (no prepare-production step).  It parses each suite's
per-test PASS/FAIL/SKIP/ERROR lines, prints one summary, and lists every
failed test.  It never needs a configured or native AMPS build.

Exit status: 0 when every selected suite ran and no test failed, errored, or
was missing; 1 otherwise; 2 for a usage error.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

SRCSEP3D = Path(__file__).resolve().parents[1]
AMPS_ROOT = SRCSEP3D.parent
# BLDL3D05 writes its fixture under this ignored, generated directory.
BLDL_OUTPUT = SRCSEP3D / "test_output" / "pd_mover_runner" / "BLDL3D05"
LIBRARY_DIR = AMPS_ROOT / "src" / "models" / "sep_common" / "parallel_diffusion"

TEST_DESCRIPTIONS = """\
suites and tests (run from any directory; no AMPS build is needed):

  library       make -f makefile verify   (in src/models/sep_common/parallel_diffusion)
                Shared parallel-diffusion library: its whole base and advanced
                suites run; the tests added for the input-section parser and
                the mover-facing pointer are
                  PD13-UNCONFIGURED     pointer returns InvalidConfiguration and
                                        no value before a model is installed
                  PD13-SECTION-GRAMMAR  ParseSection: '#'/'!' comments, blank/
                                        CRLF lines, case rules, line provenance
                  PD13-SECTION-ERRORS   12 malformed sections rejected with
                                        typed status and their line number
                  PD13-BOUND-REGISTRY   one distinct bound evaluator per model
                  PD13-BOUND-DISPATCH   section-selected pointer, kappa=v*lambda/3,
                                        model switch, failed update unchanged
                  PD13-HOST-GATE        CheckHostInputAvailability per model

  stage1-build  make -C srcSEP3D -j16 test/stage1
                Builds the AMPS-free stage-1 test executable (no tests run).

  stage1        srcSEP3D/test/stage1 --test CFG3D17 --test R3D01   (cwd srcSEP3D)
                  CFG3D17  schema-5 [parallel_diffusion]: raw lines go to the
                           library parser ('!' comments, deck line numbers in
                           errors), mover-family cross-check, restart identity,
                           Runtime installs ActiveParallelDiffusion bound to the
                           selected model, unavailable host state rejected
                  R3D01    generic PIC defaults to legacy dispatch; the deck
                           input/sep3d.input selects pointer dispatch
                           (_PIC_PARTICLE_MOVER_LEGACY_SETTINGS_ OFF), no
                           prepare-production target, audit still verifies it;
                           Parker/focused callbacks; AMPS fatal trap

  bldl3d05      python3 srcSEP3D/test/run_tests.py --test BLDL3D05 --amps-source .
                  BLDL3D05 source and copied build/main makefiles resolve the
                           same paths; strict-production delegates to make amps
                           and the audit accepts the deck-appended pointer-mode
                           definition

  stage1-all    srcSEP3D/test/stage1 --all-tests   (optional; not in default set)
                Every routine stage-1 test (regression check).

examples:
  %(prog)s                          run the default suites
  %(prog)s --suite stage1-build --suite stage1
  %(prog)s --all --verbose --json report.json
"""


@dataclass
class TestResult:
    test_id: str
    status: str            # PASS, FAIL, SKIP, ERROR, or MISSING
    suite: str
    detail: str = ""
    # Location of this result's evidence: the saved complete output of its
    # suite (LOG_DIR/<suite>.log) and the 1-based line holding the result
    # line (0 when the result was synthesized, e.g. MISSING).  Empty when no
    # log directory is in use.  Aggregators report these for failed tests.
    log: str = ""
    line: int = 0


@dataclass
class Suite:
    name: str
    command: Sequence[str]
    cwd: Path
    expected: Sequence[str] = ()
    # IDs whose FAIL line is an intended negative control inside the suite.
    ignored: Sequence[str] = ()
    timeout_s: int = 1800


@dataclass
class SuiteRun:
    suite: Suite
    exit_code: Optional[int]
    duration_s: float
    output: str
    results: List[TestResult] = field(default_factory=list)


SUITES: Dict[str, Suite] = {
    "library": Suite(
        "library", ["make", "-f", "makefile", "verify"], LIBRARY_DIR,
        expected=("PD13-UNCONFIGURED", "PD13-SECTION-GRAMMAR",
                  "PD13-SECTION-ERRORS", "PD13-BOUND-REGISTRY",
                  "PD13-BOUND-DISPATCH", "PD13-HOST-GATE")),
    "stage1-build": Suite(
        "stage1-build", ["make", "-C", str(SRCSEP3D), "-j16", "test/stage1"],
        AMPS_ROOT),
    # stage1 opens fixtures relative to srcSEP3D, so it must run there.
    "stage1": Suite(
        "stage1", [str(SRCSEP3D / "test" / "stage1"),
                   "--test", "CFG3D17", "--test", "R3D01"], SRCSEP3D,
        expected=("CFG3D17", "R3D01")),
    "bldl3d05": Suite(
        "bldl3d05", [sys.executable, str(SRCSEP3D / "test" / "run_tests.py"),
                     "--test", "BLDL3D05", "--amps-source", str(AMPS_ROOT),
                     "--output-dir", str(BLDL_OUTPUT)], AMPS_ROOT,
        expected=("BLDL3D05",)),
    "stage1-all": Suite(
        "stage1-all", [str(SRCSEP3D / "test" / "stage1"), "--all-tests"],
        SRCSEP3D),
}
# Suites run when --suite/--all is not given (stage1-all is opt-in).
DEFAULT_SUITES = ("library", "stage1-build", "stage1", "bldl3d05")

STATUSES = ("PASS", "FAIL", "SKIP", "ERROR")
# Severity used to merge repeated lines for one ID (e.g. a sub-case loop).
SEVERITY = {"PASS": 0, "SKIP": 1, "MISSING": 2, "ERROR": 3, "FAIL": 4}

# The suites use four line formats:
#   [ID] STATUS ...         sep_test_registry, run_tests.py, step-1 failures
#   STATUS ID[:] ...        parallel_diffusion library, COEF suite
#   ID STATUS               step-1 outer summary lines
#   FAIL: ID message        COEF suite failure
PATTERNS = (
    re.compile(r"^\[(?P<id>[A-Z][A-Za-z0-9_.:-]*)\] (?P<status>PASS|FAIL|SKIP|ERROR)\b"),
    re.compile(r"^(?P<status>PASS|FAIL|SKIP|ERROR) (?P<id>[A-Z][A-Za-z0-9_.:/-]*)"),
    re.compile(r"^(?P<id>[A-Z][A-Z0-9_-]*) (?P<status>PASS|FAIL|SKIP|ERROR)$"),
    re.compile(r"^FAIL: (?P<id>[A-Z][A-Z0-9_-]*\d)\b"),
)


def parse_results(suite: Suite, output: str) -> List[TestResult]:
    """Extract one merged result per test ID from a suite's combined output.

    Each result records the 1-based output line it was taken from, so a
    failed test can be located inside the saved suite log."""
    merged: Dict[str, TestResult] = {}
    for line_number, raw in enumerate(output.splitlines(), start=1):
        line = raw.strip()
        for pattern in PATTERNS:
            match = pattern.match(line)
            if not match:
                continue
            test_id = match.group("id").rstrip(":")
            status = match.group("status") if "status" in match.groupdict() else "FAIL"
            if test_id in suite.ignored:
                break
            previous = merged.get(test_id)
            if previous is None or SEVERITY[status] > SEVERITY[previous.status]:
                merged[test_id] = TestResult(test_id, status, suite.name, line,
                                             line=line_number)
            break
    return list(merged.values())


def run_suite(suite: Suite, verbose: bool,
              log_dir: Optional[Path] = None) -> SuiteRun:
    print(f"==> {suite.name}: {' '.join(suite.command)}  (cwd {suite.cwd})",
          flush=True)
    start = time.monotonic()
    try:
        completed = subprocess.run(
            list(suite.command), cwd=str(suite.cwd), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, universal_newlines=True,
            timeout=suite.timeout_s)
        exit_code: Optional[int] = completed.returncode
        output = completed.stdout
    except subprocess.TimeoutExpired as error:
        exit_code = None
        output = (error.stdout or "") if isinstance(error.stdout, str) else ""
        output += f"\n[runner] timeout after {suite.timeout_s} s\n"
    except OSError as error:
        exit_code = None
        output = f"[runner] cannot start command: {error}\n"
    duration = time.monotonic() - start
    if verbose:
        print(output, end="" if output.endswith("\n") else "\n")

    # Save the suite's complete output so every result has a log location.
    log_path = ""
    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"{suite.name}.log"
        log_file.write_text(
            f"$ {' '.join(suite.command)}\n(cwd {suite.cwd})\n{output}"
            f"[exit {exit_code}, {duration:.1f} s]\n", encoding="utf-8")
        log_path = str(log_file)
    run = SuiteRun(suite, exit_code, duration, output)
    run.results = parse_results(suite, output)
    # Line numbers refer to the output, which starts after the 2-line header.
    for result in run.results:
        result.log = log_path
        result.line = result.line + 2 if (log_path and result.line) else result.line
    seen = {result.test_id for result in run.results}
    for test_id in suite.expected:
        if test_id not in seen:
            run.results.append(TestResult(
                test_id, "MISSING", suite.name,
                "expected test did not report a result", log=log_path))
    # A failing command must never look green: if no individual FAIL/ERROR
    # explains a nonzero exit (build error, crash, timeout), record one.  A
    # MISSING entry alone does not explain it, so the exit code is kept.
    failed = any(r.status in ("FAIL", "ERROR") for r in run.results)
    if exit_code != 0 and not failed:
        tail = " | ".join(output.strip().splitlines()[-3:])
        run.results.append(TestResult(
            f"{suite.name}:command", "ERROR", suite.name,
            f"exit code {exit_code}: {tail}", log=log_path))
    counts = count(run.results)
    print(f"    exit={exit_code} time={duration:.1f}s  " + format_counts(counts),
          flush=True)
    return run


def count(results: Sequence[TestResult]) -> Dict[str, int]:
    totals = {status: 0 for status in (*STATUSES, "MISSING")}
    for result in results:
        totals[result.status] += 1
    return totals


def format_counts(totals: Dict[str, int]) -> str:
    return " ".join(f"{status}={totals[status]}" for status in
                    (*STATUSES, "MISSING"))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=TEST_DESCRIPTIONS,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--suite", action="append", choices=list(SUITES),
                        help="run only this suite (repeatable; default: "
                             + ", ".join(DEFAULT_SUITES) + ")")
    parser.add_argument("--all", action="store_true",
                        help="run every suite, including stage1-all")
    parser.add_argument("--list", action="store_true",
                        help="list suites and their commands, then exit")
    parser.add_argument("--verbose", action="store_true",
                        help="print each suite's complete output")
    parser.add_argument("--json", type=Path, metavar="FILE",
                        help="also write a machine-readable report")
    parser.add_argument("--log-dir", type=Path, metavar="DIR",
                        help="save each suite's complete output to DIR/<suite>.log "
                             "and record log/line per result (default with "
                             "--json: <json name>-logs next to the JSON file)")
    args = parser.parse_args(argv)

    if args.list:
        for name, suite in SUITES.items():
            print(f"{name:13s} {' '.join(suite.command)}  (cwd {suite.cwd})")
            if suite.expected:
                print(f"{'':13s} expected: {', '.join(suite.expected)}")
        return 0

    selected = (list(SUITES) if args.all
                else args.suite or list(DEFAULT_SUITES))
    log_dir = args.log_dir
    if log_dir is None and args.json is not None:
        log_dir = args.json.with_name(args.json.stem + "-logs")
    runs = [run_suite(SUITES[name], args.verbose, log_dir) for name in selected]
    results = [result for run in runs for result in run.results]
    totals = count(results)
    failures = [r for r in results if r.status in ("FAIL", "ERROR", "MISSING")]

    print("\n=== srcSEP3D parallel-diffusion / mover test summary ===")
    for run in runs:
        print(f"  {run.suite.name:13s} {format_counts(count(run.results))}")
    print(f"  {'TOTAL':13s} {format_counts(totals)}")
    if failures:
        print(f"\nFailed tests ({len(failures)}):")
        for result in failures:
            print(f"  [{result.suite}] {result.test_id} {result.status}: "
                  f"{result.detail}")
            if result.log:
                print(f"      log: {result.log}"
                      + (f":{result.line}" if result.line else ""))
    else:
        print("\nFailed tests: none")

    if args.json:
        report = {
            "schema": "srcsep3d-pd-mover-runner-v1",
            "passed": not failures,
            "totals": totals,
            "suites": [{
                "name": run.suite.name,
                "command": list(run.suite.command),
                "cwd": str(run.suite.cwd),
                "exit_code": run.exit_code,
                "duration_s": round(run.duration_s, 3),
                "results": [vars(r) for r in run.results],
            } for run in runs],
        }
        args.json.write_text(json.dumps(report, indent=2) + "\n",
                             encoding="utf-8")
        print(f"JSON report: {args.json}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
