#!/usr/bin/env python3
"""Run and summarize the srcSEP parallel-diffusion / particle-mover tests.

This runner executes the focused, AMPS-independent test suites added for the
shared parallel-diffusion coefficient library binding and the schema-4
``run.particle_mover`` key, parses each suite's per-test PASS/FAIL/SKIP/ERROR
lines, prints one summary, and lists every failed test.  It never needs a
configured or native AMPS build.

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

SRCSEP = Path(__file__).resolve().parents[1]
AMPS_ROOT = SRCSEP.parent
LIBRARY_DIR = AMPS_ROOT / "src" / "models" / "sep_common" / "parallel_diffusion"

TEST_DESCRIPTIONS = """\
suites and tests (run from any directory; no AMPS build is needed):

  library   make -f makefile verify   (in src/models/sep_common/parallel_diffusion)
            Shared parallel-diffusion library: its whole base and advanced
            suites run; the tests added for the input-section parser and the
            mover-facing pointer are
              PD13-UNCONFIGURED     pointer returns InvalidConfiguration and no
                                    value before any model is installed
              PD13-SECTION-GRAMMAR  ParseSection: '#'/'!' comments, blank/CRLF
                                    lines, case rules, line provenance; result
                                    equals an independently built configuration
              PD13-SECTION-ERRORS   12 malformed sections rejected with typed
                                    status and their line number, output untouched
              PD13-BOUND-REGISTRY   one distinct bound evaluator per model ID (16)
              PD13-BOUND-DISPATCH   section-selected pointer, kappa=v*lambda/3 vs
                                    independent v, model switch, failed update
              PD13-HOST-GATE        CheckHostInputAvailability accepts/rejects
                                    each model by the host's declared inputs

  binding   make -C srcSEP test-parallel-diffusion-binding-unit
            srcSEP schema-4 input path (production sources, no AMPS/MPI):
              PDB01  [parallel_diffusion] body stored verbatim with line numbers
                     and key case; schema 2/3 and duplicate sections rejected
              PDB02  section/CLI/mover cross-check and host gate reject every
                     inconsistent or unsupported selection; nothing installed
              PDB03  parallel-diffusion-library registry provider, Parker only
              PDB04  installed models through ActiveParallelDiffusion match
                     closed forms (constant lambda; rigidity/radius power law)
              PDB05  [run] particle_mover: canonical names only, schema 4 only,
                     no duplicates, optional, part of the startup fingerprint

  cli       make -C srcSEP test-cli-unit
            Production CLI parser (step-1 suite CLI01-CLI06, REGISTRY, REFINE,
            HIDDEN, REPORT).  New test:
              CLI06  mover precedence --particle-mover > run.particle_mover >
                     default fte-dmumu, and deferred mover-dependent validation
                     when --input is given without --particle-mover
            The suite deliberately prints "[HIDDEN01] FAIL" as a negative
            control; the runner does not count it as a failure.

  coefficients  make -C srcSEP test-coefficients-unit
            Coefficient registry (COEF01-COEF06, COEF-SOURCE).  Updated test:
              COEF01  three spatial providers incl. parallel-diffusion-library

examples:
  %(prog)s                          run every suite
  %(prog)s --suite binding --suite cli
  %(prog)s --verbose --json report.json
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
    "binding": Suite(
        "binding", ["make", "-C", str(SRCSEP),
                    "test-parallel-diffusion-binding-unit"], AMPS_ROOT,
        expected=("PDB01", "PDB02", "PDB03", "PDB04", "PDB05")),
    "cli": Suite(
        "cli", ["make", "-C", str(SRCSEP), "test-cli-unit"], AMPS_ROOT,
        expected=("CLI01", "CLI02", "CLI03", "CLI04", "CLI05", "CLI06"),
        ignored=("HIDDEN01",)),
    "coefficients": Suite(
        "coefficients", ["make", "-C", str(SRCSEP), "test-coefficients-unit"],
        AMPS_ROOT,
        expected=("COEF01", "COEF02", "COEF03", "COEF04", "COEF05", "COEF06")),
}

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
    parser.add_argument("--suite", action="append", choices=sorted(SUITES),
                        help="run only this suite (repeatable; default: all)")
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

    selected = args.suite or list(SUITES)
    log_dir = args.log_dir
    if log_dir is None and args.json is not None:
        log_dir = args.json.with_name(args.json.stem + "-logs")
    runs = [run_suite(SUITES[name], args.verbose, log_dir) for name in selected]
    results = [result for run in runs for result in run.results]
    totals = count(results)
    failures = [r for r in results if r.status in ("FAIL", "ERROR", "MISSING")]

    print("\n=== srcSEP parallel-diffusion / mover test summary ===")
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
            "schema": "srcsep-pd-mover-runner-v1",
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
