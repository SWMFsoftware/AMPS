#!/usr/bin/env python3
"""Run the relocated SWCME Phase-R1 common-library tests.

This runner intentionally mirrors the public selection vocabulary used by
srcSEP and srcSEP3D.  It builds the bounded R1 executable from the canonical
``src/models/swcme`` location, executes each selected ID separately, and emits
JSON/JUnit reports from the same result records printed to the terminal.

The older extended SWCME Makefile still documents a larger validation catalog,
but its source archive omitted several referenced translation units.  Missing
tests are not reported as passes here.  R1 exposes only the four executable
common-library gates registered below; the README identifies the extended
catalog as pending source restoration.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Dict, Iterable, List, Sequence, Tuple
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parent
MODEL_ROOT = ROOT.parent
BINARY = ROOT / "output" / "test_swcme_r1"


@dataclass(frozen=True)
class Definition:
    test_id: str
    group: str
    name: str


@dataclass
class Result:
    test_id: str
    group: str
    status: str
    message: str
    elapsed_seconds: float
    command: List[str]


TESTS: Tuple[Definition, ...] = (
    Definition("R1CFG01", "COMMON", "Canonical constant identity"),
    Definition("R1D01", "1D", "One-dimensional prepared-state smoke"),
    Definition("R3D01", "3D", "Three-dimensional prepared-state smoke"),
    Definition("R13D01", "1D3D", "Common ambient-state identity"),
)
BY_ID: Dict[str, Definition] = {item.test_id: item for item in TESTS}
GROUPS: Dict[str, Tuple[Definition, ...]] = {
    group: tuple(item for item in TESTS if item.group == group)
    for group in sorted({item.group for item in TESTS})
}


class UsageError(RuntimeError):
    """Invalid selection or build setup; reported with process exit code 2."""


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Run SWCME common 1-D/3-D integration tests.")
    result.add_argument("--list", action="store_true",
                        help="list test IDs and groups, then exit")
    result.add_argument("--test", dest="tests", action="append", default=[], metavar="ID",
                        help="run one test ID; repeatable")
    result.add_argument("--group", dest="groups", action="append", default=[], metavar="GROUP",
                        help="run one group; repeatable")
    result.add_argument("--routine", action="store_true",
                        help="run the complete bounded R1 set")
    result.add_argument("--all", action="store_true",
                        help="run every currently executable SWCME R1 test")
    result.add_argument("--output-dir", type=Path,
                        default=ROOT / "output" / "r1-results",
                        help="directory for logs and JSON/JUnit summaries")
    result.add_argument("--jobs", type=int,
                        default=max(1, min(16, os.cpu_count() or 1)),
                        help="parallel jobs used by GNU Make")
    result.add_argument("--rebuild", action="store_true",
                        help="remove the bounded R1 executable before building")
    result.add_argument("--no-build", action="store_true",
                        help="use an existing bounded R1 executable")
    result.add_argument("--verbose", action="store_true",
                        help="show complete build and test output")
    return result


def select(args: argparse.Namespace) -> List[Definition]:
    selectors = bool(args.tests or args.groups or args.routine or args.all)
    if args.list:
        if selectors:
            raise UsageError("--list cannot be combined with execution selectors")
        return []
    if not selectors:
        raise UsageError("select --test, --group, --routine, or --all")
    if args.routine and args.all:
        raise UsageError("--routine and --all are mutually exclusive")

    selected: List[Definition] = []
    seen = set()

    def add(items: Iterable[Definition]) -> None:
        for item in items:
            if item.test_id not in seen:
                selected.append(item)
                seen.add(item.test_id)

    if args.routine or args.all:
        add(TESTS)
    for raw in args.tests:
        key = raw.upper()
        if key not in BY_ID:
            raise UsageError(f"unknown test ID {raw!r}; use --list")
        add((BY_ID[key],))
    for raw in args.groups:
        key = raw.upper()
        if key not in GROUPS:
            raise UsageError(f"unknown group {raw!r}; use --list")
        add(GROUPS[key])
    return selected


def build(args: argparse.Namespace) -> None:
    if args.no_build:
        if not BINARY.is_file():
            raise UsageError(f"--no-build requested but {BINARY} is absent")
        return
    if args.jobs < 1:
        raise UsageError("--jobs must be at least one")
    if args.rebuild and BINARY.exists():
        BINARY.unlink()
    make = shutil.which("make")
    if make is None:
        raise UsageError("GNU Make is not available")
    command = [make, "-C", str(ROOT), f"-j{args.jobs}", "r1-smoke"]
    completed = subprocess.run(
        command, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False)
    if args.verbose and completed.stdout:
        print(completed.stdout, end="")
    if completed.returncode != 0:
        raise UsageError("SWCME R1 build failed:\n" + completed.stdout[-6000:])
    if not BINARY.is_file():
        raise UsageError(f"build completed without producing {BINARY}")


def run_one(definition: Definition, verbose: bool) -> Result:
    command = [str(BINARY), "--test", definition.test_id]
    started = time.monotonic()
    completed = subprocess.run(
        command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False)
    elapsed = time.monotonic() - started
    output = completed.stdout.strip()
    if verbose and output:
        print(output)
    if completed.returncode == 0:
        status = "PASS"
    elif completed.returncode == 1:
        status = "FAIL"
    else:
        status = "ERROR"
    return Result(definition.test_id, definition.group, status,
                  output or f"process exited {completed.returncode}",
                  elapsed, command)


def write_reports(output_dir: Path, results: Sequence[Result]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    totals = {
        status: sum(item.status == status for item in results)
        for status in ("PASS", "FAIL", "SKIP", "ERROR")
    }
    payload = {
        "schema": "swcme-r1-test-summary-v1",
        "totals": totals,
        "results": [asdict(item) for item in results],
    }
    (output_dir / "swcme-r1-tests.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    failures = totals["FAIL"]
    cases = []
    for item in results:
        body = ""
        if item.status in ("FAIL", "ERROR"):
            tag = "failure" if item.status == "FAIL" else "error"
            body = f'<{tag} message="{escape(item.message)}"/>'
        cases.append(
            f'<testcase classname="SWCME.{escape(item.group)}" '
            f'name="{escape(item.test_id)}" time="{item.elapsed_seconds:.6f}">'
            f'{body}</testcase>')
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<testsuite name="SWCME-R1" tests="{len(results)}" '
        f'failures="{failures}" errors="{totals["ERROR"]}" skipped="0">\n'
        + "\n".join(cases) + "\n</testsuite>\n")
    (output_dir / "swcme-r1-tests.xml").write_text(xml, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        selected = select(args)
        if args.list:
            for item in TESTS:
                print(f"{item.test_id:<10} {item.group:<8} {item.name}")
            return 0
        build(args)
        results = []
        for definition in selected:
            result = run_one(definition, args.verbose)
            results.append(result)
            print(f"[{result.test_id}] {result.status}")
        write_reports(args.output_dir.resolve(), results)
        failed = [item.test_id for item in results if item.status == "FAIL"]
        errors = [item.test_id for item in results if item.status == "ERROR"]
        print(f"Summary: PASS={sum(x.status == 'PASS' for x in results)} "
              f"FAIL={len(failed)} SKIP=0 ERROR={len(errors)}")
        print("FAIL: " + (", ".join(failed) if failed else "none"))
        print("ERROR: " + (", ".join(errors) if errors else "none"))
        return 2 if errors else (1 if failed else 0)
    except UsageError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
