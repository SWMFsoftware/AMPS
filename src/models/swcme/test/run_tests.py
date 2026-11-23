#!/usr/bin/env python3
"""Run the relocated SWCME Phase-R1 common-library tests.

This runner intentionally mirrors the public selection vocabulary used by
srcSEP and srcSEP3D.  It builds the bounded R1 executable from the canonical
``src/models/swcme`` location, executes each selected ID separately, and emits
JSON/JUnit reports from the same result records printed to the terminal.
The independent D10-1 material suite is a second bounded executable; its ten
checks join the routine set without altering or hiding any of the four R1 gates.

The older extended SWCME Makefile still documents a larger validation catalog,
but its source archive omitted several referenced translation units.  Missing
tests are not reported as passes here.  The four R1 integration gates remain executable; the new material checks
are independently registered below. The historical extended catalog remains
pending source restoration, never relabeled as successful routine evidence.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
import json
import math
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
MATERIAL_BINARY = ROOT / "output" / "test_material_upstream"
# Existing D10 acceptance, not a candidate-dependent performance allowance.
# This bounds execution after the build and preserves every selected ID as an
# explicit ERROR if a child consumes the remaining complete routine window.
ROUTINE_LIMIT_SECONDS = 180.0
# Explicit executable ownership, not a prefix/name heuristic: every new test
# below reaches the same canonical production Profile/Event implementation.
MATERIAL_IDS = frozenset(("MAT1D01", "MAT1D02", "MAT1D03", "MAT1D04", "MAT1D05",
                         "MAT1D06", "MAT1D07", "MAT1D08", "MAT1D09", "MAT1D10"))


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
    metrics: List[dict] = field(default_factory=list)


TESTS: Tuple[Definition, ...] = (
    Definition("R1CFG01", "COMMON", "Canonical constant identity"),
    Definition("R1D01", "1D", "One-dimensional prepared-state smoke"),
    Definition("R3D01", "3D", "Three-dimensional prepared-state smoke"),
    Definition("R13D01", "1D3D", "Common ambient-state identity"),
    Definition("MAT1D01", "MATERIAL1D", "Approved primitive/composition/EOS anchors"),
    Definition("MAT1D02", "MATERIAL1D", "Independent radial mass inventory and flux"),
    Definition("MAT1D03", "MATERIAL1D", "Full material deformation and signed Cauchy transport"),
    Definition("MAT1D04", "MATERIAL1D", "Signed magnetic flux and ideal induction"),
    Definition("MAT1D05", "MATERIAL1D", "All Cartesian gradients and second-order differentiation"),
    Definition("MAT1D06", "MATERIAL1D", "Adiabatic entropy and total-energy/work budget"),
    Definition("MAT1D07", "MATERIAL1D", "Full off-equator vector force and mechanical work"),
    Definition("MAT1D08", "MATERIAL1D", "Closed whole-step temporal/spatial coverage and large epoch"),
    Definition("MAT1D09", "MATERIAL1D", "Genuine canonical shock with independent RH fluxes"),
    Definition("MAT1D10", "MATERIAL1D", "Invalid/unsupported and preserved legacy diagnostics"),
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
        description="Run SWCME common and controlled material-profile tests.")
    result.add_argument("--list", action="store_true",
                        help="list test IDs and groups, then exit")
    result.add_argument("--test", dest="tests", action="append", default=[], metavar="ID",
                        help="run one test ID; repeatable")
    result.add_argument("--group", dest="groups", action="append", default=[], metavar="GROUP",
                        help="run one group; repeatable")
    result.add_argument("--routine", action="store_true",
                        help="run the complete bounded R1 and material set")
    result.add_argument("--all", action="store_true",
                        help="run every currently executable bounded SWCME test")
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


def build(args: argparse.Namespace, selected: Sequence[Definition]) -> None:
    # Preserve independent selection: an existing R1-only request still builds
    # only its original binary. Routine/all requests require BOTH executables;
    # an unavailable material build is an ERROR, never a filtered-out test.
    binaries = []
    targets = []
    if any(item.test_id not in MATERIAL_IDS for item in selected):
        binaries.append(BINARY)
        targets.append("r1-smoke")
    if any(item.test_id in MATERIAL_IDS for item in selected):
        binaries.append(MATERIAL_BINARY)
        targets.append("material-upstream-test-build")
    if args.no_build:
        for binary in binaries:
            if not binary.is_file():
                raise UsageError(f"--no-build requested but {binary} is absent")
        return
    if args.jobs < 1:
        raise UsageError("--jobs must be at least one")
    if args.rebuild:
        for binary in binaries:
            if binary.exists():
                binary.unlink()
    make = shutil.which("make")
    if make is None:
        raise UsageError("GNU Make is not available")
    command = [make, "-C", str(ROOT), f"-j{args.jobs}", *targets]
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, check=False)
    if args.verbose and completed.stdout:
        print(completed.stdout, end="")
    if completed.returncode != 0:
        raise UsageError("SWCME bounded build failed:\n" + completed.stdout[-6000:])
    for binary in binaries:
        if not binary.is_file():
            raise UsageError(f"build completed without producing {binary}")


def run_one(definition: Definition, verbose: bool,
            timeout_seconds: float = ROUTINE_LIMIT_SECONDS) -> Result:
    material = definition.test_id in MATERIAL_IDS
    command = [str(MATERIAL_BINARY if material else BINARY), "--test", definition.test_id]
    started = time.monotonic()
    if timeout_seconds <= 0:
        return Result(definition.test_id, definition.group, "ERROR",
                      "published180-second bounded routine execution window exhausted before this case",
                      0.0, command)
    try:
        completed = subprocess.run(
            command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        return Result(definition.test_id, definition.group, "ERROR",
                      "native case exceeded the remaining published180-second routine execution window",
                      time.monotonic()-started, command)
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
    metrics = []
    if material:
        # Read independent normalized residuals and their frozen criteria from
        # the actual native test, retaining raw stdout in message as before.
        # A success process with missing/contradictory metrics cannot qualify.
        malformed = False
        for line in output.splitlines():
            if not line.startswith("MATERIAL_METRIC "):
                continue
            fields = line.split()
            if len(fields) != 6:
                malformed = True
                continue
            _, name, actual, relation, limit, units = fields
            try:
                value, threshold = float(actual), float(limit)
            except ValueError:
                malformed = True
                continue
            finite = math.isfinite(value) and math.isfinite(threshold)
            accepted = finite and ((relation == "<=" and value <= threshold) or
                       (relation == "<" and value < threshold) or
                       (relation == ">" and value > threshold))
            metrics.append({"name": name, "value": value if finite else None,
                            "raw_value": actual, "criterion": relation,
                            "limit": threshold if math.isfinite(threshold) else None,
                            "units": units, "accepted": accepted})
            if status == "PASS" and not accepted:
                malformed = True
        if status == "PASS" and (malformed or not metrics or
                f"[{definition.test_id}] PASS" not in output):
            status = "ERROR"
            output += "\nmandatory material success record/metric evidence is missing or contradictory"
    return Result(definition.test_id, definition.group, status,
                  output or f"process exited {completed.returncode}",
                  elapsed, command, metrics)


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
        "runtime_budget": {"execution_seconds": sum(item.elapsed_seconds for item in results),
                           "limit_seconds": ROUTINE_LIMIT_SECONDS, "units": "s",
                           "scope": "bounded execution after build; selected IDs retained on timeout"},
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
        build(args, selected)
        results = []
        execution_start = time.monotonic()
        for definition in selected:
            remaining = ROUTINE_LIMIT_SECONDS-(time.monotonic()-execution_start)
            result = run_one(definition, args.verbose, remaining)
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
