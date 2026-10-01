#!/usr/bin/env python3
"""Run one SCCM test or a cumulative Stage 0--12 release gate."""

# Passion and several NASA HEC environments still provide Python 3.8.  With
# postponed annotations, expressions such as ``list[str]`` are stored as text
# instead of being evaluated during module import.  This preserves the useful
# type declarations while remaining compatible with Python 3.7/3.8, where the
# built-in collection classes are not yet runtime-subscriptable.
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]

@dataclass(frozen=True)
class Test:
    identifier: str
    stage: int
    kind: str = "cpp"

TESTS = [
    Test("DOCSCCM01", 0, "docs"), Test("ARCHSCCM01", 0, "architecture"),
    *[Test(f"CFG3D{number}", 0) for number in range(12, 21)],
    Test("RST3D04", 0), Test("RST3D05", 0),
    *[Test(f"PFSS3D{number:02d}", 1) for number in range(1, 10)],
    *[Test(f"WND3D{number:02d}", 2) for number in range(1, 21)],
    *[Test(f"CLS3D{number:02d}", 2) for number in range(1, 11)],
    Test("STR3D01", 2),
    *[Test(f"SCS3D{number:02d}", 3) for number in range(1, 10)],
    *[Test(f"HCS3D{number:02d}", 3) for number in range(1, 4)],
    *[Test(f"CPL3D{number:02d}", 3)
      for number in (*range(1, 10), 11, 12)],
    *[Test(f"PLS3D{number:02d}", 3) for number in range(1, 5)],
    Test("OFX3D01", 3), Test("LOS3D01", 3),
    *[Test(f"TUR3D{number:02d}", 4) for number in range(1, 8)],
    *[Test(f"MFP3D{number:02d}", 4) for number in range(1, 8)],
    *[Test(f"ELL3D{number:02d}", 5) for number in range(1, 11)],
    *[Test(f"RH3D{number:02d}", 6) for number in range(1, 12)],
    *[Test(f"SHK3D{number:02d}", 7) for number in range(5, 15)],
    Test("SNAP3D09", 7),
    *[Test(f"BND3D{number:02d}", 8) for number in range(1, 3)],
    Test("MESH3D01", 8),
    *[Test(f"COR3D{number:02d}", 8) for number in range(1, 4)],
    *[Test(f"TIM3D{number:02d}", 8) for number in range(1, 3)],
    *[Test(f"POP3D{number:02d}", 8) for number in range(1, 4)],
    Test("MPI3D01", 8),
    *[Test(f"INIT3D{number:02d}", 8) for number in range(1, 5)],
    Test("NAT3D13", 8), Test("RUN3D02", 8),
    *[Test(f"OBS3D{number:02d}", 8) for number in range(1, 8)],
    *[Test(f"SRC3D{number:02d}", 9) for number in range(1, 20)],
    Test("LOS3D02", 9),
    Test("HCS3D07", 10),
    *[Test(f"FLX3D{number:02d}", 10) for number in range(1, 10)],
    *[Test(f"FLX1D{number:02d}", 10) for number in range(1, 4)],
    *[Test(f"OBS1D{number:02d}", 10) for number in range(1, 3)],
    Test("NAT1D01", 10), Test("TIM1D01", 10),
    *[Test(f"POP1D{number:02d}", 10) for number in range(1, 3)],
    Test("RUN1D01", 10), Test("RST1D01", 10),
    *[Test(f"XM3D{number:02d}", 10) for number in range(1, 3)],
    *[Test(f"HCS3D{number:02d}", 11) for number in range(4, 7)],
    *[Test(f"SHEATH3D{number:02d}", 11) for number in range(1, 4)],
    *[Test(f"PROV3D{number:02d}", 12, "preprocessing") for number in range(1, 4)],
    Test("CAL3D01", 12, "preprocessing"),
]

# These are release-contract values, not counts inferred from TESTS.  Keeping
# an independent expectation is intentional: if a partial source update leaves
# a Stage 0--2 runner beside Stage 3--6 model code, a nominal Stage 6 command
# must fail instead of silently reporting the old 53-test suite as complete.
# Each value is cumulative because every stage gate re-runs all earlier stages.
EXPECTED_CUMULATIVE_COUNTS = {
    0: 13,
    1: 22,
    2: 53,
    3: 82,
    4: 96,
    5: 106,
    6: 117,
    7: 128,
    8: 153,
    9: 173,
    10: 196,
    11: 202,
    12: 206,
}

# A terminal test ID makes the diagnostic more useful than a count alone.  It
# catches a registry with the right cardinality but the wrong stage contents.
EXPECTED_STAGE_TERMINALS = {
    0: "RST3D05",
    1: "PFSS3D09",
    2: "STR3D01",
    3: "LOS3D01",
    4: "MFP3D07",
    5: "ELL3D10",
    6: "RH3D11",
    7: "SNAP3D09",
    8: "OBS3D07",
    9: "LOS3D02",
    10: "XM3D02",
    11: "SHEATH3D03",
    12: "CAL3D01",
}


def validate_registry() -> tuple[bool, str]:
    """Validate the hand-authored release registry before selecting tests.

    The runner is itself part of the release evidence.  Duplicate identifiers,
    a missing late-stage family, or an accidental stage reassignment would make
    its summary misleading, so those conditions are hard errors before any
    potentially expensive test is launched.
    """

    identifiers = [test.identifier for test in TESTS]
    if len(identifiers) != len(set(identifiers)):
        duplicates = sorted({identifier for identifier in identifiers
                             if identifiers.count(identifier) > 1})
        return False, "duplicate test ID(s): " + ", ".join(duplicates)

    for stage, expected_count in EXPECTED_CUMULATIVE_COUNTS.items():
        cumulative = [test for test in TESTS if test.stage <= stage]
        if len(cumulative) != expected_count:
            return False, (
                f"Stage {stage} registry contains {len(cumulative)} tests; "
                f"the release contract requires {expected_count}"
            )
        terminal = EXPECTED_STAGE_TERMINALS[stage]
        if not any(test.identifier == terminal and test.stage == stage
                   for test in TESTS):
            return False, (
                f"Stage {stage} registry is missing terminal test {terminal}"
            )

    return True, ""

def command_for(test: Test) -> list[str]:
    if test.kind == "cpp":
        return [str(ROOT / "build" / "sep_coronal_cme_tests"), "--test", test.identifier]
    if test.kind == "preprocessing":
        return [sys.executable, str(ROOT / "test" / "test_stage12.py"), "--test", test.identifier]
    if test.kind == "architecture":
        # Run adversarial source/artifact fixtures before the same production
        # archive/public-ABI audit. They belong to ARCHSCCM01, not DOCSCCM01.
        return [sys.executable, str(ROOT / "test" / "test_architecture.py")]
    # File-path invocation is not portable across Python versions: Python 3.8
    # converts ``test/test_model_documentation.py`` to the dotted import
    # ``test.test_model_documentation``, which fails because this test tree is
    # deliberately not an import package (and ``test`` can name the standard-
    # library test package). Discovery loads the same file directly from the
    # explicit start directory on every supported Python version.
    return [sys.executable, "-B", "-m", "unittest", "discover", "-v",
            "-s", "test", "-p", "test_model_documentation.py"]

def documentation_generated_checks() -> tuple[bool, str]:
    output = ""
    for command in ([sys.executable, "tools/generate_model.py", "--check"],
                    [sys.executable, "tools/generate_schema_registry.py", "--check"],
                    [sys.executable, "tools/generate_stage0_fixture.py", "--check"]):
        result = subprocess.run(command, cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output += result.stdout
        if result.returncode:
            return False, output
    return True, output

def write_reports(results: list[dict[str, object]], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    passed = sum(bool(item["passed"]) for item in results)
    document = {"suite": "sep_coronal_cme",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "total": len(results), "passed": passed, "failed": len(results) - passed,
        "tests": results}
    (output / "results.json").write_text(json.dumps(document, indent=2) + "\n",
                                         encoding="utf-8")
    suite = ET.Element("testsuite", name="sep_coronal_cme",
        tests=str(len(results)), failures=str(len(results) - passed),
        time=f"{sum(float(x['seconds']) for x in results):.6f}")
    for item in results:
        case = ET.SubElement(suite, "testcase", name=str(item["id"]),
            classname=f"stage{item['stage']}", time=f"{float(item['seconds']):.6f}")
        if not item["passed"]:
            failure = ET.SubElement(case, "failure", message="test failed")
            failure.text = str(item["output"])
        system = ET.SubElement(case, "system-out")
        system.text = str(item["output"])
    ET.ElementTree(suite).write(output / "junit.xml", encoding="utf-8",
                                xml_declaration=True)

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    select = parser.add_mutually_exclusive_group(required=True)
    select.add_argument("--all", action="store_true")
    select.add_argument("--stage", type=int, choices=tuple(range(0, 13)))
    select.add_argument("--test")
    select.add_argument("--list", action="store_true")
    parser.add_argument("--output-dir", default="build/test-results")
    parser.add_argument(
        "--expect-count", type=int,
        help=("fail unless the selected gate contains this many tests; Make "
              "targets use this to detect a stale or partially copied runner"),
    )
    args = parser.parse_args(argv)

    registry_ok, registry_error = validate_registry()
    if not registry_ok:
        print(f"TEST REGISTRY ERROR: {registry_error}", file=sys.stderr)
        return 2

    if args.list:
        for test in TESTS:
            print(f"{test.identifier}\tstage={test.stage}\t{test.kind}")
        return 0
    known = {test.identifier: test for test in TESTS}
    if args.test:
        if args.test not in known:
            parser.error(f"unknown test ID {args.test}")
        selected = [known[args.test]]
    elif args.all:
        selected = TESTS
    else:
        selected = [test for test in TESTS if test.stage <= args.stage]

    if args.expect_count is not None and len(selected) != args.expect_count:
        print(
            "TEST SELECTION ERROR: selected "
            f"{len(selected)} tests but --expect-count requires "
            f"{args.expect_count}",
            file=sys.stderr,
        )
        return 2

    # Print the selection before running it.  This gives batch logs an
    # unambiguous indication that a Stage 6 invocation selected all 117 tests,
    # even if a later compile, scheduler, or test failure interrupts the run.
    if args.test:
        selection_name = f"test={args.test}"
    elif args.all:
        selection_name = "all"
    else:
        selection_name = f"stage={args.stage}"
    print(f"SELECTION: {selection_name}; tests={len(selected)}; "
          f"registry-total={len(TESTS)}")

    binary = ROOT / "build" / "sep_coronal_cme_tests"
    if any(test.kind == "cpp" for test in selected) and not binary.exists():
        build = subprocess.run(["make", str(binary.relative_to(ROOT))], cwd=ROOT)
        if build.returncode:
            return build.returncode
    results: list[dict[str, object]] = []
    for test in selected:
        start = time.monotonic()
        process = subprocess.run(command_for(test), cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        passed, output = process.returncode == 0, process.stdout
        if passed and test.kind == "docs":
            passed, generated = documentation_generated_checks()
            output += generated
        elapsed = time.monotonic() - start
        print(f"[{test.identifier}] {'PASS' if passed else 'FAIL'} ({elapsed:.3f}s)")
        if not passed:
            print(output.rstrip())
        results.append({"id": test.identifier, "stage": test.stage,
            "passed": passed, "seconds": elapsed, "output": output})
    write_reports(results, ROOT / args.output_dir)
    failures = sum(not bool(item["passed"]) for item in results)
    print(f"SUMMARY: {len(results) - failures}/{len(results)} passed; evidence={ROOT / args.output_dir}")
    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
