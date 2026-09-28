#!/usr/bin/env python3
"""Run one SCCM test or a cumulative Stage 0--2 release gate."""

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
]

def command_for(test: Test) -> list[str]:
    if test.kind == "cpp":
        return [str(ROOT / "build" / "sep_coronal_cme_tests"), "--test", test.identifier]
    if test.kind == "architecture":
        return [sys.executable, str(ROOT / "tools" / "check_architecture.py")]
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
    select.add_argument("--stage", type=int, choices=(0, 1, 2))
    select.add_argument("--test")
    select.add_argument("--list", action="store_true")
    parser.add_argument("--output-dir", default="build/test-results")
    args = parser.parse_args(argv)
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
