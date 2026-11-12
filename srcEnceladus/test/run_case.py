#!/usr/bin/env python3
"""Run one Enceladus forward-model case and compare its named references."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import tempfile

import run_tests as suite


CASES = {
    "U00_exporter": None,
    "U01_source_export": None,
    "U02_core_compile": None,
    "U03_source_parser": "source",
    "U03b_canonical_source_geometry": "source_geometry",
    "U03c_catalog_source_geometry": "catalog",
    "U04_runtime_parser": "input",
    "U04b_transport_only_parser": "input",
    "U05_chemistry": "chemistry",
    "U06_source_law": "source_law",
    "U07_numerics": "numerics",
    "U08_validation": "validation",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = (args.output or
              suite.ROOT / "test_output" / "srcEnceladus" /
              "individual" / args.case).resolve()
    output.mkdir(parents=True, exist_ok=True)

    if args.case == "U00_exporter":
        result = suite.run([sys.executable, "-m", "unittest",
                            "standalone.tests.test_export_forward_sources",
                            "standalone.tests.test_fit_source_export",
                            "standalone.tests.test_export_forward_observations"])
        detail = result.stdout + result.stderr
        passed = result.returncode == 0
        values: dict[str, object] = {}
    else:
        expected_rejection = False
        with tempfile.TemporaryDirectory(prefix="enceladus_case_") as temporary:
            temp = Path(temporary)
            source = temp / "sources.in"
            exported = suite.run([
                sys.executable,
                str(suite.ROOT / "standalone" / "export_forward_sources.py"),
                "--manifest", str(suite.TEST / "data" / "shape_manifest.json"),
                "--output", str(source),
            ])
            if args.case == "U01_source_export":
                result = exported
            else:
                binary = temp / "forward_model_tests"
                compiled = suite.run([
                    "g++", "-std=c++17", "-O2", "-I..", "-I../src",
                    "ForwardModel.cpp", "test/forward_model_tests.cpp",
                    "-o", str(binary),
                ])
                if args.case == "U02_core_compile" or compiled.returncode != 0:
                    result = compiled
                elif exported.returncode != 0:
                    result = exported
                elif args.case == "U03c_catalog_source_geometry":
                    result = suite.run([sys.executable,
                        str(suite.TEST/"run_source_export_roundtrip.py"),
                        "--executable", str(binary)])
                else:
                    runtime_input = temp / "enceladus.in"
                    if args.case == "U04b_transport_only_parser":
                        suite.write_transport_only_input(runtime_input, source.name)
                    else:
                        suite.write_input(runtime_input, source.name)
                    mode = CASES[args.case]
                    if mode == "source_geometry":
                        wrong_geometry = temp / "wrong_geometry.in"
                        suite.write_wrong_geometry_source(
                            wrong_geometry,
                            source.read_text(encoding="utf-8"))
                        command = [str(binary), "source", str(wrong_geometry)]
                        expected_rejection = True
                    else:
                        command = [str(binary), str(mode)]
                    if mode in ("source", "numerics", "validation"):
                        command.append(str(source))
                    elif mode == "input":
                        command.extend((str(runtime_input), str(source)))
                    result = suite.run(command)
            detail = result.stdout + result.stderr
            passed = ((result.returncode != 0 and "must use geometry=" in detail)
                      if expected_rejection else result.returncode == 0)
            values = suite.parse_receipt(result.stdout) if passed else {}

    reference = json.loads(suite.REFERENCE.read_text(encoding="utf-8"))
    comparisons: list[tuple[str, float, float]] = []
    for key, actual in values.items():
        expected = reference.get(key)
        if not isinstance(actual, float) or not isinstance(expected, (int, float)):
            continue
        matched = math.isclose(actual, float(expected),
                               rel_tol=1e-12, abs_tol=1e-300)
        passed = passed and matched
        comparisons.append((key, actual, float(expected)))
        if not matched:
            detail += f"reference mismatch {key}: {actual} != {expected}\n"
    suite.svg_comparison(output / "reference_comparison.svg", comparisons)
    receipt = {
        "case": args.case,
        "status": "PASS" if passed else "FAIL",
        "detail": detail,
        "calculated_values": values,
        "reference": str(suite.REFERENCE.resolve()),
    }
    (output / "test_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{'PASS' if passed else 'FAIL'} {args.case}")
    if detail:
        print(detail.rstrip())
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
