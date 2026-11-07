#!/usr/bin/env python3
"""Compile and run dependency-light SI/PIC normalization tests."""

from __future__ import annotations

import re
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


EXPECTED_CASES = [
    "UNITS-FACTOR-DERIVATION-001",
    "UNITS-SCALAR-ROUNDTRIP-001",
    "UNITS-VECTOR-ROUNDTRIP-001",
    "UNITS-DIMENSIONAL-IDENTITIES-001",
    "UNITS-INVALID-SCALES-001",
]


def main() -> int:
    context = load_context("Test dependency-light PIC unit normalization")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_units.cpp")
    executable = work / "test_pic_units"
    timeout = int(context["test"].get("timeout_seconds", 120))
    expected = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    unknown = [identifier for identifier in expected
               if identifier not in EXPECTED_CASES]
    if not expected or unknown or len(set(expected)) != len(expected):
        raise RunnerError("invalid units case selection: {}".format(expected))
    try:
        compile_result = run_command(
            command(context, "cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                "-I", str(repo / "src" / "pic" / "units"),
                "-I", str(TEST_ROOT / "native"),
                str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
                "-o", str(executable),
            ], repo, timeout=timeout)
        if compile_result.returncode != 0:
            message = "compiler exit {}".format(compile_result.returncode)
            return write_result(
                context, "FAIL", "PIC units tests did not compile",
                cases=[case(identifier, "FAIL", message)
                       for identifier in expected])

        # selected_cases is already in the root runner's final order.
        run_argv = [str(executable), "--order", "normal"]
        for identifier in expected:
            run_argv.extend(["--case", identifier])
        run_result = run_command(run_argv, work, timeout=timeout)
        artifacts.mkdir(parents=True, exist_ok=True)
        transcript = artifacts / "pic-units-output.txt"
        transcript.write_text(run_result.stdout or "", encoding="utf-8")
        found = re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                           run_result.stdout or "", flags=re.MULTILINE)
        cases = [case(identifier, status, message or "ok")
                 for identifier, status, message in found]
        actual = [item["id"] for item in cases]
        if actual != expected:
            # Fail closed: dictionary reconstruction used to discard duplicate
            # or foreign FAIL records and could therefore publish a false PASS.
            # The native harness is required to emit exactly one record for
            # each requested case, in the requested order, and nothing else.
            detail = "expected exact ordered records {}; observed {}".format(
                expected, actual)
            cases = [case(identifier, "FAIL", detail)
                     for identifier in expected]
        if run_result.returncode != 0 and all(
                item["status"] == "PASS" for item in cases):
            cases[-1] = case(cases[-1]["id"], "FAIL",
                             "test executable exit {}".format(
                                 run_result.returncode))
        status = "PASS" if all(item["status"] == "PASS" for item in cases) else "FAIL"
        return write_result(
            context, status, "{} PIC units cases reported".format(len(found)),
            cases=cases, evidence=[str(transcript)])
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
