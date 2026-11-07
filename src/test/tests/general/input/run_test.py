#!/usr/bin/env python3
"""Build and run the dependency-free generic AMPS input-reader tests."""

from __future__ import annotations

import re
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    case, command, error_result, load_context, main_guard, run_command, write_result)


def main() -> int:
    context = load_context("Test the generic AMPS runtime-input parser")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source_dir = repo / "src" / "general" / "input"
    executable = work / "test_amps_input"
    timeout = int(context["test"].get("timeout_seconds", 120))
    expected = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", [])))
    try:
        compile_result = run_command(
            command(context, "cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
                "-I", str(source_dir), str(source_dir / "test_amps_input.cpp"),
                "-o", str(executable),
            ], repo, timeout=timeout)
        if compile_result.returncode != 0:
            return write_result(
                context, "FAIL", "input-reader test did not compile",
                cases=[case("INP-COMPILE", "FAIL", "compiler exit {}".format(
                    compile_result.returncode))])

        # The global runner has already materialized the requested deterministic
        # order in selected_cases.  Pass that sequence verbatim and keep the
        # child registry from applying the order transform a second time.
        run_argv = [str(executable), "--order", "normal"]
        for identifier in expected:
            run_argv.extend(["--case", identifier])
        run_result = run_command(run_argv, work, timeout=timeout)
        artifacts.mkdir(parents=True, exist_ok=True)
        transcript = artifacts / "input-reader-output.txt"
        transcript.write_text(run_result.stdout or "", encoding="utf-8")
        # The exact ordered case list is a coverage receipt, not merely output
        # decoration.  If a case is deleted, renamed, duplicated, or stops
        # executing, the adapter fails even when the executable returns zero.
        found = re.findall(r"^\[(INP-[0-9]+)\]\s+(PASS|FAIL)\s*(.*)$",
                           run_result.stdout or "", flags=re.MULTILINE)
        cases = [case(identifier, status, message or "ok")
                 for identifier, status, message in found]
        actual = [item["id"] for item in cases]
        if actual != expected:
            cases.append(case("INP-REGISTRY", "FAIL",
                              "expected {}, observed {}".format(expected, actual)))
        if run_result.returncode != 0 and not any(item["status"] == "FAIL" for item in cases):
            cases.append(case("INP-EXIT", "FAIL",
                              "test executable exit {}".format(run_result.returncode)))
        status = "PASS" if cases and all(item["status"] == "PASS" for item in cases) else "FAIL"
        return write_result(
            context, status, "{} parser cases reported".format(len(found)), cases=cases,
            evidence=[str(transcript)])
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
