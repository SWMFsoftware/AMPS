#!/usr/bin/env python3
"""Build and run dependency-light AMPS chemical-solver cases."""

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
    "CHEM-SOLVER-DECAY-001",
    "CHEM-SOLVER-REVERSIBLE-001",
    "CHEM-SOLVER-ROBERTSON-001",
    "CHEM-SOLVER-CONSERVATION-CONVERGENCE-001",
    "CHEM-SOLVER-NETWORK-VALIDATION-001",
    "CHEM-SOLVER-INVALID-001",
    "CHEM-SOLVER-SINGULAR-001",
]


def main() -> int:
    context = load_context("Test dependency-light AMPS chemical solver")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    executable = work / "test_chemical_solver"
    timeout = int(context["test"].get("timeout_seconds", 300))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get(
            "expected_cases", EXPECTED_CASES)))
    unknown = [identifier for identifier in requested
               if identifier not in EXPECTED_CASES]
    if not requested or unknown or len(set(requested)) != len(requested):
        raise RunnerError("invalid chemical-solver case selection: {}".format(
            requested))

    try:
        compiled = run_command(command(context, "cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
            "-I", str(repo / "src" / "models" / "chemical_solver"),
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"),
            str(repo / "test" / "srcChemicalSolverTests" / "main.cpp"),
            str(repo / "src" / "models" / "chemical_solver" / "chemeq2.cpp"),
            "-o", str(executable),
        ], repo, timeout=timeout)
        artifacts.mkdir(parents=True, exist_ok=True)
        build_log = artifacts / "chemical-solver-build.txt"
        build_log.write_text(compiled.stdout or "", encoding="utf-8")
        if compiled.returncode != 0:
            return write_result(
                context, "FAIL", "chemical-solver fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in requested],
                evidence=[str(build_log)])

        argv = [str(executable), "--order", "normal"]
        for identifier in requested:
            argv.extend(["--case", identifier])
        run = run_command(argv, work, timeout=timeout)
        transcript = artifacts / "chemical-solver-output.txt"
        transcript.write_text(run.stdout or "", encoding="utf-8")
        found = re.findall(
            r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
            run.stdout or "", flags=re.MULTILINE)
        observed_ids = [identifier for identifier, _, _ in found]
        exact_contract = (observed_ids == requested and
                          len(set(observed_ids)) == len(observed_ids))
        if exact_contract:
            results = [case(identifier, status, message or "ok")
                       for identifier, status, message in found]
        else:
            results = []
            for identifier in requested:
                exact = [record for record in found if record[0] == identifier]
                if len(exact) == 1:
                    results.append(case(identifier, exact[0][1],
                                        exact[0][2] or "ok"))
                else:
                    results.append(case(
                        identifier, "FAIL",
                        "expected one exact record; observed {}".format(
                            len(exact))))
            results[0] = case(
                results[0]["id"], "FAIL",
                "record contract violation: expected {} observed {}".format(
                    requested, observed_ids))
        status = "PASS" if run.returncode == 0 and all(
            item["status"] == "PASS" for item in results) else "FAIL"
        return write_result(
            context, status, "{} chemical-solver cases reported".format(
                len(results)), cases=results,
            evidence=[str(build_log), str(transcript)])
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
