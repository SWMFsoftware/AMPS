#!/usr/bin/env python3
"""Build and run dependency-light src/species radiation-table tests."""

from __future__ import annotations

import re
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard, run_command,
    write_result)


EXPECTED_CASES = [
    "SPECIES-RADIATION-REFERENCE-001",
    "SPECIES-RADIATION-SCALING-001",
    "SPECIES-RADIATION-ACCELERATION-001",
    "SPECIES-RADIATION-INVALID-001",
]


def main() -> int:
    context = load_context("Test helium/sodium radiation tables")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_species_radiation.cpp")
    executable = work / "test_species_radiation"
    timeout = int(context["test"].get("timeout_seconds", 120))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    unknown = [identifier for identifier in requested
               if identifier not in EXPECTED_CASES]
    if not requested or unknown or len(set(requested)) != len(requested):
        raise RunnerError("invalid species-radiation case selection: {}".format(
            requested))
    try:
        compiled = run_command(command(context, "cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
            "-I", str(repo / "src" / "general"),
            "-I", str(repo / "src" / "species"),
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
            str(repo / "src" / "species" / "He.cpp"),
            str(repo / "src" / "species" / "Na.cpp"),
            "-o", str(executable),
        ], repo, timeout=timeout)
        if compiled.returncode != 0:
            return write_result(
                context, "FAIL", "species radiation fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in requested])

        argv = [str(executable), "--order", "normal"]
        for identifier in requested:
            argv.extend(["--case", identifier])
        run = run_command(argv, work, timeout=timeout)
        artifacts.mkdir(parents=True, exist_ok=True)
        transcript = artifacts / "species-radiation-output.txt"
        transcript.write_text(run.stdout or "", encoding="utf-8")
        found = re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
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
            # Foreign, duplicate, or out-of-order output is a suite protocol
            # violation even when every requested case also printed PASS.
            results[0] = case(
                results[0]["id"], "FAIL",
                "record contract violation: expected {} observed {}".format(
                    requested, observed_ids))
        status = "PASS" if run.returncode == 0 and all(
            item["status"] == "PASS" for item in results) else "FAIL"
        return write_result(
            context, status, "{} species radiation cases reported".format(
                len(results)), cases=results, evidence=[str(transcript)])
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
