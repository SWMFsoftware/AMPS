#!/usr/bin/env python3
"""Build and run the dependency-light dust kernels.

The legacy ``test/srcModelDustTest`` executable exercises equilibrium and
time-dependent charging only after a complete PIC/Exosphere configuration.
The production spherical injection path is owned by the separately configured
``CORE-DUST-INJECTION-NATIVE-001`` suite.  Keeping that case out of this
dependency-light suite gives every global case ID one authoritative runner.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


CPP_CASES = [
    "CORE-DUST-SIZE-BINS-001",
    "CORE-DUST-SIZE-DISTRIBUTION-001",
    "CORE-DUST-CHARGING-POTENTIAL-001",
    "CORE-DUST-FORCE-LIMITS-001",
    "CORE-DUST-INJECTION-LAW-001",
    "CORE-DUST-MIE-INTERPOLATION-001",
    "CORE-DUST-MOMENTS-001",
    "CORE-DUST-SPLIT-MERGE-001",
]

EXPECTED_CASES = CPP_CASES


def _parse_records(output: str, requested: list[str]) -> list[dict[str, str]]:
    found = re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                       output or "", flags=re.MULTILINE)
    observed_ids = [identifier for identifier, _, _ in found]
    if observed_ids == requested and len(set(observed_ids)) == len(observed_ids):
        return [case(identifier, status, message or "ok")
                for identifier, status, message in found]

    results = []
    for identifier in requested:
        exact = [record for record in found if record[0] == identifier]
        if len(exact) == 1:
            results.append(case(identifier, exact[0][1], exact[0][2] or "ok"))
        else:
            results.append(case(
                identifier, "FAIL",
                "expected one exact record; observed {}".format(len(exact))))
    if results:
        results[0] = case(
            results[0]["id"], "FAIL",
            "record contract violation: expected {} observed {}".format(
                requested, observed_ids))
    return results


def main() -> int:
    context = load_context("Test dependency-light AMPS dust physics")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_dust_core.cpp")
    executable = work / "test_dust_core"
    timeout = int(context["test"].get("timeout_seconds", 120))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get(
            "expected_cases", EXPECTED_CASES)))
    unknown = [identifier for identifier in requested
               if identifier not in EXPECTED_CASES]
    if not requested or unknown or len(set(requested)) != len(requested):
        raise RunnerError("invalid dust case selection: {}".format(requested))

    runnable = [identifier for identifier in requested if identifier in CPP_CASES]
    by_id: dict[str, dict[str, str]] = {}
    evidence: list[str] = []
    native_returncode = 0

    try:
        artifacts.mkdir(parents=True, exist_ok=True)
        if runnable:
            compiled = run_command(command(context, "cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                "-I", str(repo / "src" / "models" / "dust"),
                "-I", str(TEST_ROOT / "native"),
                str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
                "-o", str(executable),
            ], repo, timeout=timeout)
            build_log = artifacts / "dust-core-build.txt"
            build_log.write_text(compiled.stdout or "", encoding="utf-8")
            evidence.append(str(build_log))
            if compiled.returncode != 0:
                for identifier in runnable:
                    by_id[identifier] = case(
                        identifier, "FAIL",
                        "dust fixture compiler exit {}".format(
                            compiled.returncode))
                native_returncode = compiled.returncode
            else:
                argv = [str(executable), "--order", "normal"]
                for identifier in runnable:
                    argv.extend(["--case", identifier])
                run = run_command(argv, work, timeout=timeout)
                transcript = artifacts / "dust-core-output.txt"
                transcript.write_text(run.stdout or "", encoding="utf-8")
                evidence.append(str(transcript))
                for result in _parse_records(run.stdout or "", runnable):
                    by_id[result["id"]] = result
                native_returncode = run.returncode

        results = [by_id[identifier] for identifier in requested]
        if native_returncode != 0 or any(
                result["status"] == "FAIL" for result in results):
            status = "FAIL"
        else:
            status = "PASS"
        return write_result(
            context, status,
            "{} dependency-light dust cases reported".format(len(results)),
            cases=results, evidence=evidence,
            metrics={"implemented_cases": len(runnable),
                     "unavailable_cases": 0})
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
