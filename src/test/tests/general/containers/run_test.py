#!/usr/bin/env python3
"""Compile and run tests for dependency-light utilities in src/general."""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    case, command, error_result, load_context, main_guard, run_command, write_result)


EXPECTED_CASES = [
    "GENERAL-ARRAY-001",
    "GENERAL-BITARRAY-001",
    "GENERAL-STACK-001",
    "GENERAL-LINKED-001",
    "GENERAL-STRING-001",
    "GENERAL-INDEXER4D-001",
]


def prepare_header_mirror(repo: Path, work: Path) -> Path:
    """Mirror headers so generated .general.conf is test-local.

    ``global.h`` includes ``../../.general.conf`` by design.  A clean clone may
    not have that configured file yet, and a unit test must not create or edit
    it in the repository.  Keeping the original relative layout in the work
    directory makes this test independent of the currently selected AMPS app.
    """

    source = repo / "src" / "general"
    target = work / "header-mirror" / "src" / "general"
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (item.suffix in (".h", ".dfn") or item.name.endswith(".inc")):
            shutil.copyfile(str(item), str(target / item.name))
    (work / "header-mirror" / ".general.conf").write_text("", encoding="utf-8")
    return target


def main() -> int:
    context = load_context("Test general AMPS data structures")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_general_containers.cpp")
    executable = work / "test_general_containers"
    timeout = int(context["test"].get("timeout_seconds", 120))
    expected = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    try:
        include = prepare_header_mirror(repo, work)
        compile_result = run_command(
            command(context, "mpi_cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                "-DOMPI_SKIP_MPICXX=1",
                "-I", str(include), "-I", str(TEST_ROOT / "native"),
                str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
                "-o", str(executable),
            ], repo, timeout=timeout)
        if compile_result.returncode != 0:
            return write_result(
                context, "FAIL", "general-utility tests did not compile",
                cases=[case("GENERAL-COMPILE", "FAIL",
                            "compiler exit {}".format(compile_result.returncode))])

        # selected_cases is already ordered by the global runner.  Reapplying
        # reverse/shuffle inside the C++ registry would undo or perturb it.
        run_argv = [str(executable), "--order", "normal"]
        for identifier in expected:
            run_argv.extend(["--case", identifier])
        run_result = run_command(run_argv, work, timeout=timeout)
        artifacts.mkdir(parents=True, exist_ok=True)
        transcript = artifacts / "general-containers-output.txt"
        transcript.write_text(run_result.stdout or "", encoding="utf-8")
        # Require the full ordered registry.  This catches static-linking and
        # preprocessor mistakes that can silently omit a C++ translation unit
        # while leaving the surviving tests green.
        found = re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                           run_result.stdout or "", flags=re.MULTILINE)
        cases = [case(identifier, status, message or "ok")
                 for identifier, status, message in found]
        actual = [item["id"] for item in cases]
        if actual != expected:
            cases.append(case("GENERAL-REGISTRY", "FAIL",
                              "expected {}, observed {}".format(expected, actual)))
        if run_result.returncode != 0 and not any(item["status"] == "FAIL" for item in cases):
            cases.append(case("GENERAL-EXIT", "FAIL",
                              "test executable exit {}".format(run_result.returncode)))
        status = "PASS" if cases and all(item["status"] == "PASS" for item in cases) else "FAIL"
        return write_result(
            context, status, "{} general-utility cases reported".format(len(found)),
            cases=cases, evidence=[str(transcript)])
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
