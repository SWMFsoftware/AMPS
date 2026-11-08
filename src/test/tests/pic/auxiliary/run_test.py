#!/usr/bin/env python3
"""Compile and run dependency-light T8 PIC auxiliary kernels."""

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
    "PIC-AUX-RAY-BOX-001",
    "PIC-AUX-COLUMN-CONSTANT-001",
    "PIC-AUX-COLUMN-AFFINE-001",
    "PIC-AUX-STOPPING-TABLE-001",
    "PIC-AUX-BACKGROUND-DISPATCH-001",
    "PIC-AUX-BACKGROUND-INTERPOLATION-001",
]


def main() -> int:
    context = load_context("Test PIC ray/column/table/background kernels")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 120))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    if (not requested or len(set(requested)) != len(requested) or
            any(identifier not in EXPECTED_CASES for identifier in requested)):
        raise RunnerError("invalid PIC auxiliary case selection: {}".format(
            requested))
    executable = work / "test_pic_auxiliary"
    source = Path(__file__).with_name("test_auxiliary.cpp")
    try:
        compiled = run_command(command(context, "cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
            "-I", str(repo / "src" / "pic"),
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
            "-o", str(executable),
        ], repo, timeout=timeout)
        if compiled.returncode != 0:
            return write_result(
                context, "FAIL", "PIC auxiliary fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in requested])
        argv=[str(executable),"--order","normal"]
        for identifier in requested:
            argv.extend(["--case",identifier])
        run=run_command(argv,work,timeout=timeout)
        artifacts.mkdir(parents=True,exist_ok=True)
        transcript=artifacts / "pic-auxiliary-output.txt"
        transcript.write_text(run.stdout or "",encoding="utf-8")
        found=re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                         run.stdout or "",flags=re.MULTILINE)
        observed=[identifier for identifier,_,_ in found]
        results=[]
        for identifier in requested:
            exact=[record for record in found if record[0]==identifier]
            if len(exact)==1:
                results.append(case(identifier,exact[0][1],exact[0][2] or "ok"))
            else:
                results.append(case(identifier,"FAIL",
                    "expected one exact record; observed {}".format(len(exact))))
        exact_contract=(observed==requested and len(set(observed))==len(observed))
        if not exact_contract:
            results[0]=case(results[0]["id"],"FAIL",
                "record contract violation: expected {} observed {}".format(
                    requested,observed))
        passed=run.returncode==0 and all(item["status"]=="PASS" for item in results)
        return write_result(context,"PASS" if passed else "FAIL",
            "{}/{} PIC auxiliary cases passed".format(
                sum(item["status"]=="PASS" for item in results),len(results)),
            cases=results,evidence=[str(transcript)])
    except Exception as error:
        return error_result(context,error)


if __name__=="__main__":
    main_guard(main)
