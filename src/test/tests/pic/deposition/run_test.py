#!/usr/bin/env python3
"""Compile and execute the production CIC/deposition contract tests."""

from __future__ import annotations

import re
import sys
from pathlib import Path

TEST_ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(TEST_ROOT))
from lib.runner_support import (RunnerError,case,command,error_result, # noqa: E402
    load_context,main_guard,run_command,write_result)

EXPECTED_CASES=[
    "PIC-DEPOSITION-ANALYTIC-001",
    "PIC-DEPOSITION-TOTALS-001",
    "PIC-DEPOSITION-CONTINUITY-001",
    "PIC-DEPOSITION-INVALID-001",
]

def main() -> int:
    context=load_context("Test PIC charge/current deposition kernels")
    repo=Path(context["repo_root"]); work=Path(context["work_dir"])
    artifacts=Path(context["artifact_dir"])
    timeout=int(context["test"].get("timeout_seconds",120))
    requested=list(context["test"].get("selected_cases",
        context["test"].get("expected_cases",EXPECTED_CASES)))
    if (not requested or len(set(requested))!=len(requested) or
            any(value not in EXPECTED_CASES for value in requested)):
        raise RunnerError("invalid deposition case selection: {}".format(requested))
    executable=work/"test_pic_deposition"
    try:
        built=run_command(command(context,"cxx")+[
            "-std=c++17","-O2","-Wall","-Wextra","-Wpedantic","-Werror",
            "-I",str(repo/"src"/"pic"),"-I",str(TEST_ROOT/"native"),
            str(TEST_ROOT/"native"/"test_main.cpp"),
            str(Path(__file__).with_name("test_deposition.cpp")),
            "-o",str(executable)],repo,timeout=timeout)
        if built.returncode!=0:
            return write_result(context,"FAIL","deposition fixture did not compile",
                cases=[case(value,"FAIL","compiler exit {}".format(built.returncode))
                       for value in requested])
        argv=[str(executable),"--order","normal"]
        for value in requested: argv.extend(["--case",value])
        run=run_command(argv,work,timeout=timeout)
        artifacts.mkdir(parents=True,exist_ok=True)
        transcript=artifacts/"deposition-output.txt"
        transcript.write_text(run.stdout or "",encoding="utf-8")
        records=re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                           run.stdout or "",flags=re.MULTILINE)
        observed=[value for value,_,_ in records]
        results=[]
        for value in requested:
            exact=[record for record in records if record[0]==value]
            results.append(case(value,exact[0][1] if len(exact)==1 else "FAIL",
                (exact[0][2] or "ok") if len(exact)==1 else
                "expected one exact record; observed {}".format(len(exact))))
        if observed!=requested or len(set(observed))!=len(observed):
            results[0]=case(results[0]["id"],"FAIL",
                "record contract violation: expected {} observed {}".format(
                    requested,observed))
        passed=run.returncode==0 and all(item["status"]=="PASS" for item in results)
        return write_result(context,"PASS" if passed else "FAIL",
            "{}/{} deposition cases passed".format(
                sum(item["status"]=="PASS" for item in results),len(results)),
            cases=results,evidence=[str(transcript)])
    except Exception as error:
        return error_result(context,error)

if __name__=="__main__": main_guard(main)
