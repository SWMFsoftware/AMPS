#!/usr/bin/env python3
"""Compile and execute dependency-light ECSIM boundary-row tests."""

from __future__ import annotations
import re
import sys
from pathlib import Path

TEST_ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(TEST_ROOT))
from lib.runner_support import (RunnerError,case,command,error_result,
    load_context,main_guard,run_command,write_result)  # noqa: E402

EXPECTED_CASES=[
    "ECSIM-BC-PERIODIC-001","ECSIM-BC-DIRICHLET-001",
    "ECSIM-BC-NEUMANN-001","ECSIM-BC-CONDUCTOR-UNSUPPORTED-001"]

def main()->int:
    context=load_context("Test ECSIM physical boundary-row policies")
    repo=Path(context["repo_root"]); work=Path(context["work_dir"])
    artifacts=Path(context["artifact_dir"])
    requested=list(context["test"].get("selected_cases",
        context["test"].get("expected_cases",EXPECTED_CASES)))
    if not requested or len(set(requested))!=len(requested) or any(
            identifier not in EXPECTED_CASES for identifier in requested):
        raise RunnerError("invalid ECSIM boundary case selection")
    executable=work/"test_ecsim_boundaries"
    try:
        compiled=run_command(command(context,"cxx")+[
            "-std=c++17","-O2","-Wall","-Wextra","-Wpedantic",
            "-I",str(repo/"src/pic/ecsim"),"-I",str(TEST_ROOT/"native"),
            str(TEST_ROOT/"native/test_main.cpp"),
            str(Path(__file__).with_name("test_boundaries.cpp")),
            "-o",str(executable)],repo,timeout=120)
        if compiled.returncode!=0:
            return write_result(context,"FAIL","boundary fixture did not compile",
                cases=[case(i,"FAIL","compiler failure") for i in requested])
        argv=[str(executable),"--order","normal"]
        for identifier in requested: argv.extend(["--case",identifier])
        run=run_command(argv,work,timeout=120)
        artifacts.mkdir(parents=True,exist_ok=True)
        transcript=artifacts/"ecsim-boundaries-output.txt"
        transcript.write_text(run.stdout or "",encoding="utf-8")
        found=re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                         run.stdout or "",flags=re.MULTILINE)
        observed=[x[0] for x in found]
        results=[]
        for identifier in requested:
            exact=[x for x in found if x[0]==identifier]
            results.append(case(identifier,exact[0][1] if len(exact)==1 else "FAIL",
                (exact[0][2] or "ok") if len(exact)==1 else "missing/duplicate receipt"))
        passed=(run.returncode==0 and observed==requested and
                all(x["status"]=="PASS" for x in results))
        return write_result(context,"PASS" if passed else "FAIL",
            "{}/{} ECSIM boundary cases passed".format(
                sum(x["status"]=="PASS" for x in results),len(results)),
            cases=results,evidence=[str(transcript)])
    except Exception as error:
        return error_result(context,error)

if __name__=="__main__": main_guard(main)
