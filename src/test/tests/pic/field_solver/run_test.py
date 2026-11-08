#!/usr/bin/env python3
"""Compile and execute manufactured Poisson and vacuum-wave tests."""
from __future__ import annotations
import re,sys
from pathlib import Path
TEST_ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(TEST_ROOT))
from lib.runner_support import (RunnerError,case,command,error_result,load_context, # noqa:E402
                                main_guard,run_command,write_result)
EXPECTED_CASES=["PIC-FIELD-POISSON-RESIDUAL-001","PIC-FIELD-POISSON-CONVERGENCE-001",
                "PIC-FIELD-WAVE-DISPERSION-001","PIC-FIELD-WAVE-ENERGY-001",
                "PIC-FIELD-SOLVER-INVALID-001"]
def main()->int:
  context=load_context("Test manufactured field solves and vacuum waves")
  repo=Path(context["repo_root"]); work=Path(context["work_dir"]); artifacts=Path(context["artifact_dir"])
  requested=list(context["test"].get("selected_cases",context["test"].get("expected_cases",EXPECTED_CASES)))
  if not requested or len(set(requested))!=len(requested) or any(x not in EXPECTED_CASES for x in requested):
    raise RunnerError("invalid field-solver case selection: {}".format(requested))
  timeout=int(context["test"].get("timeout_seconds",120)); executable=work/"test_field_solver"
  try:
    build=run_command(command(context,"cxx")+["-std=c++17","-O2","-Wall","-Wextra","-Wpedantic","-Werror",
      "-I",str(repo/"src"/"pic"),"-I",str(TEST_ROOT/"native"),str(TEST_ROOT/"native"/"test_main.cpp"),
      str(Path(__file__).with_name("test_field_solver.cpp")),"-o",str(executable)],repo,timeout=timeout)
    if build.returncode!=0: return write_result(context,"FAIL","field fixture did not compile",
      cases=[case(x,"FAIL","compiler exit {}".format(build.returncode)) for x in requested])
    argv=[str(executable),"--order","normal"]
    for x in requested: argv.extend(["--case",x])
    run=run_command(argv,work,timeout=timeout); artifacts.mkdir(parents=True,exist_ok=True)
    transcript=artifacts/"field-solver-output.txt"; transcript.write_text(run.stdout or "",encoding="utf-8")
    records=re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",run.stdout or "",flags=re.MULTILINE)
    observed=[x for x,_,_ in records]; results=[]
    for x in requested:
      exact=[r for r in records if r[0]==x]
      results.append(case(x,exact[0][1] if len(exact)==1 else "FAIL",
        (exact[0][2] or "ok") if len(exact)==1 else "expected one exact record; observed {}".format(len(exact))))
    if observed!=requested or len(set(observed))!=len(observed):
      results[0]=case(results[0]["id"],"FAIL","record contract violation: expected {} observed {}".format(requested,observed))
    passed=run.returncode==0 and all(x["status"]=="PASS" for x in results)
    return write_result(context,"PASS" if passed else "FAIL","{}/{} field cases passed".format(
      sum(x["status"]=="PASS" for x in results),len(results)),cases=results,evidence=[str(transcript)])
  except Exception as error: return error_result(context,error)
if __name__=="__main__": main_guard(main)
