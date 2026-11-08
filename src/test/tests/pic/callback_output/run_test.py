#!/usr/bin/env python3
from __future__ import annotations
import re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from lib.runner_support import RunnerError,case,command,error_result,load_context,main_guard,run_command,write_result # noqa:E402
EXPECTED=["PIC-CALLBACK-RUNTIME-ORDER-001","PIC-CALLBACK-PARTICLE-ONCE-001","PIC-TRAJECTORY-FORMAT-001","PIC-CALLBACK-OUTPUT-INVALID-001"]
def main()->int:
 c=load_context("Test callback dispatch and trajectory formatting"); repo=Path(c["repo_root"]); work=Path(c["work_dir"]); art=Path(c["artifact_dir"]); req=list(c["test"].get("selected_cases",c["test"].get("expected_cases",EXPECTED)))
 if not req or len(set(req))!=len(req) or any(x not in EXPECTED for x in req): raise RunnerError("invalid callback/output selection")
 timeout=int(c["test"].get("timeout_seconds",120)); exe=work/"test_callback_output"
 try:
  build=run_command(command(c,"cxx")+["-std=c++17","-O2","-Wall","-Wextra","-Wpedantic","-Werror","-I",str(repo/"src"/"pic"),"-I",str(ROOT/"native"),str(ROOT/"native"/"test_main.cpp"),str(Path(__file__).with_name("test_callback_output.cpp")),"-o",str(exe)],repo,timeout=timeout)
  if build.returncode: return write_result(c,"FAIL","callback/output fixture did not compile",cases=[case(x,"FAIL","compiler exit {}".format(build.returncode)) for x in req])
  argv=[str(exe),"--order","normal"]
  for x in req: argv += ["--case",x]
  run=run_command(argv,work,timeout=timeout); art.mkdir(parents=True,exist_ok=True); transcript=art/"callback-output.txt"; transcript.write_text(run.stdout or "",encoding="utf-8")
  records=re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",run.stdout or "",flags=re.MULTILINE); results=[]; observed=[x for x,_,_ in records]
  for x in req:
   exact=[r for r in records if r[0]==x]; results.append(case(x,exact[0][1] if len(exact)==1 else "FAIL",(exact[0][2] or "ok") if len(exact)==1 else "missing/duplicate record"))
  if observed!=req or len(set(observed))!=len(observed): results[0]=case(results[0]["id"],"FAIL","record contract violation")
  passed=run.returncode==0 and all(x["status"]=="PASS" for x in results)
  return write_result(c,"PASS" if passed else "FAIL","{}/{} callback/output cases passed".format(sum(x["status"]=="PASS" for x in results),len(results)),cases=results,evidence=[str(transcript)])
 except Exception as e: return error_result(c,e)
if __name__=="__main__": main_guard(main)
