#!/usr/bin/env python3
"""Shared adapter for the three T9 persistence/output suite owners."""
from __future__ import annotations
import json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from lib.runner_support import RunnerError,case,command,error_result,load_context,main_guard,run_command,write_result # noqa:E402
SUITES={
 "CORE-PIC-RESTART-001":["PIC-PERSISTENCE-INVENTORY-001","PIC-RESTART-ROUNDTRIP-001","PIC-RESTART-CONTINUATION-001","PIC-RESTART-EMPTY-MAX-001","PIC-RESTART-CORRUPTION-001","PIC-RESTART-ENDIAN-WIDTH-001","PIC-RESTART-ATOMIC-001","PIC-RESTART-LEGACY-POLICY-001"],
 "CORE-PIC-DATAFILE-001":["PIC-DATAFILE-BACKGROUND-001","PIC-DATAFILE-DECOMPOSITION-001","PIC-DATAFILE-COMPATIBILITY-001"],
 "CORE-PIC-OUTPUT-001":["PIC-OUTPUT-TECPLOT-LAYOUT-001","PIC-OUTPUT-NONFINITE-001"],
}
def main()->int:
 c=load_context("Test versioned persistence, restart and output integrity"); suite=c["test"]["id"]
 if suite not in SUITES: raise RunnerError("unsupported persistence suite {}".format(suite))
 repo=Path(c["repo_root"]); work=Path(c["work_dir"]); art=Path(c["artifact_dir"]); timeout=int(c["test"].get("timeout_seconds",180)); expected=SUITES[suite]
 req=list(c["test"].get("selected_cases",c["test"].get("expected_cases",expected)))
 if not req or len(set(req))!=len(req) or any(x not in expected for x in req): raise RunnerError("invalid persistence selection: {}".format(req))
 inventory=repo/"src"/"test"/"persistent_formats.json"
 try:
  document=json.loads(inventory.read_text(encoding="utf-8"))
  formats=document.get("formats")
  policy=document.get("compatibility_policy")
  if document.get("schema_version")!=1 or not isinstance(formats,list) or not formats:
   raise RunnerError("invalid persistent format inventory")
  if not isinstance(policy,dict) or set(policy)!={"enveloped_native","enveloped_canonical","legacy_native","text"}:
   raise RunnerError("persistent compatibility policy is incomplete")
  required={"id","owner","envelope","magic","version","byte_order","integer_widths","layout_identity","compatibility"}
  identifiers=[]
  for index,item in enumerate(formats):
   if not isinstance(item,dict) or not required.issubset(item):
    raise RunnerError("persistent format {} lacks required metadata".format(index))
   if not all(isinstance(item.get(key),str) and item[key].strip()
              for key in ("id","owner","envelope","compatibility")):
    raise RunnerError("persistent format {} has an empty identity/policy".format(index))
   identifiers.append(item["id"])
  if len(identifiers)!=len(set(identifiers)):
   raise RunnerError("persistent format inventory contains duplicate identifiers")
  expected_formats={"sampling_restart","particle_restart","amr_mesh",
   "amr_center_node_data","background_binary","associated_node_data",
   "trajectory_dataset","trajectory_temp_data_and_list",
   "trajectory_selected_points","stack_serialization",
   "cutcell_connectivity_scratch","tecplot_ascii"}
  if set(identifiers)!=expected_formats:
   raise RunnerError("persistent format inventory drift: {}".format(sorted(set(identifiers)^expected_formats)))
  exe=work/"test_persistence"
  define='-DAMPS_PERSISTENT_INVENTORY=\"{}\"'.format(inventory)
  build=run_command(command(c,"cxx")+["-std=c++17","-O2","-Wall","-Wextra","-Wpedantic","-Werror",define,"-I",str(repo/"src"/"pic"),"-I",str(ROOT/"native"),str(ROOT/"native"/"test_main.cpp"),str(Path(__file__).with_name("test_persistence.cpp")),"-o",str(exe)],repo,timeout=timeout)
  if build.returncode: return write_result(c,"FAIL","persistence fixture did not compile",cases=[case(x,"FAIL","compiler exit {}".format(build.returncode)) for x in req])
  argv=[str(exe),"--order","normal"]
  for x in req: argv += ["--case",x]
  run=run_command(argv,work,timeout=timeout); art.mkdir(parents=True,exist_ok=True); transcript=art/(suite+".txt"); transcript.write_text(run.stdout or "",encoding="utf-8")
  records=re.findall(r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",run.stdout or "",flags=re.MULTILINE); observed=[x for x,_,_ in records]; results=[]
  for x in req:
   exact=[r for r in records if r[0]==x]; results.append(case(x,exact[0][1] if len(exact)==1 else "FAIL",(exact[0][2] or "ok") if len(exact)==1 else "missing/duplicate record"))
  if observed!=req or len(set(observed))!=len(observed): results[0]=case(results[0]["id"],"FAIL","record contract violation: expected {} observed {}".format(req,observed))
  passed=run.returncode==0 and all(x["status"]=="PASS" for x in results)
  return write_result(c,"PASS" if passed else "FAIL","{}/{} {} cases passed".format(sum(x["status"]=="PASS" for x in results),len(results),suite),cases=results,evidence=[str(transcript),str(inventory)])
 except Exception as e: return error_result(c,e)
if __name__=="__main__": main_guard(main)
