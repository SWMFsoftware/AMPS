#!/usr/bin/env python3
"""Verify the actual PIC-independent conservative reaction primitives.

This is supporting coverage of existing chemistry contracts, not a new formal
U/I ID or a substitute for actual PIC allocation/mover checks. The C++ control
probe links the production ReactionExecution.cpp and upstream structural
validator. Its mechanical labels, masses, hazards and bath exchanges are
explicit analytical controls, never observations or approved lunar kinetics.
No Config.pl/generated/application tree is configured or modified.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO=Path(__file__).resolve().parents[2]
DEPENDENCIES=("srcMoon/species/ReactionExecution.h",
    "srcMoon/species/ReactionExecution.cpp",
    "src/models/chemical_solver/ReactionNetwork.h",
    "srcMoon/test/reaction_execution_probe.cpp",
    "srcMoon/test/M3_REACTION_EXECUTION_CONTRACT.md",
    "srcMoon/test/test_m3_reaction_execution.py")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(output,compiler):
    """Retain commands and genuine outcomes; process success alone is not PASS."""
    output=Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory must be new or empty: "+str(output))
    output.mkdir(parents=True,exist_ok=True)
    result={"schema":"srcMoon-M3-reaction-execution-controls-v1",
        "scope":"supporting-production-primitives-not-PIC-or-kinetic-qualification",
        "formal_I21_I24_status":"NOT EVALUATED",
        "fixture_notice":"Analytical mechanical/protocol controls, not observations or lunar daughter-kinematics evidence",
        "status":"ERROR","results":[],"counts":{},"commands":[]}
    manifest={"schema":"srcMoon-M3-reaction-execution-manifest-v1",
        "source_revision":subprocess.check_output(["git","rev-parse","HEAD"],cwd=REPO,text=True).strip(),
        "source_hashes":{name:sha(REPO/name) for name in DEPENDENCIES},
        "external_data":[],"random_seed":731993,"MPI":"none; standalone production primitive verification",
        "required_later_checks":"Actual PIC lists/allocator/ledger/MPI/movers and reviewed physical kinetics; formal predecessors unchanged"}
    try:
        resolved=shutil.which(compiler)
        if not resolved:
            result.update(status="SKIPPED",reason="C++17 compiler unavailable: "+compiler)
            return result
        version=subprocess.run([resolved,"--version"],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=False,timeout=15)
        (output/"compiler-version.log").write_text(version.stdout)
        manifest["compiler"]={"path":resolved,"version":version.stdout,"returncode":version.returncode}
        if version.returncode:
            raise RuntimeError("compiler identity query failed")
        executable=output/"reaction-execution-probe"
        command=[resolved,"-std=c++17","-O2","-Wall","-Wextra","-Werror","-pedantic","-I",str(REPO),
            str(REPO/"srcMoon/species/ReactionExecution.cpp"),
            str(REPO/"srcMoon/test/reaction_execution_probe.cpp"),"-o",str(executable)]
        with (output/"compile.log").open("w") as log:
            completed=subprocess.run(command,cwd=REPO,stdout=log,stderr=subprocess.STDOUT,check=False,timeout=120)
        result["commands"].append({"command":command,"returncode":completed.returncode,"log":str(output/"compile.log")})
        if completed.returncode:
            raise RuntimeError("production primitive compilation failed; see compile.log")
        manifest["executable"]={"path":str(executable),"sha256":sha(executable)}
        command=[str(executable)]
        completed=subprocess.run(command,cwd=output,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=False,timeout=60)
        (output/"stdout.log").write_text(completed.stdout)
        (output/"stderr.log").write_text(completed.stderr)
        result["commands"].append({"command":command,"returncode":completed.returncode,
            "stdout":str(output/"stdout.log"),"stderr":str(output/"stderr.log")})
        rows=[{"name":name,"status":status} for status,name in re.findall(r"^(PASS|FAIL) (.+)$",completed.stdout,re.M)]
        summary=re.findall(r"^SUMMARY checks=(\d+) failed=(\d+) first_channel=(\d+) seed=(\d+)$",completed.stdout,re.M)
        if len(summary)!=1 or not rows or len({r["name"] for r in rows})!=len(rows):
            raise RuntimeError("probe emitted missing/duplicate/malformed control receipts")
        total,failed,first,seed=map(int,summary[0])
        if total!=len(rows) or failed!=sum(r["status"]=="FAIL" for r in rows) or seed!=731993:
            raise RuntimeError("probe summary contradicts actual control rows/seed")
        expected_exit=1 if failed else 0
        if completed.returncode!=expected_exit:
            raise RuntimeError("probe exit disagrees with actual evaluated control outcomes")
        after={name:sha(REPO/name) for name in DEPENDENCIES}
        if after!=manifest["source_hashes"]:
            raise RuntimeError("production/test source changed during compilation or execution")
        # A reproducible scalar comparison artifact. The expected counts and
        # predeclared Hoeffding bounds are an independent sampler oracle, not
        # a frozen current production output or physical channel probability.
        bound=math.sqrt(.5*100000*math.log(4/.01))
        with (output/"branch_control.csv").open("w",newline="") as stream:
            writer=csv.writer(stream)
            writer.writerow(["control_channel","trials","count","expected_count","simultaneous_99pct_lower_count","simultaneous_99pct_upper_count","scope"])
            for name,count,probability in (("first",first,.25),("second",100000-first,.75)):
                expected=100000*probability
                writer.writerow([name,100000,count,expected,expected-bound,expected+bound,"ANALYTICAL CONTROL ONLY"])
        result.update(status="FAIL" if failed else "PASS",results=rows,
            counts={"PASS":total-failed,"FAIL":failed,"ERROR":0,"SKIPPED":0},
            metrics={"evaluated_controls":total,"failed_controls":failed,
                "channel_trials":100000,"first_channel_count":first,"seed":seed},
            artifacts={"compile_log":str(output/"compile.log"),"stdout":str(output/"stdout.log"),
                "stderr":str(output/"stderr.log"),"branch_comparison":str(output/"branch_control.csv")})
    except (OSError,ValueError,RuntimeError,subprocess.TimeoutExpired) as error:
        result.update(status="ERROR",reason=str(error))
    finally:
        result["run_manifest"]=str(output/"run_manifest.json")
        (output/"run_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
        (output/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path)
    parser.add_argument("--cxx",default="g++",help="C++17 compiler executable; no application build is performed")
    args=parser.parse_args()
    output=args.output_dir or REPO/"test_output/srcMoon"/("m3-reaction-execution-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ"))
    try:
        result=run(output,args.cxx)
    except (OSError,ValueError,subprocess.SubprocessError) as error:
        print(str(error),file=sys.stderr);return 2
    print(json.dumps(result,indent=2))
    return {"PASS":0,"FAIL":1,"ERROR":2,"SKIPPED":77}[result["status"]]


if __name__=="__main__":raise SystemExit(main())
