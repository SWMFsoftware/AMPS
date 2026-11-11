#!/usr/bin/env python3
"""Qualify dependency-light performance kernels without configuring AMPS.

The same production MPI helper is exercised at 1/2/4 ranks in private work
directories. Each result must appear once per launch, with a successful
launcher exit and no foreign/missing IDs; rank-local failures are reduced by
the fixture. A fast kernel with an incomplete result is never accepted.
"""
from pathlib import Path
import re
import sys

TEST_ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(TEST_ROOT))
from lib.runner_support import (case,command,error_result,load_context,
                                main_guard,run_command,write_result)

EXPECTED_CASES=["PERF-PERIODIC-LOCAL-001","PERF-PERIODIC-REMOTE-001",
    "PERF-PERIODIC-LAYOUT-001","PERF-PERIODIC-CHUNK-001",
    "PERF-PERIODIC-ALIAS-001","PERF-PERIODIC-EMPTY-001",
    "PERF-PERIODIC-COMMUNICATOR-001"]

UNIT_CASES=["PERF-STENCIL-EXACT-001","PERF-STENCIL-FALLBACK-001","PERF-STENCIL-REFRESH-001",
    "PERF-SPMV-ORDER-001","PERF-SPMV-REFRESH-001","PERF-SPMV-EMPTY-INVALID-001"]
EXPECTED_CASES+=UNIT_CASES

def main():
    context=load_context("Qualify order-preserving performance kernels")
    repo=Path(context["repo_root"]); work=Path(context["work_dir"])
    artifacts=Path(context["artifact_dir"])
    requested=list(context["test"].get("selected_cases",EXPECTED_CASES))
    executable=work/"ordered_exchange"
    timeout=int(context["test"].get("timeout_seconds",180))
    try:
        build=run_command(command(context,"mpi_cxx")+[
            "-std=c++17","-O2","-Wall","-Wextra","-Werror","-DOMPI_SKIP_MPICXX=1",
            "-I",str(repo/"src/pic"),str(Path(__file__).with_name("test_ordered_exchange.cpp")),
            "-o",str(executable)],repo,timeout=timeout)
        if build.returncode: raise RuntimeError("performance fixture compilation failed")
        periodic=[identifier for identifier in EXPECTED_CASES if identifier not in UNIT_CASES]
        observed={identifier:[] for identifier in EXPECTED_CASES}
        transcripts=[]
        for ranks in (1,2,4):
            for identifier in requested:
                if identifier in periodic: print("[CASE] RUN {} ranks={}".format(identifier,ranks),flush=True)
            run=run_command(command(context,"mpiexec")+[
                str(context["execution"].get("mpi_rank_flag","-n")),str(ranks),str(executable)],
                work,timeout=timeout)
            transcripts.append("ranks={} exit={}\n{}".format(ranks,run.returncode,run.stdout or ""))
            records=re.findall(r"^\[([A-Z0-9-]+)\] (PASS|FAIL) .*$",run.stdout or "",re.M)
            for identifier in periodic:
                rows=[status for name,status in records if name==identifier]
                observed[identifier].append(run.returncode==0 and rows==["PASS"] and
                    all(name in periodic for name,status in records))
        unit=work/"performance_unit"
        compiled=run_command(command(context,"cxx")+[
            "-std=c++17","-O2","-Wall","-Wextra","-Werror",
            "-I",str(repo/"src/pic"),"-I",str(repo/"src/general"),"-I",str(TEST_ROOT/"native"),
            str(TEST_ROOT/"native/test_main.cpp"),
            str(Path(__file__).with_name("test_complete_stencil.cpp")),
            str(Path(__file__).with_name("test_sparse_matvec.cpp")),"-o",str(unit)],
            repo,timeout=timeout)
        if compiled.returncode: raise RuntimeError("stencil fixture compilation failed")
        for identifier in UNIT_CASES: print("[CASE] RUN "+identifier,flush=True)
        run=run_command([str(unit)],work,timeout=timeout)
        transcripts.append("unit kernels exit={}\n{}".format(run.returncode,run.stdout or ""))
        records=re.findall(r"^\[([A-Z0-9-]+)\] (PASS|FAIL) .*$",run.stdout or "",re.M)
        for identifier in UNIT_CASES:
            observed[identifier].append(run.returncode==0 and
                [status for name,status in records if name==identifier]==["PASS"])
        artifacts.mkdir(parents=True,exist_ok=True)
        transcript=artifacts/"performance-kernels.log"
        transcript.write_text("\n".join(transcripts))
        cases=[case(identifier,"PASS" if all(observed[identifier]) else "FAIL",
                    "independent stencil/sparse reference" if identifier in UNIT_CASES else
                    "independent payload/order matrix at 1/2/4 ranks") for identifier in requested]
        status="PASS" if all(row["status"]=="PASS" for row in cases) else "FAIL"
        return write_result(context,status,"ordered transport qualification",cases=cases,
                            evidence=[str(transcript)],metrics={"rank_matrix":[1,2,4]})
    except Exception as error: return error_result(context,error)

if __name__=="__main__": main_guard(main)
