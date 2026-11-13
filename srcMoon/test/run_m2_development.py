#!/usr/bin/env python3
"""Bounded supporting M2 verification before M1 closure.

Formal I14-I20 retain their existing predecessor checks and SKIPPED results.
Supporting production checks are a separate result namespace: their PASS
cannot establish numerical convergence, observational validation or M2 PASS.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

from integration_common import (REPO,check,manifest,parse_receipt,
    prepare_na_ar_input,replace_moon_setting,run,sha256)
from m1_execution import execute_levels


def command_log(command,path):
    """Keep exact commands and complete child streams without shell quoting."""
    with path.open("a") as stream:
        stream.write("COMMAND: "+" ".join(command)+"\n");stream.flush()
        return subprocess.run(command,cwd=REPO,stdout=stream,
                              stderr=subprocess.STDOUT,check=False).returncode


def run_campaign(a):
    output=a.output_dir.resolve()
    if output.exists() and any(output.iterdir()):raise ValueError("output must be new or empty")
    output.mkdir(parents=True,exist_ok=True)
    records=[];supporting=[];built=False
    executable=a.amps.resolve() if a.amps else output/"build/amps-moon-m2-na-ar"
    log=output/"build.log"
    try:
        code=command_log([sys.executable,"-B",str(REPO/"srcMoon/test/integration-tests/I19/test_budget.py")],output/"budget-infrastructure.log")
        if code:raise RuntimeError("budget infrastructure tests failed; see budget-infrastructure.log")
        if a.amps is None:
            built=True
            code=command_log([sys.executable,"-B",str(REPO/"srcMoon/test/build_m2_qualification.py"),
                "--output-dir",str(output/"build"),"--spice-path",str(a.spice_path),
                "--spice-kernels",str(a.spice_kernels)],log)
            if code:raise RuntimeError(f"M2 build failed; see {log}")
        if not executable.is_file():raise FileNotFoundError(str(executable))
        # Every ordinary U runner remains registered and owns its own oracle.
        local=output/"local"
        command=[sys.executable,"-B",str(REPO/"srcMoon/test/run_all_test.py"),
                 "--skip-build","--output-dir",str(local)]
        for number in range(8,14):command.extend(["--only",f"U{number:02d}"])
        command_log(command,output/"local.log")
        records.extend(json.loads((local/"run_all_summary.json").read_text())["results"])
        data=output/"data-audit"
        command_log([sys.executable,"-B",str(REPO/"srcMoon/test/data_processing/audit_m2_packages.py"),
                     "--output-dir",str(data)],output/"data-audit.log")
        audit=json.loads((data/"result.json").read_text())
        # Preserve the actual delegated runner/log even for a failed package
        # check. Data readiness and linked physics stay in separate scopes.
        for record in audit["records"]:
            record.update(runner=str(REPO/"srcMoon/test/data_processing/audit_m2_packages.py"),
                          log=str(output/"data-audit.log"))
            supporting.append(record)
        for test_id in ("I14","I18","I19"):
            directory=output/"supporting"/test_id
            directory.mkdir(parents=True);(directory/"runtime").mkdir()
            receipt=directory/"receipt.txt";application=directory/"application.in"
            prepare_na_ar_input(a.input.resolve(),application)
            application.write_text(replace_moon_setting(application.read_text(),
                "surface_inventory_checkpoint_file","none",application))
            cmd=["mpiexec","-n",str(a.mpi_ranks),str(executable),
                 "--moon-integration-test",test_id,"--integration-output",str(receipt)]
            if test_id=="I19":cmd.append("--m2-budget-check")
            cmd.extend(["-input",str(application)])
            identity=manifest(test_id,executable,cmd,731993,application,a.mpi_ranks)
            identity.update(execution_scope="supporting-M2-development",
                            species_fixture="Na,Ar",
                            formal_test_status="SKIPPED",formal_reason="I13 not qualified by this campaign")
            (directory/"run_manifest.json").write_text(json.dumps(identity,indent=2,sort_keys=True)+"\n")
            print(f"RUNNING supporting {test_id}; log={directory/'stdout.log'}",flush=True)
            monitor=directory/"monitor";monitor.mkdir()
            codes=execute_levels({test_id:{"directory":directory,"command":cmd}},
                monitor,1,a.test_runtime_seconds,render_previews=False)
            result={"check":test_id+"-production-support","id":test_id,
                    "scope":"supporting-not-formal-qualification","status":"ERROR",
                    "runner":str(Path(__file__).resolve()),"log":str(directory/"stdout.log")}
            try:
                if codes.get(test_id)!=0:raise RuntimeError(f"MPI exit {codes.get(test_id)}; see stdout/stderr/monitor")
                parsed=parse_receipt(receipt)
                if test_id=="I19":
                    spec=importlib.util.spec_from_file_location("m2_budget_oracle",
                        REPO/"srcMoon/test/integration-tests/I19/analyze_budget.py")
                    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
                    result.update(module.evaluate(directory/"m2_budget.csv"))
                    module.render(directory/"m2_budget.csv",directory)
                    if parsed.get("budget_valid")!=1:result["status"]="FAIL"
                    if (parsed.get("positive_energy_exit_classified_escape")!=1 or
                        parsed.get("bound_exit_classified_numerical")!=1):
                        result.update(status="FAIL",reason="I09-defined exit classification did not match controlled energy signs")
                else:
                    status,metrics,reason=check(test_id,parsed,directory)
                    result.update(status=status,metrics=metrics,reason=reason)
            except (OSError,ValueError,KeyError,RuntimeError) as error:
                result.update(status="ERROR",reason=str(error))
            (directory/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
            supporting.append(result)
            print(f"{result['status']} supporting {test_id}",flush=True)
        # Run the existing per-ID policy to preserve REAL predecessor-blocked
        # artifacts. No synthetic I13 PASS is created and no flag bypasses it.
        for number in range(14,21):
            test_id=f"I{number:02d}"
            run(test_id,executable,output/"formal",application_input=a.input.resolve(),mpi_ranks=a.mpi_ranks)
            records.append(json.loads((output/"formal"/test_id/"result.json").read_text()))
    except (OSError,ValueError,KeyError,RuntimeError) as error:
        supporting.append({"check":"development-infrastructure","status":"ERROR",
            "scope":"supporting-infrastructure","reason":str(error),
            "runner":str(Path(__file__).resolve()),"log":str(output)})
    finally:
        if built:
            # Restore only when this runner changed configuration. Preserve
            # the dedicated M2 executable/log before replacing root build/.
            build=REPO/"build"
            if build.exists():shutil.rmtree(build)
            # Restoring the application name does not reset persisted nightly
            # TESTMODE. The delivered normal executable must explicitly leave
            # that mode, including its default 100-step/debug/sorting policy.
            code=command_log([str(REPO/"Config.pl"),"-application=moon","-amps-test=off",
                f"-spice-path={a.spice_path}",f"-spice-kernels={a.spice_kernels}"],output/"restore.log")
            # Restore through the same bounded native build policy as the
            # dedicated M2 build; the retained physics input is unchanged.
            if not code:code=command_log(["make","-j16","OPT3=-O3 -g"],output/"restore.log")
            if code:supporting.append({"check":"restore-build","status":"ERROR",
                "scope":"supporting-infrastructure","reason":"canonical restore failed; see restore.log",
                "runner":str(Path(__file__).resolve()),"log":str(output/"restore.log")})
    counts=lambda rows:{s:sum(r.get("status")==s for r in rows) for s in ("PASS","FAIL","SKIPPED","ERROR")}
    result={"schema":"srcMoon-M2-development-campaign-v1","formal_m2_status":"NOT QUALIFIED",
            "counts":counts(records),"results":records,"supporting_counts":counts(supporting),
            "supporting_checks":supporting,"roadmap_sha256":sha256(REPO/"srcMoon/development_and_validation.md")}
    (output/"development_summary.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("results","supporting_checks")},indent=2))
    all_counts=counts(records+supporting)
    return 2 if all_counts["ERROR"] else 1 if all_counts["FAIL"] else 0


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--amps",type=Path,help="reuse an explicitly supplied M2 executable and matching current generated build")
    p.add_argument("--input",type=Path,default=REPO/"srcMoon/examples/lola_surface.in")
    p.add_argument("--mpi-ranks",type=int,default=8)
    p.add_argument("--test-runtime-seconds",type=float,default=3600)
    p.add_argument("--spice-path",type=Path,default=Path("/home/vtenishe/SPICE/cspice"))
    p.add_argument("--spice-kernels",type=Path,default=Path("/home/vtenishe/SPICE/Kernels"))
    a=p.parse_args()
    if a.mpi_ranks<=0 or not math.isfinite(a.test_runtime_seconds) or a.test_runtime_seconds<=0:p.error("ranks and runtime must be positive and finite")
    raise SystemExit(run_campaign(a))
