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
import time

from campaign_budget import CampaignBudget
import m2_protocol
from build_m2_qualification import verify_retained_executable

from integration_common import (REPO,check,manifest,parse_receipt,
    prepare_na_ar_input,replace_moon_setting,run,sha256,moon_setting)
from m1_execution import execute_levels

# The existing I18 acceptance contract freezes this base seed. Numerical
# replicate seeds belong only to the separately declared 36/288-case design.
# Never let a generated numerical template silently change the baseline or
# select a different seed after a stochastic comparison has been examined.
BASELINE_SEED=731993


class SupportingGateStopped(Exception):
    """An already recorded FAIL/SKIPPED halts dependent work, not an ERROR."""


def remaining_seconds():
    """Every phase consumes the same enclosing monotonic work deadline."""
    deadline=os.environ.get("AMPS_MOON_WORK_DEADLINE")
    return max(0.01,float(deadline)-time.monotonic()) if deadline else 3000.0


def command_log(command,path):
    """Keep exact commands and complete child streams without shell quoting."""
    with path.open("a") as stream:
        stream.write("COMMAND: "+" ".join(command)+"\n");stream.flush()
        return subprocess.run(command,cwd=REPO,stdout=stream,
                              stderr=subprocess.STDOUT,check=False,timeout=remaining_seconds()).returncode


def run_campaign(a):
    output=a.output_dir.resolve()
    if output.exists() and any(output.iterdir()):raise ValueError("output must be new or empty")
    output.mkdir(parents=True,exist_ok=True)
    records=[];supporting=[];built=False
    def persist():
        """Publish completed evidence after every phase for a forced shutdown."""
        counts=lambda rows:{v:sum(r.get("status")==v for r in rows) for v in ("PASS","FAIL","SKIPPED","ERROR")}
        value={"schema":"srcMoon-M2-development-campaign-v2","formal_m2_status":"NOT QUALIFIED",
               "counts":counts(records),"results":records,"supporting_counts":counts(supporting),
               "supporting_checks":supporting,"roadmap_sha256":sha256(REPO/"srcMoon/development_and_validation.md")}
        temporary=output/"development_summary.tmp"
        temporary.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")
        temporary.replace(output/"development_summary.json")
    persist()
    executable=a.amps.resolve() if a.amps else output/"build/amps-moon-m2-na-ar"
    log=output/"build.log"
    try:
        code=command_log([sys.executable,"-B",str(REPO/"srcMoon/test/integration-tests/I19/test_budget.py")],output/"budget-infrastructure.log")
        if code:raise RuntimeError("budget infrastructure tests failed; see budget-infrastructure.log")
        code=command_log([sys.executable,"-B",str(REPO/"srcMoon/test/test_m2_protocol.py")],output/"protocol-infrastructure.log")
        if code:raise RuntimeError("M2 protocol/deadline checks failed; inspect protocol-infrastructure.log")
        # Completion-control fixtures exercise real orchestration/metadata
        # refusal paths, never supply observations or bypass formal gates.
        # Keep them in the SINGLE global development entry so new backend
        # code cannot silently escape routine regression coverage.
        for label,relative in (
            ("formal-contract-set","test_m2_contracts.py"),
            ("retained-build-identity","test_m2_build_identity.py"),
            ("illumination-contract","test_m2_illumination.py"),
            ("illumination-spatial-processor","data_processing/test_i15_spatial_support.py"),
            ("i16-grey-boundary-control","test_i16_grey_boundary.py"),
            ("i17-holdout-control","integration-tests/I17/test_holdout.py"),
            ("i17-launcher-control","integration-tests/I17/test_infrastructure.py"),
            ("i19-formal-control","test_m2_formal_campaign.py"),
            ("i20-polar-control","test_m2_polar.py"),
            ("native-operator-control","test_m2_operator_checks.py")):
            runner=REPO/"srcMoon/test"/relative;test_log=output/(label+".log")
            code=command_log([sys.executable,"-B",str(runner)],test_log)
            supporting.append({"check":label,"status":"PASS" if code==0 else "ERROR",
                "scope":"supporting-infrastructure-not-observational-validation",
                "runner":str(runner),"log":str(test_log)})
            persist()
            if code:raise RuntimeError(label+" infrastructure checks failed; inspect "+str(test_log))
        if a.amps is None:
            # Baseline requalification uses its own retained configurations and
            # binaries. It is inside the SAME command budget as M2, including
            # all clean native variants, rather than a hidden pre-test build.
            # This coordinator owns configuration restoration from the first
            # native mutation, including a failed M0 child. Avoid building an
            # intermediate normal Moon binary that the following clean M2
            # configuration would immediately discard. All M0 variants/gates
            # still run; the final ordinary build is restored in our finally.
            built=True
            code=command_log([sys.executable,"-B",str(REPO/"srcMoon/test/run_integration_tests.py"),
                "--no-restore-build","--output-dir",str(output/"m0")],output/"m0.log")
            if code:raise RuntimeError("required M0 requalification failed; inspect m0.log")
            baseline=json.loads((output/"m0/baseline_summary.json").read_text())
            supporting.append({"check":"M0-prerequisite","status":baseline["m0_gate_status"],
                "scope":"required-linked-baseline","runner":str(REPO/"srcMoon/test/run_integration_tests.py"),
                "log":str(output/"m0.log"),"counts":baseline["counts"]})
            persist()
            code=command_log([sys.executable,"-B",str(REPO/"srcMoon/test/build_m2_qualification.py"),
                "--output-dir",str(output/"build"),"--spice-path",str(a.spice_path),
                "--spice-kernels",str(a.spice_kernels)],log)
            if code:raise RuntimeError(f"M2 build failed; see {log}")
        if not executable.is_file():raise FileNotFoundError(str(executable))
        verify_retained_executable(executable)
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
        # Retain a separate physical receipt for the original supporting and
        # formal adapters. The numerical template's complete I19-only group
        # must not dispatch I14/I18 or replace the old positive-coverage I19
        # fixture. The actual sweep below still receives the unmodified group.
        production_template=output/"production-template.in"
        prepare_na_ar_input(a.input.resolve(),production_template)
        production_template.write_text(m2_protocol.without_verification_controls(
            production_template.read_text()))
        production_template.write_text(replace_moon_setting(
            production_template.read_text(),"random_seed",str(BASELINE_SEED),production_template))
        for test_id in ("I14","I18","I19"):
            directory=output/"supporting"/test_id
            directory.mkdir(parents=True);(directory/"runtime").mkdir()
            receipt=directory/"receipt.txt";application=directory/"application.in"
            prepare_na_ar_input(production_template,application)
            application.write_text(replace_moon_setting(application.read_text(),
                "surface_inventory_checkpoint_file","none",application))
            cmd=["mpiexec","-n",str(a.mpi_ranks),str(executable),
                 "--moon-integration-test",test_id,"--integration-output",str(receipt)]
            if test_id=="I19":cmd.append("--m2-budget-check")
            cmd.extend(["-input",str(application)])
            identity=manifest(test_id,executable,cmd,
                int(moon_setting(application,"random_seed")),application,a.mpi_ranks)
            identity.update(execution_scope="supporting-M2-development",
                            species_fixture="Na,Ar",
                            formal_test_status="SKIPPED",formal_reason="I13 not qualified by this campaign")
            (directory/"run_manifest.json").write_text(json.dumps(identity,indent=2,sort_keys=True)+"\n")
            print(f"RUNNING supporting {test_id}; log={directory/'stdout.log'}",flush=True)
            monitor=directory/"monitor";monitor.mkdir()
            codes=execute_levels({test_id:{"directory":directory,"command":cmd}},
                monitor,1,min(a.test_runtime_seconds,remaining_seconds()),
                memory_gib=200000000000/1024**3,render_previews=False)
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
            persist()
            print(f"{result['status']} supporting {test_id}",flush=True)
        def numerical_result(result):
            supporting.append(result);persist()
        # Exercise the new native parser/export paths through their compiled
        # production entry points. Their explicit tiny forcing-replay fixture
        # verifies wiring/schema only; it does not qualify lunar spin-up.
        native_controls=output/"native-input-controls"
        native_code=command_log([sys.executable,"-B",str(REPO/"srcMoon/test/test_m2_native_controls.py"),
            "--build-dir",str(REPO/"build"),"--input",str(a.input.resolve()),
            "--output-dir",str(native_controls)],output/"native-input-controls.log")
        native_report=json.loads((native_controls/"result.json").read_text())
        native_status=native_report.get("status")
        if native_status not in ("PASS","FAIL","SKIPPED","ERROR") or native_code!={
                "PASS":0,"FAIL":1,"SKIPPED":77,"ERROR":2}[native_status]:
            raise RuntimeError("native input-control exit/result status disagree")
        supporting.append({"check":"native-m2-input-controls","status":native_status,
            "scope":"supporting-production-parser-not-validation",
            "runner":str(REPO/"srcMoon/test/test_m2_native_controls.py"),"log":str(output/"native-input-controls.log")})
        persist()
        if native_status in ("FAIL","SKIPPED"):raise SupportingGateStopped()
        if native_status=="ERROR":raise RuntimeError("native M2 input-control checks errored")
        from run_m2_operator_checks import run_checks as operator_checks
        operator_results=operator_checks(executable,a.input.resolve(),output/"native-operators",
            mpi_ranks=a.mpi_ranks,wall_seconds=remaining_seconds())
        supporting.extend(operator_results["results"]);persist()
        if operator_results["status"] in ("FAIL","SKIPPED"):raise SupportingGateStopped()
        if operator_results["status"]=="ERROR":
            raise RuntimeError("native M2 operator checks errored; preserve native-operators outputs")
        m2_protocol.execute_design(executable,a.input.resolve(),output/"numerical",
            tier="development",mpi_ranks=a.mpi_ranks,workers=2,
            remaining_seconds=remaining_seconds,on_result=numerical_result)
        # A deterministic source-only reference is cheap enough to include
        # in the same development deadline. It advances no model particles.
        # Request the actual I19 full initial-column CSV before the source
        # integral's early return, so the global runner also exercises this
        # production export branch. Its explicitly unqualified candidate is
        # diagnostic evidence; it cannot create initial-state approval or
        # promote the normalization reference to formal I19/I20 PASS.
        from prepare_m2_normalization import prepare as prepare_normalization
        normalization=prepare_normalization(executable,a.input.resolve(),output/"source-normalization",
            horizon_s=2550240.0,mpi_ranks=a.mpi_ranks,wall_s=remaining_seconds(),
            record_initial_state=True)
        supporting.append({"check":"source-normalization-reference","status":"PASS" if normalization["status"]=="READY" else "ERROR",
            "scope":"deterministic-model-input-reference-not-qualification",
            "runner":str(REPO/"srcMoon/test/prepare_m2_normalization.py"),
            "log":str(output/"source-normalization/stdout.log")})
        persist()
        supporting.append({"check":"cluster-qualification","status":"SKIPPED",
            "scope":"deferred-long-qualification","reason":"cluster tier not selected; short diagnostics cannot close I19/I20",
            "runner":str(REPO/"srcMoon/test/run_m2_cluster.py"),"log":str(output)})
        persist()
        # Run the existing per-ID policy to preserve REAL predecessor-blocked
        # artifacts. No synthetic I13 PASS is created and no flag bypasses it.
        for number in range(14,21):
            test_id=f"I{number:02d}"
            run(test_id,executable,output/"formal",application_input=production_template,mpi_ranks=a.mpi_ranks)
            records.append(json.loads((output/"formal"/test_id/"result.json").read_text()))
    except SupportingGateStopped:
        # The failed/unavailable child and its exact log were already saved.
        # Preserve that scientific disposition while still restoring the
        # ordinary build in finally. Do not manufacture an infrastructure
        # ERROR merely because an evaluated gate correctly stopped execution.
        persist()
    except (OSError,ValueError,KeyError,RuntimeError,subprocess.TimeoutExpired) as error:
        supporting.append({"check":"development-infrastructure","status":"ERROR",
            "scope":"supporting-infrastructure","reason":str(error),
            "runner":str(Path(__file__).resolve()),"log":str(output)})
        persist() # Preserve the fault before a potentially interrupted restore.
    finally:
        if built and remaining_seconds()>1:
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
            else:
                # Actual compiled evidence, not only a metadata fixture:
                # ordinary Na-only restoration must leave the retained Na/Ar
                # executable verifiable from its OWN immutable ABI/RNG bytes.
                try:
                    retained=verify_retained_executable(executable,require_current_generated_tree=False)
                    try:
                        verify_retained_executable(executable)
                    except ValueError:
                        live_rejected=True
                    else:live_rejected=False
                    if not live_rejected:raise ValueError("ordinary restored tree unexpectedly passed the dedicated M2 probe-ABI check")
                    (output/"retained-build-after-restoration.json").write_text(json.dumps({
                        "status":"PASS","retained_executable_sha256":retained["executable"]["sha256"],
                        "strict_live_tree_rejected":live_rejected,"scope":"build-identity-only"},indent=2)+"\n")
                    supporting.append({"check":"retained-build-after-restoration","status":"PASS",
                        "scope":"supporting-build-identity-not-science","runner":str(REPO/"srcMoon/test/build_m2_qualification.py"),
                        "log":str(output/"retained-build-after-restoration.json")})
                except (OSError,ValueError,KeyError) as error:
                    supporting.append({"check":"retained-build-after-restoration","status":"ERROR",
                        "scope":"supporting-build-identity-not-science","reason":str(error),
                        "runner":str(REPO/"srcMoon/test/build_m2_qualification.py"),"log":str(output/"restore.log")})
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
    p.add_argument("--campaign-wall-seconds",type=float,default=3600)
    p.add_argument("--cleanup-reserve-seconds",type=float,default=600)
    p.add_argument("--budget-worker",action="store_true",help=argparse.SUPPRESS)
    p.add_argument("--spice-path",type=Path,default=Path("/home/vtenishe/SPICE/cspice"))
    p.add_argument("--spice-kernels",type=Path,default=Path("/home/vtenishe/SPICE/Kernels"))
    a=p.parse_args()
    if not 1<=a.mpi_ranks<=1009 or not math.isfinite(a.test_runtime_seconds) or a.test_runtime_seconds<=0:
        p.error("1..1009 ranks and positive finite runtime required by frozen seed stride")
    if a.budget_worker:
        raise SystemExit(run_campaign(a))
    if not 0<a.campaign_wall_seconds<=3600 or not 0<=a.cleanup_reserve_seconds<a.campaign_wall_seconds:
        p.error("cleanup reserve must be nonnegative and smaller than whole-command budget")
    output=a.output_dir.resolve()
    if output.exists() and any(output.iterdir()):p.error("output must be new or empty")
    output.mkdir(parents=True,exist_ok=True)
    budget=CampaignBudget(a.campaign_wall_seconds,a.cleanup_reserve_seconds)
    command=[sys.executable,"-B",str(Path(__file__).resolve()),"--budget-worker",
        "--output-dir",str(output/"campaign"),"--input",str(a.input.resolve()),
        "--mpi-ranks",str(a.mpi_ranks),"--test-runtime-seconds",str(a.test_runtime_seconds),
        "--spice-path",str(a.spice_path),"--spice-kernels",str(a.spice_kernels)]
    if a.amps:command.extend(["--amps",str(a.amps.resolve())])
    environment=dict(os.environ,AMPS_MOON_WORK_DEADLINE=str(budget.started+budget.wall_s-budget.reserve_s))
    print(f"M2 whole-command budget {a.campaign_wall_seconds:g} s, cleanup reserve {a.cleanup_reserve_seconds:g} s; log={output/'campaign.log'}",flush=True)
    supervision=budget.run(command,REPO,output/"campaign.log",output=output,environment=environment)
    child=output/"campaign/development_summary.json"
    result=json.loads(child.read_text()) if child.is_file() else {
        "schema":"srcMoon-M2-development-campaign-v2","results":[],"supporting_checks":[]}
    inconsistent_exit=(supervision["returncode"] not in (0,) and
        not any(r.get("status") in ("FAIL","ERROR") for r in result["results"]+result["supporting_checks"]))
    if supervision["reason"] or not child.is_file() or inconsistent_exit:
        result["supporting_checks"].append({"check":"whole-campaign-execution","status":"ERROR",
            "scope":"infrastructure","reason":supervision["reason"] or f"worker exit {supervision['returncode']} without a matching final disposition",
            "runner":str(Path(__file__).resolve()),"log":str(output/"campaign.log")})
    counts=lambda rows:{v:sum(r.get("status")==v for r in rows) for v in ("PASS","FAIL","SKIPPED","ERROR")}
    result.update(counts=counts(result["results"]),supporting_counts=counts(result["supporting_checks"]),
        execution_budget=supervision,formal_m2_status="NOT QUALIFIED")
    (output/"development_summary.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    totals=counts(result["results"]+result["supporting_checks"])
    print(json.dumps({"counts":result["counts"],"supporting_counts":result["supporting_counts"],
                      "wall_seconds":supervision["wall_seconds"],"summary":str(output/"development_summary.json")},indent=2))
    raise SystemExit(2 if totals["ERROR"] else 1 if totals["FAIL"] else 0)
