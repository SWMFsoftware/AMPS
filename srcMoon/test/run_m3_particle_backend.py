#!/usr/bin/env python3
"""Clean-build and verify the real conservative PIC reaction backend at1/8MPI.

Run in an isolated checkout with --isolated-checkout. The selected checkout's
generated build/configuration changes; the caller's active checkout need not.
Analytical nominal-mass/zero-Q controls are NOT physical lunar CEX/ionization
kinematics, observations, or formal I21--I24 closure. No empirical law is added.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity, replace_once, sha256
from integration_common import REPO
sys.path.insert(0, str(Path(__file__).resolve().parent / "common"))
from moon_testlib import _probe_compile_command


def derive_input(source: str, application: str) -> str:
    """Select declared numerical controls only; no active kinetic process.

    The eight actual species retain their repository masses/charges. Local
    reference weights and local timesteps are required to expose accidental
    W/dt birth resampling; individual corrections preserve represented number.
    Unique parser/kernel include names prevent collisions with other builders.
    """
    text=replace_once(source,r"^(\s*SpeciesList\s*=\s*)Na(\s*!.*)$",
        r"\1Na,NA_PLUS,He,HE_PLUS,Ne,NE_PLUS,Ar,AR_PLUS\2","canonical Na species list")
    for key,old,new in (("ParticleWeightCorrectionMode","off","on"),
            ("ParticleWeightMode","SpeciesGlobalParticleWeight","SpeciesLocalParticleWeight"),
            ("TimeStepMode","SpeciesGlobalTimeStep","SpeciesLocalTimeStep"),
            ("VelocityDistributionSampling","on","off")):
        text=replace_once(text,rf"^(\s*{key}\s*=\s*){old}\b",rf"\g<1>{new}",key)
    # Canonical Moon velocity-sampler coordinates lie outside the analytical
    # unit box. Disable that unrelated diagnostic, retaining ordinary native
    # Init_AfterParser and all allocator/mesh/weight/mover production code.
    text=replace_once(text,r"^#block[ \t]+\.\./moon\.input\.pl[ \t]*[^\r\n]*$",
        "#block ../"+application+".input.pl ! derived analytical PIC backend fixture","Moon parser block")
    return replace_once(text,r"^\s*#include\s+moon\.input\.spicekernels\s*$",
        "#include "+application+".input.spicekernels","Moon kernel include")


def _version(command):
    value=subprocess.run(command,cwd=REPO,capture_output=True,text=True,check=False)
    return {"argv":command,"returncode":value.returncode,"stdout":value.stdout,"stderr":value.stderr}


def read_layout(directory: Path, ranks: int):
    """Require completed real native evidence, never infer PASS from exit0.

    Each rank must evaluate the complete named contract. The reported native
    count independently agrees with CSV row count. Duplicate/repeated product
    checks are retained and matched across MPI layouts, not silently collapsed.
    """
    required={"actual_eight_species_and_registered_storage","actual_atomic_mass_charge_and_species_map","real_owned_block_and_empty_cell_list",
        "actual_unequal_species_local_timesteps","complete_packet_capacity_deferral",
        "deferral_preserves_parent_list_count_and_ledger","same_plan_retry_after_real_capacity_release",
        "atomic_parent_delete_product_list_and_ledger","committed_exact_integer_and_kinematic_ledger",
        "double_commit_refused_without_publication","destructor_retains_committed_products",
        "real_private_slots_reserved_before_publish","staged_records_unpublished_parent_unchanged",
        "explicit_rollback_releases_all_private_slots","rollback_destructor_idempotent",
        "changed_opaque_parent_bytes_refuse_commit","failed_commit_rollback_preserves_restored_parent",
        "same_plan_retry_after_parent_invariant_restored","stale_catalog_state_identity_refused",
        "stale_reference_weight_refused","stale_parent_position_refused",
        "detached_retain_parent_preserves_mover_ownership","detached_products_allocated_but_unlinked",
        "detached_destructor_does_not_delete_parent_or_products","mover_owner_deletes_retained_parent_exactly_once",
        "detached_continuation_owner_can_release_products","retain_parent_requires_detached_placement",
        "rate_only_cannot_publish_pic_particles","all_native_slots_and_cell_lists_clean",
        "interleaved_ledger_allocation_failure_before_physical_publish","interleaved_shared_ledger_same_plan_retry",
        "native_nonzero_analytical_mechanical_commit"}
    required.update({
        "collector_rejects_actual_pic_mass_mismatch_before_install",
        "collector_installs_explicit_nominal_state_resolver","collector_reinstall_refused",
        "actual_sphere_dispatch_records_collector_and_returns_delete",
        "collector_callback_preserves_mover_owned_particle_and_remaining_interval",
        "collector_uses_actual_local_weight_and_individual_correction",
        "collector_actual_charge_and_nuclear_inventory",
        "collector_distinguishes_actual_impact_and_frozen_geometry_epochs",
        "collector_actual_sphere_frame_receipt",
        "collector_separates_kinetic_and_formation_energy_without_heat_claim",
        "collector_caller_deletes_mover_owned_incident_exactly_once",
        "collector_neutral_bypass_preserves_bytes_and_ledger",
        "collector_native_slots_and_lists_clean"})
    required_products={"complete_product_packet","actual_compiled_and_catalog_species",
        "original_local_interval_not_daughter_dt","fresh_product_clock_and_internal_state",
        "represented_number_without_dt_birth_factor","mover_epoch_matches_physical_event",
        "opaque_lineage_attributes_preserved","actual_product_position_and_velocity"}
    all_checks=[];all_products=[];all_species=[];rank_signatures=[]
    transcript=(directory/"stdout.log").read_text()
    import re
    summaries=re.findall(r"M3 native particle backend rank=(\d+) checks=(\d+) status=(PASS|FAIL|ERROR)",transcript)
    if len(summaries)!=ranks or {int(row[0]) for row in summaries}!=set(range(ranks)):
        raise ValueError("missing/duplicate/incomplete native per-rank completion summaries")
    summary={int(row[0]):(int(row[1]),row[2]) for row in summaries}
    for rank in range(ranks):
        prefix=directory/("backend.rank"+str(rank))
        with Path(str(prefix)+".checks.csv").open(newline="") as stream:
            reader=csv.DictReader(stream)
            if reader.fieldnames!=["rank","check","status"]:raise ValueError("native check CSV header mismatch")
            checks=list(reader)
        if not checks or any(row["rank"]!=str(rank) or row["status"] not in ("PASS","FAIL") for row in checks):
            raise ValueError("malformed/incomplete native checks")
        names=collections.Counter(row["check"] for row in checks)
        if not (required|required_products)<=names.keys() or summary[rank][0]!=len(checks):
            raise ValueError("native evidence omitted a required control/count")
        if any(names[name]!=1 for name in required) or names["complete_product_packet"]!=4 or any(
                names[name]!=8 for name in required_products-{"complete_product_packet"}):
            raise ValueError("native contract/repeated product coverage mismatch")
        if summary[rank][1]=="ERROR":raise ValueError("native rank reports runtime ERROR")
        if (summary[rank][1]=="FAIL")!=any(row["status"]=="FAIL" for row in checks):
            raise ValueError("native status contradicts completed check evidence")
        rank_signatures.append(dict(names));all_checks.extend(checks)
        with Path(str(prefix)+".products.csv").open(newline="") as stream:
            reader=csv.DictReader(stream)
            expected=["rank","context","product","compiled_species","catalog_species","mass_kg","charge_c",
                "reference_weight","individual_correction","represented_number","local_dt_s","start_s",
                "duration_s","elapsed_s","remaining_s","optical_depth","x_m","y_m","z_m","vx_m_s","vy_m_s","vz_m_s"]
            if reader.fieldnames!=expected:raise ValueError("native product CSV header mismatch")
            products=list(reader)
        if len(products)!=8 or {(row["context"],row["product"]) for row in products}!={
                (name,str(i)) for name in ("capacity-retry","opaque-retry","detached","elastic-control") for i in range(2)}:
            raise ValueError("native product evidence does not cover complete packets/contexts")
        for row in products:
            if row["rank"]!=str(rank):raise ValueError("native product rank mismatch")
            values={key:float(row[key]) for key in expected[5:]}
            if not all(math.isfinite(value) for value in values.values()) or values["mass_kg"]<=0:
                raise ValueError("nonfinite/invalid native product values")
        all_products.extend(products)
        with Path(str(prefix)+".species.csv").open(newline="") as stream:
            reader=csv.DictReader(stream)
            if reader.fieldnames!=["rank","compiled_species","mass_kg","charge_c"]:raise ValueError("species CSV header mismatch")
            species=list(reader)
        if len(species)!=8 or {int(row["compiled_species"]) for row in species}!=set(range(8)) or any(
                row["rank"]!=str(rank) or not math.isfinite(float(row["mass_kg"])) or float(row["mass_kg"])<=0 or
                not math.isfinite(float(row["charge_c"])) for row in species):
            raise ValueError("native species inventory malformed/incomplete")
        all_species.extend(species)
    if any(signature!=rank_signatures[0] for signature in rank_signatures):
        raise ValueError("different ranks evaluated different native controls")
    return {"status":"FAIL" if any(row["status"]=="FAIL" for row in all_checks) else "PASS",
        "mpi_ranks":ranks,"checks":all_checks,"check_count":len(all_checks),"products":all_products,
        "species":all_species,"named_control_counts":rank_signatures[0]}


def run(output,*,isolated_checkout=False,wall_seconds=3600,build_jobs=8,
        spice_path="/home/vtenishe/SPICE/cspice",spice_kernels="/home/vtenishe/SPICE/Kernels",retained_build=None):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();budget=CampaignBudget(wall_seconds,min(180,wall_seconds*.1))
    result={"schema":"srcMoon-M3-native-particle-backend-v1","status":"ERROR",
        "scope":"supporting real PIC allocator/list adapter; analytical nominal-mass topology only",
        "formal_I21_I24_status":"NOT EVALUATED","observational_scoring":False,"commands":[]}
    identity={"schema":"srcMoon-M3-native-particle-backend-manifest-v1","execution":"not-run",
        "source_root":str(REPO),"source_revision":subprocess.run(["git","rev-parse","HEAD"],cwd=REPO,
            capture_output=True,text=True,check=False).stdout.strip(),
        "fixture":"ANALYTICAL-NOMINAL-MASS-STORAGE-NOT-CEX","random_seed":None,
        "random_control":"explicit fixed open uniform clock inputs; no transaction RNG",
        "external_validation_data":[],"uses_spice_kernels":False,
        "spice_role":"build dependency only; geometry/ephemerides not initialized by probe",
        "mpi_layouts":[1,8],"omp_threads":1,"created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    created=[];application="moon_m3_backend_"+uuid.uuid4().hex[:12]
    def execute(command,label,cwd=REPO,*,native=False):
        execution=budget.run(command,cwd,output/(label+".log"),output=output,
            allocated_cores=8 if label.startswith("np8/") else 1 if label.startswith("np1/") else build_jobs,
            environment=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1"))
        result["commands"].append(execution)
        if execution["returncode"]!=0 and not (native and execution["returncode"]==1 and not execution.get("reason")):
            raise RuntimeError(label+" failed: "+str(execution.get("reason") or execution["returncode"]))
        return execution
    try:
        if not isolated_checkout:
            result.update(status="SKIPPED",reason="explicit --isolated-checkout required before changing selected build/configuration")
            return result
        identity["production_source_hashes"]=production_source_identity()
        identity["builder"]={"path":str(Path(__file__).resolve()),"sha256":sha256(Path(__file__).resolve())}
        identity["compiler"]=_version([os.environ.get("MPICXX","mpicxx"),"--version"])
        identity["mpi"]=_version(["mpiexec","--version"])
        if retained_build is not None:
            # A successful APPLICATION build may be reused after a failed thin
            # probe/launcher attempt. This does not promote that failed native
            # result: all actual controls run again below. Require recorded
            # make success plus unchanged source/library/generated/input bytes.
            receipt=Path(retained_build).resolve();prior=json.loads(receipt.read_text())
            prior_result=json.loads(receipt.with_name("result.json").read_text())
            phases=[row for row in prior_result["commands"] if row["command"][:1]==["make"]]
            if len(phases)!=1 or phases[0]["returncode"]!=0 or phases[0].get("reason"):
                raise ValueError("retained producer has no successful normal application build")
            if Path(prior["source_root"]).resolve()!=REPO or prior["production_source_hashes"]!=production_source_identity():
                raise ValueError("retained producer source identity changed")
            library=Path(prior["library"]["path"]);build=library.parent
            if sha256(library)!=prior["library"]["sha256"]:raise ValueError("retained library bytes changed")
            for path,digest in prior["derived_compile_inputs"].items():
                if sha256(Path(path))!=digest:raise ValueError("retained compile-input snapshot corrupted")
            for path,digest in prior["generated_snapshots"].items():
                relative=Path(path).relative_to(receipt.parent/"generated")
                if sha256(Path(path))!=digest or sha256(build/relative)!=digest:
                    raise ValueError("retained generated snapshot/current ABI mismatch")
            identity["derived_compile_inputs"]=prior["derived_compile_inputs"]
            identity["retained_successful_application_build"]={"manifest":str(receipt),"sha256":sha256(receipt),
                "result":str(receipt.with_name("result.json")),"result_sha256":sha256(receipt.with_name("result.json")),
                "normal_build_phase":phases[0],"prior_native_status":prior_result["status"]}
        else:
            source=(REPO/"input/moon.input").read_text()
            staged={".input":derive_input(source,application),".input.pl":(REPO/"input/moon.input.pl").read_text(),
                ".input.spicekernels":(REPO/"input/moon.input.spicekernels").read_text()}
            for suffix,content in staged.items():
                path=REPO/"input"/(application+suffix);path.write_text(content);created.append(path)
                snapshot=output/"compile-inputs"/(application+suffix);snapshot.parent.mkdir(exist_ok=True);snapshot.write_text(content)
            identity["derived_compile_inputs"]={str(path):sha256(path) for path in (output/"compile-inputs").iterdir()}
            build=REPO/"build"
            if build.exists():shutil.rmtree(build)
            execute([str(REPO/"Config.pl"),"-application="+application,"-amps-test=off",
                "-spice-path="+str(spice_path),"-spice-kernels="+str(spice_kernels)],"configure")
            execute(["make","-j"+str(build_jobs),"OPT3=-O3 -g"],"build")
        executable=output/"native-particle-backend-probe"
        command=_probe_compile_command(build,executable)
        # PIC mesh/CMPI_channel templates are instantiated in this probe.
        # MPI_ON is an actual production preprocessing flag, not merely a
        # linker choice; omitting it silently compiles serial channel logic.
        command.insert(1,"-DMPI_ON")
        old=str(Path(__file__).resolve().parent/"common/production_kernel_probe.cpp")
        command[command.index(old)]=str(Path(__file__).resolve().parent/"m3_native_reaction_adapter_probe.cpp")
        # Normal srcMoon/makefile owns all new services. Extract exactly those
        # generated production objects from libAMPS, sharing the actual PIC
        # ABI and flags; no second source copy or alternative executor is built.
        execute(command,"compile-probe")
        snapshots={}
        for relative in ("pic/pic.h","pic/picGlobal.dfn","pic/picSpeciesMacro.dfn","pic/picParticleDataMacro.h","general/rnd.h"):
            destination=output/"generated"/relative;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(build/relative,destination);snapshots[str(destination)]=sha256(destination)
        identity.update(generated_snapshots=snapshots,executable={"path":str(executable),"sha256":sha256(executable)},
            library={"path":str(build/"libAMPS.a"),"sha256":sha256(build/"libAMPS.a")})
        layouts=[]
        for ranks in (1,8):
            directory=output/("np"+str(ranks));directory.mkdir();runtime=directory/"runtime";runtime.mkdir();(runtime/"PT").mkdir()
            command=["mpiexec","-n",str(ranks),str(executable),str(directory/"backend")]
            identity["execution"]="native PIC adapter controls"
            (output/"run_manifest.json").write_text(json.dumps(identity,indent=2,sort_keys=True)+"\n")
            execution=execute(command,"np"+str(ranks)+"/stdout",runtime,native=True)
            (directory/"stderr.log").write_text("Combined native stdout/stderr is in stdout.log.\n")
            layout=read_layout(directory,ranks)
            if (execution["returncode"]==1)!=(layout["status"]=="FAIL"):
                raise ValueError("native exit and retained evaluated outcome disagree")
            (directory/"result.json").write_text(json.dumps(layout,indent=2,sort_keys=True)+"\n");layouts.append(layout)
        if layouts[0]["named_control_counts"]!=layouts[1]["named_control_counts"]:
            raise ValueError("MPI1/MPI8 evaluated different adapter contracts")
        if identity["production_source_hashes"]!=production_source_identity():
            raise ValueError("production source changed during native build/run; retain evidence and rebuild from stable source")
        result.update(status="FAIL" if any(layout["status"]=="FAIL" for layout in layouts) else "PASS",
            layouts=[{"mpi_ranks":layout["mpi_ranks"],"status":layout["status"],"check_count":layout["check_count"]} for layout in layouts],
            executable=identity["executable"],generated_snapshots=snapshots,
            counts={key:sum(row["status"]==key for layout in layouts for row in layout["checks"]) for key in ("PASS","FAIL","ERROR","SKIPPED")})
    except (OSError,ValueError,RuntimeError,KeyError,TypeError) as error:
        result.update(status="ERROR",reason=str(error))
    finally:
        # Only unique fixture inputs and Config.pl's corresponding root copies
        # are cleaned. No wildcard touches an unrelated user/test artifact.
        for path in created:path.unlink(missing_ok=True)
        for suffix in (".input",".input.pl",".input.spicekernels"):
            (REPO/(application+suffix)).unlink(missing_ok=True)
        result["elapsed_wall_seconds"]=time.monotonic()-started
        (output/"run_manifest.json").write_text(json.dumps(identity,indent=2,sort_keys=True)+"\n")
        (output/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--isolated-checkout",action="store_true",help="explicitly permit configuration/build changes in this isolated checkout")
    parser.add_argument("--wall-seconds",type=float,default=3600)
    parser.add_argument("--build-jobs",type=int,default=8)
    parser.add_argument("--retained-build",type=Path,help="reuse only a hash-verified successful normal build; rerun all native controls")
    parser.add_argument("--spice-path",default="/home/vtenishe/SPICE/cspice")
    parser.add_argument("--spice-kernels",default="/home/vtenishe/SPICE/Kernels")
    args=parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):parser.error("output must be new or empty")
    if not 0<args.wall_seconds<=3600 or args.build_jobs<1:parser.error("wall budget must be(0,3600], build jobs positive")
    result=run(args.output_dir,isolated_checkout=args.isolated_checkout,wall_seconds=args.wall_seconds,
        build_jobs=args.build_jobs,spice_path=args.spice_path,spice_kernels=args.spice_kernels,retained_build=args.retained_build)
    print(json.dumps(result,indent=2,sort_keys=True))
    raise SystemExit({"PASS":0,"FAIL":1,"ERROR":2,"SKIPPED":77}[result["status"]])
