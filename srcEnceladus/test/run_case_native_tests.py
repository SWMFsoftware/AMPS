#!/usr/bin/env python3
"""Native case CLI/MPI/transport campaign, using explicitly synthetic sources.

The executable must already be built with a model-data root containing the
fixture source set. This runner never configures or builds AMPS. The global
native owner creates the data first, serializes the build, then invokes this
runner. One/four-rank input checks allocate no particle buffer; short transport
runs check exact global budgets, not agreement with Cassini observations.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from case_fixtures import CASES

GROUPS=("N06_case_local","N07_case_mpi","N08_case_transport","N09_case_uvis","N10_runtime_mesh","N11_case_errors")


class MissingLauncher(Exception):
    """An unavailable MPI executable is SKIP, never a physics PASS."""


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable",type=Path,required=True)
    parser.add_argument("--fixture-deck",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--launcher-json",default='["mpiexec"]')
    parser.add_argument("--rank-flag",default="-n")
    parser.add_argument("--case",action="append",choices=GROUPS)
    args=parser.parse_args()
    executable=args.executable.resolve();deck=args.fixture_deck.resolve()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    launcher=json.loads(args.launcher_json)
    if not isinstance(launcher,list) or not all(isinstance(x,str) and x for x in launcher):parser.error("launcher must be an argv list")
    results=[]

    def command(label,argv,cwd=None,ok=True,timeout=180):
        print("[native case] RUN "+label,flush=True)
        start=time.monotonic()
        try:
            p=subprocess.run(list(map(str,argv)),cwd=cwd,text=True,capture_output=True,timeout=timeout)
        except subprocess.TimeoutExpired as error:
            raise AssertionError(label+" timed out; MPI failure path may be hanging") from error
        text=p.stdout+p.stderr
        (output/(label+".log")).write_text(text)
        if (p.returncode==0)!=ok:raise AssertionError(label+" exit="+str(p.returncode)+"; "+text[-1800:])
        return text,time.monotonic()-start

    def mpi(label,case,ranks=1,run_deck=deck,validate=False,extra_launcher=None,ok=True):
        if not launcher or shutil.which(launcher[0]) is None:
            raise MissingLauncher("MPI launcher unavailable; configure tools.mpiexec for native runs")
        directory=output/label;directory.mkdir();(directory/"PT/plots").mkdir(parents=True)
        argv=[*launcher,args.rank_flag,str(ranks)]
        if extra_launcher:argv+=extra_launcher
        argv += [str(executable),"-i",str(run_deck),"-case",case]
        if validate:argv += ["--validate-runtime-input"]
        text,elapsed=command(label,argv,directory,ok=ok,timeout=300)
        if not validate and ok:
            if "completed_iterations=4" not in text:raise AssertionError("native run did not reach all iterations")
            ledger=re.findall(r'final_species=(\S+) injected_macroparticles=(\d+) expected_injected_macroparticles=(\d+).*ledger_residual=(-?\d+)',text)
            if not ledger or any((actual,expected,residual)!=("40","40","0") for _,actual,expected,residual in ledger):raise AssertionError("native global particle budget differs")
            manifest=json.loads((directory/"run_manifest.json").read_text())
            if manifest["phase"]!="complete" or manifest["completed_iterations"]!=4 or manifest["mpi_ranks"]!=ranks:raise AssertionError("incomplete run snapshot")
            if not manifest.get("executable_sha256") or not (directory/"case_inputs/sources.in").is_file():raise AssertionError("native snapshot identity missing")
            # Verify retained bytes independently with Python's SHA-256 rather
            # than comparing two values emitted by the C++ digest helper.
            for asset in manifest["assets"]:
                if hashlib.sha256((directory/asset["snapshot"]).read_bytes()).hexdigest()!=asset["sha256"]:raise AssertionError("snapshot bytes differ from recorded identity")
            if not any(a["logical_path"].startswith("runtime_input/") for a in manifest["assets"]):raise AssertionError("raw runtime input snapshot missing")
            for values in manifest["species"].values():
                expected_weight=values["Q_s-1"]*manifest["global_dt_s"]/10
                if abs(values["weight"]-expected_weight)>1e-12*expected_weight:raise AssertionError("statistical weight does not match Q*dt/N")
            if manifest["minimum_cell_size_m"]<=0:raise AssertionError("realized cell size missing from run identity")
        return text,directory,elapsed

    def local():
        listing,_=command("list",[executable,"-i",deck,"--list-cases"])
        ids=re.findall(r'^READY (\S+)',listing,re.M)
        if len(ids)!=28:raise AssertionError("native list did not report 28 ready cases")
        for index,case in enumerate(ids):
            text,_=command("local_"+str(index),[executable,"--input="+str(deck),"--case="+case,"--validate-input"])
            if "input_validation=PASS" not in text:raise AssertionError("missing local case receipt")
        text,_=command("describe",[executable,"-i",deck,"-case",ids[0],"--describe-case"])
        if "pending_realized_AMR" not in text:raise AssertionError("inspection fabricated final mesh numerics")
        return "28 source/observation validations, case inventory and deferred-mesh inspection"

    def collective():
        for ranks in (1,4):
            text,_,_=mpi("validate_"+str(ranks),"E3--"+CASES[0],ranks,validate=True)
            if "mpi_input_validation=PASS" not in text or "allocated_particle_buffer=false" not in text:raise AssertionError("MPI inspection allocated particles or missed receipt")
        # Different rank-local runtime bytes must fail collectively. The
        # wrapper changes only rank 1's input path; the executable is identical.
        alternate=output/"divergent.in"
        alternate.write_text(deck.read_text().replace("random_seed=20261113","random_seed=20261114"))
        wrapper=output/"rank_wrapper.py"
        wrapper.write_text("import os,sys\nrank=int(os.environ.get('OMPI_COMM_WORLD_RANK',os.environ.get('PMI_RANK','0')))\na=sys.argv[1:]\nif rank==1:a[a.index('-i')+1]="+repr(str(alternate))+"\nos.execv(a[0],a)\n")
        prefix=[sys.executable,str(wrapper)]
        mpi("divergent_input","E3--"+CASES[0],4,validate=True,extra_launcher=prefix,ok=False)
        alternate.write_text(deck.read_text().replace("source_set=enceladus_case_validation_20261113_synthetic","source_set=missing_data").replace("source_set=synthetic_v1","source_set=missing_data"))
        mpi("missing_rank_data","E3--"+CASES[0],4,validate=True,extra_launcher=prefix,ok=False)
        # Failure must precede writing any replacement receipt or building a
        # mesh. Use byte sentinels rather than a prior production run so this
        # regression cannot accidentally damage another run's artifacts.
        for artifact in ("run_manifest.json","mesh.msh","PT/plots/old.dat"):
            directory=output/("occupied_"+artifact.split('/')[-1]);directory.mkdir()
            sentinel=directory/artifact;sentinel.parent.mkdir(parents=True,exist_ok=True)
            sentinel.write_bytes(b"preserve existing run\n")
            text,_=command(directory.name,[*launcher,args.rank_flag,"1",str(executable),"-i",str(deck),"-case","E3--"+CASES[0]],directory,ok=False)
            if "fresh output directory" not in text or sentinel.read_bytes()!=b"preserve existing run\n":raise AssertionError("case run overwrote or accepted existing output")
        return "one/four-rank agreement; divergent/missing inputs and occupied outputs fail collectively without overwriting"

    def transport():
        for i,basis in enumerate(CASES):
            for event in ("E3","E14"):mpi("transport_"+str(i)+"_"+event,event+"--"+basis)
        mpi("transport_four_ranks","E3--"+CASES[0],4)
        return "nine native runs; all bases, both primary species, exact global 40-particle ledgers"

    def uvis():
        for event in ("uvis_2007_zeta_orionis","uvis_2010_e10_solar"):
            _,directory,_=mpi("run_"+event,event+"--"+CASES[1])
            files=list((directory/"PT/plots").glob("observation.*.dat"))
            if len(files)!=1:raise AssertionError("UVIS case did not write its selected observation")
            data=files[0].read_text()
            if 'AMPS UVIS comparison '+event not in data or ('event_id='+event) not in data or data.count('ZONE ')!=2:raise AssertionError("UVIS output cadence/identity differs")
            from plot_flyby_density import read_tecplot
            variables,zones=read_tecplot(files[0])
            path=variables.index("integrated_path_in_domain_m")
            if any(len(zone.rows)!=3 or any(row[path]<=0 for row in zone.rows) for zone in zones):raise AssertionError("UVIS ray integration omitted a fixture row or domain path")
        return "stellar chord and finite-observer UVIS cases reach two sampling outputs"

    def mesh():
        # A noncubic domain plus off-center fine box exercises the aspect
        # conversion and conservative intersection criterion, not just probes.
        altered=output/"rectangular_mesh.in"
        text=deck.read_text().replace("mesh.domain_max_m=1260500,1260500,1260500","mesh.domain_max_m=1521000,1260500,1260500")
        additions="""#begin refinement: observer_box
geometry=box
min_m=350000,-30000,-230000
max_m=450000,30000,-170000
cell_size_m=80000
#end refinement
#begin refinement: near_shell
geometry=shell
center_m=0,0,0
inner_radius_m=300000
outer_radius_m=350000
cell_size_m=100000
#end refinement
#begin refinement: small_sphere
geometry=sphere
center_m=-700000,0,0
radius_m=10000
cell_size_m=100000
#end refinement
"""
        altered.write_text(text.replace("#end application",additions+"#end application"))
        log,directory,_=mpi("rectangular_mesh","E3--"+CASES[0],run_deck=altered)
        if "domain_max_m=1.521e+06" not in log and "domain_max_m=1521000" not in log:raise AssertionError("runtime domain was not applied")
        from plot_flyby_density import read_tecplot
        variables,zones=read_tecplot(next((directory/"PT/plots").glob("observation.*.dat")))
        edge=variables.index("cell_max_edge_m");inside=variables.index("in_domain")
        widths=[row[edge] for row in zones[-1].rows if row[inside]>0.5]
        if not widths or max(widths)>80000:raise AssertionError("observer box did not reach requested actual cell edges")
        manifest=json.loads((directory/"run_manifest.json").read_text())
        # Both the observer box and shell request level 3 on this root.
        minimum=min(2781500/(5*8),2521000/(5*8))
        if not abs(manifest["global_dt_s"]-0.5*minimum/4000)<1e-9:raise AssertionError("runtime mesh timestep differs from realized refinement")
        return "rectangular bounds, box/sphere/shell refinement and realized-cell Courant timestep"

    def errors():
        for label,change in (("inside_domain",("-1260500,-1260500,-1260500","-100,-100,-100")),
                             ("too_fine",("mesh.background_cell_size_m=600000","mesh.background_cell_size_m=0.01")),
                             ("unknown_mesh",("mesh.background_cell_size_m=600000","mesh.backgroun_cell_size_m=600000"))):
            invalid=output/(label+".in");invalid.write_text(deck.read_text().replace(*change))
            command(label,[executable,"-i",invalid,"-case","E3--"+CASES[0],"--validate-input"],ok=False)
        return "invalid domain, unattainable target and unknown mesh key rejected before allocation"

    functions=dict(zip(GROUPS,(local,collective,transport,uvis,mesh,errors)))
    for name in args.case or GROUPS:
        print("[CASE] RUN "+name,flush=True)
        started=time.monotonic()
        try:detail=functions[name]();status="PASS"
        except MissingLauncher as error:detail=str(error);status="SKIP"
        except Exception as error:detail=str(error);status="FAIL"
        results.append({"test":name,"status":status,"detail":detail,"seconds":time.monotonic()-started})
        print(status+" "+name+" "+detail,flush=True)
    payload={"status":"FAIL" if any(r["status"]=="FAIL" for r in results) else "PASS","results":results}
    (output/"test_receipt.json").write_text(json.dumps(payload,indent=2)+"\n")
    return 0 if payload["status"]=="PASS" else 1

if __name__=="__main__":raise SystemExit(main())
