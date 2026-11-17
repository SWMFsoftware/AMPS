#!/usr/bin/env python3
"""Prepare, explicitly run, or import a reproducible controlled warm study.

The default action only prepares inputs/provenance. --run is an explicit slow
calculation; --import-studies verifies existing evidence without solving roots.
Every physical plasma input remains mandatory in the supplied JSON. Native
centers/edges are obtained through the production value API, including their
exact double representations. Nominal-grid and native-grid coordinates remain
distinct inputs, with their differences recorded independently of root error.
"""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import tempfile
from typing import Any, Dict, List

sys.dont_write_bytecode=True
import mpmath as mp
import warm_parallel_maxwellian as warm
from plot_warm_parallel_study import verify_coordinate_record

DIRECTORY=Path(__file__).resolve().parent
ROOT=DIRECTORY.parents[2]
SCHEMA="srcsep-controlled-warm-protocol-v1"


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def atomic_bytes(path: Path,value: bytes,immutable: bool=False) -> None:
    """Publish a complete file; immutable evidence is never replaced silently."""
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and immutable:
        if path.read_bytes()!=value:
            raise ValueError("existing immutable artifact differs: "+str(path))
        return
    descriptor,temporary=tempfile.mkstemp(prefix=path.name+".",dir=str(path.parent))
    try:
        with os.fdopen(descriptor,"wb") as stream:
            stream.write(value);stream.flush();os.fsync(stream.fileno())
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def atomic_json(path: Path,value: Any,immutable: bool=False) -> None:
    atomic_bytes(path,(json.dumps(value,indent=2,allow_nan=False)+"\n").encode(),immutable)


def archive(path: Path,output: Path) -> Dict[str,str]:
    value=path.read_bytes();name="sources/"+digest(value)+"-"+path.name
    atomic_bytes(output/name,value,True)
    return {"sha256":digest(value),"archive":name}


def capture_grid(output: Path,cxx: str,overrides: List[str]) -> Dict[str,Any]:
    sources=[DIRECTORY/"export_wave_grid.cpp",ROOT/"srcSEP/util/sep_polarized_wave_state.cpp",
             ROOT/"srcSEP/util/sep_resonant_collision.cpp",ROOT/"src/models/sep_common/sep_transport_common.cpp"]
    headers=[ROOT/"srcSEP/util/sep_polarized_wave_state.h",ROOT/"srcSEP/util/sep_resonant_collision.h",
             ROOT/"srcSEP/util/sep_turbulence_core.h",ROOT/"srcSEP/util/sep_common_header_path.h",
             ROOT/"src/models/sep_common/sep_transport_common.h"]
    pinned={str(path.relative_to(ROOT)):archive(path,output) for path in sources+headers}
    compiler=shlex.split(cxx)
    if not compiler:raise ValueError("empty C++ compiler command")
    version=subprocess.run(compiler+["--version"],check=True,capture_output=True,text=True).stdout
    # The executable is a temporary tool, never a committed/reference result.
    with tempfile.TemporaryDirectory(prefix="srcsep-wave-grid-") as temporary:
        executable=Path(temporary)/"export_wave_grid"
        command=compiler+["-std=c++11","-O1","-I"+str(ROOT),"-I"+str(ROOT/"srcSEP/util")]+[
            str(path) for path in sources]+["-o",str(executable)]
        completed=subprocess.run(command,capture_output=True,text=True)
        if completed.returncode:
            raise RuntimeError("native grid exporter compilation failed:\n"+completed.stderr)
        binary_hash=digest(executable.read_bytes())
        extracted=subprocess.run([str(executable)]+overrides,check=True,capture_output=True,text=True)
    for path in sources+headers:
        if digest(path.read_bytes())!=pinned[str(path.relative_to(ROOT))]["sha256"]:
            raise RuntimeError("source changed during native-grid extraction: "+str(path))
    grid=json.loads(extracted.stdout)
    if grid.get("schema")!="srcsep-native-centered-grid-v1":raise ValueError("unknown native grid schema")
    if len(grid["centers_hex"])!=grid["bins"] or len(grid["edges_hex"])!=grid["bins"]+1:
        raise ValueError("native grid shape differs from its recorded bin count")
    for key in ("centers_hex","edges_hex"):
        values=[float.fromhex(value) for value in grid[key]]
        if any(not math.isfinite(value) or value<=0 for value in values) or any(a>=b for a,b in zip(values,values[1:])):
            raise ValueError("native grid contains invalid/unordered physical wavenumbers")
    return {"grid":grid,"source_files":pinned,"compiler_command":compiler,"compiler_version":version,
            "compiler_flags":["-std=c++11","-O1"],"executable_sha256":binary_hash,
            "capture_scope":"Production value API compiled with the recorded toolchain/configuration; not a dump of a running PIC application's arrays.",
            "grid_overrides":overrides,"platform":platform.platform()}


def coordinates(grid: Dict[str,Any],plasma: warm.Plasma,sampling: str,raw_points: Any,
                digits: int) -> Dict[str,Any]:
    with mp.workdps(digits+25):
        if raw_points is not None:
            if not isinstance(raw_points,list) or not raw_points or any(not isinstance(p,str) for p in raw_points):
                raise ValueError("raw requested coordinates must be a nonempty JSON decimal-string list")
            for point in raw_points:
                if warm.scalar(point,"requested k*d_p")<=0:raise ValueError("requested coordinates must be positive")
        length=plasma.scales()["ion_inertial_length_m"]
        actual=[grid["edges_hex"][0]]+grid["centers_hex"]+[grid["edges_hex"][-1]]
        native=[mp.mpf(float.fromhex(value))*length for value in actual]
        first=mp.mpf(grid["configured_first_center_per_m"])
        last=mp.mpf(grid["configured_last_center_per_m"])
        width=mp.log(last/first)/(grid["bins"]-1)
        nominal=[first*mp.exp(-width/2)]+[first*mp.exp(j*width) for j in range(grid["bins"])]+[last*mp.exp(width/2)]
        nominal=[point*length for point in nominal]
        if sampling=="explicit-kdi":
            if raw_points is None:raise ValueError("explicit-kdi sampling requires --points-json")
            requested=raw_points
        else:
            reference=native if sampling=="native-grid" else nominal
            requested=[mp.nstr(point,digits) for point in reference] if raw_points is None else raw_points
            if len(requested)!=len(reference):raise ValueError("requested point count differs from declared grid sampling")
            # Existing records may use fewer decimal digits. Accept only the
            # same mathematical coordinate rounded to that declared precision,
            # never a freely chosen physical comparison tolerance.
            for text,point in zip(requested,reference):
                significant=len(Decimal(text).as_tuple().digits)
                if Decimal(text)!=Decimal(mp.nstr(point,significant)):
                    raise ValueError("requested coordinate is not the rounded declared grid point: "+text)
        if not isinstance(requested,list) or not requested:raise ValueError("requested coordinates must be a nonempty JSON list")
        values=[warm.scalar(point,"requested k*d_p") for point in requested]
        if any(point<=0 for point in values) or any(a>=b for a,b in zip(values,values[1:])):
            raise ValueError("requested coordinates must be positive and strictly increasing")
        mapping=[]
        for j,(binary,point,nominal_point) in enumerate(zip(actual,native,nominal)):
            row={"index":j,"kind":"outer-edge" if j in (0,len(actual)-1) else "center",
                 "actual_k_per_m_hex":binary,"actual_k_di":mp.nstr(point,digits),
                 "nominal_k_di":mp.nstr(nominal_point,digits),
                 "native_vs_nominal_relative_difference":mp.nstr(abs(point-nominal_point)/point,digits)}
            if sampling!="explicit-kdi":
                row["requested_k_di"]=requested[j]
                row["native_vs_requested_relative_difference"]=mp.nstr(abs(point-values[j])/point,digits)
            mapping.append(row)
        return {"sampling":sampling,"requested_k_di":requested,"native_grid_map":mapping,
                "coordinate_scope":"Root numerical errors are conditional on the exact requested decimal strings. Native/nominal coordinate differences are separate inputs, not root errors or accepted physical tolerances."}


def validate_study(study: Dict[str,Any],protocol: Dict[str,Any]) -> None:
    for key,expected in (("schema","srcsep-controlled-warm-parallel-study-v1"),
                         ("input_sha256",protocol["input"]["sha256"]),
                         ("implementation_sha256",protocol["solver"]["sha256"]),
                         ("base_digits",protocol["base_digits"]),
                         ("verification_digits",protocol["verification_digits"]),
                         ("requested_k_di",protocol["requested_k_di"])):
        if study.get(key)!=expected:raise ValueError("study provenance mismatch: "+key)
    if study.get("physical_production_qualified") is not False or study.get("domain_limits_selected") is not False:
        raise ValueError("this unqualified protocol cannot import a production/domain qualification claim")
    labels=study.get("requested_helicities",[])
    if not labels or len(labels)!=len(set(labels)) or any(h not in protocol["requested_helicities"] for h in labels):
        raise ValueError("invalid or duplicate study helicity coverage")
    expected={(h,k) for h in labels for k in protocol["requested_k_di"]};seen=set()
    for key in ("python_version","mpmath_version","algorithm","root_isolation"):
        if not isinstance(study.get(key),str) or not study[key]:raise ValueError("study lacks numerical provenance: "+key)
    for row in study["rows"]:
        key=(row["spatial_helicity"],row.get("requested_k_di",row["k_di"]))
        if key not in expected or key in seen:raise ValueError("duplicate/unrequested study row")
        seen.add(key)
        verify_coordinate_record(row,key[1],protocol["verification_digits"])
        if row["status"]=="UNRESOLVED":
            if not isinstance(row.get("reason"),str) or not row["reason"]:raise ValueError("unresolved row lacks a reason")
        elif row["status"]=="RESOLVED":
            for field in ("k_di","omega_real_radps","omega_imaginary_raw_radps","frequency_numerical_error_radps",
                          "phase_speed_relative_mhd_error","group_speed_relative_mhd_error",
                          "damping_over_real_frequency_raw","weak_damping_hermitian_energy_over_magnetic_energy"):
                if not Decimal(row[field]).is_finite():raise ValueError("nonfinite resolved study field: "+field)
            if Decimal(row["omega_real_radps"])<=0 or Decimal(row["frequency_numerical_error_radps"])<0:
                raise ValueError("invalid positive-frequency/numerical-error record")
            if not isinstance(row.get("damping_resolved"),bool):raise ValueError("damping availability must be explicit")
        else:raise ValueError("unknown study-row status")
    if seen!=expected:raise ValueError("study artifact does not cover its requested points")


def combine_studies(paths: List[Path],protocol: Dict[str,Any],output: Path) -> Dict[str,Any]:
    parts=[];by_key={};evidence=[]
    case=json.loads((output/protocol["input"]["archive"]).read_text())
    for path in paths:
        value=path.read_bytes();part=json.loads(value);validate_study(part,protocol)
        if parts and any(part[key]!=parts[0][key] for key in ("python_version","mpmath_version","algorithm","root_isolation")):
            raise ValueError("cannot label mixed numerical runtime/algorithm versions as one homogeneous campaign")
        if part["input"]!=case:raise ValueError("study input document differs from archived raw input")
        source=path.parent/part["implementation_archive"]
        if digest(source.read_bytes())!=protocol["solver"]["sha256"]:
            raise ValueError("study implementation archive differs from its declared digest")
        stored="raw/studies/"+digest(value)+".json";atomic_bytes(output/stored,value,True)
        evidence.append({"archive":stored,"sha256":digest(value)})
        for row in part["rows"]:
            key=(row["spatial_helicity"],row.get("requested_k_di",row["k_di"]))
            if key in by_key:raise ValueError("overlapping imported study artifacts")
            by_key[key]=row
        parts.append(part)
    expected={(h,k) for h in protocol["requested_helicities"] for k in protocol["requested_k_di"]}
    if set(by_key)!=expected:raise ValueError("study artifacts do not cover the prepared helicities/points")
    combined=dict(parts[0]);combined["requested_helicities"]=protocol["requested_helicities"]
    combined["implementation_archive"]=protocol["solver"]["archive"]
    combined["rows"]=[by_key[(h,k)] for h in protocol["requested_helicities"] for k in protocol["requested_k_di"]]
    combined["source_artifacts"]=evidence;combined["sampling"]=protocol["sampling"]
    combined["native_grid_map"]="native-grid-map.json";combined["coverage"]={}
    combined["reference_percentage_sampling_brackets"]=[]
    for h in protocol["requested_helicities"]:
        rows=[by_key[(h,k)] for k in protocol["requested_k_di"]]
        combined["coverage"][h]={status:sum(row["status"]==status for row in rows) for status in ("RESOLVED","UNRESOLVED")}
        for field in ("phase_speed_relative_mhd_error","group_speed_relative_mhd_error"):
            for percentage in (1,5,10):
                for low,high in zip(rows,rows[1:]):
                    if low["status"]!="RESOLVED" or high["status"]!="RESOLVED":continue
                    with localcontext() as context:
                        context.prec=protocol["verification_digits"]+25
                        crosses=Decimal(low[field])*100<percentage<=Decimal(high[field])*100
                    if crosses:
                        combined["reference_percentage_sampling_brackets"].append({"spatial_helicity":h,
                            "metric":field,"reference_percentage":percentage,"lower_k_di":low["requested_k_di"],
                            "upper_k_di":high["requested_k_di"],"accepted_cutoff":False})
                        break
    return combined


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",type=Path,required=True)
    parser.add_argument("--output-directory",type=Path,required=True)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--prepare-only",action="store_true")
    mode.add_argument("--run",action="store_true",help="explicit slow root calculation")
    mode.add_argument("--import-studies",type=Path,nargs="+",help="verify/archive existing studies without solving roots")
    parser.add_argument("--sampling",choices=("native-grid","nominal-grid","explicit-kdi"),default="native-grid")
    parser.add_argument("--points-json",type=Path,help="raw JSON decimal-string list; grid modes verify exact rounding")
    parser.add_argument("--helicities",choices=("Right","Left"),nargs="+",default=["Right","Left"])
    parser.add_argument("--digits",type=int,default=40)
    parser.add_argument("--verification-digits",type=int,default=70)
    parser.add_argument("--grid",nargs=3,metavar=("N","FIRST_CENTER","LAST_CENTER"))
    parser.add_argument("--cxx",default=os.environ.get("CXX","c++"))
    parser.add_argument("--plot",action="store_true",help="render verified results with repository plotting source")
    args=parser.parse_args()
    if args.digits<25 or args.verification_digits<args.digits+15:parser.error("need >=25 base digits and >=15 refinement digits")
    if len(args.helicities)!=len(set(args.helicities)):parser.error("duplicate helicities")
    if args.plot and not (args.run or args.import_studies):parser.error("--plot requires --run or --import-studies; preparation calculates no roots")
    output=args.output_directory.resolve();output.mkdir(parents=True,exist_ok=True)
    if (output/"protocol.json").exists():raise ValueError("campaign directory already has a protocol; choose a new output directory")
    raw=args.case.read_bytes();case=json.loads(raw);plasma=warm.Plasma(case)
    with mp.workdps(args.verification_digits+25):plasma.numbers()
    atomic_bytes(output/"raw/input.json",raw,True)
    raw_points=None
    if args.points_json:
        point_bytes=args.points_json.read_bytes();raw_points=json.loads(point_bytes)
        atomic_bytes(output/"raw/requested-points.json",point_bytes,True)
    extracted=capture_grid(output,args.cxx,args.grid or [])
    mapped=coordinates(extracted["grid"],plasma,args.sampling,raw_points,args.verification_digits+15)
    atomic_json(output/"native-grid.json",extracted,True);atomic_json(output/"native-grid-map.json",mapped,True)
    atomic_json(output/"requested-k-di.json",mapped["requested_k_di"],True)
    protocol={"schema":SCHEMA,"status":"PREPARED","input":{"archive":"raw/input.json","sha256":digest(raw)},
              "solver":archive(DIRECTORY/"warm_parallel_maxwellian.py",output),
              "runner":archive(Path(__file__),output),"plotter":archive(DIRECTORY/"plot_warm_parallel_study.py",output),
              "sampling":args.sampling,"requested_k_di":mapped["requested_k_di"],"requested_helicities":args.helicities,
              "base_digits":args.digits,"verification_digits":args.verification_digits,
              "python_version":platform.python_version(),"mpmath_version":mp.__version__,
              "physical_production_qualified":False,"domain_limits_selected":False,"physical_error_acceptance_budget":None,
              "numerical_scope":"Convergence and numerical winding evidence; no interval proof or complete global mode-family certification."}
    command=[sys.executable,str(output/protocol["solver"]["archive"]),"--case",str(output/"raw/input.json"),
             "--output",str(output/"calculated-study.json"),"--digits",str(args.digits),"--verification-digits",str(args.verification_digits),
             "--helicities"]+args.helicities+["--kdi"]+mapped["requested_k_di"]
    protocol["calculation_command"]=command
    protocol["bundle_working_directory"]="."
    protocol["calculation_command_relative"]=["python3",protocol["solver"]["archive"],"--case","raw/input.json",
        "--output","calculated-study.json","--digits",str(args.digits),"--verification-digits",str(args.verification_digits),
        "--helicities"]+args.helicities+["--kdi"]+mapped["requested_k_di"]
    protocol["plot_command_relative"]=["python3",protocol["plotter"]["archive"],"--study","study.json","--output-prefix","study"]
    protocol["uncertainty_scope"]="The solver reports a frequency convergence-error estimate. Group-derivative and Hermitian-response values have no separately certified per-point error bars."
    atomic_json(output/"protocol.json",protocol)
    paths=args.import_studies
    if args.run:
        with (output/"calculation.log").open("wb") as stream:
            completed=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
        if completed.returncode not in (0,1) or not (output/"calculated-study.json").is_file():
            raise RuntimeError("warm solver failed before producing a complete study; inspect calculation.log")
        paths=[output/"calculated-study.json"]
    if not paths:
        print("PREPARED "+str(output/"protocol.json")+"; no roots calculated")
        return 0
    combined=combine_studies(paths,protocol,output);atomic_json(output/"study.json",combined)
    unresolved=sum(count["UNRESOLVED"] for count in combined["coverage"].values())
    protocol["status"]="COMPLETED_WITH_UNRESOLVED" if unresolved else "COMPLETED"
    protocol["coverage"]=combined["coverage"];protocol["study_sha256"]=digest((output/"study.json").read_bytes())
    if args.plot:
        subprocess.run([sys.executable,str(output/protocol["plotter"]["archive"]),"--study",str(output/"study.json"),
                        "--output-prefix",str(output/"study")],check=True)
        plotted=json.loads((output/"study.plot.json").read_bytes())
        if plotted["plotter_sha256"]!=protocol["plotter"]["sha256"] or plotted["study_sha256"]!=protocol["study_sha256"]:
            raise ValueError("rendered figure provenance differs from the pinned plotter/study")
        protocol["figures"]={name:{"path":"study."+name,"sha256":digest((output/("study."+name)).read_bytes())}
                             for name in ("png","pdf","plot.json")}
    atomic_json(output/"protocol.json",protocol)
    print(protocol["status"]+" "+str(output/"protocol.json"))
    return 1 if unresolved else 0


if __name__=="__main__":
    try:raise SystemExit(main())
    except (ValueError,KeyError,OSError,RuntimeError,subprocess.SubprocessError) as error:
        print("warm-study protocol ERROR: "+str(error),file=sys.stderr);raise SystemExit(2)
