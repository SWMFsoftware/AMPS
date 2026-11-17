#!/usr/bin/env python3
"""Plot archived controlled-study evidence, preserving unresolved gaps.

All plasma labels come from the input document. Percentages are reference
lines, never selected physical cutoffs. A real derivative of a complex mode
is a dispersion diagnostic; for strong damping it need not be an energy
transport speed. The Hermitian response is likewise weak-damping-only.
"""
import argparse
from decimal import Decimal,InvalidOperation,localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import sys
from typing import Any,Dict

sys.dont_write_bytecode=True

def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def finite_decimal(value: Any,name: str) -> Decimal:
    if not isinstance(value,str):raise ValueError(name+" must retain its decimal-string representation")
    try:result=Decimal(value)
    except InvalidOperation as error:raise ValueError(name+" is not a decimal number") from error
    if not result.is_finite():raise ValueError(name+" is nonfinite")
    return result


def verify_coordinate_record(row: Dict[str,Any],requested: str,verification_digits: int) -> None:
    """Bind a row to its raw request, allowing only declared formatting error.

    The solver parses k at verification_digits+15 decimal working digits,
    then serializes at verification_digits+5 significant digits. The bound
    below accounts for one-half output decimal ulp plus a conservative input
    conversion ulp. It is a numeric representation bound, not a physical
    wavenumber comparison tolerance or permission to move a study point.
    """
    wanted=finite_decimal(requested,"requested k_di")
    actual=finite_decimal(row["k_di"],"reported k_di")
    if wanted<=0 or actual<=0:raise ValueError("study coordinates must be positive")
    if row["status"]=="UNRESOLVED":
        if row["k_di"]!=requested:raise ValueError("unresolved coordinate differs from its exact raw request")
        return
    if not isinstance(verification_digits,int) or verification_digits<25:
        raise ValueError("invalid coordinate serialization precision")
    digits=verification_digits+5
    if row.get("serialized_significant_digits")!=digits:
        raise ValueError("row coordinate serialization precision differs from the solver protocol")
    with localcontext() as context:
        context.prec=max(digits+30,len(wanted.as_tuple().digits)+10,len(actual.as_tuple().digits)+10)
        decimal_ulp=Decimal(10)**(wanted.adjusted()-digits+1)
        input_bound=Decimal(10)**(wanted.adjusted()-(verification_digits+15)+1)
        if abs(actual-wanted)>decimal_ulp/2+input_bound:
            raise ValueError("reported coordinate differs from requested point beyond serialization rounding")


def plot_number(value: Any,name: str) -> float:
    decimal=finite_decimal(value,name);result=float(decimal)
    if not math.isfinite(result) or (decimal!=0 and result==0):
        raise ValueError(name+" is outside the plotting float range; no replacement or floor supplied")
    return result


def atomic_json(path: Path,value: Any) -> None:
    # Keep the content-addressed plotter standalone: replaying the archived
    # source needs only its recorded Python/Matplotlib environment, not imports
    # from a mutable repository checkout or unhashed sibling modules.
    path.parent.mkdir(parents=True,exist_ok=True)
    descriptor,temporary=tempfile.mkstemp(prefix=path.name+".",dir=str(path.parent))
    try:
        with os.fdopen(descriptor,"w") as stream:
            json.dump(value,stream,indent=2,allow_nan=False);stream.write("\n");stream.flush();os.fsync(stream.fileno())
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def equilibrium_caption(document: Dict[str,Any]) -> str:
    def quantity(name: str,scale: str="1") -> str:
        value=finite_decimal(document[name],name)/Decimal(scale)
        if value<=0:raise ValueError("controlled plasma label requires positive "+name)
        return format(value,".5g")
    return ("n_p="+quantity("proton_density_m3","1000000")+", n_e="+quantity("electron_density_m3","1000000")
            +" cm⁻³; B="+quantity("magnetic_field_t","1e-9")+" nT; Tp="+quantity("proton_temperature_k")
            +" K; Te="+quantity("electron_temperature_k")+" K")


def load_study(path: Path) -> Dict[str,Any]:
    result=json.loads(path.read_bytes())
    if result.get("schema")!="srcsep-controlled-warm-parallel-study-v1":raise ValueError("unknown controlled-study schema")
    if result.get("physical_production_qualified") is not False or result.get("domain_limits_selected") is not False:
        raise ValueError("this figure protocol is for explicitly unqualified, no-cutoff controlled studies")
    labels=result.get("requested_helicities",[]);points=result.get("requested_k_di",[])
    if not labels or len(set(labels))!=len(labels) or any(h not in ("Right","Left") for h in labels):
        raise ValueError("invalid/duplicate spatial-helicity labels")
    if not points or len(set(points))!=len(points):raise ValueError("missing/duplicate requested coordinates")
    expected={(h,k) for h in labels for k in points};seen=set()
    for row in result["rows"]:
        key=(row["spatial_helicity"],row.get("requested_k_di",row["k_di"]))
        if key not in expected or key in seen:raise ValueError("unrequested/duplicate plot row")
        seen.add(key)
        verify_coordinate_record(row,key[1],result["verification_digits"])
        if finite_decimal(row["k_di"],"k_di")<=0:raise ValueError("plot wavenumber must be positive")
        if row["status"]=="UNRESOLVED":
            if not row.get("reason"):raise ValueError("unresolved plot row lacks its reason")
        elif row["status"]=="RESOLVED":
            for field in ("phase_speed_relative_mhd_error","group_speed_relative_mhd_error",
                          "damping_over_real_frequency_raw","weak_damping_hermitian_energy_over_magnetic_energy"):
                value=finite_decimal(row[field],field)
                if field!="weak_damping_hermitian_energy_over_magnetic_energy" and value<0:
                    raise ValueError("negative absolute diagnostic: "+field)
            if not isinstance(row.get("damping_resolved"),bool):raise ValueError("missing damping availability flag")
        else:raise ValueError("unknown plot row status")
    if seen!=expected:raise ValueError("plot source lacks requested points; gaps require explicit UNRESOLVED records")
    if result["input"].get("distribution")!="isotropic-nonrelativistic-Maxwellian" or result["input"].get("propagation")!="parallel":
        raise ValueError("figure labels require the explicitly declared parallel Maxwellian model")
    equilibrium_caption(result["input"])
    return result


def atomic_figure(figure: Any,path: Path,format_name: str,metadata: Dict[str,Any]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    descriptor,temporary=tempfile.mkstemp(prefix=path.name+".",suffix="."+format_name,dir=str(path.parent))
    os.close(descriptor)
    try:
        figure.savefig(temporary,format=format_name,dpi=180,metadata=metadata)
        with open(temporary,"rb") as stream:os.fsync(stream.fileno())
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def render(source: Path,prefix: Path) -> Dict[str,Any]:
    # Matplotlib is needed only by an explicit plot request, never preparation.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data=load_study(source);figure,axes=plt.subplots(2,2,figsize=(11,8))
    coverage={};zero_errors=0;unresolved=0;any_damping=False
    for helicity,color in (("Right","#2166ac"),("Left","#b2182b")):
        rows=sorted((row for row in data["rows"] if row["spatial_helicity"]==helicity),
                    key=lambda row:Decimal(row["k_di"]))
        if not rows:continue
        resolved=sum(row["status"]=="RESOLVED" for row in rows)
        missing=len(rows)-resolved;coverage[helicity]={"RESOLVED":resolved,"UNRESOLVED":missing};unresolved+=missing
        label=helicity+" ("+str(resolved)+"/"+str(len(rows))+" resolved)"
        x=[plot_number(row["k_di"],"k_di") for row in rows]
        if any(not math.isfinite(value) or value<=0 for value in x):raise ValueError("coordinate is outside plotting float range")
        # NaN here is a plotting mask, never an inserted physical datum. It
        # breaks the line at an explicitly unavailable row instead of joining
        # samples across an unresolved interval or adding a plotting floor.
        gap=float("nan")
        for axis,field in ((axes[0,0],"phase_speed_relative_mhd_error"),(axes[0,1],"group_speed_relative_mhd_error")):
            values=[]
            for row in rows:
                if row["status"]!="RESOLVED":values.append(gap);continue
                decimal=finite_decimal(row[field],field)
                if decimal==0:zero_errors+=1;values.append(gap)
                else:values.append(plot_number(str(decimal*100),field+" percentage"))
            axis.loglog(x,values,"o-",markersize=3,color=color,label=label)
        damping=[plot_number(row["damping_over_real_frequency_raw"],"resolved damping ratio") if row["status"]=="RESOLVED" and
                 row["damping_resolved"] and Decimal(row["damping_over_real_frequency_raw"])>0 else gap for row in rows]
        energy=[plot_number(row["weak_damping_hermitian_energy_over_magnetic_energy"],"Hermitian response") if row["status"]=="RESOLVED" else gap for row in rows]
        if any(not math.isfinite(value) for value in damping+energy if not math.isnan(value)):
            raise ValueError("damping/energy diagnostic is outside plotting float range")
        if any(not math.isnan(value) for value in damping):
            axes[1,0].loglog(x,damping,"o-",markersize=3,color=color,label=label);any_damping=True
        else:axes[1,0].plot([],[],"o-",markersize=3,color=color,label=label)
        axes[1,1].semilogx(x,energy,"o-",markersize=3,color=color,label=label)
    for axis in axes.flat:
        axis.set_xlabel(r"$k d_p$ (ion inertial-length normalization)")
        axis.grid(True,which="both",alpha=.2);axis.legend(loc="best",fontsize=8)
    for axis in axes[0,:]:
        axis.set_ylabel(r"Absolute relative error against $v_{A,p}$ [%]")
        for percentage in (1,5,10):
            axis.axhline(percentage,color=".5",linestyle=":",linewidth=.8)
            axis.text(.03,percentage*1.04,str(percentage)+"% reference",fontsize=8,color=".4",transform=axis.get_yaxis_transform())
    axes[0,0].set_title("Real phase-speed error");axes[0,1].set_title("Real group-derivative error")
    axes[1,0].set_title("Resolved damping\nValues below numerical error omitted")
    if not any_damping:
        axes[1,0].set_xscale("log")
        axes[1,0].text(.5,.5,"No damping values numerically resolved",ha="center",transform=axes[1,0].transAxes)
    axes[1,0].set_ylabel(r"$|\mathrm{Im}\,\omega|/|\mathrm{Re}\,\omega|$")
    axes[1,1].set_title("Hermitian energy response\nWeak-damping diagnostic only")
    axes[1,1].set_ylabel(r"Diagnostic $W/E_B$");axes[1,1].axhline(2,color=".5",linestyle=":",linewidth=.8)
    figure.suptitle("Controlled parallel proton-electron Maxwellian equilibrium\n"+equilibrium_caption(data["input"]),fontsize=12)
    figure.text(.5,.018,"Reference percentages are not accepted cutoffs. No production qualification or automatic spectral truncation.\n"
                "Strong damping prevents a conserved-energy interpretation of the Hermitian response or an energy-speed interpretation of the real group derivative.\n"
                +str(unresolved)+"/"+str(len(data["rows"]))+" roots unresolved; curves retain gaps. Full inputs, numerical estimates and provenance remain in JSON.",
                ha="center",fontsize=7.5)
    figure.tight_layout(rect=(0,.065,1,.93))
    plot_hash=digest(Path(__file__).read_bytes())
    png=Path(str(prefix)+".png");pdf=Path(str(prefix)+".pdf")
    atomic_figure(figure,png,"png",{"Software":"srcSEP controlled-study plot "+plot_hash})
    atomic_figure(figure,pdf,"pdf",{"Creator":"srcSEP controlled-study plot "+plot_hash,"CreationDate":None,"ModDate":None})
    plt.close(figure)
    metadata={"schema":"srcsep-controlled-warm-plot-v1","study_sha256":digest(source.read_bytes()),
              "plotter_sha256":plot_hash,"matplotlib_version":matplotlib.__version__,"coverage":coverage,
              "zero_error_samples_omitted_from_log_panels":zero_errors,"reference_percentages":[1,5,10],
              "physical_error_acceptance_budget":None,"physical_production_qualified":False,
              "curve_gaps":"Explicit UNRESOLVED rows remain gaps; no interpolation or damping floor is supplied.",
              "plasma_inputs":data["input"],"figures":{str(path.name):digest(path.read_bytes()) for path in (png,pdf)}}
    atomic_json(Path(str(prefix)+".plot.json"),metadata)
    return metadata


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",type=Path,required=True);parser.add_argument("--output-prefix",type=Path,required=True)
    args=parser.parse_args();result=render(args.study,args.output_prefix)
    print(json.dumps({"coverage":result["coverage"],"output_prefix":str(args.output_prefix)}));return 0


if __name__=="__main__":
    try:raise SystemExit(main())
    except (ValueError,KeyError,OSError) as error:
        print("warm-study plot ERROR: "+str(error),file=__import__("sys").stderr);raise SystemExit(2)
