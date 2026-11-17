#!/usr/bin/env python3
"""Plot pinned finite-time response evidence, retaining unresolved gaps.

Standalone archive: only standard-library validation and Matplotlib are used;
no mutable sibling imports, orbit calculations or physical cutoff selection.
"""
import argparse
from decimal import Decimal, InvalidOperation, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile
from typing import Any, Dict


def digest(value: bytes) -> str:return hashlib.sha256(value).hexdigest()


def number(value: Any,name: str) -> float:
    if isinstance(value,bool):raise ValueError(name+" is a Boolean, not a numerical result")
    try:exact=Decimal(str(value));result=float(exact)
    except (InvalidOperation,ValueError,OverflowError) as error:raise ValueError("invalid "+name) from error
    if not exact.is_finite() or not math.isfinite(result) or (exact!=0 and result==0):
        raise ValueError(name+" cannot be plotted without nonfinite values or underflow")
    return result


def load_evidence(protocol_path: Path) -> Any:
    root=protocol_path.resolve().parent;raw=protocol_path.read_bytes();protocol=json.loads(raw)
    if protocol.get("schema")!="srcsep-finite-amplitude-protocol-v1" or protocol.get("physical_limits_selected") is not False:
        raise ValueError("unknown/falsely qualified protocol")
    sources={}
    for key,record in protocol["sources"].items():
        path=(root/record["archive"]).resolve()
        if root not in path.parents or digest(path.read_bytes())!=record["sha256"]:raise ValueError("pinned source/input identity changed: "+key)
        sources[key]=path
    controls=json.loads(sources["controls"].read_text());study_path=root/"study.json";study_bytes=study_path.read_bytes();study=json.loads(study_bytes)
    if study.get("schema")!="srcsep-finite-amplitude-study-v1" or study["protocol_sha256"]!=digest(raw) or \
        study["resolved_inputs"]!=protocol["resolved_inputs"] or study.get("physical_production_qualified") is not False:
        raise ValueError("study does not bind the pinned protocol/physical coordinates")
    cases={case["label"]:case for case in protocol["resolved_inputs"]["mode_experiments"]}
    times=dict(zip(controls["observation_gyroperiods"],protocol["resolved_inputs"]["actual_tau_hex"]))
    expected={(a,label,t) for a in controls["amplitudes_delta_b_over_b0"] for label in cases for t in times};seen=set()
    for row in study["rows"]:
        key=(row["amplitude_delta_b_over_b0"],row["mode_experiment"],row["observation_gyroperiods"])
        if key not in expected or key in seen:raise ValueError("study contains duplicate or unrequested coordinates")
        seen.add(key)
        if row.get("status") not in ("RESOLVED","UNRESOLVED"):raise ValueError("unknown response status")
        if row.get("physical_production_qualified") is not False or row["actual_amplitude_hex"]!=number(key[0],"amplitude").hex():
            raise ValueError("row amplitude or qualification metadata differs from pinned request")
        if "actual_tau_hex" in row and row["actual_tau_hex"]!=times[key[2]]:raise ValueError("row time differs from pinned coordinate")
        finest=cases[key[1]]["levels"][-1]
        if "actual_k_per_m_hex" in row and (row["actual_k_per_m_hex"]!=finest["k_per_m_hex"] or row["actual_probability_weights_hex"]!=finest["probability_weight_hex"]):
            raise ValueError("row k/weights differ from pinned wave coordinates")
        if row["status"]=="RESOLVED":
            for value,name in ((row["finite_time_full_lorentz"]["pitch_change_variance"],"full response"),(row["linear_phase_variance"],"linear response")):
                if number(value,name)<0:raise ValueError("negative variance")
    if seen!=expected:raise ValueError("study has missing requested response coordinates")
    if study["resolved_rows"]!=sum(r["status"]=="RESOLVED" for r in study["rows"]) or study["unresolved_rows"]!=sum(r["status"]=="UNRESOLVED" for r in study["rows"]):
        raise ValueError("study coverage counts differ from actual row statuses")
    return protocol,controls,study,study_path,raw,study_bytes


def atomic_json(path: Path,value: Dict[str,Any]) -> None:
    descriptor,temporary=tempfile.mkstemp(prefix=path.name+".",dir=str(path.parent))
    try:
        with os.fdopen(descriptor,"w") as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write("\n");stream.flush();os.fsync(stream.fileno())
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def render(protocol_path: Path,prefix: Path) -> Dict[str,Any]:
    source_bytes=Path(__file__).read_bytes()
    protocol,controls,study,study_path,protocol_bytes,study_bytes=load_evidence(protocol_path)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if protocol["plot_dependency_version"]!=matplotlib.__version__:raise ValueError("Matplotlib version differs from pinned plotting dependency")
    figure,axes=plt.subplots(1,3,figsize=(16,4.8),constrained_layout=True)
    index={(r["amplitude_delta_b_over_b0"],r["mode_experiment"],r["observation_gyroperiods"]):r for r in study["rows"]}
    unavailable=0;zero_amplitudes=0;zero_or_undefined_values=0;unidentified=0;marker_label_used=False
    for case in protocol["resolved_inputs"]["mode_experiments"]:
        for time in controls["observation_gyroperiods"]:
            amp=[];full=[];linear=[];difference=[];uncertainty=[];resolution_markers=[]
            for value in controls["amplitudes_delta_b_over_b0"]:
                a=number(value,"amplitude")
                if a==0:zero_amplitudes+=1;continue  # Exact zero has no log coordinate.
                row=index[(value,case["label"],time)];amp.append(a)
                if row["status"]!="RESOLVED":
                    full.append(math.nan);linear.append(math.nan);difference.append(math.nan);uncertainty.append(math.nan);resolution_markers.append(math.nan);unavailable+=1;continue
                vf=number(row["finite_time_full_lorentz"]["pitch_change_variance"],"full variance")
                vl=number(row["linear_phase_variance"],"linear variance")
                full.append(vf if vf>0 else math.nan);linear.append(vl if vl>0 else math.nan)
                relative=row["relative_variance_difference"]
                difference.append(abs(number(relative,"relative difference")) if row["finite_amplitude_difference_resolved"] and relative is not None and relative!=0 else math.nan)
                if vl>0:
                    with localcontext() as context:
                        context.prec=80
                        ratio=Decimal(str(row["numerical_uncertainty_estimates"]["sum"]))/Decimal(str(vl))
                    uncertainty.append(number(ratio,"relative uncertainty") if ratio!=0 else math.nan)
                else:uncertainty.append(math.nan)
                # A computed difference below its numerical estimate is not
                # an identified physical discrepancy. Show its resolution
                # estimate as a separate marker, not as a nonlinear-error point
                # or a certified upper bound. No orbit data are recalculated.
                if not row["finite_amplitude_difference_resolved"] and math.isfinite(uncertainty[-1]):
                    resolution_markers.append(uncertainty[-1]);unidentified+=1
                else:resolution_markers.append(math.nan)
                zero_or_undefined_values+=int(vf==0)+int(vl==0)+int(relative is None or relative==0)+int(not math.isfinite(uncertainty[-1]))
            label=case["label"]+", T="+time+" gyroperiods"
            line=axes[0].plot(amp,full,"o-",label=label)[0];axes[0].plot(amp,linear,"--",color=line.get_color(),alpha=.55)
            axes[1].plot(amp,difference,"o-");axes[2].plot(amp,uncertainty,"o-",label=label)
            if any(math.isfinite(value) for value in resolution_markers):
                axes[1].plot(amp,resolution_markers,"v",color="grey",markerfacecolor="none",
                    label="unidentified difference: numerical estimate" if not marker_label_used else None)
                marker_label_used=True
    # Every nonpositive/undefined datum was explicitly replaced by NaN above,
    # before Matplotlib sees it. No library-version-specific clipping keyword
    # can invent a finite log point or silently change the evidence.
    for ax in axes:ax.set_xscale("log");ax.set_yscale("log");ax.set_xlabel("declared deltaB/B0");ax.grid(True,which="both",alpha=.2)
    requested_positive=[number(value,"amplitude") for value in controls["amplitudes_delta_b_over_b0"] if number(value,"amplitude")>0]
    if len(requested_positive)>1:
        for ax in axes:ax.set_xlim(min(requested_positive),max(requested_positive))
    for value in controls["amplitudes_delta_b_over_b0"]:
        selected=[r for r in study["rows"] if r["amplitude_delta_b_over_b0"]==value]
        if number(value,"amplitude")>0 and all(r["status"]=="UNRESOLVED" for r in selected):
            for ax in axes:
                ax.axvline(number(value,"amplitude"),color="grey",linestyle=":",alpha=.3)
                ax.text(number(value,"amplitude"),.02,"unresolved",rotation=90,ha="right",va="bottom",transform=ax.get_xaxis_transform(),fontsize=6,color="grey")
    axes[0].set_ylabel("finite-time phase variance of delta mu");axes[0].set_title("Lorentz solid; linear dashed")
    axes[1].set_ylabel("absolute relative variance difference");axes[1].set_title("difference above numerical estimate; grey shows resolution")
    axes[2].set_ylabel("numerical uncertainty estimate / linear variance");axes[2].set_title("separate numerical resolution")
    for level in (.01,.05,.1):axes[1].axhline(level,color="grey",linestyle=":",alpha=.4)
    handles,labels=axes[0].get_legend_handles_labels()
    if handles:axes[0].legend(handles,labels,fontsize=6)
    if marker_label_used:axes[1].legend(fontsize=6)
    baseline=json.loads((protocol_path.parent/protocol["sources"]["equilibrium"]["archive"]).read_text())
    actual_u=float.fromhex(protocol["resolved_inputs"]["actual_probe"]["wave_speed_mps"])
    context=(controls["probe"]["species"]+" K="+controls["probe"]["kinetic_energy_ev"]+" eV; mu0="+
        controls["probe"]["mean_pitch"]+"; u="+format(actual_u,".7g")+" m/s; B0="+baseline["magnetic_field_t"]+
        " T; n_p="+baseline["proton_density_m3"]+" /m3; Tp/Te="+baseline["proton_temperature_k"]+"/"+baseline["electron_temperature_k"]+" K")
    figure.suptitle("Prescribed single-wave finite-time response; no diffusivity or selected physical cutoff\n"+context+"\n"+
        str(study["resolved_rows"])+" resolved, "+str(study["unresolved_rows"])+" unresolved; gaps retained; 1/5/10% are reporting references",fontsize=9)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    outputs={}
    for suffix in ("png","pdf"):
        target=Path(str(prefix)+"."+suffix);descriptor,temporary=tempfile.mkstemp(prefix=target.name+".",suffix="."+suffix,dir=str(target.parent));os.close(descriptor)
        try:
            metadata={"CreationDate":None,"ModDate":None} if suffix=="pdf" else {"Software":"srcSEP pinned finite-time study plotter"}
            figure.savefig(temporary,format=suffix,dpi=180,metadata=metadata);os.replace(temporary,target)
        finally:
            if os.path.exists(temporary):os.unlink(temporary)
        outputs[suffix]={"path":str(target),"sha256":digest(target.read_bytes())}
    plt.close(figure)
    if protocol_path.read_bytes()!=protocol_bytes or study_path.read_bytes()!=study_bytes or Path(__file__).read_bytes()!=source_bytes:
        raise ValueError("plot source/protocol/study changed during rendering; output verdict discarded")
    # A replot may use a newer explicitly invoked plotter without relabeling
    # the old calculation as new. Archive the ACTUAL rendering implementation
    # and record both identities; default protocol commands still reproduce
    # the original archived plotter, while this derived artifact is auditable.
    archive_directory=Path(str(prefix)+".sources");archive_directory.mkdir(exist_ok=True)
    archive_path=archive_directory/(digest(source_bytes)+"-plot_finite_amplitude_study.py")
    if archive_path.exists() and archive_path.read_bytes()!=source_bytes:raise ValueError("immutable plot-source archive differs")
    if not archive_path.exists():
        descriptor,temporary=tempfile.mkstemp(prefix=archive_path.name+".",dir=str(archive_directory))
        try:
            with os.fdopen(descriptor,"wb") as stream:stream.write(source_bytes);stream.flush();os.fsync(stream.fileno())
            os.replace(temporary,archive_path)
        finally:
            if os.path.exists(temporary):os.unlink(temporary)
    result={"schema":"srcsep-finite-amplitude-plot-v1","study_sha256":digest(study_bytes),
        "protocol_sha256":digest(protocol_bytes),"plotter_sha256":digest(source_bytes),
        "original_protocol_plotter_sha256":protocol["sources"]["plotter"]["sha256"],
        "uses_original_protocol_plotter":digest(source_bytes)==protocol["sources"]["plotter"]["sha256"],
        "actual_plotter_archive":os.path.relpath(archive_path,prefix.parent),
        "matplotlib":matplotlib.__version__,"resolved":study["resolved_rows"],"unresolved":study["unresolved_rows"],
        "unresolved_log_points_omitted":unavailable,"exact_zero_amplitudes_omitted_from_log":zero_amplitudes,
        "zero_or_undefined_log_values_omitted":zero_or_undefined_values,
        "unidentified_difference_resolution_markers":unidentified,
        "equilibrium_and_probe_caption":context,"wave_speed_source":controls["alfven_speed_source"],
        "physical_limits_selected":False,"outputs":outputs}
    atomic_json(Path(str(prefix)+".plot.json"),result);return result


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--protocol",type=Path,required=True);parser.add_argument("--output-prefix",type=Path,required=True)
    options=parser.parse_args()
    try:render(options.protocol.resolve(),options.output_prefix.resolve());return 0
    except (ValueError,KeyError,OSError,TypeError,RuntimeError,ImportError) as error:print("ERROR "+str(error),file=sys.stderr);return 2


if __name__=="__main__":raise SystemExit(main())
