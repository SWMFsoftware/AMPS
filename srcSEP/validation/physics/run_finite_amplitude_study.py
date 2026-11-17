#!/usr/bin/env python3
"""Prepare or explicitly execute a source-pinned finite-time Lorentz study.

Default preparation computes no orbits. --run executes the archived runner,
core and equilibrium helper; raw JSON/decimal requests and actual binary64
coordinates remain distinct provenance. This study returns finite-time
phase-ensemble responses, never a diffusivity, domain cutoff or production
qualification. Unresolved numerical rows remain in the artifact and exit1.
"""
import argparse
import hashlib
import importlib.util
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from typing import Any, Dict, List

import mpmath as mp
import numpy as np
import scipy

sys.dont_write_bytecode=True
DIRECTORY=Path(__file__).resolve().parent
SCHEMA="srcsep-finite-amplitude-protocol-v1"


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def load_module(name: str,path: Path) -> Any:
    specification=importlib.util.spec_from_file_location(name,path)
    if specification is None or specification.loader is None:
        raise ValueError("cannot load pinned study dependency: "+str(path))
    module=importlib.util.module_from_spec(specification)
    sys.modules[name]=module  # Required for dataclass annotations in archives.
    specification.loader.exec_module(module)
    return module


def atomic_json(path: Path,value: Any) -> None:
    data=(json.dumps(value,indent=2,allow_nan=False)+"\n").encode()
    descriptor,temporary=tempfile.mkstemp(prefix=path.name+".",dir=str(path.parent))
    try:
        with os.fdopen(descriptor,"wb") as stream:
            stream.write(data);stream.flush();os.fsync(stream.fileno())
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def validate_control(control: Dict[str,Any],core: Any) -> None:
    required={"schema","provenance","alfven_speed_source","probe","amplitudes_delta_b_over_b0",
              "observation_gyroperiods","ensemble","orbit_steps_per_gyroperiod","phase_counts","numerical_controls"}
    if set(control)!=required or control["schema"]!="srcsep-finite-amplitude-controls-v1" or not control["provenance"]:
        raise ValueError("incomplete/unknown finite-amplitude study controls")
    if set(control["probe"])!={"species","kinetic_energy_ev","mean_pitch","branch","helicity"}:
        raise ValueError("probe inputs must be explicit and recognized")
    amplitudes=[core.explicit_number(x,"amplitude") for x in control["amplitudes_delta_b_over_b0"]]
    if not amplitudes or amplitudes[0]<0 or any(a>=b for a,b in zip(amplitudes,amplitudes[1:])):
        raise ValueError("nonnegative strictly increasing amplitude study coordinates required")
    for name in ("orbit_steps_per_gyroperiod","phase_counts"):
        values=control[name]
        if len(values)!=3 or any(type(x) is not int or x<4 for x in values) or values[1]!=2*values[0] or values[2]!=2*values[1]:
            raise ValueError(name+" requires three explicit factor-two numerical refinements")
    ensemble=control["ensemble"]
    if ensemble.get("kind")=="single-mode-probes":
        if set(ensemble)!={"kind","resonant_k_ratios"}:raise ValueError("unknown single-mode inputs")
        ratios=[core.explicit_number(x,"k ratio") for x in ensemble["resonant_k_ratios"]]
        if not ratios or ratios[0]<=0 or any(a>=b for a,b in zip(ratios,ratios[1:])):raise ValueError("positive ordered single-mode k ratios required")
    elif ensemble.get("kind")=="uniform-k-mixture":
        if set(ensemble)!={"kind","lower_resonant_k_ratio","upper_resonant_k_ratio","mode_counts"}:raise ValueError("unknown mode-mixture inputs")
        low=core.explicit_number(ensemble["lower_resonant_k_ratio"],"lower k ratio")
        high=core.explicit_number(ensemble["upper_resonant_k_ratio"],"upper k ratio")
        values=ensemble["mode_counts"]
        if not 0<low<high or len(values)!=3 or any(type(x) is not int or x<2 for x in values) or values[1]!=2*values[0] or values[2]!=2*values[1]:
            raise ValueError("finite positive k support and three explicit mode refinements required")
    else:raise ValueError("field ensemble must explicitly select separate single-wave probes or uniform-k mixture")
    numerical=control["numerical_controls"]
    names={"relative_response_target","absolute_pitch_target","sampled_invariant_over_initial_kinetic_target",
           "dop853_rtol","dop853_atol","dop853_phase_samples","linear_reference_digits",
           "maximum_steps_per_trajectory","maximum_total_trajectory_updates"}
    if set(numerical)!=names:raise ValueError("explicit complete numerical-control table required")
    for name in ("relative_response_target","absolute_pitch_target","sampled_invariant_over_initial_kinetic_target","dop853_rtol","dop853_atol"):
        if core.explicit_number(numerical[name],name)<=0:raise ValueError("numerical targets must be positive")
    for name in ("dop853_phase_samples","maximum_steps_per_trajectory","maximum_total_trajectory_updates"):
        if type(numerical[name]) is not int or numerical[name]<=0:raise ValueError("positive integer numerical work controls required")
    digits=numerical["linear_reference_digits"]
    if len(digits)!=2 or any(type(x) is not int for x in digits) or digits[0]<25 or digits[1]<digits[0]+15:
        raise ValueError("analytic-reference precision refinement requires>=25 digits and>=15 added digits")


def mode_sets(control: Dict[str,Any],probe: Any,core: Any) -> List[Dict[str,Any]]:
    scale=core.resonance_scale_per_m(probe);ensemble=control["ensemble"]
    if ensemble["kind"]=="single-mode-probes":
        cases=[{"label":"single-mode-k-ratio="+ratio,"request":ratio,
                 "levels":[(np.array([float(ratio)*scale]),np.array([1.0]))]*3,"mixture":None}
                for ratio in ensemble["resonant_k_ratios"]]
        if any(not np.all(np.isfinite(k)) or np.any(k<=0) for case in cases for k,_ in case["levels"]):raise core.Unresolved("requested single-mode k is not representable")
        return cases
    low=float(ensemble["lower_resonant_k_ratio"])*scale;high=float(ensemble["upper_resonant_k_ratio"])*scale
    levels=[]
    for count in ensemble["mode_counts"]:
        nodes,weights=np.polynomial.legendre.leggauss(count)
        levels.append(((low+high)/2+(high-low)*nodes/2,weights/2))
    if not 0<low<high or any(not np.all(np.isfinite(k)) for k,_ in levels):raise core.Unresolved("finite-k mixture support is not representable")
    return [{"label":"uniform-k-mixture","request":dict(ensemble),"levels":levels,"mixture":(low,high)}]


def resolve_inputs(baseline: Dict[str,Any],control: Dict[str,Any],core: Any,helper: Any) -> Dict[str,Any]:
    validate_control(control,core)
    with mp.workdps(80):
        plasma=helper.Plasma(baseline);numbers=plasma.numbers();scales=plasma.scales()
        probe=core.resolve_probe(numbers,scales,control)
        times=core.observation_times(probe,control["observation_gyroperiods"])
        for cycles in control["observation_gyroperiods"]:
            for steps in control["orbit_steps_per_gyroperiod"]:
                count=float(cycles)*steps
                if count!=round(count):raise ValueError("raw observation coordinates do not align with declared numerical step grids")
        records=[]
        for case in mode_sets(control,probe,core):
            records.append({"label":case["label"],"request":case["request"],
                "levels":[{"k_per_m_hex":[float(x).hex() for x in k],"probability_weight_hex":[float(x).hex() for x in w]}
                          for k,w in case["levels"]],
                "mixture_support_per_m_hex":None if case["mixture"] is None else [x.hex() for x in case["mixture"]]})
        actual={name:float(value).hex() for name,value in vars(probe).items() if name not in ("momentum","helicity")}
        actual["momentum_over_mc_hex"]=[x.hex() for x in probe.momentum];actual["helicity"]=probe.helicity
        return {"actual_probe":actual,"actual_tau_hex":[float(x).hex() for x in times],"mode_experiments":records,
                "actual_amplitudes_hex":[float(x).hex() for x in control["amplitudes_delta_b_over_b0"]],
                "equilibrium_scales_decimal":{key:mp.nstr(value,75) for key,value in scales.items()},
                "frame":"mean plasma rest; pitch relative to mean B0, not instantaneous total B",
                "mode_mixture_scope":"one circular wave per trajectory, not a coherent multimode field",
                "physical_production_qualified":False,"physical_error_acceptance_budget":None}


def dependencies() -> Dict[str,str]:
    return {"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"mpmath":mp.__version__}


def prepare(equilibrium: Path,control_path: Path,output: Path) -> Path:
    if output.exists() and any(output.iterdir()):raise ValueError("prepare requires a fresh empty output directory")
    output.mkdir(parents=True,exist_ok=True);(output/"sources").mkdir()
    files={"runner":Path(__file__).resolve(),"core":DIRECTORY/"finite_amplitude_response.py","plotter":DIRECTORY/"plot_finite_amplitude_study.py",
           "equilibrium_helper":DIRECTORY/"warm_parallel_maxwellian.py","equilibrium":equilibrium.resolve(),"controls":control_path.resolve()}
    archived={}
    for name,path in files.items():
        value=path.read_bytes();relative="sources/"+digest(value)+"-"+path.name;(output/relative).write_bytes(value)
        archived[name]={"archive":relative,"sha256":digest(value)}
    core=load_module("finite_amplitude_prepare_core",files["core"])
    helper=load_module("finite_amplitude_prepare_equilibrium",files["equilibrium_helper"])
    baseline=json.loads(files["equilibrium"].read_text());control=json.loads(files["controls"].read_text())
    try:plot_version=importlib.metadata.version("matplotlib")
    except importlib.metadata.PackageNotFoundError:plot_version=None
    protocol={"schema":SCHEMA,"sources":archived,"dependencies":dependencies(),"platform":platform.platform(),"plot_dependency_version":plot_version,
        "resolved_inputs":resolve_inputs(baseline,control,core,helper),
        "calculation_command_relative":[sys.executable,archived["runner"]["archive"],"--execute-protocol","protocol.json"],
        "plot_command_relative":[sys.executable,archived["plotter"]["archive"],"--protocol","protocol.json","--output-prefix","finite-amplitude-response"],
        "physical_production_qualified":False,"physical_limits_selected":False,
        "interpretation":"finite-time prescribed-circular-wave response only; no diffusivity or physical acceptance cutoff"}
    atomic_json(output/"protocol.json",protocol)
    # Do not accept a live source edit between initial capture and preparation.
    for name,path in files.items():
        if digest(path.read_bytes())!=archived[name]["sha256"]:raise ValueError("source/input changed during preparation: "+str(path))
    return output/"protocol.json"


def read_protocol(path: Path) -> Any:
    protocol=json.loads(path.read_text());root=path.resolve().parent
    if protocol.get("schema")!=SCHEMA or protocol.get("physical_production_qualified") is not False or protocol.get("physical_limits_selected") is not False:
        raise ValueError("unknown or falsely qualified finite-amplitude protocol")
    paths={}
    for name,record in protocol["sources"].items():
        target=(root/record["archive"]).resolve()
        if root not in target.parents or digest(target.read_bytes())!=record["sha256"]:
            raise ValueError("pinned study source/input identity mismatch: "+name)
        paths[name]=target
    if protocol["dependencies"]!=dependencies():raise ValueError("replay dependency versions differ from pinned numerical protocol")
    core=load_module("finite_amplitude_execution_core",paths["core"])
    helper=load_module("finite_amplitude_execution_equilibrium",paths["equilibrium_helper"])
    baseline=json.loads(paths["equilibrium"].read_text());control=json.loads(paths["controls"].read_text())
    if resolve_inputs(baseline,control,core,helper)!=protocol["resolved_inputs"]:
        raise ValueError("resolved physical/numerical coordinates differ from pinned raw inputs")
    return protocol,paths,baseline,control,core,helper


def observed_order(coarse: float,middle: float,fine: float) -> Any:
    first=abs(coarse-middle);second=abs(middle-fine)
    return math.log(first/second,2) if first>0 and second>0 else None


def execute(path: Path) -> int:
    protocol,paths,baseline,control,core,helper=read_protocol(path)
    if digest(Path(__file__).read_bytes())!=protocol["sources"]["runner"]["sha256"]:
        raise ValueError("execute the actual archived runner, not a different live implementation")
    with mp.workdps(80):
        plasma=helper.Plasma(baseline);scales=plasma.scales();probe=core.resolve_probe(plasma.numbers(),scales,control)
    numerical=control["numerical_controls"];rows=[];updates=0
    for case in mode_sets(control,probe,core):
        for amplitude_text in control["amplitudes_delta_b_over_b0"]:
            amplitude=float(amplitude_text);fine_k,fine_weights=case["levels"][-1]
            common={"mode_experiment":case["label"],"amplitude_delta_b_over_b0":amplitude_text,"actual_amplitude_hex":amplitude.hex()}
            case_rows=[]
            try:
                def calculate(k: np.ndarray,w: np.ndarray,phases: int,steps: int) -> Any:
                    nonlocal updates
                    requested=int(float(control["observation_gyroperiods"][-1])*steps)*len(k)*phases
                    if updates+requested>numerical["maximum_total_trajectory_updates"]:
                        raise core.Unresolved("declared total Lorentz trajectory-update resource bound exceeded")
                    orbit=core.integrate_rk4(probe,k,amplitude,phases,steps,control["observation_gyroperiods"],numerical["maximum_steps_per_trajectory"])
                    updates+=orbit["trajectory_updates"]
                    return orbit,core.moments(probe,orbit,k,w,amplitude)
                orbit_levels=[calculate(fine_k,fine_weights,control["phase_counts"][-1],steps)
                              for steps in control["orbit_steps_per_gyroperiod"]]
                finest_orbit,finest=orbit_levels[-1]
                phase_levels=[calculate(fine_k,fine_weights,count,control["orbit_steps_per_gyroperiod"][-1])[1]
                              for count in control["phase_counts"][:-1]]+[finest]
                if case["mixture"] is None:mode_levels=[finest]*3
                else:
                    mode_levels=[calculate(k,w,control["phase_counts"][-1],control["orbit_steps_per_gyroperiod"][-1])[1]
                                 for k,w in case["levels"][:-1]]+[finest]
                cross=core.dop853_comparison(probe,fine_k,amplitude,finest_orbit,float(numerical["dop853_rtol"]),
                                            float(numerical["dop853_atol"]),numerical["dop853_phase_samples"])
                for sample,(cycles,tau) in enumerate(zip(control["observation_gyroperiods"],finest_orbit["times"])):
                    row=dict(common);row["observation_gyroperiods"]=cycles;row["actual_tau_hex"]=float(tau).hex();row["time_s"]=float(tau)/probe.gyro_rest_per_s
                    row["actual_k_per_m_hex"]=[float(k).hex() for k in fine_k];row["actual_probability_weights_hex"]=[float(w).hex() for w in fine_weights]
                    row["controlled_dimensionless_wave_coordinates"]={
                        "k_thermal_proton_gyroradius":[float(k)*float(scales["proton_thermal_gyroradius_m"]) for k in fine_k],
                        "k_ion_inertial_length":[float(k)*float(scales["ion_inertial_length_m"]) for k in fine_k],
                        "wave_frequency_over_proton_gyro_magnitude":[abs(float(k)*probe.wave_speed_mps)/float(scales["ion_gyrofrequency_radps"]) for k in fine_k]}
                    row["finite_time_full_lorentz"]=finest[sample];row["independent_dop853_subset"]=cross
                    low_digits,high_digits=numerical["linear_reference_digits"]
                    low=core.analytic_variance(probe,fine_k,fine_weights,amplitude,float(tau),low_digits,case["mixture"])
                    high=core.analytic_variance(probe,fine_k,fine_weights,amplitude,float(tau),high_digits,case["mixture"])
                    with mp.workdps(high_digits):
                        linear=float(high);linear_error=float(abs(high-low)+abs(high-core.exact_mp(linear)))
                        row["linear_phase_variance_high_precision"]=mp.nstr(high,high_digits)
                        row["linear_precision_serialization_error_high_precision"]=mp.nstr(abs(high-low)+abs(high-core.exact_mp(linear)),high_digits)
                        if high!=0 and linear==0:
                            row.update(status="UNRESOLVED",unresolved_reason="nonzero analytic finite-time response underflows binary64; high-precision value retained",physical_production_qualified=False)
                            row["linear_phase_variance"]=None;row["relative_variance_difference"]=None
                            case_rows.append(row);continue
                        if abs(high-low)+abs(high-core.exact_mp(linear))!=0 and linear_error==0:
                            # Round a diagnostic error bound UP; this never
                            # changes wave power, a physical response or a bin.
                            linear_error=float(np.nextafter(0.0,1.0))
                    full=finest[sample]["pitch_change_variance"]
                    orbit_values=[values[sample]["pitch_change_variance"] for _,values in orbit_levels]
                    phase_values=[values[sample]["pitch_change_variance"] for values in phase_levels]
                    mode_values=[values[sample]["pitch_change_variance"] for values in mode_levels]
                    orbit_error=abs(orbit_values[-1]-orbit_values[-2]);phase_error=abs(phase_values[-1]-phase_values[-2])
                    mode_error=0 if case["mixture"] is None else abs(mode_values[-1]-mode_values[-2])
                    pitch_cross=cross["maximum_absolute_pitch_difference"]+finest[sample]["momentum_encoding_pitch_resolution_estimate"]
                    cross_variance_error=2*math.sqrt(full)*pitch_cross+pitch_cross*pitch_cross
                    total_error=orbit_error+phase_error+mode_error+linear_error+cross_variance_error
                    pitch_allowance=float(numerical["absolute_pitch_target"])+float(numerical["relative_response_target"])*math.sqrt(max(full,linear))
                    # The declared relative target applies separately to the
                    # variance observable and to the DOP pitch crosscheck.
                    # The absolute pitch target contributes its square to the
                    # variance criterion; it is numerical accuracy, not power.
                    variance_allowance=float(numerical["absolute_pitch_target"])**2+float(numerical["relative_response_target"])*max(full,linear)
                    # This target controls numerical reference resolution only.
                    # It never accepts a physical model discrepancy or sets a
                    # production amplitude/k cutoff. Raw differences remain.
                    resolved=bool(max(orbit_error,phase_error,mode_error,linear_error)<=variance_allowance and \
                        finest[sample]["maximum_sampled_invariant_residual_over_initial_kinetic"]<=float(numerical["sampled_invariant_over_initial_kinetic_target"]) and \
                        cross["maximum_absolute_pitch_difference"]<=pitch_allowance)
                    row["status"]="RESOLVED" if resolved else "UNRESOLVED"
                    row["unresolved_reason"]=None if resolved else "declared orbit/phase/mode/linear/invariant numerical controls did not resolve this finite-time response"
                    row["linear_phase_variance"]=linear
                    row["raw_numerical_refinement"]={"orbit_variance":orbit_values,"phase_variance":phase_values,
                        "mode_variance":None if case["mixture"] is None else mode_values,
                        "observed_orbit_variance_order":observed_order(*orbit_values),
                        "actual_dt_tau":[orbit["dt_tau"] for orbit,_ in orbit_levels],
                        "actual_steps_per_observation":[orbit["steps"] for orbit,_ in orbit_levels]}
                    row["numerical_uncertainty_estimates"]={"orbit":orbit_error,"phase":phase_error,"mode":mode_error,
                        "linear_precision_and_serialization":linear_error,"DOP_subset_and_momentum_encoding":cross_variance_error,
                        "sum":total_error,"variance_target":variance_allowance,
                        "are_convergence_estimates_not_certified_bounds":True}
                    row["variance_difference"]=full-linear
                    row["relative_variance_difference"]=(full-linear)/linear if linear>total_error and linear>0 else None
                    row["relative_difference_status"]="defined" if row["relative_variance_difference"] is not None else "undefined: linear response zero or below numerical resolution"
                    row["finite_amplitude_difference_resolved"]=bool(resolved and abs(full-linear)>total_error)
                    row["physical_production_qualified"]=False
                    case_rows.append(row)
                rows.extend(case_rows)
            except (core.Unresolved,FloatingPointError,OverflowError,ZeroDivisionError) as error:
                for cycles in control["observation_gyroperiods"]:
                    rows.append(dict(common,observation_gyroperiods=cycles,status="UNRESOLVED",unresolved_reason=str(error),physical_production_qualified=False))
    # Source identity must remain stable even though execution uses archives.
    read_protocol(path)
    artifact={"schema":"srcsep-finite-amplitude-study-v1","protocol_sha256":digest(path.read_bytes()),
        "resolved_inputs":protocol["resolved_inputs"],"trajectory_updates":updates,"rows":rows,
        "resolved_rows":sum(r["status"]=="RESOLVED" for r in rows),"unresolved_rows":sum(r["status"]=="UNRESOLVED" for r in rows),
        "observable":"finite-time phase-ensemble response; no diffusivity inferred",
        "finite_band_policy":"only declared positive k; nonresonant finite-time response is not resonance-gap broadening",
        "physical_production_qualified":False,"physical_limits_selected":False}
    atomic_json(path.parent/"study.json",artifact)
    print(json.dumps({"resolved":artifact["resolved_rows"],"unresolved":artifact["unresolved_rows"],"artifact":str(path.parent/"study.json")}))
    return 1 if artifact["unresolved_rows"] else 0


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equilibrium",type=Path,default=DIRECTORY/"controlled_ep_maxwellian_baseline.json")
    parser.add_argument("--controls",type=Path,default=DIRECTORY/"controlled_finite_amplitude_pilot.json")
    parser.add_argument("--output-directory",type=Path)
    parser.add_argument("--run",action="store_true")
    parser.add_argument("--plot",action="store_true",help="plot the pinned study after an explicit run")
    parser.add_argument("--execute-protocol",type=Path)
    options=parser.parse_args()
    try:
        if options.execute_protocol:return execute(options.execute_protocol.resolve())
        if not options.output_directory:raise ValueError("--output-directory is required")
        path=prepare(options.equilibrium,options.controls,options.output_directory.resolve())
        print("Prepared "+str(path))
        if options.run:
            protocol=json.loads(path.read_text())
            result=subprocess.run(protocol["calculation_command_relative"],cwd=str(path.parent)).returncode
            if options.plot and result in (0,1):
                plotted=subprocess.run(protocol["plot_command_relative"],cwd=str(path.parent)).returncode
                if plotted:return plotted
            return result
        if options.plot:raise ValueError("--plot requires --run; existing artifacts use their archived plot_command_relative")
        return 0
    except (ValueError,KeyError,OSError,RuntimeError,TypeError) as error:
        print("ERROR "+str(error),file=sys.stderr);return 2


if __name__=="__main__":raise SystemExit(main())
