"""Linked-application runner and acceptance scoring for IV01-IV06.

Numerical rows must originate from the selected AMPS executable.  References
run as separate case-local Python processes, and every metric is recomputed
from the saved CSVs.  This preserves the same fail-closed evidence boundary as
CV01-CV12 while allowing each integrated case to own different observables.
"""
from __future__ import annotations
import math, subprocess, sys, time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence
from linked_case_common import atomic_json, finish_result, load_input, metric, read_csv, run_linked_model, sha256, write_csv

def _list(v:Iterable[Any])->str:return ",".join(str(x) for x in v)
def _mean(v:Sequence[float])->float:return sum(v)/len(v)
def _order(errors:Dict[float,float])->float:
 items=sorted((h,e) for h,e in errors.items() if h>0 and e>0)
 if len(items)<2:return 99.0
 return math.log(items[-1][1]/items[0][1])/math.log(items[-1][0]/items[0][0])

def _plot(cid:str,out:Path,title:str,x:Sequence[float],model:Sequence[float],ref:Sequence[float],xlabel:str,ylabel:str)->List[Path]:
 """Create review figures from the exact same series saved for scoring."""
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 fig,ax=plt.subplots(2,1,figsize=(8,7),sharex=True,gridspec_kw={"height_ratios":[3,1]})
 ax[0].plot(x,ref,"k-",label="independent reference");ax[0].plot(x,model,"o",label="linked srcSEP/AMPS");ax[0].set_title(f"{cid}: {title}");ax[0].set_ylabel(ylabel);ax[0].legend();ax[0].grid(alpha=.25)
 ax[1].axhline(0,color="black",linewidth=1);ax[1].plot(x,[a-b for a,b in zip(model,ref)],"o-");ax[1].set_xlabel(xlabel);ax[1].set_ylabel("residual");ax[1].grid(alpha=.25);fig.tight_layout();paths=[]
 for ext in ("png","eps"):
  p=out/f"{cid}_comparison.{ext}";fig.savefig(p,dpi=180 if ext=="png" else None);paths.append(p)
 plt.close(fig);return paths

def _solution(cid:str,out:Path,x,model,ref,units,title,xlabel,ylabel):
 p=out/f"{cid}_solution.csv";write_csv(p,("coordinate","numerical","analytical","residual","units"),({"coordinate":a,"numerical":b,"analytical":c,"residual":b-c,"units":units}for a,b,c in zip(x,model,ref)));return[p]+_plot(cid,out,title,x,model,ref,xlabel,ylabel)

def _args(cid:str,c:Dict[str,Any])->List[str]:
 p,n=c["physics"],c["numerics"];seed=["--campaign-seed",str(n["campaign_seed"])]
 if cid=="IV01":return["--inner-radius-m",str(p["inner_radius_m"]),"--outer-radius-m",str(p["outer_radius_m"]),"--particle-speed-m-per-s",str(p["particle_speed_m_per_s"]),"--wind-speeds-m-per-s",_list(p["wind_speeds_m_per_s"]),"--weak-d0-per-s",str(p["weak_d0_per_s"]),"--time-steps-s",_list(n["time_steps_s"]),"--particles",str(n["particles"])]+seed
 if cid=="IV02":return["--speed-m-per-s",str(p["speed_m_per_s"]),"--d0-per-s",str(p["d0_per_s"]),"--focusing-ratios",_list(p["focusing_ratios"]),"--duration-s",str(n["duration_s"]),"--time-steps-s",_list(n["time_steps_s"]),"--particles",str(n["particles"])]+seed
 if cid=="IV03":return["--steps",_list(n["steps"])]+seed
 if cid=="IV04":return["--resolutions",_list(n["resolutions"])]+seed
 if cid=="IV05":return["--shock-start-m",str(p["shock_start_m"]),"--node-m",str(p["node_m"]),"--shock-speeds-m-per-s",_list(p["shock_speeds_m_per_s"]),"--upstream-speed-m-per-s",str(p["upstream_speed_m_per_s"]),"--particle-speed-m-per-s",str(p["particle_speed_m_per_s"]),"--time-steps-s",_list(n["time_steps_s"]),"--particles",str(n["particles"])]+seed
 if cid=="IV06":return["--base-d0-per-s",str(p["base_d0_per_s"]),"--growth-per-s",str(p["growth_per_s"]),"--initial-wave-j",str(p["initial_wave_j"]),"--duration-s",str(n["duration_s"]),"--dt-s",str(n["dt_s"]),"--particles",str(n["particles"])]+seed
 raise ValueError(cid)

def _iv01(c,rows,refs,out):
 a=c["acceptance"];finest=min(map(float,c["numerics"]["time_steps_s"]));lookup={float(r["wind_m_per_s"]):r for r in refs};arrival=muerr=orient=0.;errors={};x=[];model=[];exact=[]
 for dt in map(float,c["numerics"]["time_steps_s"]):
  local=0
  for wind in map(float,c["physics"]["wind_speeds_m_per_s"]):
   g=[r for r in rows if r["mode"]=="ballistic" and float(r["dt_s"])==dt and float(r["wind_m_per_s"])==wind];r=lookup[wind]
   local=max(local,abs(float(g[0]["arrival_time_s"])/float(r["arrival_time_s"])-1));muerr=max(muerr,abs(float(g[0]["final_mu"])-float(r["final_mu"])));orient=max(orient,abs(float(g[0]["arrival_time_s"])/float(g[1]["arrival_time_s"])-1))
   if dt==finest:x.append(wind);model.append(float(g[0]["arrival_time_s"]));exact.append(float(r["arrival_time_s"]))
  errors[dt]=local
 arrival=errors[finest]
 # Reversed storage order uses identical keyed streams, so weak-scattering
 # ensemble means must also be bitwise-equivalent up to output roundoff.
 weak=[]
 for wind in map(float,c["physics"]["wind_speeds_m_per_s"]):
  for dt in (finest,):
   means=[]
   for orientation in (-1,1):means.append(_mean([float(r["final_mu"]) for r in rows if r["mode"]=="weak" and float(r["wind_m_per_s"])==wind and float(r["dt_s"])==dt and int(r["orientation"])==orientation]))
   weak.append(abs(means[0]-means[1]))
 metrics=[metric("ballistic_arrival_relative_error",arrival,a["ballistic_arrival_relative_error_max"],"<=","dimensionless"),metric("ballistic_mu_absolute_error",muerr,a["ballistic_mu_absolute_error_max"],"<=","dimensionless"),metric("orientation_relative_difference",max(orient,max(weak)),a["orientation_relative_difference_max"],"<=","dimensionless"),metric("arrival_refinement_order",_order(errors),a["refinement_order_min"],">=","dimensionless")]
 return metrics,_solution("IV01",out,x,model,exact,"s","Parker-spiral arrival","solar-wind speed [m/s]","arrival time [s]")

def _iv02(c,rows,refs,out):
 a=c["acceptance"];finest=min(map(float,c["numerics"]["time_steps_s"]));l1=initdiff=meanerr=norm=0.;x=[];model=[];exact=[]
 for ratio in map(float,c["physics"]["focusing_ratios"]):
  ref=[r for r in refs if float(r["ratio"])==ratio]
  groups=[]
  for initial in ("isotropic","beam"):
   g=[r for r in rows if float(r["ratio"])==ratio and r["initial"]==initial and float(r["dt_s"])==finest];groups.append(g)
   err=sum(abs(float(r["probability"])-float(q["probability"])) for r,q in zip(g,ref));l1=max(l1,err);norm=max(norm,abs(sum(float(r["probability"]) for r in g)-1));meanerr=max(meanerr,abs(float(g[0]["mean_mu"])-float(ref[0]["mean_mu"])))
  initdiff=max(initdiff,sum(abs(float(a1["probability"])-float(a2["probability"])) for a1,a2 in zip(*groups)))
  if ratio==c["physics"]["focusing_ratios"][-1]:
   x=[.5*(float(r["mu_left"])+float(r["mu_right"])) for r in groups[0]];model=[float(r["probability"]) for r in groups[0]];exact=[float(r["probability"]) for r in ref]
 metrics=[metric("stationary_L1_error",l1,a["stationary_l1_max"],"<=","probability"),metric("initial_condition_L1_difference",initdiff,a["initial_condition_l1_max"],"<=","probability"),metric("mean_mu_absolute_error",meanerr,a["mean_mu_error_max"],"<=","dimensionless"),metric("normalization_error",norm,a["normalization_error_max"],"<=","probability")]
 return metrics,_solution("IV02",out,x,model,exact,"bin_probability","scattering-focusing equilibrium","pitch-angle cosine","probability")

def _iv03(c,rows,refs,out):
 a=c["acceptance"];errors={};positive=0;x=[];model=[];exact=[]
 for h in map(float,c["numerics"]["steps"]):
  g=[r for r in rows if float(r["step"])==h];num=math.sqrt(sum(float(r["error"])**2 for r in g));den=math.sqrt(sum(float(r["exact_residual"])**2 for r in g));errors[h]=num/den;positive+=sum(int(r["positive"])==0 for r in g)
  if h==min(map(float,c["numerics"]["steps"])):x=[float(r["s"]) for r in g];model=[float(r["numerical_residual"]) for r in g];exact=[float(r["exact_residual"]) for r in g]
 metrics=[metric("finest_relative_L2_error",errors[min(errors)],a["relative_l2_error_max"],"<=","dimensionless"),metric("manufactured_refinement_order",_order(errors),a["refinement_order_min"],">=","dimensionless"),metric("positivity_failures",positive,a["positivity_failure_max"],"<=","count")]
 return metrics,_solution("IV03",out,x,model,exact,"residual","manufactured transport residual","s","transport residual")

def _iv04(c,rows,refs,out):
 a=c["acceptance"];uniform=cons=residual=0.;minimum=1.;finest=max(map(int,c["numerics"]["resolutions"]));x=[];model=[];exact=[]
 for motion in c["physics"]["motions"]:
  for n in map(int,c["numerics"]["resolutions"]):
   g=[r for r in rows if r["motion"]==motion and int(r["resolution"])==n];uniform=max(uniform,max(abs(float(r["energy"])-float(r["expected"])) for r in g));cons=max(cons,abs(float(g[0]["total"])-1));residual=max(residual,max(abs(float(r["remap_residual"])) for r in g));minimum=min(minimum,min(float(r["min_length"]) for r in g))
   if motion=="sinusoidal" and n==finest:x=[float(r["left"]) for r in g];model=[float(r["energy"]) for r in g];exact=[float(r["expected"]) for r in g]
 metrics=[metric("uniform_state_deviation",uniform,a["uniform_state_deviation_max"],"<=","J"),metric("conservation_relative_error",cons,a["conservation_relative_error_max"],"<=","dimensionless"),metric("remap_residual",residual,a["remap_residual_max"],"<=","J"),metric("minimum_segment_length",minimum,a["minimum_segment_length_min"],">=","m")]
 return metrics,_solution("IV04",out,x,model,exact,"J","moving-grid free-stream preservation","physical coordinate","cell integral [J]")

def _iv05(c,rows,refs,out):
 a=c["acceptance"];cross=frame=slope=node=0.;x=[];model=[];exact=[]
 for speed in map(float,c["physics"]["shock_speeds_m_per_s"]):
  for dt in map(float,c["numerics"]["time_steps_s"]):
   gs=[r for r in rows if float(r["shock_speed"])==speed and float(r["dt"])==dt];cross=max(cross,max(abs(float(r["crossing_time"])/float(r["exact_crossing_time"])-1) for r in gs));node+=sum(int(r["node_case"])!=1 for r in gs)
   moving=[r for r in gs if r["frame"]=="moving"];stationary=[r for r in gs if r["frame"]=="stationary"]
   frame=max(frame,max(abs(float(q["momentum"])/float(r["momentum"])-1) for q,r in zip(moving,stationary)))
  g=[r for r in rows if float(r["shock_speed"])==speed];slope=max(slope,abs(float(g[0]["q_sample"])-float(g[0]["q_exact"])));x.append(speed);model.append(float(g[0]["q_sample"]));exact.append(float(g[0]["q_exact"]))
 metrics=[metric("crossing_time_relative_error",cross,a["crossing_time_relative_error_max"],"<=","dimensionless"),metric("frame_momentum_relative_difference",frame,a["frame_momentum_relative_difference_max"],"<=","dimensionless"),metric("spectral_slope_absolute_error",slope,a["spectral_slope_absolute_error_max"],"<=","index"),metric("node_edge_case_failures",node,a["node_case_failure_max"],"<=","count")]
 return metrics,_solution("IV05",out,x,model,exact,"index","moving-shock frame equivalence","shock speed [m/s]","DSA q")

def _iv06(c,rows,refs,out):
 a=c["acceptance"];groups={name:[r for r in rows if r["control"]==name] for name in ("frozen","one-way","two-way")};frozen=abs(float(groups["frozen"][-1]["wave_energy"])/float(groups["frozen"][0]["wave_energy"])-1);reduction=1-float(groups["two-way"][-1]["streaming"])/float(groups["one-way"][-1]["streaming"]);energy=max(abs(float(r["ledger_residual"]))/float(groups[r["control"]][0]["total_energy"]) for r in rows);resbin=max(abs(int(r["resonant_bin"])-3) for r in rows);r0,r1=groups["one-way"][:2];observed=math.log(float(r1["wave_energy"])/float(r0["wave_energy"]))/(2*float(c["numerics"]["dt_s"])*float(r0["anisotropy"]));growth=abs(observed/float(c["physics"]["growth_per_s"])-1);x=[float(r["time"]) for r in groups["two-way"]];model=[float(r["streaming"]) for r in groups["two-way"]];exact=[float(r["streaming"]) for r in groups["one-way"]]
 metrics=[metric("early_growth_relative_error",growth,a["early_growth_relative_error_max"],"<=","dimensionless"),metric("two_way_streaming_reduction",reduction,a["two_way_streaming_reduction_min"],">=","fraction"),metric("resonant_bin_error",resbin,a["resonant_bin_error_max"],"<=","bins"),metric("total_energy_relative_error",energy,a["total_energy_relative_error_max"],"<=","dimensionless"),metric("frozen_wave_relative_change",frozen,a["frozen_wave_relative_change_max"],"<=","dimensionless")]
 return metrics,_solution("IV06",out,x,model,exact,"particle_sum_mu","self-generated feedback","time [s]","streaming proxy")

SCORERS={"IV01":_iv01,"IV02":_iv02,"IV03":_iv03,"IV04":_iv04,"IV05":_iv05,"IV06":_iv06}
def run_integrated_case(cid:str,*,source_root:Path,input_path:Path,output_dir:Path,executable:Path,timeout:Optional[float]):
 """Execute a complete linked IV case and return registry-compatible evidence."""
 started=time.monotonic();c=load_input(input_path,cid);output_dir.mkdir(parents=True,exist_ok=True);resolved=output_dir/"resolved_input.json";atomic_json(resolved,c);native=run_linked_model(case_id=cid,arguments=_args(cid,c),source_root=source_root,output_dir=output_dir,executable=executable,timeout=timeout);reference=output_dir/f"{cid}_reference.csv";script=Path(__file__).parent/cid/"reference_solution.py";q=subprocess.run([sys.executable,str(script),"--input",str(resolved),"--output",str(reference)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=timeout); 
 if q.returncode:raise RuntimeError(f"{cid} reference failed: {q.stdout}")
 metrics,derived=SCORERS[cid](c,read_csv(native["model"]),read_csv(reference),output_dir);provenance=output_dir/"provenance.json";atomic_json(provenance,{"schema":"srcsep-validation-provenance-v1","case_id":cid,"execution":"linked-srcsep-amps","sha256":{"executable":sha256(executable),"input":sha256(resolved),"model":sha256(native["model"]),"reference":sha256(reference)}});artifacts=[resolved,native["manifest"],native["model"],native["report"],native["junit"],native["log"],reference,provenance]+derived;return finish_result(case_id=cid,started=started,seed=int(c["numerics"]["campaign_seed"]),input_path=input_path,executable=executable,metrics=metrics,artifacts=artifacts,message=f"{cid} linked integrated validation passed")
