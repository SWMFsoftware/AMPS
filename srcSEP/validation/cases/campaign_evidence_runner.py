"""Campaign-level validation for EV01 and EV02 using observed CCMC SEP events.

The observational reference is never synthesized.  The committed table is a
machine-readable transcription of the GOES event quantities published by NASA
CCMC for the SHINE/ISWAT SEP model-validation challenge.  srcSEP creates only
model predictions.  EV01 fits a single CME-speed amplitude law on a sealed
training partition; EV02 repeats that fit without reading held-out targets and
then scores the predeclared hold-out events.
"""
from __future__ import annotations
import csv, datetime as dt, json, math, random, shutil, struct, time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from linked_case_common import atomic_json, finish_result, load_input, metric, read_csv, run_linked_model, sha256, write_csv

UTC = dt.timezone.utc


def _next_float_toward_positive_infinity(value: float) -> float:
    """Return the next IEEE-754 binary64 value greater than *value*.

    This is a small compatibility replacement for ``math.nextafter(value,
    math.inf)``.  ``math.nextafter`` was added to Python in 3.9, while srcSEP
    is also run on NASA/HPC systems whose system Python can be older.  The
    validation input only needs a strictly larger representable timestamp;
    operating directly on the binary64 bit pattern preserves the exact
    next-representable-value semantics without introducing an arbitrary time
    increment.

    The implementation handles finite positive and negative values, signed
    zero, and +infinity.  NaN is rejected because a source timestamp must be
    an ordered finite number.
    """
    value = float(value)
    if math.isnan(value):
        raise ValueError("cannot advance a NaN source timestamp")
    if value == math.inf:
        return value
    if value == 0.0:
        # Smallest positive subnormal IEEE-754 binary64 number.
        return struct.unpack(">d", struct.pack(">Q", 1))[0]

    bits = struct.unpack(">Q", struct.pack(">d", value))[0]
    # IEEE-754 bit ordering is monotone for positive values.  For negative
    # values, moving toward +infinity means decrementing the magnitude bits.
    bits = bits + 1 if value > 0.0 else bits - 1
    return struct.unpack(">d", struct.pack(">Q", bits))[0]

def _parse_time(text: str) -> Optional[dt.datetime]:
    text = (text or "").strip()
    if not text:
        return None
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00")).astimezone(UTC)

def _read_reference(path: Path) -> List[Dict[str,str]]:
    with path.open(newline="", encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    required={"event_id","threshold_mev","threshold_pfu","peak_flux_pfu","cme_speed_km_s"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("campaign observation table has an invalid schema")
    return rows

def _event_rows(rows: Sequence[Dict[str,str]]) -> Dict[str,List[Dict[str,str]]]:
    out: Dict[str,List[Dict[str,str]]]={}
    for r in rows: out.setdefault(r["event_id"],[]).append(r)
    return out

def _write_source(rows: Sequence[Dict[str,str]], path: Path) -> dt.datetime:
    """Create a flare-timing source history without altering the observations.

    The history is a *model input*, never an observational reference solution.
    For ordinary events it rises linearly from the reported flare onset to the
    reported peak and falls to zero at the reported flare end.  Some CCMC rows
    legitimately report the flare peak and end at the same timestamp (notably
    2012-07-12).  The native source reader requires strictly increasing sample
    times, so that zero-duration right-hand drop is encoded with the immediately
    next representable floating-point hour.  This introduces no resolved source
    duration and, unlike an arbitrary one-second/minute tail, does not invent an
    observational timescale.

    When the CCMC challenge provides no flare timing (2014-01-06), the CME
    21.5-Rs time remains the explicitly documented source-limited fallback and
    the caller marks the event accordingly.
    """
    r=rows[0]
    onset, peak, end = map(_parse_time, (r["flare_onset_utc"],r["flare_peak_utc"],r["flare_end_utc"]))
    if onset and peak and end:
        if not (onset <= peak <= end):
            raise ValueError("flare timing must satisfy onset <= peak <= end")
        epoch=onset
        peak_h=(peak-epoch).total_seconds()/3600.0
        end_h=(end-epoch).total_seconds()/3600.0
        points=[(0.0,0.0)]
        # A peak exactly at onset cannot coexist with a zero-rate sample at the
        # same time in the native strictly-monotone source format.  Represent
        # the observed right-hand peak at the immediately following float.
        if peak_h == 0.0:
            peak_h=_next_float_toward_positive_infinity(0.0)
        points.append((peak_h,1.0))
        if end_h > peak_h:
            shutdown_h=end_h
        else:
            # peak==end is an observed zero-duration decay, not a missing datum.
            # _next_float_toward_positive_infinity() preserves that limiting interpretation while meeting
            # the native parser's strict monotonic-time contract.
            shutdown_h=_next_float_toward_positive_infinity(peak_h)
        points.append((shutdown_h,0.0))
    else:
        epoch=_parse_time(r["cme_21p5rs_utc"])
        if epoch is None: raise ValueError("event has neither flare nor CME epoch")
        points=[(0.0,0.0),(0.25,1.0),(0.5,0.0)]
    if not points[-1][0] < 168.0:
        raise ValueError("observed source timing extends beyond the 168 h campaign window")
    points.append((168.0,0.0))
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["elapsed_hours","relative_source_rate"]); w.writerows(points)
    return epoch

def _args(case: Dict[str,Any], source: Path, event_id: str, seed_offset:int) -> List[str]:
    p,n=case["physics"],case["numerics"]
    return [
      "--observer",event_id,"--output-mode","time-profile","--source-radius-mode","inner-boundary",
      "--energies-mev",",".join(str(x) for x in p["energy_grid_mev"]),"--source-history-csv",str(source),
      "--injection-radius-solar-radii",str(p["injection_radius_solar_radii"]),"--observer-radius-au","1.0",
      "--solar-wind-speed-m-per-s",str(p["solar_wind_speed_m_per_s"]),"--shock-speed-m-per-s",str(p["shock_speed_m_per_s"]),
      "--solar-rotation-rate-rad-per-s",str(p["solar_rotation_rate_rad_per_s"]),"--launch-offset-s","0","--connection-delay-s","0",
      "--mfp-normalization-au",str(p["mean_free_path_normalization_au"]),"--mfp-radial-exponent",str(p["mean_free_path_radial_exponent"]),
      "--mfp-rigidity-exponent",str(p["mean_free_path_rigidity_exponent"]),"--injection-momentum-index",str(p["injection_momentum_index"]),
      "--particles-per-energy",str(n["particles_per_energy"]),"--time-step-s",str(n["time_step_s"]),
      "--duration-s",str(3600*float(n["duration_hours"])),"--output-cadence-s",str(3600*float(n["output_cadence_hours"])),
      "--campaign-seed",str(int(n["campaign_seed"])+seed_offset)]

def _integral_profiles(model: Sequence[Dict[str,str]], thresholds=(10.,100.)) -> Dict[float,List[Tuple[float,float]]]:
    by_t: Dict[float,List[Tuple[float,float]]]={}
    for r in model:
        by_t.setdefault(float(r["elapsed_hours"]),[]).append((float(r["energy_mev"]),float(r["relative_intensity"])))
    result={q:[] for q in thresholds}
    for t,spectrum in sorted(by_t.items()):
        spectrum=sorted(spectrum)
        for q in thresholds:
            pts=[p for p in spectrum if p[0]>=q]
            total=0.0
            for (e0,j0),(e1,j1) in zip(pts[:-1],pts[1:]): total += 0.5*(j0+j1)*(e1-e0)
            result[q].append((t,total))
    return result

def _fit_amplitude(training: Sequence[Dict[str,Any]]) -> Tuple[float,float]:
    """Fit one global log-amplitude covariate law: logA=a+b log(Vcme/1000)."""
    xs=[]; ys=[]
    for r in training:
        if r["model_peak"]>0 and r["obs_peak"]>0:
            xs.append(math.log10(r["cme_speed"]/1000.0)); ys.append(math.log10(r["obs_peak"]/r["model_peak"]))
    if len(xs)<2: raise RuntimeError("insufficient positive training rows for amplitude calibration")
    xm=sum(xs)/len(xs); ym=sum(ys)/len(ys); den=sum((x-xm)**2 for x in xs)
    b=0.0 if den==0 else sum((x-xm)*(y-ym) for x,y in zip(xs,ys))/den
    return ym-b*xm,b

def _bootstrap_fit(training: Sequence[Dict[str,Any]], *, seed: int, samples: int = 1000) -> Dict[str,Any]:
    """Estimate calibration-parameter stability by event-level bootstrap.

    Resampling is performed by *event*, not by threshold/channel row, so the
    >10 and >100 MeV measurements from one physical SEP event cannot leak into
    separate pseudo-independent samples.  This is the EV01 identifiability
    diagnostic requested by the validation plan; it never changes the
    observational reference values.
    """
    by_event: Dict[str,List[Dict[str,Any]]]={}
    for row in training:
        by_event.setdefault(str(row["event_id"]),[]).append(row)
    ids=sorted(by_event)
    if len(ids)<2:
        raise RuntimeError("at least two training events are required for bootstrap identifiability")
    rng=random.Random(seed)
    fits=[]
    for i in range(samples):
        picked=[rng.choice(ids) for _ in ids]
        draw=[row for eid in picked for row in by_event[eid]]
        try:
            a,b=_fit_amplitude(draw)
        except RuntimeError:
            continue
        fits.append({"replicate":i,"log10_amplitude_intercept":a,"log10_cme_speed_slope":b})
    if not fits:
        raise RuntimeError("bootstrap produced no valid calibration replicates")
    av=[x["log10_amplitude_intercept"] for x in fits]
    bv=[x["log10_cme_speed_slope"] for x in fits]
    def mean(v): return sum(v)/len(v)
    def std(v):
        m=mean(v); return math.sqrt(sum((x-m)**2 for x in v)/max(1,len(v)-1))
    am,bm=mean(av),mean(bv); sa,sb=std(av),std(bv)
    cov=sum((x-am)*(y-bm) for x,y in zip(av,bv))/max(1,len(av)-1)
    corr=0.0 if sa==0 or sb==0 else cov/(sa*sb)
    return {
        "replicates_requested":samples,
        "replicates_valid":len(fits),
        "intercept_mean":am,
        "intercept_std":sa,
        "slope_mean":bm,
        "slope_std":sb,
        "intercept_slope_correlation":corr,
        "replicates":fits,
    }

def _write_campaign_figures(scored: Sequence[Dict[str,Any]], output_dir: Path, case_id: str) -> List[Path]:
    """Create reproducible PNG/EPS campaign summaries from saved score rows.

    Plotting is intentionally downstream of the model calculation.  The figure
    contains only observed NASA CCMC quantities and predictions already written
    to the campaign score table, so regenerating it cannot alter scoring.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise RuntimeError("matplotlib is required for EV campaign figures") from exc

    obs=[float(r["obs_peak"]) for r in scored if float(r["obs_peak"])>0 and float(r["predicted_peak_pfu"])>0]
    pred=[float(r["predicted_peak_pfu"]) for r in scored if float(r["obs_peak"])>0 and float(r["predicted_peak_pfu"])>0]
    terr=[float(r["peak_time_error_h"]) for r in scored if r["peak_time_error_h"] != ""]
    labels=[f"{r['event_id']} >{int(float(r['threshold_mev']))}" for r in scored if r["peak_time_error_h"] != ""]

    fig,axes=plt.subplots(1,2,figsize=(12,5))
    axes[0].loglog(obs,pred,"o")
    if obs:
        lo=min(obs+pred); hi=max(obs+pred); axes[0].loglog([lo,hi],[lo,hi],"--")
    axes[0].set_xlabel("Observed peak flux [pfu]")
    axes[0].set_ylabel("Predicted peak flux [pfu]")
    axes[0].set_title(f"{case_id}: campaign peak-flux comparison")
    axes[0].grid(True,which="both",alpha=0.25)
    x=list(range(len(terr)))
    axes[1].axhline(0.0,linewidth=1.0)
    axes[1].plot(x,terr,"o")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels,rotation=75,ha="right",fontsize=7)
    axes[1].set_ylabel("Model - observed peak time [h]")
    axes[1].set_title(f"{case_id}: peak-time residuals")
    axes[1].grid(True,axis="y",alpha=0.25)
    fig.tight_layout()
    paths=[]
    for suffix in ("png","eps"):
        path=output_dir/f"{case_id}_campaign_summary.{suffix}"
        fig.savefig(path,dpi=180 if suffix=="png" else None,bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths

def _scores(rows: Sequence[Dict[str,Any]], a:float,b:float) -> Tuple[List[Dict[str,Any]],Dict[str,float]]:
    out=[]; timing=[]; logs=[]; tp=tn=fp=fn=0
    for r in rows:
        amp=10**(a+b*math.log10(r["cme_speed"]/1000.0)); pred=amp*r["model_peak"]
        observed_event=bool(r["obs_crossing"]); predicted_event=pred>=r["threshold_pfu"]
        if observed_event and predicted_event: tp+=1
        elif observed_event: fn+=1
        elif predicted_event: fp+=1
        else: tn+=1
        terr=None
        if r["obs_peak_time_h"] is not None and r["model_peak_time_h"] is not None:
            terr=r["model_peak_time_h"]-r["obs_peak_time_h"]; timing.append(abs(terr))
        if r["obs_peak"]>0 and pred>0: logs.append(math.log10(pred/r["obs_peak"]))
        out.append({**r,"amplitude_factor":amp,"predicted_peak_pfu":pred,"predicted_event":int(predicted_event),"observed_event":int(observed_event),"peak_time_error_h":"" if terr is None else terr,"log10_peak_ratio":"" if r["obs_peak"]<=0 or pred<=0 else math.log10(pred/r["obs_peak"])})
    rmse=math.sqrt(sum(x*x for x in logs)/len(logs)) if logs else float("inf")
    mae=sum(timing)/len(timing) if timing else float("inf")
    tss=(tp/(tp+fn) if tp+fn else 0)-(fp/(fp+tn) if fp+tn else 0)
    return out,{"tp":tp,"tn":tn,"fp":fp,"fn":fn,"log_rmse":rmse,"timing_mae_h":mae,"tss":tss}

def run_campaign_case(case_id:str, *, source_root:Path,input_path:Path,output_dir:Path,executable:Path,timeout:Optional[float]) -> Dict[str,Any]:
    started=time.monotonic(); case=load_input(input_path,case_id); output_dir.mkdir(parents=True,exist_ok=True); atomic_json(output_dir/"resolved_input.json",case)
    ref=(input_path.parent/case["reference"]["csv"]).resolve(); shutil.copyfile(ref,output_dir/"ccmc_observations.csv")
    ref_prov=(input_path.parent/case["reference"].get("provenance", "reference/provenance.json")).resolve()
    if ref_prov.exists(): shutil.copyfile(ref_prov,output_dir/"reference_provenance.json")
    rows=_read_reference(ref); grouped=_event_rows(rows); split=case["split"]
    selected=list(split["training"])+list(split.get("validation",[]))+(list(split.get("holdout",[])) if case_id=="EV02" else [])
    summaries=[]; artifacts=[output_dir/"resolved_input.json",output_dir/"ccmc_observations.csv"]
    if (output_dir/"reference_provenance.json").exists(): artifacts.append(output_dir/"reference_provenance.json")
    for idx,eid in enumerate(selected):
        if eid not in grouped: raise ValueError(f"split event {eid} absent from observational reference")
        edir=output_dir/"events"/eid; edir.mkdir(parents=True,exist_ok=True); source=edir/"source.csv"; epoch=_write_source(grouped[eid],source)
        native=run_linked_model(case_id=case_id,arguments=_args(case,source,eid,10000*idx),source_root=source_root,output_dir=edir,executable=executable,timeout=timeout)
        profiles=_integral_profiles(read_csv(native["model"]))
        for obs in grouped[eid]:
            q=float(obs["threshold_mev"]); prof=profiles[q]; mt,mv=max(prof,key=lambda x:x[1]); opt=_parse_time(obs["peak_utc"]); octime=_parse_time(obs["crossing_utc"])
            summaries.append({"event_id":eid,"partition":"training" if eid in split["training"] else ("validation" if eid in split.get("validation",[]) else "holdout"),"threshold_mev":q,"threshold_pfu":float(obs["threshold_pfu"]),"cme_speed":float(obs["cme_speed_km_s"]),"model_peak":mv,"model_peak_time_h":mt,"obs_peak":float(obs["peak_flux_pfu"] or 0.0),"obs_peak_time_h":None if opt is None else (opt-epoch).total_seconds()/3600,"obs_crossing":octime is not None,"source_limited":not bool(obs["flare_onset_utc"])})
        artifacts += [source,native["model"],native["report"],native["log"]]
    train=[r for r in summaries if r["partition"]=="training" and r["obs_peak"]>0]; a,b=_fit_amplitude(train)
    bootstrap=_bootstrap_fit(train,seed=int(case["numerics"]["campaign_seed"])+7919,samples=int(case.get("identifiability",{}).get("bootstrap_replicates",1000)))
    target=[r for r in summaries if r["partition"] in (("holdout",) if case_id=="EV02" else ("training","validation"))]
    scored,stats=_scores(target,a,b)
    comparison=output_dir/f"{case_id}_campaign_scores.csv"; fields=list(scored[0].keys()); write_csv(comparison,fields,scored); artifacts.append(comparison)
    bootstrap_csv=output_dir/f"{case_id}_calibration_bootstrap.csv"
    write_csv(bootstrap_csv,["replicate","log10_amplitude_intercept","log10_cme_speed_slope"],bootstrap["replicates"]); artifacts.append(bootstrap_csv)
    calibration={
        "schema":"srcsep-ev-calibration-v2",
        "fit_partition":"training-only",
        "log10_amplitude_intercept":a,
        "log10_cme_speed_slope":b,
        "bootstrap":{k:v for k,v in bootstrap.items() if k!="replicates"},
        "split":split,
        "reference_sha256":sha256(ref),
    }
    atomic_json(output_dir/"calibration.json",calibration); artifacts.append(output_dir/"calibration.json")
    artifacts.extend(_write_campaign_figures(scored,output_dir,case_id))
    # Campaign completeness is visible but is not conflated with model skill: the
    # committed 2021 challenge is a real-observation pilot; the full EV campaign
    # should replace it with the current CLEAR event/non-event package.
    metrics=[metric("observational_event_count",len(grouped),case["acceptance"]["pilot_event_count_min"],">=","events"),metric("peak_log10_rmse",stats["log_rmse"],case["acceptance"]["peak_log10_rmse_max"],"<=","dex"),metric("peak_time_mae",stats["timing_mae_h"],case["acceptance"]["peak_time_mae_hours_max"],"<=","h"),metric("true_skill_statistic",stats["tss"],case["acceptance"]["true_skill_statistic_min"],">=","score")]
    return finish_result(case_id=case_id,started=started,seed=int(case["numerics"]["campaign_seed"]),input_path=input_path,executable=executable,metrics=metrics,artifacts=artifacts,message=f"{case_id} real-observation campaign evidence passed")
