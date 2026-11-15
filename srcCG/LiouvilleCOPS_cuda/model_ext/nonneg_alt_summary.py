#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_alt_summary.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_alt_summary.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: bddab2947c4fffff8d17bd962143e66156d60a5cd6e852d53bb6c14d19b0a4d8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Summary of one finished non-negative-model fit (lead 2026-11-14): costs and penalties, termination,
#   rms ln NG/RG overall / before the gap / late window, parameters at bounds, Q_H2O at 10:33:10, active seen area,
#   seen area holding 50 / 90 % of the seen production, NG/RG shares from the active area.
# USED BY: jobs/nonneg_alt_summary.sbatch (snapshot model_ext_frozen/nonneg_summary_20261114T084417Z; test 125030).
# 2026-11-14: run copy updated (cost labels: fit_cops cost = data + spread penalty; snapshot nonneg_summary_20261114T090853Z)
# 2026-11-14: run copy updated (u0, u_bar, (u_bar - u0)^2 and the two penalty parts; the start record; snapshot nonneg_summary_20261115T020314Z)
# ==== end of header ====
"""Summary of ONE finished fit of the non-negative activity models (model_ext/nonneg_*.py; written for the stage-1 fits
of options A, B, D, user 2026-11-14, and usable for softplus / squared fits): the numbers the lead asked for
(2026-11-14 ~02:35 CST), from the fit's own products only (nothing is refitted).

    nonneg_alt_summary.py --fit <fit dir: fit/final> --base <base> --geometry <validated nonneg_maps geometry>
        --out <summary.txt> [--procs N]
Writes <out> (text) and <out>.json (every number, the definitions, the inputs with sha256). Never overwrites.

DEFINITIONS
  costs        fit_manifest.json optimizer_results[0]: data_cost (0.5 r.r of the data rows) and cost (fit_cops: data +
               the option-D spread penalty P = f C_ref Var_w(u), without the crossing rows); the NG = RG crossing
               penalty and the objective with it (model_extension.final); P from the fit's records (model_check.
               alternative.spread.final); termination = scipy's message, status, nfev
  rms ln       sqrt(mean(ln(p_model / p_obs)^2)) over the fit's used rows (fit_restart.npz pred_*, ref_*, use_*; p_model
               as the residuals see it: W10x5 for the late observations, + grain for NG), per gauge, overall / before
               the gap (UTC < 2016-09-30T10:33:10) / late window (>= 10:33:10); the overall values are checked against
               fit_result.txt's rms_ln_residual lines
  at bounds    a parameter with x - lb <= 1e-8 max(1, |lb|) (or ub - x <= ...) in the optimizer's coordinates (fit_restart
               x, lo, hi); scipy TRF keeps iterates ~1e-10 inside a bound
  Q_H2O        sum_e F_H2O,e area_e at 2016-09-30T10:33:10 UTC, F_H2O = q_H2O(cos SZA, night/shadow, r_h) A_e
               (nonneg_maps.source_rates, the fitter's law; geometry label t11), A = nonneg_facets.npz A (the fit's
               activity incl. multipliers and the held fine-tuning term); for all facets and the seen facets
  active area  facets with A_e > 1 % of the median of A_orig (srcCG map, elements.bin AH2O) over the seen facets
               (plain median over facets); reported as the fraction of the seen area
  production   seen production at 10:33:10: P_e = F_H2O,e area_e over the seen facets; the seen area holding 50 % /
               90 % of it = the smallest area of seen facets, taken in decreasing order of F_H2O,e (production per
               area; ties: lower facet index first), whose production reaches 50 % / 90 % of the seen total; also
               the same with the time-independent activity (P_e = A_e area_e, order by A_e)
  NG/RG shares from the active area: per base row the surface pressure of the active facets / that of all facets
               (fit_cops.Model term by term with the fit's FP32 records, nonneg_contrib.record_contributions: H2O + CO2
               per facet, the fitted T and kappa; late observations W10x5 of the 1-s grid rows); the NG grain term is
               not a facet contribution; reported over the used rows: mean, min, and the ratio of sums, overall /
               before the gap / late window. Verification: the surface sums (+ grain for NG) against the fit's
               predictions (fit_restart.npz pred_*), max rel.
"""

import argparse
import datetime as dt
import json
import math
import multiprocessing as mp
import os
import re
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
for _p in (str(HERE), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import nonneg_maps as NM            # noqa: E402  readers, source law, geometry
import nonneg_contrib as NC         # noqa: E402  per-record contributions (the fitted forward model)
import positive_problem as pp       # noqa: E402  time keys, late layout (W10x5), point tables

SCHEMA = "nonneg_alt_summary_v1"
LATE_UTC = "2016-09-30T10:33:10"
AT_BOUND_RTOL = 1e-8
ACTIVE_FRACTION_OF_MEDIAN = 0.01
_W = {}


def say(text):
    print("nonneg_alt_summary: [%s] %s" % (dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


def row_work(task):
    row, ip = task
    kp, h, rec = NM.read_kernel_full(Path(_W["base"]) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
    face, ((ng, rg),), _ = NC.record_contributions(rec, h, (_W["A"],), "fp32", row)
    act = _W["active"][face]
    return (row, float(np.sum(ng)), float(np.sum(rg)), float(np.sum(ng[act])), float(np.sum(rg[act])),
            int(np.sum(ng < 0) + np.sum(rg < 0)), int(np.sum(~np.isfinite(ng)) + np.sum(~np.isfinite(rg))), h["et"])


def area_holding(order_key, prod, area, levels=(0.5, 0.9)):
    """Smallest area (facets in decreasing order_key, ties lower index first) whose prod reaches each level."""
    idx = np.lexsort((np.arange(len(order_key)), -order_key))
    cs = np.cumsum(prod[idx])
    tot = float(cs[-1]) if len(cs) else 0.0
    ca = np.cumsum(area[idx])
    out = {}
    for lv in levels:
        k = int(np.searchsorted(cs, lv * tot, side="left"))
        out["%d%%" % round(lv * 100)] = {"area_m2": float(ca[k]), "facets": k + 1}
    return out, tot


def stats_share(sh, m):
    if not np.any(m):
        return None
    return {"rows": int(np.sum(m)), "mean": float(np.mean(sh[m])), "min": float(np.min(sh[m])),
            "median": float(np.median(sh[m]))}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fit", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--geometry", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    args = ap.parse_args(argv)
    t0 = time.time()
    out = Path(args.out)
    if out.exists() or Path(str(out) + ".json").exists():
        raise SystemExit("%s exists (never overwritten)" % out)
    fit, base = Path(args.fit).resolve(), Path(args.base).resolve()
    man = json.load(open(fit / "fit_manifest.json"))
    rec = man["model_extension"]["nonneg"]
    model = rec["config"]["model"]
    opt = man["optimizer_results"][0]
    with np.load(fit / "fit_restart.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        names = list(meta["names"])
        x, lo, hi = (np.asarray(z[k], float) for k in ("x", "lo", "hi"))
        point_ids = np.asarray(z["point_ids"], np.int64)
        utc_r = [str(u) for u in z["utc"]]
        ref = {"NG": np.asarray(z["ref_ng"], float), "RG": np.asarray(z["ref_rg"], float)}
        pred = {"NG": np.asarray(z["pred_ng"], float), "RG": np.asarray(z["pred_rg"], float)}
        use = {"NG": np.asarray(z["use_ng"]) > 0, "RG": np.asarray(z["use_rg"]) > 0}
        T = tuple(float(v) for v in np.asarray(z["T"]))
        kappa = float(np.asarray(z["kappa"]).reshape(-1)[0])
    if meta.get("status") != "complete":
        raise SystemExit("%s: restart status %r, not a complete fit" % (fit, meta.get("status")))
    S = {"schema": SCHEMA, "created_utc": NM.utc_now(), "host": socket.gethostname(),
         "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv, "fit": str(fit), "model": model,
         "config": {k: rec["config"].get(k) for k in ("model", "sos_rank", "bg_min", "spread_f", "spread_cref",
                                                       "spread_u0", "slope_mult_min", "start_choice")
                    if k in rec["config"]},
         "start": {k: (rec.get("start_restart") or {}).get(k) for k in ("start_choice", "path", "source_fit_directory",
                                                                        "source_cost", "source_nfev")}}
    # ---------------------------------------------------------------- costs and termination
    alt = (rec.get("model_check") or {}).get("alternative") or {}
    sp = alt.get("spread")
    final = (man["model_extension"].get("final") or {})
    S["costs"] = {"cost_fitter": float(opt["cost"]), "data_cost": float(opt["data_cost"]),
                  "crossing_penalty": final.get("penalty_cost"),
                  "objective_with_crossing": final.get("objective_cost_with_penalty"),
                  "spread_penalty": (sp.get("final") or {}).get("penalty") if sp else None,
                  "spread": ({k: sp.get(k) for k in ("f", "C_ref")} | {"final": sp.get("final")}) if sp else None}
    S["termination"] = {"message": (meta.get("optimizer") or {}).get("message", opt.get("message")),
                        "status": (meta.get("optimizer") or {}).get("status"),
                        "nfev": sum(int(s.get("nfev", 0)) for s in meta.get("segments", [])),
                        "segments": len(meta.get("segments", [])),
                        "optimality": (meta.get("optimizer") or {}).get("optimality")}
    # ---------------------------------------------------------------- rms ln
    late_key = pp.time_key(LATE_UTC)
    keys = np.array([pp.time_key(u) for u in utc_r], dtype=np.int64)
    early = keys < late_key
    rms = {}
    for g in ("NG", "RG"):
        m = use[g] & (ref[g] > 0) & (pred[g] > 0)
        d = np.log(pred[g] / np.where(m, ref[g], 1.0))
        rms[g] = {lab: (float(np.sqrt(np.mean(d[mm] ** 2))) if np.any(mm) else None, int(np.sum(mm)))
                  for lab, mm in (("all", m), ("before_gap", m & early), ("late_window", m & ~early))}
    fres = (fit / "fit_result.txt").read_text()
    chk = {}
    for g in ("NG", "RG"):
        mm = re.search(r"^rms_ln_residual_%s (\S+)" % g, fres, re.M)
        chk[g] = {"fit_result_txt": float(mm.group(1)) if mm else None,
                  "agrees_4_digits": bool(mm) and abs(float(mm.group(1)) - rms[g]["all"][0]) <= 5e-5}
    S["rms_ln"] = {"definition": "sqrt(mean(ln(p_model / p_obs)^2)) over the used rows; before the gap: UTC < %s; "
                                 "late window: >= %s" % (LATE_UTC, LATE_UTC), "values": rms, "check_fit_result": chk}
    # ---------------------------------------------------------------- bounds
    def at_lo(k):
        return bool(np.isfinite(lo[k]) and x[k] - lo[k] <= AT_BOUND_RTOL * max(1.0, abs(lo[k])))

    def at_hi(k):
        return bool(np.isfinite(hi[k]) and hi[k] - x[k] <= AT_BOUND_RTOL * max(1.0, abs(hi[k])))
    sel = [k for k, n in enumerate(names) if n.startswith("c_slope_") or n in ("A_bg", "kappa", "grain_A_uPa",
                                                                               "grain_B_uPa_per_100s")]
    S["parameters"] = {names[k]: {"x": float(x[k]), "lower": float(lo[k]), "upper": float(hi[k]),
                                  "at_lower": at_lo(k), "at_upper": at_hi(k)} for k in sel}
    S["at_bounds_all"] = [names[k] for k in range(len(names)) if at_lo(k) or at_hi(k)]
    S["parameters_note"] = "optimizer coordinates (fit_restart.npz x, lo, hi): multipliers and A_bg physical; kappa, grain as fitted"
    # ---------------------------------------------------------------- activity, production, Q
    el = NM.read_elements(base / "elements.bin")
    nf = el["nfaces"]
    area, aorig = el["area"], el["AH2O"]
    with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
        A = np.asarray(z["A"], np.float64)
        seen = np.asarray(z["seen"]).astype(bool)
    if len(A) != nf:
        raise SystemExit("nonneg_facets.npz has %d facets, the base %d" % (len(A), nf))
    if not (np.all(np.isfinite(A)) and np.all(A >= 0)):
        raise SystemExit("A is negative or not finite")
    med = float(np.median(aorig[seen]))
    thr = ACTIVE_FRACTION_OF_MEDIAN * med
    active = A > thr
    a_seen = float(area[seen].sum())
    S["active_area"] = {"threshold": thr, "srcCG_median_A_orig_seen": med, "fraction_of_median": ACTIVE_FRACTION_OF_MEDIAN,
                        "seen_area_m2": a_seen, "seen_facets": int(seen.sum()),
                        "active_seen_area_m2": float(area[seen & active].sum()),
                        "active_seen_area_fraction": float(area[seen & active].sum()) / a_seen,
                        "zero_seen_area_fraction": float(area[seen & (A == 0)].sum()) / a_seen,
                        "active_facets_all": int(active.sum())}
    gval = json.load(open(Path(args.geometry) / "geometry_validation.json"))
    if not gval.get("passed") or Path(gval["base"]) != base:
        raise SystemExit("geometry %s not validated against %s" % (args.geometry, base))
    fc, header, cols, rows = NM.base_table(base)
    helio, sza = fc.header_value(header, "HeliocentricModel"), fc.header_value(header, "SZAMode")
    G_ = NM.read_geometry(args.geometry)
    Tt = G_["times"]["t11"]
    if Tt["utc"] != LATE_UTC:
        raise SystemExit("geometry t11 is %s, expected %s" % (Tt["utc"], LATE_UTC))
    cos, sh = NM.geometry_time(args.geometry, "t11", nf)
    F_h, F_c, night, qf = NM.source_rates(cos, sh, Tt["rh_m"], helio, sza, A, el["ACO2"], kappa, fc.AU)
    P_e = F_h * area
    hold_t, tot_t = area_holding(F_h[seen], P_e[seen], area[seen])
    hold_a, tot_a = area_holding(A[seen], (A * area)[seen], area[seen])
    for d_ in (hold_t, hold_a):
        for v in d_.values():
            v["fraction_of_seen_area"] = v["area_m2"] / a_seen
    S["Q_H2O_103310"] = {"utc": LATE_UTC, "Q_all_s-1": float(P_e.sum()), "Q_seen_s-1": float(P_e[seen].sum()),
                         "Q_srcCG_base_s-1": Tt["Q_H2O_base"], "Q_CO2_all_s-1": float((F_c * area).sum())}
    S["seen_production"] = {"at_103310": {"holding": hold_t, "total_s-1": tot_t,
                                          "definition": "P_e = F_H2O area at 10:33:10, ordered by F_H2O"},
                            "activity": {"holding": hold_a, "total": tot_a,
                                         "definition": "P_e = A area (time independent), ordered by A"}}
    # ---------------------------------------------------------------- NG / RG shares from the active area
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("restart rows differ from the base rows")
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    bm, cNG, cRG, gconst = NC.gauge_constants(fc, header, TRG)
    ptab = rec["config"]["point_table"]
    layout = pp.LateLayout(pp.read_point_table(ptab), pp.read_run_points(base), ipoint)
    W = layout.W.tocsr()
    from cops_cuda import device_store
    NC._G.update({"A": A, "ACO2": el["ACO2"], "T": T, "kappa": kappa, "bm": bm, "cNG": cNG, "cRG": cRG,
                  "helio": helio, "AU": float(fc.AU), "COS75": float(fc.COS75), "KT": float(fc.KT),
                  "fp32_record_fields": device_store.fp32_record_fields})
    _W.update({"base": str(base), "A": A, "active": active})
    nrow = len(ipoint)
    tasks = sorted([(r, int(ipoint[r])) for r in range(nrow)], key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    say("contributions: %d rows, %d processes" % (nrow, args.procs))
    tot = {k: np.zeros(nrow) for k in ("ng", "rg", "ng_act", "rg_act")}
    et = np.full(nrow, np.nan)
    bad = 0
    with mp.get_context("fork").Pool(args.procs) as pool:
        for r, ng, rg, nga, rga, neg, nonf, et_r in pool.imap_unordered(row_work, tasks, chunksize=1):
            tot["ng"][r], tot["rg"][r], tot["ng_act"][r], tot["rg_act"][r] = ng, rg, nga, rga
            et[r] = et_r
            bad += neg + nonf
    if bad:
        raise SystemExit("%d negative or non-finite record contributions" % bad)
    lr, gr = layout.late_rows, layout.grid_rows
    for k in tot:
        tot[k][lr] = W @ tot[k][gr]
    grain_et = None
    gt = json.load(open(fit / "grain_term.json")) if (fit / "grain_term.json").exists() else None
    t_end = float(gt["t_end_et"]) if gt else None
    gA = float(x[names.index("grain_A_uPa")]) if "grain_A_uPa" in names else 0.0
    gB = float(x[names.index("grain_B_uPa_per_100s")]) if "grain_B_uPa_per_100s" in names else 0.0
    if t_end is not None and np.all(np.isfinite(et)):
        grain_et = np.where(et >= t_end, np.maximum(0.0, gA * 1e-6 - gB * 1e-8 * (et - t_end)), 0.0)
    ver = {}
    for g, k in (("NG", "ng"), ("RG", "rg")):
        m = use[g]
        mod = tot[k] + (grain_et if (g == "NG" and grain_et is not None) else 0.0)
        ver[g] = float(np.max(np.abs(mod[m] / pred[g][m] - 1.0)))
    shares = {}
    for g, k in (("NG", "ng"), ("RG", "rg")):
        m = use[g] & (tot[k] > 0)
        sh_ = np.where(m, tot[k + "_act"] / np.where(m, tot[k], 1.0), np.nan)
        shares[g] = {lab: stats_share(sh_, mm) for lab, mm in (("all", m), ("before_gap", m & early),
                                                              ("late_window", m & ~early))}
        for lab, mm in (("all", m), ("before_gap", m & early), ("late_window", m & ~early)):
            if shares[g][lab] is not None:
                shares[g][lab]["ratio_of_sums"] = float(tot[k + "_act"][mm].sum() / tot[k][mm].sum())
        if g == "NG" and grain_et is not None:
            gm = m & ~early
            shares[g]["grain_share_of_NG_late_window_mean"] = (float(np.mean(grain_et[gm] / (tot[k][gm] + grain_et[gm])))
                                                               if np.any(gm) else None)
    S["shares_from_active_area"] = {"values": shares, "verification_max_rel_vs_fit_predictions": ver,
                                    "grain_t_end_et": t_end, "note": "surface pressures (H2O + CO2 per facet); the "
                                    "NG grain term is not a facet contribution"}
    S["inputs_sha256"] = {n: NM.sha256_file(fit / n) for n in ("fit_manifest.json", "fit_restart.npz", "nonneg_facets.npz",
                                                               "fit_result.txt")}
    S["code_sha256"] = {p.name: NM.sha256_file(p) for p in (Path(__file__).resolve(), Path(NM.__file__).resolve(),
                                                            Path(NC.__file__).resolve(), Path(pp.__file__).resolve())}
    S["geometry"] = str(Path(args.geometry).resolve())
    S["base"] = str(base)
    S["wall_seconds"] = time.time() - t0
    write_text(out, S)
    with open(str(out) + ".json", "x") as f:
        json.dump(S, f, indent=1, sort_keys=True, default=float)
    say("wrote %s (+ .json), %.0f s" % (out, time.time() - t0))
    return 0


def fmt(v, f="%.4g"):
    return "-" if v is None else (f % v)


def write_text(out, S):
    c, t, r = S["costs"], S["termination"], S["rms_ln"]["values"]
    L = ["# %s -- %s, written %s by model_ext/nonneg_alt_summary.py (job %s)" % (Path(out).name, S["fit"], S["created_utc"],
                                                                              S["slurm_job_id"]),
         "model %s; config %s" % (S["model"], json.dumps(S["config"], sort_keys=True)),
         "start: %s" % json.dumps(S["start"], sort_keys=True),
         "",
         "COSTS: data %.10e, spread penalty %s, NG = RG crossing penalty %s; fit_cops cost (data + spread penalty) %.10e, "
         "objective with the crossing rows %s%s" % (
             c["data_cost"], fmt(c["spread_penalty"], "%.10e"), fmt(c["crossing_penalty"], "%.6e"), c["cost_fitter"],
             fmt(c["objective_with_crossing"], "%.10e"),
             ("; spread f %r, C_ref %r, Var_w(u) %s (rms %s neper)%s" % (
                 c["spread"]["f"], c["spread"]["C_ref"], fmt((c["spread"]["final"] or {}).get("var_w_u"), "%.6e"),
                 fmt((c["spread"]["final"] or {}).get("rms_spread"), "%.4f"),
                 ("; fixed level u0 %r: u_bar %s, (u_bar - u0)^2 %s; penalty = variance part %s + level part %s" % (
                     (c["spread"]["final"] or {}).get("u0"), fmt((c["spread"]["final"] or {}).get("mean_u"), "%.6f"),
                     fmt((c["spread"]["final"] or {}).get("dev2_mean_u0"), "%.6e"),
                     fmt((c["spread"]["final"] or {}).get("penalty_variance_part"), "%.6e"),
                     fmt((c["spread"]["final"] or {}).get("penalty_level_part"), "%.6e")))
                 if (c["spread"]["final"] or {}).get("u0") is not None else "")) if c["spread"] else ""),
         "TERMINATION: %s (status %s), nfev %d in %d segment(s), first-order optimality %s" % (
             t["message"], t["status"], t["nfev"], t["segments"], fmt(t["optimality"], "%.3e")),
         "RMS LN (%s):" % S["rms_ln"]["definition"]]
    for g in ("NG", "RG"):
        L.append("  %s: all %s (%d rows), before the gap %s (%d), late window %s (%d); fit_result.txt %s" % (
            g, fmt(r[g]["all"][0], "%.5f"), r[g]["all"][1], fmt(r[g]["before_gap"][0], "%.5f"), r[g]["before_gap"][1],
            fmt(r[g]["late_window"][0], "%.5f"), r[g]["late_window"][1], S["rms_ln"]["check_fit_result"][g]))
    L.append("PARAMETERS (multipliers, A_bg, kappa, grain; optimizer coordinates):")
    for n, v in S["parameters"].items():
        L.append("  %-22s %.10g  [%.6g, %.6g]%s" % (n, v["x"], v["lower"], v["upper"],
                                                  " AT LOWER BOUND" if v["at_lower"] else (" AT UPPER BOUND" if v["at_upper"] else "")))
    L.append("  every parameter at a bound: %s" % (", ".join(S["at_bounds_all"]) or "none"))
    q = S["Q_H2O_103310"]
    L.append("Q_H2O at %s UTC: %.6e s^-1 (seen facets %.6e; srcCG base %.6e)" % (q["utc"], q["Q_all_s-1"], q["Q_seen_s-1"],
                                                                                q["Q_srcCG_base_s-1"]))
    a = S["active_area"]
    L.append("ACTIVE AREA (A > %.0f %% of the srcCG median A_orig over the seen facets = %.6e): %.4f of the seen area "
             "(%.6e of %.6e m^2); A = 0 on %.4f of the seen area" % (
                 100 * a["fraction_of_median"], a["threshold"], a["active_seen_area_fraction"], a["active_seen_area_m2"],
                 a["seen_area_m2"], a["zero_seen_area_fraction"]))
    for key, lab in (("at_103310", "production at 10:33:10 (F_H2O area)"), ("activity", "activity (A area)")):
        h = S["seen_production"][key]["holding"]
        L.append("SEEN AREA HOLDING 50 %% / 90 %% of the seen %s: %.4f / %.4f of the seen area (%.4e / %.4e m^2)" % (
            lab, h["50%"]["fraction_of_seen_area"], h["90%"]["fraction_of_seen_area"], h["50%"]["area_m2"],
            h["90%"]["area_m2"]))
    sh = S["shares_from_active_area"]["values"]
    L.append("SHARES FROM THE ACTIVE AREA (surface pressure of the active facets / all facets, used rows; verification "
             "vs the fit's predictions max rel %s):" % json.dumps(S["shares_from_active_area"]["verification_max_rel_vs_fit_predictions"]))
    for g in ("NG", "RG"):
        parts = []
        for lab in ("all", "before_gap", "late_window"):
            v = sh[g][lab]
            parts.append("%s mean %.6f min %.6f ratio of sums %.6f (%d rows)" % (lab, v["mean"], v["min"], v["ratio_of_sums"],
                                                                                v["rows"]) if v else "%s -" % lab)
        L.append("  %s: %s" % (g, "; ".join(parts)))
    if sh["NG"].get("grain_share_of_NG_late_window_mean") is not None:
        L.append("  NG grain share of the NG model in the late window (mean): %.4f" % sh["NG"]["grain_share_of_NG_late_window_mean"])
    L.append("")
    L.append("Definitions: model_ext/nonneg_alt_summary.py docstring; every number with its inputs: %s.json" % Path(out).name)
    with open(out, "x") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
