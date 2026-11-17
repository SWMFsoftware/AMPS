#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_report.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_report.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: dfb9ba03ed59c13d7beb0eef54e184d523b8562837b9ccf63c243e367a2b6e45
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (option C): the per-fit numbers of the lead's item 5 (costs incl. P_m, rms ln overall / before the gap /
#   late, bounds, Q_H2O at 10:33:10, active seen area, seen area holding 50 / 90 % of the seen production, NG / RG shares
#   of the active area; D2 and c9 for comparison), nonneg_bins.txt in fit/final, run-root README.md and HOW-TO-RERUN.md
#   (--docs-only after the GPU check); lambda and P_m described from the recorded settings.
# RUN: jobs/localnn_report.sbatch (localnn_anch / localnn_post2 / model_ext_frozen/localnn_mp_20261115T020233Z).
# 2026-11-15: stage 2: the Problem rebuilt with problem_from_meta (f_m, the held term).
# RUN 2026-11-15: stage-2 chain post jobs 132018, 132027, 132036, 132045.
# ==== end of header ====
"""Option C (local non-negative basis): the summary numbers of one stage-1 fit (lead 2026-11-14, item 5), the bins
table and the run root's README.md / HOW-TO-RERUN.md. CPU (one pass over the base's kernel records for the active-area
shares). Nothing in the fit's segment is changed; one new file is added there (nonneg_bins.txt).

Numbers (all from the final segment, <run root>/fit/final):
* data cost, crossing cost, penalty lambda_C x_c^T Lg x_c, objective, lambda_C; termination, nfev (fit_manifest)
* rms ln NG / RG overall (the used rows), before the gap (et <= the last NG point before the gap), late window (the 95
  observations from 10:33:10 UTC, W10x5)
* parameters at bounds (within 1e-8 max(1, |lb|)): the six multipliers at 1, kappa at 0, grain A, B at 0; c_j at 0
  (fraction of all 2745 nodes and of the nodes with seen support)
* Q_H2O and Q_CO2 at 10:33:10 UTC [s^-1] = sum over facets of F area (nonneg_maps.source_rates with the validated
  geometry's label t11 = 10:33:10: cos SZA, cast shadow, r_h; the law of fit_cops), all facets and the seen facets
* active area: A_e > 1 % of the srcCG median activity of the seen facets (median over the seen facets, one value per
  facet, of elements.bin AH2O = the run's srcCG map); the fraction of the seen area that is active (also for D2's and
  the softplus c9 final activity, for comparison)
* the smallest seen area holding 50 % and 90 % of the seen H2O production at 10:33:10 UTC (facets in decreasing order
  of F area); also for D2 and c9
* NG / RG shares from the active area: the active facets' contribution to the modelled surface pressure (frozen
  nonneg_contrib arithmetic per record; late observations W10x5 of the grid rows) divided by the total surface model
  (H2O + kappa CO2), over the 95 late observations and over every used observation row (min / mean / max)
* bins table (nonneg_problem.write_bin_table's columns; sigma = the local Gauss-Newton sigma with the smoothness rows,
  as fit_cops' posterior sigma with its Tikhonov prior)

CLI: localnn_report.py --runroot RR --base BASE --d2 D2 --point-table PT --design DESIGN --grid GRID --geometry GEOM
     --contrib-code FROZEN --c9-ft C9FT --c9-fit C9FIT [--procs N] -> RR/analysis/report.json, report.txt; RR/fit/final/nonneg_bins.txt;
     RR/README.md, RR/HOW-TO-RERUN.md (an existing one is moved aside as *.before_<UTC>); --docs-only: only the two
     documents, from RR/analysis/report.json (after the GPU evaluation)
"""

import argparse
import datetime as dt
import glob
import json
import multiprocessing as mp
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_design as ld        # noqa: E402
import localnn_problem as lp       # noqa: E402

GEOM_LABEL = "t11"                  # nonneg_maps MAP_TIMES[11] = 2016-09-30T10:33:10
ACTIVE_FRACTION = 0.01
_S = {}


def utcstamp():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def work(task):
    r, ip = task
    S = _S
    nc, NM = S["CS"]["nc"], S["CS"]["NM"]
    kp, h, rec = NM.read_kernel_full(Path(S["base"]) / "kernels" / ("point=%04d.bin" % ip))
    nc._G["kappa"] = 0.0
    _, ((ng, rg), (nga, rga)), _ = nc.record_contributions(rec, h, (S["A"], S["A_act"]), "fp32", r)
    return r, float(np.sum(ng)), float(np.sum(rg)), float(np.sum(nga)), float(np.sum(rga))


def area_holding(P_e, area, fracs):
    o = np.argsort(-P_e, kind="stable")
    cs = np.cumsum(P_e[o])
    tot = cs[-1]
    out = {}
    for f in fracs:
        k = int(np.searchsorted(cs, f * tot, side="left"))
        out["%d%%" % round(100 * f)] = float(np.sum(area[o[:k + 1]]))
    return out


def production_stats(A, kappa, el, seen, geom, NM, helio, sza, AU, thr):
    G = NM.read_geometry(geom)
    Tt = G["times"][GEOM_LABEL]
    cos, sh = NM.geometry_time(geom, GEOM_LABEL, el["nfaces"])
    F_h, F_c, _, _ = NM.source_rates(cos, sh, Tt["rh_m"], helio, sza, A, el["ACO2"], kappa, AU)
    area = el["area"]
    P_e = F_h * area
    sa = float(np.sum(area[seen]))
    hold = area_holding(P_e[seen], area[seen], (0.5, 0.9))
    return {"utc": Tt["utc"], "Q_H2O_all_s-1": float(np.sum(P_e)), "Q_H2O_seen_s-1": float(np.sum(P_e[seen])),
            "Q_CO2_all_s-1": float(np.sum(F_c * area)),
            "seen_area_m2": sa, "active_seen_area_fraction": float(np.sum(area[seen & (A > thr)]) / sa),
            "seen_area_holding_production_m2": hold,
            "seen_area_holding_production_fraction": {k: v / sa for k, v in hold.items()}}


def gn_sigma(P, x):
    """Local Gauss-Newton sigma (fit_cops.fit_one: H = Jd^T Jd + R^T R, v = |r|^2 / (N - p_eff)); R = smoothness rows."""
    Jd = P.data_jacobian(x)
    rd = P.data_residuals(x)
    S = np.zeros((P.R.shape[0], len(x)))
    S[:, P.ic] = np.sqrt(2.0 * P.lam) * P.R
    H = Jd.T @ Jd + S.T @ S
    Hi = np.linalg.pinv(H, rcond=1e-14)
    peff = float(np.trace(Hi @ (Jd.T @ Jd)))
    v = float(rd @ rd) / max(len(rd) - peff, 1.0)
    return np.sqrt(np.maximum(np.diag(v * Hi), 0.0)), peff


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("runroot", "base", "d2", "point-table", "design", "grid", "geometry", "contrib-code", "c9-ft", "c9-fit"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    ap.add_argument("--docs-only", action="store_true", help="rewrite README.md / HOW-TO-RERUN.md from analysis/report.json "
                    "(e.g. after the GPU evaluation); the old ones are moved aside")
    args = ap.parse_args(argv)
    t0 = time.time()
    RR = Path(args.runroot).resolve()
    seg = (RR / "fit" / "final").resolve()
    adir = RR / "analysis"
    if args.docs_only:
        write_docs(RR, json.load(open(adir / "report.json")))
        lp.say("docs rewritten in %s" % RR)
        return 0
    if adir.exists():
        raise SystemExit("%s exists" % adir)
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        x = np.asarray(z["x"], float)
        meta = json.loads(str(z["meta_json"]))
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        A = np.asarray(z["A"], np.float64)
        seen = np.asarray(z["seen"], bool)
        bins = np.asarray(z["bin"]).astype(np.intp)
    design = lp.load_design(args.design)
    grid = lp.load_grid(args.grid)
    D = lp.Data(args.base, args.d2, args.point_table, design)
    P = lp.problem_from_meta(meta, D, design, grid)
    c, m, kappa, gA, gB = P.unpack(x)
    cst, rows = P.evaluate(x)
    r = cst["r"]
    rn, rr = r[:D.NNG] / D.wNG, r[D.NNG:] / D.wRG
    bef_n = D.et[D.rows_ng] <= D.et_before
    bef_r = D.et[D.rows_rg] <= D.et_before
    late = set(int(v) for v in D.late_rows)
    lat_n = np.array([int(v) in late for v in D.rows_ng])
    lat_r = np.array([int(v) in late for v in D.rows_rg])
    rms = lambda a: float(np.sqrt(np.mean(a ** 2)))      # noqa: E731
    rmsd = {"overall": [rms(rn), rms(rr)], "before_gap": [rms(rn[bef_n]), rms(rr[bef_r])],
            "late_window": [rms(rn[lat_n]), rms(rr[lat_r])],
            "points": {"overall": [int(len(rn)), int(len(rr))], "before_gap": [int(bef_n.sum()), int(bef_r.sum())],
                       "late_window": [int(lat_n.sum()), int(lat_r.sum())]}}
    tol = 1e-8
    node_seen = np.asarray(grid["node_seen"], bool)
    c_zero = x[P.ic] - P.lo[P.ic] <= tol
    bounds = {"multipliers_at_1": [lp.BIN_NAMES[1 + k] for k in range(6) if x[P.im][k] - 1.0 <= tol],
              "kappa_at_0": bool(x[P.ik] <= tol), "grain_A_at_0": bool(x[P.iA] <= tol),
              "grain_B_at_0": bool(x[P.iB] <= tol),
              "c_at_0_fraction_all_nodes": float(c_zero.mean()), "c_at_0_nodes": int(c_zero.sum()),
              "c_at_0_fraction_seen_nodes": float(c_zero[node_seen].mean()) if node_seen.size else None}
    # srcCG median and the active mask
    CS = ld.contrib_setup(args.base, args.contrib_code, lp.USER_T0)
    el, NM, fc = CS["el"], CS["NM"], CS["fc"]
    helio, sza = fc.header_value(CS["header"], "HeliocentricModel"), fc.header_value(CS["header"], "SZAMode")
    med = float(np.median(el["AH2O"][seen]))
    thr = ACTIVE_FRACTION * med
    active = A > thr
    prod = {"localnn": production_stats(A, kappa, el, seen, args.geometry, NM, helio, sza, float(fc.AU), thr)}
    with np.load(Path(args.d2) / "fit_restart.npz", allow_pickle=False) as z:
        A2 = np.maximum(np.asarray(z["activity_linear"], np.float64), 0.0)
        k2 = float(np.asarray(z["kappa"]).reshape(-1)[0])
    prod["D2"] = production_stats(A2, k2, el, seen, args.geometry, NM, helio, sza, float(fc.AU), thr)
    del A2
    with np.load(Path(args.c9_ft) / "finetune_state.npz", allow_pickle=False) as z:
        A9 = np.asarray(z["A_final"], np.float64)
    with np.load(Path(args.c9_fit) / "fit_restart.npz", allow_pickle=False) as z:
        k9 = float(np.asarray(z["kappa"]).reshape(-1)[0])
    prod["c9"] = production_stats(A9, k9, el, seen, args.geometry, NM, helio, sza, float(fc.AU), thr)
    del A9
    lp.say("production: " + json.dumps(prod))
    # active-area shares (one pass over the records)
    _S.update({"CS": CS, "base": args.base, "A": A, "A_act": np.where(active, A, 0.0)})
    tot = np.zeros((D.nrow, 4))
    with mp.get_context("fork").Pool(args.procs) as pool:
        for rw, a1, a2, a3, a4 in pool.imap_unordered(work, [(i, int(D.ipoint[i])) for i in range(D.nrow)], chunksize=1):
            tot[rw] = (a1, a2, a3, a4)
    sN = D.E @ (tot[:, 0] + kappa * design["C_NG"])
    sR = D.E @ (tot[:, 1] + kappa * design["C_RG"])
    aN, aR = D.E @ tot[:, 2], D.E @ tot[:, 3]
    shN, shR = aN / sN, aR / sR
    lr = D.late_rows
    used_n, used_r = D.rows_ng, D.rows_rg
    st = lambda a: {"min": float(np.min(a)), "mean": float(np.mean(a)), "max": float(np.max(a))}   # noqa: E731
    shares = {"late_window_95": {"NG": st(shN[lr]), "RG": st(shR[lr])},
              "all_used_observations": {"NG": st(shN[used_n]), "RG": st(shR[used_r])},
              "direct_vs_design_total_max_rel": [float(np.max(np.abs(sN / rows["surface_NG"] - 1))),
                                                 float(np.max(np.abs(sR / rows["surface_RG"] - 1)))]}
    lp.say("active-area shares: " + json.dumps(shares))
    sig, peff = gn_sigma(P, x)
    rep = {"created_utc": lp.utc_now(), "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "command": sys.argv, "runroot": str(RR), "segment": str(seg),
           "lambda_C": meta["lambda_C"], "lambda_source": meta["lambda_source"], "templates": meta["templates"],
           "costs": {k: cst[k] for k in ("data_cost", "crossing_cost", "penalty", "objective", "wrong_side",
                                         "mult_penalty") if k in cst},
           "mult_penalty_f": float((meta.get("mult_penalty") or {}).get("f_m", 0.0)),
           "start_file": meta.get("start_file"),
           "costs_segment_manifest": meta["costs"],
           "termination": {"stop": meta["stop"], "message": meta["message"], "status": meta["status"],
                           "max_nfev": meta.get("max_nfev"),
                           "nfev": meta["nfev"], "nfev_total": meta["nfev_total"], "njev": meta["njev"],
                           "optimality": meta["optimality"], "segment": meta["segment"]},
           "rms_ln": rmsd, "bounds": bounds, "multipliers": list(m), "kappa": kappa, "grain_A_uPa": gA,
           "grain_B_uPa_per_100s": gB, "srcCG_median_activity_seen": med, "active_threshold": thr,
           "active_definition": "A_e > %g x median over the seen facets of elements.bin AH2O" % ACTIVE_FRACTION,
           "production_1033": prod, "active_area_shares": shares,
           "gn_sigma_multipliers": {lp.BIN_NAMES[1 + k]: float(sig[P.im][k]) for k in range(6)},
           "gn_sigma_kappa_grain": [float(s) for s in sig[P.ik:]], "p_eff": peff,
           "code_sha256": {p.name: lp.sha256_file(p) for p in sorted(HERE.glob("localnn_*.py"))},
           "wall_seconds": time.time() - t0}
    stage = adir.parent / ("analysis.incomplete_%s" % utcstamp())
    stage.mkdir()
    with open(stage / "report.json", "x") as f:
        json.dump(rep, f, indent=1, sort_keys=True, default=float)
    with open(stage / "report.txt", "x") as f:
        f.write(report_text(rep) + "\n")
    stage.rename(adir)
    write_bins(seg, el, seen, A, bins, m, sig, P, meta)
    write_docs(RR, rep)
    lp.say("report done: %s (%.0f s)" % (adir, time.time() - t0))
    return 0


def lambda_text(src):
    """One line on how lambda_C was set (the start rule, a rule file such as the solution-anchored rule, or explicit)."""
    rec = src.get("rule_record")
    if rec:
        return "%s; f = %s (rule file %s)" % (rec.get("description", rec.get("rule")), src.get("factor"), src.get("rule_file"))
    if "factor" in src:
        return "%s; f = %s" % (src.get("rule"), src["factor"])
    return "explicit"


def report_text(r):
    p = r["production_1033"]
    L = ["Option C stage-1 fit %s" % r["runroot"],
         "lambda_C %.9e (%s), templates %s" % (r["lambda_C"], lambda_text(r["lambda_source"]), r["templates"]),
         "data cost %.10e  crossing %.3e  penalty %.6e  objective %.10e%s" % (
             r["costs"]["data_cost"], r["costs"]["crossing_cost"], r["costs"]["penalty"], r["costs"]["objective"],
             "  multiplier penalty P_m %.6e (f_m %g)" % (r["costs"]["mult_penalty"], r.get("mult_penalty_f", 0.0))
             if "mult_penalty" in r["costs"] else ""),
         "termination: %s (%s), nfev %s (total %s), optimality %s" % (
             r["termination"]["stop"], r["termination"]["message"], r["termination"]["nfev"],
             r["termination"]["nfev_total"], r["termination"]["optimality"]),
         "rms ln NG / RG: overall %.4f / %.4f, before the gap %.4f / %.4f, late window %.4f / %.4f" % (
             tuple(r["rms_ln"]["overall"]) + tuple(r["rms_ln"]["before_gap"]) + tuple(r["rms_ln"]["late_window"])),
         "multipliers %s; at 1: %s; kappa %.3e (at 0 %s); grain A %.6g uPa (at 0 %s), B %.3e (at 0 %s)" % (
             ["%.4g" % v for v in r["multipliers"]], r["bounds"]["multipliers_at_1"], r["kappa"], r["bounds"]["kappa_at_0"],
             r["grain_A_uPa"], r["bounds"]["grain_A_at_0"], r["grain_B_uPa_per_100s"], r["bounds"]["grain_B_at_0"]),
         "c_j at 0: %d nodes (%.3f of all, %.3f of the seen-support nodes)" % (
             r["bounds"]["c_at_0_nodes"], r["bounds"]["c_at_0_fraction_all_nodes"], r["bounds"]["c_at_0_fraction_seen_nodes"]),
         "active threshold %.4e (%s)" % (r["active_threshold"], r["active_definition"])]
    for k in ("localnn", "D2", "c9"):
        q = p[k]
        L.append("%-8s Q_H2O(10:33:10) %.4e s^-1 (seen %.4e); active seen area %.4f; seen area holding 50%% / 90%% of the "
                 "seen production %.4e / %.4e m^2 (%.4f / %.4f of the seen area)" % (
                     k, q["Q_H2O_all_s-1"], q["Q_H2O_seen_s-1"], q["active_seen_area_fraction"],
                     q["seen_area_holding_production_m2"]["50%"], q["seen_area_holding_production_m2"]["90%"],
                     q["seen_area_holding_production_fraction"]["50%"], q["seen_area_holding_production_fraction"]["90%"]))
    s = r["active_area_shares"]
    L.append("share of the modelled surface NG / RG from the active area: late window mean %.4f / %.4f (min %.4f / %.4f); "
             "all used observations mean %.4f / %.4f" % (
                 s["late_window_95"]["NG"]["mean"], s["late_window_95"]["RG"]["mean"], s["late_window_95"]["NG"]["min"],
                 s["late_window_95"]["RG"]["min"], s["all_used_observations"]["NG"]["mean"],
                 s["all_used_observations"]["RG"]["mean"]))
    return "\n".join(L)


def write_bins(seg, el, seen, A, bins, m, sig, P, meta):
    area = el["area"]
    tot = float(np.sum(area))
    prod = area * A
    ptot, pseen = float(np.sum(prod)), float(np.sum(prod[seen]))
    z = A <= 0.0
    labels = ("[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]")
    L = ["# nonneg_bins.txt -- written by model_ext/localnn_report.py (%s) for this localnn segment" % lp.utc_now(),
         "# %s: A = m_b (sum_j c_j phi_j + TS) (templates %s), local hat-function basis, slope-bin multipliers" % (
             lp.PARAMETERIZATION, meta["templates"]),
         "# zero = A = 0; seen = facet contributes to a fit point; production share = sum_bin area*A / sum_all area*A "
         "(activity x area; q(SZA) not applied); *_seen over seen facets. sigma: local Gauss-Newton sigma with the "
         "smoothness rows (as fit_cops' posterior sigma with its Tikhonov prior).",
         "#bin interval_deg multiplier sigma facets area_fraction seen_area_fraction_of_bin zero_facets_frac "
         "zero_area_frac seen_zero_area_frac production_share production_share_seen"]
    for k, name in enumerate(lp.BIN_NAMES):
        sel = bins == k
        ss = sel & seen
        a_bin, a_seen = float(area[sel].sum()), float(area[ss].sum())
        sg = "fixed" if k == 0 else "%.3e" % sig[P.im][k - 1]
        L.append("%s %s %.9e %s %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e" % (
            name, labels[k], m[k], sg, int(sel.sum()), a_bin / tot, a_seen / a_bin, float(np.sum(sel & z)) / sel.sum(),
            float(area[sel & z].sum()) / a_bin, float(area[ss & z].sum()) / a_seen, float(prod[sel].sum()) / ptot,
            float(prod[ss].sum()) / pseen))
    sa = float(area[seen].sum())
    L.append("all [0,180] - - %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e" % (
        len(area), 1.0, sa / tot, float(np.sum(z)) / len(area), float(area[z].sum()) / tot,
        float(area[seen & z].sum()) / sa, 1.0, 1.0))
    L.append("# total production sum(area*A) = %.9e (seen %.9e)" % (ptot, pseen))
    L.append("# multiplier lower bound 1 (user 2026-11-12 / lead 2026-11-14), no upper bound; at the lower bound "
             "(m - 1 <= 1e-8): %s" % [lp.BIN_NAMES[k] for k in range(1, 7) if m[k] - 1.0 <= 1e-8])
    with open(seg / "nonneg_bins.txt", "x") as f:
        f.write("\n".join(L) + "\n")


def write_docs(RR, rep):
    for name in ("README.md", "HOW-TO-RERUN.md"):
        p = RR / name
        if p.exists():
            p.rename(RR / ("%s.before_%s" % (name, utcstamp())))
    seg = Path(rep["segment"])
    gpu = RR / "gpu_eval" / "gpu_check.json"
    gtxt = "not run yet (jobs/localnn_gpu_eval.sbatch SEG=%s/fit/final)" % RR
    if gpu.exists():
        g = json.load(open(gpu))
        gtxt = ("%s: GPU data cost %.13e vs the CPU state at the GPU's grain %.13e (rel %.1e; the D2 path cannot hold "
                "the grain: GPU grain A %.9g uPa, B %.3g vs the fit's A %.9g, B %.3g); COPS_model.dat pNG / pRG printed "
                "differences %d / %d (max rel %.1e / %.1e); late window surface NG / RG max rel %.1e / %.1e" % (
                    "PASSED" if g["passed"] else "FAILED", g["cost_gpu_manifest_data_cost"], g["cost_cpu_with_gpu_grain"],
                    g["cost_rel_diff"], g["grain_gpu"]["A_uPa"], g["grain_gpu"]["B_uPa_per_100s"],
                    g["grain_fit"]["A_uPa"], g["grain_fit"]["B_uPa_per_100s"], g["COPS_model_pNG"]["printed_differences"],
                    g["COPS_model_pRG"]["printed_differences"], g["COPS_model_pNG"]["max_rel"],
                    g["COPS_model_pRG"]["max_rel"], g["late_surface_pNG"]["max_rel"], g["late_surface_pRG"]["max_rel"]))
    dc = sorted(glob.glob(str(RR / "fit" / "direct_check_seg*" / "direct_check.json")))
    dtxt = "-"
    if dc:
        d = json.load(open(dc[-1]))
        dtxt = "%s: data cost direct %.13e vs segment %.13e (rel %.1e); design vs direct max rel %s" % (
            "PASSED" if d["passed"] else "FAILED", d["data_cost_direct"], d["data_cost_segment"],
            d["data_cost_rel_diff"], d["design_vs_direct_instantaneous_max_rel"])
    L = ["# %s" % RR.name, "",
         "Option C of the non-negative activity models (local non-negative basis), stage-1 fit; user / lead decisions "
         "2026-11-14 (~02:35 CST). Model and code: `srcCG/LiouvilleCOPS_cuda/model_ext/README.md`, section \"Option C: "
         "local non-negative basis\"; `localnn_problem.py` docstring.", "",
         "* A_e = m_b(e) (sum_j c_j phi_j(e) + A_s sum_k ft_ck(e)), c_j >= 0, 2745 hat functions of the 150 m coarse grid "
         "`runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/grid_150m`; D2's ft_c1..ft_c6 held at 1, added as is (their sum "
         "is >= 0 on every facet: 265 positive, 0 negative; `template_signs_20261114.json`);",
         "* multipliers: 7 slope bins, 0-40 deg held at 1, six >= 1 (start 1); kappa >= 0, grain A, B >= 0 (start D2's); "
         "W10x5 from 10:33:10 UTC, grain at the time stamp, crossing weight 10, original T law held; lambda 0 Tikhonov;",
         "* smoothness lambda_C x_c^T Lg x_c: lambda_C = %.9e (%s); start = the NNLS projection of D2's activity "
         "minus the template part, m = 1, kappa and grain from D2;" % (rep["lambda_C"], lambda_text(rep["lambda_source"])),
         "* scipy least_squares TRF (bounds, x_scale jac, exact), ftol = xtol = gtol = 1e-10, max_nfev %s, CPU." % (
             rep["termination"].get("max_nfev", 4000))] + ([
         "* multiplier penalty (user 2026-11-14 ~17:30 CST): P_m = f_m C_ref sum_b (ln m_b)^2, f_m = %g, C_ref = 4.268806657481e-03; "
         "P_m = %.6e at the end; start: %s" % (rep["mult_penalty_f"], rep["costs"]["mult_penalty"], rep.get("start_file"))]
         if "mult_penalty" in rep["costs"] else []) + ["",
         "## Result (`analysis/report.txt`, `analysis/report.json`)", "", "```", report_text(rep), "```", "",
         "Reaching max_nfev is not convergence. Costs are the data cost 0.5 sum r^2 of fit_cops; the penalty and the "
         "crossing cost are separate.", "",
         "## Checks", "",
         "* direct facet sums of the final state at every row (`fit/direct_check_*`): %s" % dtxt,
         "* GPU standard evaluation (`gpu_eval/`): %s" % gtxt, "",
         "## Files", "",
         "* `fit/final` -> the last segment: `fit_result.txt`, `localnn_state.npz`, `nonneg_facets.npz` (A, P, u, TS, "
         "bin, seen, held = 0, area), `COPS_model.dat`, `late_window_model.dat` (CPU design route), `nonneg_bins.txt`, "
         "`fit_manifest.json`, `README.txt`, `HOW-TO-RERUN.md` (segment command)",
         "* `gpu_eval/fit`: the standard fitter products of the GPU evaluation (D2's path, clip_refit.py --variant D2 with "
         "this activity as a held template; grain A, B = the fitter's start, not the fitted ones: D2's code cannot hold "
         "them); `gpu_eval/gpu_check.json`",
         "* `analysis/`: this report; `provenance-<job>/`: job script, commands.sh, logs, git, modules, toolchain, pip list"]
    (RR / "README.md").write_text("\n".join(L) + "\n")
    cmds = []
    for p in sorted(glob.glob(str(RR / "provenance-*" / "commands.sh"))):
        cmds += ["### %s" % Path(p).parent.name, "```bash", Path(p).read_text().rstrip(), "```", ""]
    gp = sorted(glob.glob(str(RR / "gpu_eval" / "provenance-*" / "commands.sh")))
    for p in gp:
        cmds += ["### gpu_eval/%s" % Path(p).parent.name, "```bash", Path(p).read_text().rstrip(), "```", ""]
    H = ["# HOW TO RERUN %s" % RR.name, "",
         "AMPS commit and tree state: `provenance-<job>/git.txt`; code: the frozen snapshot named in each "
         "`provenance-<job>/env.txt` / commands.sh (MODEL_EXT; SHA256SUMS there; version-controlled copies in "
         "srcCG/LiouvilleCOPS_cuda/model_ext, RUN_COPIES.sha256).", "",
         "## From the existing inputs (grid, design, start)", "",
         "Job submission (as run):", "```bash", "cd /data001/highresolu_vt/cops_output/jobs",
         "LX=%s" % os.environ.get("MODEL_EXT", "<the frozen localnn snapshot>"),
         "sbatch --time=12:00:00 --export=ALL,MODEL_EXT=$LX,TEMPLATES=%s,%s,NFEV=%s,RUN=%s localnn_fit.sbatch"
         % (rep["templates"], ("SMOOTH_LAMBDA=%r,SMOOTH_RULE_FILE=%s" % (rep["lambda_C"], rep["lambda_source"]["rule_file"])
                               if rep["lambda_source"].get("rule_file") else
                               "SMOOTH_FACTOR=%s" % rep["lambda_source"].get("factor")) +
            (",MULT_PENALTY_F=%r,START=%s" % (rep["mult_penalty_f"], rep.get("start_file"))
             if rep.get("mult_penalty_f") else ""),
            rep["termination"].get("max_nfev", 4000), RR.name.rsplit("_", 1)[0]),
         "# after an exit 75 (wall-time stop): the same with RUNROOT=%s" % RR,
         "sbatch --export=ALL,MODEL_EXT=$LX,SEG=%s/fit/final localnn_gpu_eval.sbatch     # GPU standard evaluation" % RR,
         "sbatch --export=ALL,MODEL_EXT=$LX,RUNROOT=%s localnn_report.sbatch            # this report" % RR, "```", "",
         "Exact commands with every argument (each job's provenance commands.sh; `heavy` is the job scripts' shell "
         "function, define it first):", "```bash",
         "heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec \"$@\"' heavy \"$@\"; }", "```", ""] + cmds + [
         "## From scratch", "",
         "1. base, D2, slope bins: model_ext/README.md, 'Softplus and squared-harmonics models', steps 1-3;",
         "2. grid, design, verification, starts, tests: `sbatch --export=ALL,MODEL_EXT=$LX localnn_build.sbatch` (job 125004; "
         "exact commands in `runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/grid_150m/provenance-125004/commands.sh`);",
         "3. the fit, GPU evaluation and report as above."]
    (RR / "HOW-TO-RERUN.md").write_text("\n".join(H) + "\n")


if __name__ == "__main__":
    sys.exit(main())
