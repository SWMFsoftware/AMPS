#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nightfix_diag.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nightfix_diag.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 57a35f0068ec968bfb63c143f55285cbbef5a8821403018831067b2a3a466c9f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (2026-11-18, H2O night-side fix). Re-evaluates ONE final state with the activity held fixed under the old
#   (q_min/q_max 0.2) and the new (0.02) H2O law: data cost, rms ln NG/RG (all / before the gap / late window), Q_H2O
#   and the seen production at 10:33:10, night/shadow shares, the post-hoc uniform A-scale bound, per-row arrays;
#   --base-run: fit_cops.main's base consistency check for both laws. No fitting.
# USED BY: jobs/nightfix_diag.sbatch (v1 frozen model_ext_frozen/nightfix_diag_20261118T073237Z, jobs 136436-136444;
#   v2 nightfix_diag_20261118T074014Z, jobs 136457-136464); results tests/nightfix_diagnostic_20261118/.
# RUN: see tests/nightfix_diagnostic_20261118/HOW-TO-RERUN.md
# ==== end of header ====
"""H2O night-side fix diagnostic (2026-11-18): re-evaluate ONE final state with the activity A HELD FIXED under the old
(q_min/q_max = 0.2) and the new (0.02) H2O production law. No fitting: A (incl. slope multipliers, the held term and the
late fine-tuning factors), T, kappa and the grain term stay as the state has them; only q_min of H2O changes.

User decision 2026-11-18 (relayed by the lead): "(b) -- use 0.02; x10 was not intended; use the same law but with
decreased activity". Lead's task 4: for each final state the data cost, rms ln NG / RG (overall and late window),
Q_H2O total and seen at 10:33:10 UTC, the share of Q_H2O and of the seen production from night / cast-shadow facets,
for the OLD law (reproducing the recorded cost to rounding first) and the NEW law. This tells how much a refit is needed.

THE MODEL (the arithmetic of the contribution tool, reused unchanged): nonneg_contrib.record_contributions of the given
contribution code (a nonneg_contrib snapshot whose nonneg_maps.py has H2O_QMIN_RATIO; q from nonneg_maps.q_factors,
which reads the constant at call time). Per kernel record: FP32 record fields as the fits ran (--precision mixed =
cops_cuda.device_store.fp32_record_fields; shadow folded into cos = -2), FP64 arithmetic; night = cos32 <= 0;
q_H2O = qmin if night else qmin + (qmax - qmin) min(cos, 1), qmin = r qmax, r = 0.2 (old) or 0.02 (new). The temperature
law depends on cos only (unchanged). Surface model per base row = sum of the records; late observations = W10x5 of the
1-s grid rows (positive_problem.LateLayout of the state's own positive_problem.py and point table); + the grain term
max(0, A_g - B_g (t - t_end)) for t >= t_end at the late NG observations (the state's finetune_state.npz grain_A_uPa,
grain_B_uPa_per_100s; t_end of grain_term.json). The OLD-law per-row model is checked against the state's own tables
(late_window_model.dat surface_p*_Pa, pressure_comparison.txt p*_tuned_Pa; 10 significant digits) as nonneg_contrib does.

DEFINITIONS
  data cost   0.5 sum_g (1/N_g) sum_{used rows} ln(p_model / p_obs)^2 (the fitter's data cost; no crossing, no
              regularization), used = fit_restart.npz use_ng / use_rg, p_obs = the base's COPS_model.dat refNG / refRG
  rms ln      sqrt(mean ln(p_model / p_obs)^2) over the used rows, all / before the gap (et < t_end) / late window
              (et >= t_end = 2016-09-30T10:33:10) (nonneg_finetune's late_window_rms definition)
  Q_H2O       sum_e F_H2O,e area_e at 10:33:10 UTC (geometry label t11 of the state's validated nonneg_maps geometry),
              F = q_H2O A (nonneg_maps.source_rates with the given ratio), all facets and the seen facets (the state's
              seen mask); night / shadow = cos SZA <= 0 or cast shadow (nonneg_maps' rule); shares of Q and of the seen
              production from those facets; also the share of the modelled late-window surface pressure from night /
              shadow records (mean over the used late rows, per gauge)
  --base-run  instead of a state: fit_cops.main's consistency check of the base (the run's own activity elements.bin
              AH2O, ACO2, T0 = the header's analytical temperatures, kappa = 1, FP64 records, instantaneous rows) vs the
              base's COPS_model.dat pNG / pRG written by the C++ (Source.cpp) when the base was built: max rel diff for
              both laws (fit_cops refuses to fit above --check-tolerance, default 1e-6)
  uniform A scale (v2, 2026-11-18; a post-hoc bound, NOT a refit): the single factor s on the H2O activity (A -> s A;
              CO2, T, kappa, grain, fine-tuning factors unchanged) that minimizes the data cost with each law
              (scipy minimize_scalar on ln s, bounded [-6, 6]); tells how much of the change is a global level
  rows file   (v2) <out>.rows.npz: per base row and law the surface NG / RG, their night / shadow parts, the grain, the
              model, p_obs, use flags, et, UTC, iPoint

usage (cwd $APP, .venv/bin/python):
  nightfix_diag.py --key KEY --ft-dir FT --fit-dir FIT --base BASE --geometry GEOM --contrib-code CONTRIB_SNAPSHOT
                   [--adapter localnn_tools.py|d2_tools.py] [--recorded-cost X] --out OUT.json [--procs N]
  nightfix_diag.py --base-run --base BASE --contrib-code CONTRIB_SNAPSHOT --out OUT.json [--procs N]
"""

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import re
import socket
import sys
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
if str(APP) not in sys.path:
    sys.path.insert(0, str(APP))
RATIOS = (("old", 0.2), ("new", 0.02))
LATE_UTC = "2016-09-30T10:33:10"
T0 = time.time()
NM = NC = None
_W = {}


def say(t):
    print("nightfix_diag [%6.0f s]: %s" % (time.time() - T0, t), flush=True)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_code(cdir, adapter):
    """nonneg_maps / nonneg_contrib of the contribution snapshot (as localnn_tools.load_frozen), optional state adapter."""
    global NM, NC
    cdir = Path(cdir)
    if str(cdir) not in sys.path:
        sys.path.insert(0, str(cdir))
    spec = importlib.util.spec_from_file_location("nonneg_maps", str(cdir / "nonneg_maps.py"))
    NM = importlib.util.module_from_spec(spec)
    sys.modules["nonneg_maps"] = NM
    spec.loader.exec_module(NM)
    spec = importlib.util.spec_from_file_location("nonneg_contrib", str(cdir / "nonneg_contrib.py"))
    NC = importlib.util.module_from_spec(spec)
    sys.modules["nonneg_contrib"] = NC
    spec.loader.exec_module(NC)
    if NC.NM is not NM or not hasattr(NM, "H2O_QMIN_RATIO"):
        raise SystemExit("contribution code %s: not a snapshot with nonneg_maps.H2O_QMIN_RATIO" % cdir)
    info = {"contrib_code": str(cdir.resolve()), "nonneg_maps.py": sha256_file(cdir / "nonneg_maps.py"),
            "nonneg_contrib.py": sha256_file(cdir / "nonneg_contrib.py"), "H2O_QMIN_RATIO_of_code": NM.H2O_QMIN_RATIO}
    if adapter:
        ap = Path(adapter).resolve()
        if str(ap.parent) not in sys.path:
            sys.path.insert(1, str(ap.parent))
        spec = importlib.util.spec_from_file_location("state_adapter", str(ap))
        ad = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ad)
        NM.state_activity = ad.make_state_activity(NM)
        info["adapter"] = {"path": str(ap), "sha256": sha256_file(ap)}
    return info


def row_work(task):
    """One base row: per law the surface NG / RG (sum of records) and the part from night / shadow records."""
    row, ip = task
    G = NC._G
    kp, h, rec = NM.read_kernel_full(Path(G["base"]) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
    prec = _W["prec"]
    if prec == "fp32":
        f32 = G["fp32_record_fields"](rec, G["COS75"], {})
        night = f32["cos"].astype(np.float64) <= 0.0
    else:
        night = (np.asarray(rec["cos"], np.float64) <= 0.0) | (rec["shadow"] != 0)
    out = {"row": row, "ip": ip, "records": int(len(rec)), "et": h["et"]}
    for lab, r in RATIOS:
        NM.H2O_QMIN_RATIO = r
        _, ((ng, rg),), _ = NC.record_contributions(rec, h, (G["A"],), prec, row)
        out[lab] = (float(np.sum(ng)), float(np.sum(rg)), float(np.sum(ng[night])), float(np.sum(rg[night])))
    NM.H2O_QMIN_RATIO = _W["code_ratio"]
    return out


def rows_pass(tasks, procs):
    R = {}
    with mp.get_context("fork").Pool(procs) as pool:
        for n, res in enumerate(pool.imap_unordered(row_work, tasks, chunksize=1), start=1):
            R[res["row"]] = res
            if n % 200 == 0 or n == len(tasks):
                say("  %d / %d rows" % (n, len(tasks)))
    return R


def read_table(path, keys):
    names, out = None, {}
    for line in open(path):
        if line.startswith("#iPoint"):
            names = line[1:].split()
            continue
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        c = {n: i for i, n in enumerate(names)}
        out[int(t[c["iPoint"]])] = {k: (t[c[k]] if k in ("UTC", "kind") else float(t[c[k]])) for k in keys}
    return out


def data_cost(model, ref, use):
    cost = 0.0
    for g in ("NG", "RG"):
        u = use[g] & (ref[g] > 0)
        lr = np.log(model[g][u] / ref[g][u])
        cost += 0.5 * float(np.sum(lr ** 2)) / int(u.sum())
    return cost


def uniform_scale(surf, grain, ref, use):
    """v2: the factor s on the H2O activity minimizing the data cost (the CO2 part is kappa ~ 0 in every state and is
    scaled with it here; reported). Returns s, the cost at s and at s = 1."""
    from scipy.optimize import minimize_scalar
    f = lambda ls: data_cost({"NG": np.exp(ls) * surf["NG"] + grain, "RG": np.exp(ls) * surf["RG"]}, ref, use)  # noqa
    r = minimize_scalar(f, bounds=(-6.0, 6.0), method="bounded", options={"xatol": 1e-10})
    return {"s": float(np.exp(r.x)), "ln_s": float(r.x), "data_cost_at_s": float(r.fun), "data_cost_at_1": f(0.0),
            "note": "post-hoc bound, not a refit: A -> s A for every facet (surface model x s; grain unchanged)"}


def cost_stats(model, ref, use, et, t_end):
    out, cost = {}, 0.0
    for g in ("NG", "RG"):
        u = use[g] & (ref[g] > 0)
        lr = np.log(model[g][u] / ref[g][u])
        cost += 0.5 * float(np.sum(lr ** 2)) / int(u.sum())
        late = et[u] >= t_end
        out["mean_ln_%s" % g] = float(np.mean(lr))
        out["rms_ln_%s" % g] = float(np.sqrt(np.mean(lr ** 2)))
        out["rms_ln_%s_before_gap" % g] = float(np.sqrt(np.mean(lr[~late] ** 2)))
        out["rms_ln_%s_late_window" % g] = float(np.sqrt(np.mean(lr[late] ** 2)))
        out["mean_ln_%s_late_window" % g] = float(np.mean(lr[late]))
        out["used_rows_%s" % g], out["late_rows_%s" % g] = int(u.sum()), int(late.sum())
    out["data_cost"] = cost
    return out


def main_base_run(a, info):
    base = Path(a.base)
    fc, header, cols, rows = NM.base_table(base)
    el = NM.read_elements(base / "elements.bin")
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    T = tuple(float(fc.header_value(header, "Analytical" + k)) for k in fc.TNAMES)
    bm, cNG, cRG, gconst = NC.gauge_constants(fc, header, TRG)
    from cops_cuda import device_store
    NC._G.update({"base": str(base), "A": el["AH2O"], "ACO2": el["ACO2"], "T": T, "kappa": 1.0, "bm": bm, "cNG": cNG,
                  "cRG": cRG, "helio": fc.header_value(header, "HeliocentricModel"), "AU": float(fc.AU),
                  "COS75": float(fc.COS75), "KT": float(fc.KT), "fp32_record_fields": device_store.fp32_record_fields})
    _W.update({"prec": "fp64", "code_ratio": NM.H2O_QMIN_RATIO})
    tasks = sorted([(r, int(ipoint[r])) for r in range(len(ipoint))],
                   key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    say("base run: %d rows, FP64 records, T0 %s, kappa 1" % (len(tasks), T))
    R = rows_pass(tasks, a.procs)
    tab = {"NG": np.array([float(r[cols["pNG"]]) for r in rows]), "RG": np.array([float(r[cols["pRG"]]) for r in rows])}
    ok = (tab["NG"] > 0) & (tab["RG"] > 0)
    res = {"kind": "base_run", "base": str(base.resolve()), "rows": len(rows), "rows_compared": int(ok.sum()), "T0": T,
           "gauge_constants": gconst, "code": info, "laws": {}}
    for lab, r in RATIOS:
        m = {"NG": np.array([R[k][lab][0] for k in range(len(rows))]), "RG": np.array([R[k][lab][1] for k in range(len(rows))])}
        d = {g: np.abs(m[g][ok] / tab[g][ok] - 1.0) for g in ("NG", "RG")}
        res["laws"][lab] = {"ratio": r, "max_rel_NG": float(d["NG"].max()), "max_rel_RG": float(d["RG"].max()),
                            "median_rel_NG": float(np.median(d["NG"])), "median_rel_RG": float(np.median(d["RG"])),
                            "night_share_NG_median": float(np.median([R[k][lab][2] / R[k][lab][0] for k in range(len(rows))
                                                                      if R[k][lab][0] > 0])),
                            "passes_check_tolerance_1e-6": bool(max(d["NG"].max(), d["RG"].max()) <= 1e-6)}
        say("%s law (r = %g): consistency max rel NG %.3e RG %.3e" % (lab, r, d["NG"].max(), d["RG"].max()))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--key")
    ap.add_argument("--ft-dir")
    ap.add_argument("--fit-dir")
    ap.add_argument("--base", required=True)
    ap.add_argument("--geometry")
    ap.add_argument("--contrib-code", required=True)
    ap.add_argument("--adapter")
    ap.add_argument("--recorded-cost", type=float)
    ap.add_argument("--base-run", action="store_true")
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists():
        raise SystemExit("%s exists (never overwritten)" % out)
    info = load_code(a.contrib_code, a.adapter)
    head = {"schema": "nightfix_diag_v2", "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
            "script": str(Path(__file__).resolve()), "script_sha256": sha256_file(__file__), "procs": a.procs,
            "ratios": dict(RATIOS)}
    if a.base_run:
        res = main_base_run(a, info)
        res.update(head)
        res["wall_seconds"] = time.time() - T0
        json.dump(res, open(out, "x"), indent=1, default=float)
        say("done: %s" % out)
        return 0
    base, ft_dir, fit_dir = Path(a.base), Path(a.ft_dir), Path(a.fit_dir)
    gval = json.load(open(Path(a.geometry) / "geometry_validation.json"))
    if not gval.get("passed") or Path(gval["base"]) != base.resolve():
        raise SystemExit("geometry %s not validated against %s" % (a.geometry, base))
    fc, header, cols, rows = NM.base_table(base)
    law = NM.law_guard(fc)                                  # the fitter's law = the new constant (0.02)
    helio, sza = fc.header_value(header, "HeliocentricModel"), fc.header_value(header, "SZAMode")
    if sza != "continuous":
        raise SystemExit("SZAMode %s" % sza)
    el = NM.read_elements(base / "elements.bin")
    nf = el["nfaces"]
    S = NM.state_activity(ft_dir, fit_dir, el, a.procs)
    if not S["checks"]["passed"]:
        raise SystemExit("state checks failed")
    with np.load(fit_dir / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
        point_ids = np.asarray(z["point_ids"], np.int64)
        use = {"NG": np.asarray(z["use_ng"]) > 0, "RG": np.asarray(z["use_rg"]) > 0}
    with np.load(ft_dir / "finetune_state.npz", allow_pickle=False) as z:
        gA, gB = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    t_end = float(json.load(open(fit_dir / "grain_term.json"))["t_end_et"])
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("restart rows differ from the base rows")
    ref = {"NG": np.array([float(r[cols["refNG"]]) for r in rows]), "RG": np.array([float(r[cols["refRG"]]) for r in rows])}
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    et = np.array([float(r[cols["et"]]) for r in rows])
    utc = [r[cols["UTC"]] for r in rows]
    bm, cNG, cRG, gconst = NC.gauge_constants(fc, header, TRG)
    ftj = S["ft_json"]
    pprec = ftj["modules"]["positive_problem.py"]
    pp = NM.load_module(pprec["path"], "positive_problem_of_state", pprec["sha256"])
    ptab = ftj["fit"]["point_table"]
    if NM.sha256_file(ptab["path"]) != ptab["sha256"]:
        raise SystemExit("point table %s changed" % ptab["path"])
    layout = pp.LateLayout(pp.read_point_table(ptab["path"]), pp.read_run_points(base), ipoint)
    W = layout.W.tocsr()
    W.sort_indices()
    late_rows, gr = layout.late_rows, layout.grid_rows
    if len(late_rows) != 95:
        raise SystemExit("%d late observations, expected 95" % len(late_rows))
    from cops_cuda import device_store
    NC._G.update({"base": str(base), "A": S["A"], "ACO2": el["ACO2"], "T": T, "kappa": S["kappa"], "bm": bm, "cNG": cNG,
                  "cRG": cRG, "helio": helio, "AU": float(fc.AU), "COS75": float(fc.COS75), "KT": float(fc.KT),
                  "fp32_record_fields": device_store.fp32_record_fields, "nfaces": nf})
    _W.update({"prec": "fp32", "code_ratio": NM.H2O_QMIN_RATIO})
    tasks = sorted([(r, int(ipoint[r])) for r in range(len(ipoint))],
                   key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    say("%s: %d rows, %d procs; A from %s (model %s), kappa %.4g, T %s, grain A %.6g uPa" % (
        a.key, len(tasks), a.procs, ft_dir, S["model"], S["kappa"], T, gA))
    R = rows_pass(tasks, a.procs)
    nrow = len(ipoint)
    late_obs = np.zeros(nrow, bool)
    late_obs[late_rows] = True
    late_et = layout.late_et
    grain = np.zeros(nrow)
    grain[late_rows] = np.maximum(0.0, gA * 1e-6 - gB * 1e-8 * (late_et - t_end)) * (late_et >= t_end)
    lw = read_table(ft_dir / "late_window_model.dat", ("UTC", "kind", "surface_pNG_Pa", "surface_pRG_Pa", "grain_NG_Pa"))
    pc = read_table(ft_dir / "pressure_comparison.txt", ("UTC", "kind", "pNG_tuned_Pa", "pRG_tuned_Pa"))
    rows_out = {"iPoint": ipoint, "et": et, "utc": np.array(utc), "ref_NG": ref["NG"], "ref_RG": ref["RG"],
                "use_NG": use["NG"], "use_RG": use["RG"], "grain_NG": grain, "late_obs": late_obs}
    res = {"kind": "state", "key": a.key, "model": S["model"], "ft_dir": str(ft_dir), "fit_dir": str(fit_dir),
           "ft_dir_resolved": str(ft_dir.resolve()), "base": str(base.resolve()), "geometry": str(Path(a.geometry).resolve()),
           "code": info, "law_guard": law["checked"], "fit_cops_sha256": law["fit_cops_sha256"], "T": T,
           "kappa": S["kappa"], "grain_A_uPa": gA, "grain_B_uPa_per_100s": gB, "t_end_et": t_end,
           "state_checks_passed": True, "records": int(sum(R[r]["records"] for r in R)), "laws": {}}
    for lab, r in RATIOS:
        tot = {g: np.array([R[k][lab][i] for k in range(nrow)]) for i, g in ((0, "NG"), (1, "RG"))}
        dark = {g: np.array([R[k][lab][i] for k in range(nrow)]) for i, g in ((2, "NG"), (3, "RG"))}
        surf = {g: tot[g].copy() for g in tot}
        sdark = {g: dark[g].copy() for g in dark}
        for g in tot:
            surf[g][late_rows] = W @ tot[g][gr]
            sdark[g][late_rows] = W @ dark[g][gr]
        model = {"NG": surf["NG"] + grain, "RG": surf["RG"]}
        st = cost_stats(model, ref, use, et, t_end)
        st["uniform_A_scale"] = uniform_scale(surf, grain, ref, use)          # v2
        for g in ("NG", "RG"):
            rows_out["%s_surface_%s" % (lab, g)], rows_out["%s_night_shadow_%s" % (lab, g)] = surf[g], sdark[g]
            rows_out["%s_model_%s" % (lab, g)] = model[g]
        lu = {g: use[g] & (et >= t_end) & (ref[g] > 0) for g in ("NG", "RG")}
        st["late_window_night_shadow_share_of_surface_mean"] = {g: float(np.mean(sdark[g][lu[g]] / surf[g][lu[g]]))
                                                                for g in ("NG", "RG")}
        st["all_used_night_shadow_share_of_surface_mean"] = {
            g: float(np.mean(sdark[g][use[g] & (ref[g] > 0)] / surf[g][use[g] & (ref[g] > 0)])) for g in ("NG", "RG")}
        st["late_window_grain_share_of_NG_mean"] = float(np.mean(grain[lu["NG"]] / model["NG"][lu["NG"]]))
        if lab == "old":
            # verification against the state's own tables (10 significant digits; nonneg_contrib's criterion)
            tabs = {"NG": np.empty(nrow), "RG": np.empty(nrow)}
            for k in range(nrow):
                ip = int(ipoint[k])
                if ip in lw:
                    tabs["NG"][k], tabs["RG"][k] = lw[ip]["surface_pNG_Pa"], lw[ip]["surface_pRG_Pa"]
                else:
                    tabs["NG"][k], tabs["RG"][k] = pc[ip]["pNG_tuned_Pa"], pc[ip]["pRG_tuned_Pa"]
            ver = {g: float(np.max(np.abs(surf[g] / tabs[g] - 1.0))) for g in ("NG", "RG")}
            gt = np.array([lw[int(i)]["grain_NG_Pa"] for i in ipoint[late_rows]])
            with np.errstate(divide="ignore", invalid="ignore"):
                ver["grain_vs_table_max_rel"] = float(np.nanmax(np.where(gt > 0, np.abs(grain[late_rows] / gt - 1.0),
                                                                         np.abs(grain[late_rows] - gt))))
            ver["passed"] = bool(ver["NG"] <= 5.0e-10 * (1 + 1e-6) + 1e-14 and ver["RG"] <= 5.0e-10 * (1 + 1e-6) + 1e-14)
            res["verification_old_law_vs_tables"] = ver
            if a.recorded_cost is not None:
                res["recorded_data_cost"] = a.recorded_cost
                res["old_law_cost_vs_recorded_rel"] = abs(st["data_cost"] / a.recorded_cost - 1.0)
            say("old law vs tables: %s; data cost %.10e (recorded %s)" % (json.dumps(ver), st["data_cost"],
                                                                           a.recorded_cost))
        res["laws"][lab] = {"ratio": r, **st}
        say("%s law (r = %g): data cost %.6e, rms ln NG %.5f RG %.5f, late NG %.5f RG %.5f" % (
            lab, r, st["data_cost"], st["rms_ln_NG"], st["rms_ln_RG"], st["rms_ln_NG_late_window"],
            st["rms_ln_RG_late_window"]))
    # production at 10:33:10
    G_ = NM.read_geometry(a.geometry)
    Tt = G_["times"]["t11"]
    if Tt["utc"] != LATE_UTC:
        raise SystemExit("geometry t11 is %s" % Tt["utc"])
    cos, sh = NM.geometry_time(a.geometry, "t11", nf)
    area, seen = el["area"], S["seen"].astype(bool)
    for lab, r in RATIOS:
        NM.H2O_QMIN_RATIO = r
        F_h, F_c, night, qf = NM.source_rates(cos, sh, Tt["rh_m"], helio, sza, S["A"], el["ACO2"], S["kappa"], fc.AU)
        P = F_h * area
        Q, Qs = float(P.sum()), float(P[seen].sum())
        res["laws"][lab]["production_103310"] = {
            "utc": Tt["utc"], "q_min_over_q_max": float(qf["H2O"][0] / qf["H2O"][1]), "Q_H2O_s-1": Q,
            "Q_H2O_seen_s-1": Qs, "seen_share_of_Q_H2O": Qs / Q,
            "night_shadow_share_of_Q_H2O": float(P[night].sum() / Q),
            "night_shadow_share_of_seen_production": float(P[night & seen].sum() / Qs),
            "Q_CO2_s-1": float(np.sum(F_c * area)), "night_shadow_area_fraction": float(area[night].sum() / area.sum()),
            "night_shadow_area_fraction_seen": float(area[night & seen].sum() / area[seen].sum())}
    NM.H2O_QMIN_RATIO = _W["code_ratio"]
    o, n = res["laws"]["old"], res["laws"]["new"]
    res["new_over_old"] = {"data_cost": n["data_cost"] / o["data_cost"],
                           "Q_H2O": n["production_103310"]["Q_H2O_s-1"] / o["production_103310"]["Q_H2O_s-1"],
                           "Q_H2O_seen": n["production_103310"]["Q_H2O_seen_s-1"] / o["production_103310"]["Q_H2O_seen_s-1"]}
    res.update(head)
    res["wall_seconds"] = time.time() - T0
    np.savez(str(out) + ".rows.npz", **rows_out)                            # v2
    res["rows_file"] = str(out) + ".rows.npz"
    json.dump(res, open(out, "x"), indent=1, default=float)
    say("done: %s (new/old data cost %.4f; uniform A scale new law s = %.4f -> cost %.4e)" % (
        out, res["new_over_old"]["data_cost"], n["uniform_A_scale"]["s"], n["uniform_A_scale"]["data_cost_at_s"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
