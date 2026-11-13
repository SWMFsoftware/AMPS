#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_contrib.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_contrib.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 7c75c5abf1779d2470d71981b696fafc7f166073f66f0028a2e62bf0ef898895
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: analysis module. Per kernel record, the contribution of every facet to the modelled NG and RG surface
#   signal of a chain's final state, as the fit computes it (FP32 record fields, FP64 arithmetic, F = q A with the
#   night/shadow rule, the fit's T and kappa, W10x5 for the 95 late observations); P% source areas S_P (after Schmid
#   1994), source-area occurrence N_P, slope-bin fractions. Frozen nonneg_contrib_20261112T210012Z.
# RUN: jobs/nonneg_contrib.sbatch (CPU, normal128cores; CHAIN, CYCLE, GEOM, MODEL_EXT[, CONTRIB_NAME]); jobs 124549,
#   124668 (bounded c5), 124776 (c9).
# ==== end of header ====
"""Per-facet contributions of the fitted model to the COPS NG and RG surface pressures, the P% source areas of the
post-gap (late) observations, and the slope-bin fractions of every fit point, for ONE final state of the non-negative
H2O activity models (model_ext/nonneg_*.py).

User request (2026-11-12, relayed by the lead): "for each of the observations after the gap determine the area that
contributed more than 50% to the RG measurement; start counting from the surface elements that have higher
contribution to the observations; combine all those plots into a single one ...; then repeat the same for 60, 70, 80,
90%; then create another map with this data where the colour would indicate how many times the surface element was in
the most important area ...; repeat the same for NG and find an appropriate name for this variable"; "as the fitting
shows importance of the cliffs, also add a plot that illustrates that". This tool computes the DATA only (no plots).

Local, untracked model extension: nothing in the git tree is edited; the reference fitter (fit_cops) and the GPU
application's FP32 record conversion (cops_cuda.device_store.fp32_record_fields) are imported unchanged.

CONTRIBUTION OF A FACET (exactly the fitted forward model; the fits ran --precision mixed = cops_cuda
  fp32_records_fp64_arithmetic_v1, so Kn, KPhi, KV, cos are the FP32 values of fp32_record_fields, shadow folded into
  cos = -2, everything else FP64; cops_cuda/kernels/pressure.cu static_factors_m / fused_body_m = fit_cops.Model
  load_kernels / prepare / evaluate). Per kernel record r (point p, facet e):
    night = cos32 <= 0; cz = min(cos32, 1); q_s = qmin_s if night else qmin_s + (qmax_s - qmin_s) cz (r_h of point p)
    wnH = qH Kn, wvH = qH VdotB KV, wpH = qH KPhi; qC = q_CO2 ACO2_e, wnC = qC Kn, wvC = qC VdotB KV, phiC = qC KPhi
    Te = max(Tmin, Tmax + KT (Tmax - T75)(1/c - 1)) where lit (c >= cos 75 deg), else Tmin; u = Te^-1/2
    NG_r = cNG_H bm_H wnH A_e u + kappa cNG_C bm_C wnC u
    RG_r = cRG_H,p (wpH A_e + bm_H wvH A_e u) + kappa cRG_C,p (phiC + bm_C wvC u)
  A_e = the final state's activity (finetune_state.npz A_final); T and kappa = the fit's (held in the fine tuning).
  Model surface pressure of a point = sum_r. A late observation (>= 10:33:10 UTC) is modelled by W10x5 of the 1-s grid
  rows (positive_problem.LateLayout, the chain's own copy): its facet contribution is sum_j W_kj c_j,e. The NG grain
  term (grain(t) = max(0, A - B (t - t_end)), unaveraged) is not a facet contribution; its share is reported.

P% SOURCE AREA (proposed name; after the footprint / source-area definition of Schmid 1994): for an observation k and a
  gauge g, S_P(k, g) = the smallest set of facets, taken in decreasing order of their contribution c_k,g,e (ties: lower
  facet index first), whose cumulative contribution reaches P % of the observation's total surface contribution
  sum_e c_k,g,e. SOURCE-AREA OCCURRENCE N_P,g(e) = the number of the 95 late observations whose S_P(k, g) contains e.

Subcommand-free CLI: nonneg_contrib.py --ft-dir c5/ft/final --fit-dir c5/fit/final --base BASE --geometry GEOM
  --out <chain>/contributions [--procs N]
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
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
for _p in (str(HERE), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import nonneg_maps as NM            # noqa: E402  (readers, state checks, law guard, source law, geometry)

SCHEMA = "nonneg_contrib_v1"
P_LEVELS = (50, 60, 70, 80, 90)
GAUGES = ("NG", "RG")
NBINS = 7
_G = {}                              # worker globals (set before the pools fork)


def say(text):
    print("nonneg_contrib: [%s] %s" % (dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


# ------------------------------------------------------------------------------------------------ constants
def gauge_constants(fc, header, TRG):
    """bm_s, cNG_s (scalars) and cRG_s per row: fit_cops.Model.prepare / crg (the functions themselves)."""
    TNG = float(fc.header_value(header, "TemperatureNG"))
    S = fc.sensitivities(header)
    bm = {s: np.sqrt(fc.MASS[s] / (2.0 * fc.KB)) * fc.SQRTPI / 2.0 for s in fc.SPECIES}
    cNG = {s: fc.KB * TNG / Ss for s, Ss in zip(fc.SPECIES, S)}
    ns = types.SimpleNamespace(TRG=np.asarray(TRG, float), rows=np.arange(len(TRG)),
                               RGModel=fc.header_value(header, "RGPressureModel"))
    cRG = {s: fc.Model.crg(ns, s) / Ss for s, Ss in zip(fc.SPECIES, S)}
    return bm, cNG, cRG, {"TemperatureNG": TNG, "sensitivities": list(S), "RGPressureModel": ns.RGModel}


def record_contributions(rec, h, A_list, prec, row):
    """Per-record NG and RG contributions [Pa] of one point for each activity in A_list (fp32: the fit's mixed-precision
    records; fp64: the file's values, diagnostic). Returns (face int64, [(ng, rg) per A])."""
    G = _G
    face = rec["face"].astype(np.int64)
    if prec == "fp32":
        conv = {}
        f32 = G["fp32_record_fields"](rec, G["COS75"], conv)
        kn, kphi, kv, c0 = (f32[k].astype(np.float64) for k in ("kn", "kphi", "kv", "cos"))
        night = c0 <= 0.0
    else:
        conv = None
        kn, kphi, kv, c0 = (np.asarray(rec[k], np.float64) for k in ("Kn", "KPhi", "KV", "cos"))
        night = (c0 <= 0.0) | (rec["shadow"] != 0)
    qf = NM.q_factors(h["rh"], G["helio"], G["AU"])
    (qmin_h, qmax_h), (qmin_c, qmax_c) = qf["H2O"], qf["CO2"]
    vb = h["VdotB"]
    cz = np.where(c0 < 1.0, c0, 1.0)
    qh = np.where(night, qmin_h, qmin_h + (qmax_h - qmin_h) * cz)
    qc = np.where(night, qmin_c, qmin_c + (qmax_c - qmin_c) * cz) * G["ACO2"][face]
    wnH, wvH, wpH = qh * kn, qh * vb * kv, qh * kphi
    wnC, wvC, phiC = qc * kn, qc * vb * kv, qc * kphi
    clipped = np.clip(c0, -1.0, 1.0)
    c = np.where(night, -1.0, clipped)
    lit = c >= G["COS75"]
    x = np.where(lit, 1.0 / np.where(lit, c, 1.0) - 1.0, 0.0)
    Tmax, T75, Tmin = G["T"]
    tl = Tmax + G["KT"] * (Tmax - T75) * x
    te = np.where(lit & (tl > Tmin), tl, Tmin)
    u = 1.0 / np.sqrt(te)
    kappa, bm, cNG = G["kappa"], G["bm"], G["cNG"]
    cRGH, cRGC = G["cRG"]["H2O"][row], G["cRG"]["CO2"][row]
    out = []
    for A in A_list:
        af = A[face]
        ng = cNG["H2O"] * bm["H2O"] * (wnH * af * u) + (kappa * cNG["CO2"] * bm["CO2"]) * (wnC * u)
        rg = cRGH * (wpH * af + bm["H2O"] * (wvH * af * u)) + kappa * cRGC * (phiC + bm["CO2"] * (wvC * u))
        out.append((ng, rg))
    return face, out, conv


def row_work(task):
    """One base row: totals (final A fp32 / fp64 records, fit A fp32), slope-bin sums, negatives; per-record facet
    contributions for the W10x5 grid rows."""
    row, ip, want = task
    G = _G
    kp, h, rec = NM.read_kernel_full(Path(G["base"]) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
    face, ((ng, rg), (ng_fit, rg_fit)), conv = record_contributions(rec, h, (G["A"], G["A_fit"]), "fp32", row)
    _, ((ng64, rg64),), _ = record_contributions(rec, h, (G["A"],), "fp64", row)
    b = G["bins"][face]
    res = {"row": row, "ip": ip, "records": int(len(rec)), "et": h["et"], "rh": h["rh"], "VdotB": h["VdotB"],
           "ng": float(np.sum(ng)), "rg": float(np.sum(rg)), "ng_fp64": float(np.sum(ng64)), "rg_fp64": float(np.sum(rg64)),
           "ng_fit": float(np.sum(ng_fit)), "rg_fit": float(np.sum(rg_fit)),
           "bins_ng": np.bincount(b, weights=ng, minlength=NBINS), "bins_rg": np.bincount(b, weights=rg, minlength=NBINS),
           "neg_ng": int(np.sum(ng < 0)), "neg_rg": int(np.sum(rg < 0)),
           "nonfinite": int(np.sum(~np.isfinite(ng)) + np.sum(~np.isfinite(rg))),
           "lit_threshold_adjusted": int((conv or {}).get("lit_threshold_adjusted", 0))}
    if want:
        res["face"] = face.astype(np.int32)
        res["c_ng"], res["c_rg"] = ng, rg
    return res


def obs_work(k):
    """One late observation: dense facet contributions sum_j W_kj c_j (grid rows in ascending order), the ranking and
    the P% source areas of both gauges."""
    G = _G
    W = G["W"]
    lo, hi = W.indptr[k], W.indptr[k + 1]
    cols, wts = W.indices[lo:hi], W.data[lo:hi]
    o = np.argsort(cols, kind="stable")
    cols, wts = cols[o], wts[o]
    nf = G["nfaces"]
    C = {"NG": np.zeros(nf), "RG": np.zeros(nf)}
    touched = np.zeros(nf, dtype=bool)
    wsum_ng = wsum_rg = 0.0
    for j, w in zip(cols.tolist(), wts.tolist()):
        F = G["FAC"][int(G["grid_rows"][j])]
        face = F["face"]
        C["NG"][face] += w * F["c_ng"]
        C["RG"][face] += w * F["c_rg"]
        touched[face] = True
        wsum_ng += w * F["ng"]
        wsum_rg += w * F["rg"]
    idx = np.flatnonzero(touched)
    out = {"k": k, "W_columns": int(len(cols)), "W_rowsum": float(np.sum(wts)), "W_rows_total": {"NG": wsum_ng,
                                                                                              "RG": wsum_rg}}
    area = G["area"]
    for g in GAUGES:
        v = C[g][idx]
        neg = int(np.sum(v < 0))
        pos = v > 0
        f, vv = idx[pos], v[pos]
        order = np.lexsort((f, -vv))
        fs, vs = f[order].astype(np.int32), vv[order]
        cs = np.cumsum(vs)
        tot = float(cs[-1]) if len(cs) else 0.0
        kp = [int(np.searchsorted(cs, p / 100.0 * tot, side="left")) + 1 for p in P_LEVELS]
        out[g] = {"total_face_order": float(np.sum(v)), "total_sorted": tot, "negative": neg,
                  "negative_sum": float(np.sum(v[v < 0])), "n_contrib": int(len(fs)),
                  "area_contrib": float(np.sum(area[fs])), "k": kp,
                  "area": [float(np.sum(area[fs[:n]])) for n in kp],
                  "c_min": [float(vs[n - 1]) for n in kp],
                  "cum_at_k": [float(cs[n - 1] / tot) for n in kp],
                  "seen_area": [float(np.sum(area[fs[:n]][G["seen"][fs[:n]]])) for n in kp],
                  "faces": fs, "contrib": vs}
    return out


# ------------------------------------------------------------------------------------------------ tables
def read_table(path, key_cols):
    names, rows = None, {}
    for line in open(path):
        if line.startswith("#iPoint"):
            names = line[1:].split()
            continue
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        c = {n: i for i, n in enumerate(names)}
        rows[int(t[c["iPoint"]])] = {n: (t[c[n]] if n in ("UTC", "kind") else float(t[c[n]])) for n in key_cols}
    return rows


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a / b - 1.0)


# ------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ft-dir", required=True)
    ap.add_argument("--fit-dir", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--geometry", required=True, help="a validated nonneg_maps geometry (labels t11, t12)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--expect-model", choices=("softplus", "squared"))
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    args = ap.parse_args(argv)
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists (never overwritten)" % out)
    stage = out.parent / ("%s.incomplete_%s" % (out.name, dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    base, ft_dir, fit_dir = Path(args.base), Path(args.ft_dir), Path(args.fit_dir)
    gval = json.load(open(Path(args.geometry) / "geometry_validation.json"))
    if not gval.get("passed") or Path(gval["base"]) != base.resolve():
        raise SystemExit("geometry %s not validated against %s" % (args.geometry, base))
    fc, header, cols, rows = NM.base_table(base)
    law = NM.law_guard(fc)
    helio, sza = fc.header_value(header, "HeliocentricModel"), fc.header_value(header, "SZAMode")
    if sza != "continuous":
        raise SystemExit("SZAMode %s" % sza)
    mrun = re.search(r"^# run: (\S+)$", (fit_dir / "fit_result.txt").read_text(), re.M)
    if mrun is None or Path(mrun.group(1)).resolve() != base.resolve():
        raise SystemExit("the fit did not run on %s" % base)
    el = NM.read_elements(base / "elements.bin")
    nf = el["nfaces"]
    S = NM.state_activity(ft_dir, fit_dir, el, args.procs)
    if not S["checks"]["passed"]:
        raise SystemExit("state checks failed: %s" % json.dumps({k: v for k, v in S["checks"].items()
                                                                  if k not in ("nonneg_mult", "slope_mult")}))
    if args.expect_model and S["model"] != args.expect_model:
        raise SystemExit("model %s, expected %s" % (S["model"], args.expect_model))
    with np.load(fit_dir / "nonneg_facets.npz", allow_pickle=False) as z:
        A_fit = np.asarray(z["A"], np.float64)
    with np.load(fit_dir / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
        point_ids = np.asarray(z["point_ids"], np.int64)
        use_ng, use_rg = np.asarray(z["use_ng"]), np.asarray(z["use_rg"])
        tv = np.asarray(z["template_values"], float)
    with np.load(ft_dir / "finetune_state.npz", allow_pickle=False) as z:
        grain_A_uPa, grain_B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("restart rows differ from the base rows")
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    utc = [r[cols["UTC"]] for r in rows]
    bm, cNG, cRG, gconst = gauge_constants(fc, header, TRG)
    # the late layout of the chain's own positive_problem
    ftj = S["ft_json"]
    pprec = ftj["modules"]["positive_problem.py"]
    pp = NM.load_module(pprec["path"], "positive_problem_of_state", pprec["sha256"])
    ptab = ftj["fit"]["point_table"]
    if NM.sha256_file(ptab["path"]) != ptab["sha256"]:
        raise SystemExit("point table %s changed" % ptab["path"])
    layout = pp.LateLayout(pp.read_point_table(ptab["path"]), pp.read_run_points(base), ipoint)
    W = layout.W.tocsr()
    W.sort_indices()
    nlate = len(layout.late_rows)
    say("late layout: %s" % json.dumps(layout.summary, sort_keys=True))
    if nlate != 95:
        raise SystemExit("%d late observations, expected 95" % nlate)
    if np.any(W.data < 0):
        raise SystemExit("W10x5 has negative weights")
    from cops_cuda import device_store
    _G.update({"base": str(base), "A": S["A"], "A_fit": A_fit, "ACO2": el["ACO2"], "bins": S["bins"].astype(np.intp),
               "T": T, "kappa": S["kappa"], "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio, "AU": float(fc.AU),
               "COS75": float(fc.COS75), "KT": float(fc.KT), "fp32_record_fields": device_store.fp32_record_fields,
               "nfaces": nf, "area": el["area"], "seen": S["seen"]})
    grid_rows = set(int(r) for r in layout.grid_rows[np.unique(W.indices)])
    tasks = [(r, int(ipoint[r]), r in grid_rows) for r in range(len(ipoint))]
    tasks.sort(key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))      # largest first (load balance)
    say("phase A: %d rows (%d W10x5 grid rows with facet contributions), %d processes" % (len(tasks), len(grid_rows),
                                                                                       args.procs))
    R, FAC = {}, {}
    with mp.get_context("fork").Pool(args.procs) as pool:
        for n, res in enumerate(pool.imap_unordered(row_work, tasks, chunksize=1), start=1):
            r = res["row"]
            if "face" in res:
                FAC[r] = {k: res.pop(k) for k in ("face", "c_ng", "c_rg")}
                FAC[r]["ng"], FAC[r]["rg"] = res["ng"], res["rg"]
            R[r] = res
            if n % 200 == 0 or n == len(tasks):
                say("  %d / %d rows (%.0f s)" % (n, len(tasks), time.time() - t0))
    nrow = len(ipoint)
    neg = {"NG": sum(R[r]["neg_ng"] for r in R), "RG": sum(R[r]["neg_rg"] for r in R)}
    nonfinite = sum(R[r]["nonfinite"] for r in R)
    say("records %d; negative record contributions NG %d RG %d; non-finite %d; FP32 lit-threshold adjustments %d" % (
        sum(R[r]["records"] for r in R), neg["NG"], neg["RG"], nonfinite, sum(R[r]["lit_threshold_adjusted"] for r in R)))
    if neg["NG"] or neg["RG"] or nonfinite:
        raise SystemExit("negative or non-finite contributions: stop and ask before defining the ranking")
    # per-row arrays
    tot = {k: np.array([R[r][k] for r in range(nrow)]) for k in ("ng", "rg", "ng_fp64", "rg_fp64", "ng_fit", "rg_fit")}
    bins = {"NG": np.array([R[r]["bins_ng"] for r in range(nrow)]), "RG": np.array([R[r]["bins_rg"] for r in range(nrow)])}
    # model value per row: instantaneous; late observations W10x5 of the grid rows
    late_rows, gr = layout.late_rows, layout.grid_rows
    model = {k: v.copy() for k, v in tot.items()}
    for k in model:
        model[k][late_rows] = W @ tot[k][gr]
    mbins = {g: bins[g].copy() for g in GAUGES}
    for g in GAUGES:
        mbins[g][late_rows] = W @ bins[g][gr]
    # verification against the written model tables (10 significant digits)
    lw_ft = read_table(ft_dir / "late_window_model.dat", ("UTC", "kind", "surface_pNG_Pa", "surface_pRG_Pa", "grain_NG_Pa"))
    lw_fit = read_table(fit_dir / "late_window_model.dat", ("UTC", "kind", "surface_pNG_Pa", "surface_pRG_Pa",
                                                            "grain_NG_Pa"))
    pc = read_table(ft_dir / "pressure_comparison.txt", ("UTC", "kind", "pNG_fit_Pa", "pNG_tuned_Pa", "pRG_fit_Pa",
                                                         "pRG_tuned_Pa"))
    ref = {k: np.empty(nrow) for k in ("ng", "rg", "ng_fit", "rg_fit")}
    src = []
    for r in range(nrow):
        ip = int(ipoint[r])
        if ip in lw_ft:
            ref["ng"][r], ref["rg"][r] = lw_ft[ip]["surface_pNG_Pa"], lw_ft[ip]["surface_pRG_Pa"]
            ref["ng_fit"][r], ref["rg_fit"][r] = lw_fit[ip]["surface_pNG_Pa"], lw_fit[ip]["surface_pRG_Pa"]
            src.append("late_window_model")
        else:
            if pc[ip]["kind"] == "late_obs":
                raise SystemExit("late observation %d not in late_window_model.dat" % ip)
            ref["ng"][r], ref["rg"][r] = pc[ip]["pNG_tuned_Pa"], pc[ip]["pRG_tuned_Pa"]       # grain 0 before t_end
            ref["ng_fit"][r], ref["rg_fit"][r] = pc[ip]["pNG_fit_Pa"], pc[ip]["pRG_fit_Pa"]
            src.append("pressure_comparison")
    ver = {"table_precision": "the tables print 10 significant digits: agreement to rounding is max rel <= 5e-10"}
    for key, label in (("ng", "NG final"), ("rg", "RG final"), ("ng_fit", "NG fit state"), ("rg_fit", "RG fit state")):
        d = rel(model[key], ref[key])
        ver[label] = {"max_rel": float(np.max(d)), "late_obs_max_rel": float(np.max(d[late_rows])),
                      "argmax_iPoint": int(ipoint[int(np.argmax(d))])}
    for key, label in (("ng_fp64", "NG final, FP64 records (diagnostic)"), ("rg_fp64", "RG final, FP64 records (diagnostic)")):
        d = rel(model[key], ref[key[:2]])
        ver[label] = {"max_rel": float(np.max(d)), "late_obs_max_rel": float(np.max(d[late_rows]))}
    ver["bins_sum_vs_total_max_rel"] = {g: float(np.max(rel(mbins[g].sum(axis=1), model[g.lower()]))) for g in GAUGES}
    ver["passed_tables"] = bool(all(ver[l]["max_rel"] <= 5.0e-10 * (1 + 1e-6) + 1e-14 for l in
                                    ("NG final", "RG final", "NG fit state", "RG fit state")))
    say("verification vs tables: " + json.dumps(ver, sort_keys=True))
    # grain of the late observations: table and formula
    t_end = float(json.load(open(fit_dir / "grain_term.json"))["t_end_et"]) if (fit_dir / "grain_term.json").exists() else None
    late_ip = ipoint[late_rows]
    grain_tab = np.array([lw_ft[int(i)]["grain_NG_Pa"] for i in late_ip])
    late_et = layout.late_et
    grain_formula = np.maximum(0.0, grain_A_uPa * 1e-6 - grain_B * 1e-8 * (late_et - t_end)) * (late_et >= t_end)
    ver["grain_formula_vs_table_max_rel"] = float(np.max(rel(grain_formula, grain_tab)))
    # ---------------------------------------------------------------- phase B: the late observations
    _G.update({"W": W, "grid_rows": gr, "FAC": FAC})
    say("phase B: %d late observations" % nlate)
    OB = [None] * nlate
    with mp.get_context("fork").Pool(min(args.procs, 48)) as pool:
        for res in pool.imap_unordered(obs_work, range(nlate), chunksize=1):
            OB[res["k"]] = res
    FAC.clear()
    obs_ng = model["ng"][late_rows]
    obs_rg = model["rg"][late_rows]
    vb = {}
    for g, m in (("NG", obs_ng), ("RG", obs_rg)):
        tf = np.array([o[g]["total_face_order"] for o in OB])
        ts = np.array([o[g]["total_sorted"] for o in OB])
        wr = np.array([o["W_rows_total"][g] for o in OB])
        tab = np.array([lw_ft[int(i)]["surface_p%s_Pa" % g] for i in late_ip])
        vb[g] = {"sum_facets_vs_table_max_rel": float(np.max(rel(tf, tab))),
                 "sum_facets_vs_W_rowtotals_max_rel": float(np.max(rel(tf, wr))),
                 "sorted_vs_face_order_max_rel": float(np.max(rel(ts, tf))),
                 "negative_facet_contributions": int(sum(o[g]["negative"] for o in OB))}
    vb["W_rowsum_range"] = [float(min(o["W_rowsum"] for o in OB)), float(max(o["W_rowsum"] for o in OB))]
    ver["late_facet_sums"] = vb
    ver["passed_late"] = bool(all(vb[g]["sum_facets_vs_table_max_rel"] <= 5.0e-10 * (1 + 1e-6) + 1e-14 and
                                  vb[g]["sum_facets_vs_W_rowtotals_max_rel"] < 1e-12 and
                                  vb[g]["negative_facet_contributions"] == 0 for g in GAUGES))
    say("late facet sums: " + json.dumps(vb, sort_keys=True))
    if not (ver["passed_tables"] and ver["passed_late"]):
        raise SystemExit("verification failed")
    # occurrence counts
    stage.mkdir(parents=True)
    occ = {g: np.zeros((len(P_LEVELS), nf), dtype=np.uint8) for g in GAUGES}
    for o in OB:
        for g in GAUGES:
            for i, n in enumerate(o[g]["k"]):
                occ[g][i][o[g]["faces"][:n]] += 1
    area = el["area"]
    late_utc = [utc[int(r)] for r in late_rows]
    for g in GAUGES:
        offs = np.concatenate([[0], np.cumsum([len(o[g]["faces"]) for o in OB])]).astype(np.int64)
        np.savez(stage / ("source_areas_%s.npz" % g),
                 p_levels=np.array(P_LEVELS), iPoint=late_ip, utc=np.array(late_utc), et=late_et,
                 stop_et=layout.stop_et, model_surface_Pa=(obs_ng if g == "NG" else obs_rg),
                 table_surface_Pa=np.array([lw_ft[int(i)]["surface_p%s_Pa" % g] for i in late_ip]),
                 grain_Pa=(grain_tab if g == "NG" else np.zeros(nlate)),
                 total_facets_Pa=np.array([o[g]["total_sorted"] for o in OB]),
                 k=np.array([o[g]["k"] for o in OB]), area_m2=np.array([o[g]["area"] for o in OB]),
                 seen_area_m2=np.array([o[g]["seen_area"] for o in OB]),
                 c_min_Pa=np.array([o[g]["c_min"] for o in OB]), cum_fraction_at_k=np.array([o[g]["cum_at_k"] for o in OB]),
                 n_contrib=np.array([o[g]["n_contrib"] for o in OB]), area_contrib_m2=np.array([o[g]["area_contrib"] for o in OB]),
                 offsets=offs, faces_sorted=np.concatenate([o[g]["faces"] for o in OB]),
                 contrib_sorted_Pa=np.concatenate([o[g]["contrib"] for o in OB]),
                 occurrence=occ[g], union=(occ[g] > 0))
        say("wrote source_areas_%s.npz" % g)
    # ---------------------------------------------------------------- slope bins
    bins_e = S["bins"].astype(np.intp)
    seen = S["seen"]
    G_ = NM.read_geometry(args.geometry)
    prod = {}
    for lab in ("t11", "t12"):
        Tt = G_["times"][lab]
        cos, sh = NM.geometry_time(args.geometry, lab, nf)
        F_h, F_c, _, _ = NM.source_rates(cos, sh, Tt["rh_m"], helio, sza, S["A"], el["ACO2"], S["kappa"], fc.AU)
        P_e = F_h * area
        prod[Tt["utc"]] = {"seen": np.bincount(bins_e[seen], weights=P_e[seen], minlength=NBINS) / np.sum(P_e[seen]),
                           "all": np.bincount(bins_e, weights=P_e, minlength=NBINS) / np.sum(P_e),
                           "Q_seen_s-1": float(np.sum(P_e[seen])), "Q_all_s-1": float(np.sum(P_e))}
    bin_area = np.bincount(bins_e, weights=area, minlength=NBINS)
    bin_seen = np.bincount(bins_e[seen], weights=area[seen], minlength=NBINS)
    frac = {g: mbins[g] / mbins[g].sum(axis=1, keepdims=True) for g in GAUGES}
    kind = np.array([pc[int(i)]["kind"] for i in ipoint])
    np.savez(stage / "slope_bins.npz", bin_names=np.array(S["checks"]["bin_names"]),
             bin_labels_deg=np.array(S["checks"]["bin_labels_deg"]), multiplier=tv[:NBINS],
             bin_area_m2=bin_area, bin_area_fraction=bin_area / bin_area.sum(), bin_seen_area_m2=bin_seen,
             seen_area_fraction_of_bin=bin_seen / bin_area, bin_share_of_seen_area=bin_seen / bin_seen.sum(),
             production_share_seen_103310=prod["2016-09-30T10:33:10"]["seen"],
             production_share_seen_103926=prod["2016-09-30T10:39:26"]["seen"],
             production_share_all_103310=prod["2016-09-30T10:33:10"]["all"],
             production_share_all_103926=prod["2016-09-30T10:39:26"]["all"],
             iPoint=ipoint, utc=np.array(utc), et=np.array([R[r]["et"] for r in range(nrow)]), kind=kind,
             use_ng=use_ng, use_rg=use_rg, model_surface_NG_Pa=model["ng"], model_surface_RG_Pa=model["rg"],
             table_surface_NG_Pa=ref["ng"], table_surface_RG_Pa=ref["rg"], bins_NG_Pa=mbins["NG"], bins_RG_Pa=mbins["RG"],
             fraction_NG=frac["NG"], fraction_RG=frac["RG"])
    say("wrote slope_bins.npz")
    write_texts(stage, OB, late_ip, late_utc, obs_ng, obs_rg, grain_tab, S, tv, bin_area, bin_seen, prod, ipoint, utc,
                kind, use_ng, use_rg, frac, area)
    man = {"schema": SCHEMA, "created_utc": NM.utc_now(), "host": socket.gethostname(),
           "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv, "model": S["model"],
           "state": {"ft_dir": str(ft_dir), "ft_dir_resolved": str(ft_dir.resolve()), "fit_dir": str(fit_dir),
                     "fit_dir_resolved": str(fit_dir.resolve()),
                     "inputs_sha256": {n: NM.sha256_file(Path(d) / n) for d, n in (
                         (ft_dir, "finetune_state.npz"), (ft_dir, "late_window_model.dat"),
                         (ft_dir, "pressure_comparison.txt"), (fit_dir, "nonneg_facets.npz"),
                         (fit_dir, "fit_restart.npz"), (fit_dir, "late_window_model.dat"))}},
           "base": str(base.resolve()), "geometry": str(Path(args.geometry).resolve()), "law": law,
           "gauge_constants": gconst, "T": T, "kappa": S["kappa"], "grain_A_uPa": grain_A_uPa,
           "grain_B_uPa_per_100s": grain_B, "t_end_et": t_end, "state_checks": S["checks"],
           "late_layout": layout.summary, "positive_problem": pprec, "point_table": ptab,
           "precision": "records as the fit: cops_cuda.device_store.fp32_record_fields (FP32 Kn, KPhi, KV, cos; shadow "
                        "-> cos -2), FP64 arithmetic",
           "negative_record_contributions": neg, "records": int(sum(R[r]["records"] for r in R)),
           "verification": ver, "p_levels": list(P_LEVELS), "procs": args.procs,
           "code_sha256": {p.name: NM.sha256_file(p) for p in sorted(HERE.glob("*nonneg_contrib*")) if p.suffix == ".py"},
           "code_dir": str(HERE), "wall_seconds": time.time() - t0}
    with open(stage / "contrib_manifest.json", "x") as f:
        json.dump(man, f, indent=1, sort_keys=True, default=float)
    write_readme(stage, man)
    stage.rename(out)
    say("done: %s (%.0f s)" % (out, time.time() - t0))
    return 0


def write_texts(stage, OB, late_ip, late_utc, obs_ng, obs_rg, grain, S, tv, bin_area, bin_seen, prod, ipoint, utc, kind,
                use_ng, use_rg, frac, area):
    tot_area = float(np.sum(area))
    L = ["# per late observation (W10x5, 10:33:10 .. 10:39:26 UTC): model surface pressure, NG grain share, and the P% "
         "source areas (k = number of facets, A = area [m^2]) of NG and RG",
         "# grain_share_NG = grain / (surface + grain); area fractions are of the total nucleus area %.6e m^2" % tot_area,
         "#iPoint UTC surf_NG_Pa grain_NG_Pa grain_share_NG surf_RG_Pa " +
         " ".join("k%d_NG A%d_NG" % (p, p) for p in P_LEVELS) + " " + " ".join("k%d_RG A%d_RG" % (p, p) for p in P_LEVELS)]
    for k, o in enumerate(OB):
        L.append("%d %s %.9e %.9e %.6f %.9e %s %s" % (
            late_ip[k], late_utc[k], obs_ng[k], grain[k], grain[k] / (obs_ng[k] + grain[k]), obs_rg[k],
            " ".join("%d %.6e" % (n, a) for n, a in zip(o["NG"]["k"], o["NG"]["area"])),
            " ".join("%d %.6e" % (n, a) for n, a in zip(o["RG"]["k"], o["RG"]["area"]))))
    with open(stage / "source_areas_per_observation.txt", "x") as f:
        f.write("\n".join(L) + "\n")
    B = ["# per effective-gravitational-slope bin: fitted multiplier, area, seen area, H2O production shares (F = q A, "
         "sum over facets of F area) at 10:33:10 and 10:39:26 UTC, of the seen facets and of all facets",
         "#bin interval_deg multiplier area_fraction seen_area_fraction_of_bin share_of_seen_area prod_seen_103310 "
         "prod_seen_103926 prod_all_103310 prod_all_103926"]
    p1, p2 = prod["2016-09-30T10:33:10"], prod["2016-09-30T10:39:26"]
    for b in range(NBINS):
        B.append("%s %s %.9g %.6e %.6e %.6e %.6e %.6e %.6e %.6e" % (
            S["checks"]["bin_names"][b], S["checks"]["bin_labels_deg"][b], tv[b], bin_area[b] / bin_area.sum(),
            bin_seen[b] / bin_area[b], bin_seen[b] / bin_seen.sum(), p1["seen"][b], p2["seen"][b], p1["all"][b],
            p2["all"][b]))
    B += ["", "# per fit point (all base rows; late_obs = W10x5 of the grid rows): fraction of the model surface pressure "
              "from each slope bin, NG then RG",
          "#iPoint UTC kind useNG useRG " + " ".join("fNG_%s" % n for n in S["checks"]["bin_names"]) + " " +
          " ".join("fRG_%s" % n for n in S["checks"]["bin_names"])]
    for r in range(len(ipoint)):
        B.append("%d %s %s %d %d %s %s" % (ipoint[r], utc[r], kind[r], use_ng[r], use_rg[r],
                                           " ".join("%.6e" % v for v in frac["NG"][r]),
                                           " ".join("%.6e" % v for v in frac["RG"][r])))
    with open(stage / "slope_bins.txt", "x") as f:
        f.write("\n".join(B) + "\n")


def write_readme(stage, man):
    v = man["verification"]
    L = ["Facet contributions, P%% source areas and slope-bin fractions -- %s model, FINAL state" % man["model"],
         "Written %s by nonneg_contrib.py (%s) on %s (Slurm job %s)." % (man["created_utc"], man["schema"], man["host"],
                                                                      man["slurm_job_id"]),
         "State: ft %s, fit %s; base %s." % (man["state"]["ft_dir"], man["state"]["fit_dir"], man["base"]),
         "", "DEFINITIONS",
         "Contribution c_k,g,e [Pa] of facet e to the MODEL surface pressure of point k at gauge g (NG, RG): the fitted",
         "  forward model term by term (fit_cops.Model load_kernels / prepare / evaluate = cops_cuda pressure.cu mixed",
         "  kernel): records in the fit's mixed precision (FP32 Kn, KPhi, KV, cos; shadow folded into cos = -2;",
         "  cops_cuda.device_store.fp32_record_fields), flux F = q A (night/shadow rule, r_h of the point), A = the final",
         "  state's activity (finetune_state.npz A_final), T = %s K and kappa = %.6g held; H2O + CO2 terms." % (
             list(man["T"]), man["kappa"]),
         "  Late observations (95, 10:33:10 .. 10:39:26 UTC): c_k,g,e = sum_j W_kj c_j,g,e over the 1-s grid rows j (W10x5,",
         "  positive_problem.LateLayout of the chain). The NG grain term is not a facet contribution (grain_share_NG).",
         "P% SOURCE AREA S_P(k, g) (proposed name, after the footprint 'source area' of Schmid 1994): the smallest set of",
         "  facets, taken in decreasing order of c_k,g,e (ties: lower facet index first), whose cumulative contribution",
         "  reaches P % of sum_e c_k,g,e (the observation's total surface contribution), P = 50, 60, 70, 80, 90.",
         "SOURCE-AREA OCCURRENCE N_P,g(e): the number of the 95 late observations whose S_P(k, g) contains facet e (0..95);",
         "  the union mask is N_P,g > 0.",
         "", "VERIFICATION (tables print 10 significant digits: rounding <= 5e-10)",
         "  per row sum vs the model tables (late_window_model.dat, pressure_comparison.txt), all 1551 rows:",
         "    final state NG %.2e, RG %.2e; fit state NG %.2e, RG %.2e (max rel)" % (
             v["NG final"]["max_rel"], v["RG final"]["max_rel"], v["NG fit state"]["max_rel"], v["RG fit state"]["max_rel"]),
         "    diagnostic, FP64 records instead of the fit's FP32: NG %.2e, RG %.2e" % (
             v["NG final, FP64 records (diagnostic)"]["max_rel"], v["RG final, FP64 records (diagnostic)"]["max_rel"]),
         "  late observations, sum over facets vs late_window_model.dat: NG %.2e, RG %.2e; vs W x row totals: NG %.2e, RG %.2e"
         % (v["late_facet_sums"]["NG"]["sum_facets_vs_table_max_rel"], v["late_facet_sums"]["RG"]["sum_facets_vs_table_max_rel"],
            v["late_facet_sums"]["NG"]["sum_facets_vs_W_rowtotals_max_rel"],
            v["late_facet_sums"]["RG"]["sum_facets_vs_W_rowtotals_max_rel"]),
         "  negative record contributions: NG %d, RG %d (of %d records); grain formula vs table: %.2e" % (
             man["negative_record_contributions"]["NG"], man["negative_record_contributions"]["RG"], man["records"],
             v["grain_formula_vs_table_max_rel"]),
         "  slope-bin sums vs row totals: NG %.2e, RG %.2e" % (v["bins_sum_vs_total_max_rel"]["NG"],
                                                             v["bins_sum_vs_total_max_rel"]["RG"]),
         "", "FILES",
         "  source_areas_NG.npz, source_areas_RG.npz: per late observation (index k, 95): iPoint, utc, et, stop_et,",
         "    model_surface_Pa (W10x5 of the per-row sums), table_surface_Pa, grain_Pa (NG), total_facets_Pa, and per P",
         "    (p_levels, axis 1): k (facets in S_P), area_m2, seen_area_m2, c_min_Pa (smallest contribution in S_P),",
         "    cum_fraction_at_k; n_contrib, area_contrib_m2 (facets with c > 0); the ranking: faces_sorted (int32) and",
         "    contrib_sorted_Pa for observation k in [offsets[k], offsets[k+1]) in decreasing order, so S_P(k) =",
         "    faces_sorted[offsets[k] : offsets[k] + k[k, P index]]; occurrence (5 x nfaces uint8) = N_P,g; union.",
         "  slope_bins.npz / slope_bins.txt: per bin multiplier, area fraction, seen area fraction of the bin, share of",
         "    the seen area, H2O production shares (seen / all facets) at 10:33:10 and 10:39:26; per base row (all 1551;",
         "    kind obs / grid / late_obs = W10x5; use_ng, use_rg = the fit's selection) bins_*_Pa and fraction_* (7 bins).",
         "  source_areas_per_observation.txt: the per-observation table (k, area per P and gauge, grain share).",
         "  contrib_manifest.json: inputs (sha256), constants, checks, verification, code sha256.",
         "", "COMMAND", "  " + " ".join(man["command"])]
    with open(stage / "README.txt", "x") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
