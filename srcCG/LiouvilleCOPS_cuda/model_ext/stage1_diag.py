# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/checks/stage1_diag_20261114/stage1_diag.py (on Voyager: /data001/highresolu_vt/cops_output/checks/stage1_diag_20261114/stage1_diag.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 809dbdd0ab0289cfb45d2f99a01fc9602fa5be136669637a0ed3abb55992e797
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: check. Convergence diagnosis of a stage-1 fit on CPU (lead 2026-11-14 16:20 CST): exact data cost and
#   Jacobian (fit's FP32 records), gradient, bound-constrained Gauss-Newton step, actual objective along it and along
#   scaled steepest descent, column norms at the start / end (scipy x_scale='jac' history).
# USED BY: checks/stage1_diag_20261114/run.sbatch (jobs 125197-125203, 125211-125215, 125216-125229).
# ==== end of header ====
"""Convergence diagnosis of a finished stage-1 fit of the alternative non-negative models (lead 2026-11-14 ~16:20 CST):
is the end point a stationary point of the fitter's objective, or an artefact of the optimizer?

Everything is evaluated on the CPU from the fit's own products, with the fit's frozen model code:
  activity and its analytic column factors: nonneg_mult.host_fields (exactly the fit's host fields; checked bitwise
    against nonneg_facets.npz A);
  pressures: per-record H2O coefficients and the CO2 unit term of the fitted forward model with the fit's FP32 records
    (nonneg_contrib.record_contributions, frozen; kappa = 0 for the H2O-linear part, kappa = 1 with A = 0 for the CO2
    unit), p = H2O(A) + kappa CO2; W10x5 for the late observations; + grain for NG; floors;
  residuals r_d = w (ln p - ln obs) over the fit's used rows (w = 1/sqrt(N_gauge)), the option-D penalty rows
    sqrt(2 f C_ref) S_x x_a; the crossing rows are inactive at every end point (penalty 0, no wrong side) and left out.
Checks: 0.5 r_d.r_d vs the manifest data cost; singular values of J_d vs the manifest's data_jacobian_singular_values.
Diagnosis at the end point x*:
  gradient g = J^T r (data and penalty parts), scaled by the column norms of J (scipy x_scale='jac'), with the
    bound-active components (x at its lower bound and g > 0) reported apart;
  the bound-constrained Gauss-Newton step (scipy.optimize.lsq_linear) and its predicted decrease;
  the ACTUAL objective along that step (t = 1, 0.5, 0.25, 0.1) and along the projected scaled steepest-descent
    direction (t = 0.5, 1, 2 x the Cauchy step).
    stage1_diag.py <summary snapshot> <nonneg snapshot> <fit dir (final)> <out.json>
"""
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

SX, NX, FIT, OUTJ = sys.argv[1:5]
sys.path.insert(0, "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, NX)
sys.path.insert(0, SX)                      # nonneg_maps / nonneg_contrib / positive_problem of the summary snapshot
import nonneg_maps as NM                    # noqa: E402
import nonneg_contrib as NC                 # noqa: E402
import positive_problem as pp               # noqa: E402
import nonneg_mult as nm                    # noqa: E402  (the nonneg snapshot)
import nonneg_alt as na                     # noqa: E402
import nonneg_problem as npb                # noqa: E402
R = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5")
BASE = R / "base"
T0 = (182.09922, 136.11149, 133.0)
CH = 1 << 17
_W = {}


def say(t):
    print("stage1_diag: [%s] %s" % (time.strftime("%H:%M:%S", time.gmtime()), t), flush=True)


def work(task):
    """Row r: H2O per-record coefficients (kappa 0, A = 1) dotted with every column of X, and the CO2 unit sum."""
    r, ip = task
    kp, h, rec = NM.read_kernel_full(BASE / "kernels" / ("point=%04d.bin" % ip))
    NC._G["kappa"] = 0.0
    face, ((an, ar),), _ = NC.record_contributions(rec, h, (_W["ones"],), "fp32", r)
    NC._G["kappa"] = 1.0
    _, ((cn, cr),), _ = NC.record_contributions(rec, h, (_W["zeros"],), "fp32", r)
    X = _W["X"]
    sn = np.zeros(X.shape[1])
    sr = np.zeros(X.shape[1])
    for c0 in range(0, len(face), CH):
        Xc = X[face[c0:c0 + CH]]
        sn += an[c0:c0 + CH] @ Xc
        sr += ar[c0:c0 + CH] @ Xc
    return r, sn, sr, float(np.sum(cn)), float(np.sum(cr)), h["et"]


def record_pass(X, rows, cols, ipoint, procs):
    _W["X"] = np.ascontiguousarray(X)
    n = len(ipoint)
    SN, SR = np.zeros((n, X.shape[1])), np.zeros((n, X.shape[1]))
    CN, CR, ET = np.zeros(n), np.zeros(n), np.zeros(n)
    tasks = sorted([(r, int(ipoint[r])) for r in range(n)], key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    with mp.get_context("fork").Pool(procs) as pool:
        for r, sn, sr, cn, cr, et in pool.imap_unordered(work, tasks, chunksize=1):
            SN[r], SR[r], CN[r], CR[r], ET[r] = sn, sr, cn, cr, et
    return SN, SR, CN, CR, ET


def main():
    t0 = time.time()
    procs = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
    fit = Path(FIT).resolve()
    man = json.load(open(fit / "fit_manifest.json"))
    rec = man["model_extension"]["nonneg"]
    cfg = rec["config"]
    model = cfg["model"]
    model_cfg = {"model": model, "sos_rank": cfg.get("sos_rank")}
    opt = man["optimizer_results"][0]
    with np.load(fit / "fit_restart.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        names = list(meta["names"])
        x, lo, hi = (np.asarray(z[k], float) for k in ("x", "lo", "hi"))
        ref = {"NG": np.asarray(z["ref_ng"], float), "RG": np.asarray(z["ref_rg"], float)}
        use = {"NG": (np.asarray(z["use_ng"]) > 0) & (np.asarray(z["ref_ng"]) > 0),
               "RG": (np.asarray(z["use_rg"]) > 0) & (np.asarray(z["ref_rg"]) > 0)}
        point_ids = np.asarray(z["point_ids"], np.int64)
    ascale, A_s = float(meta["ascale"]), float(meta["enhancement_scale"])
    DIAG_X = os.environ.get("DIAG_X", "final")
    if DIAG_X == "start":            # the fit's start point (column norms at the start: scipy's x_scale='jac' history)
        x = np.asarray(rec["start_restart"]["x_start"], float)
        if list(rec["start_restart"].get("names") or names) != names:
            raise SystemExit("start layout differs")
    if cfg.get("held"):
        raise SystemExit("a held fine-tuning term: stage-1 fits only")
    from cops_cuda import reference, device_store
    reference.bind()
    fc = reference.load("fit_cops")
    at = reference.load("activity_templates")
    el = fc.read_elements(str(BASE / "elements.bin"))
    Y = np.ascontiguousarray(el["Y"], np.float64)
    area, aorig = np.asarray(el["area"]), np.asarray(el["AH2O"], np.float64)
    nf = len(area)
    with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
        A_fit, seen, bins = np.asarray(z["A"]), np.asarray(z["seen"]).astype(bool), np.asarray(z["bin"])
    fts, _ = npb.ft_templates_from(cfg["ft_templates_from"])
    Bft = np.ascontiguousarray(np.column_stack([at.select_templates(at.load_bundle(t["bundle"]["path"], el["centroid"],
                                                                                  area), [t["name"]]) for t in fts]))
    cft = np.ones(Bft.shape[1])
    s_w = cfg["s_ratio"] * aorig if model == "background" else None
    ia = [names.index("a%02d" % j) for j in range(25)]
    im = [names.index("c_" + b) for b in nm.BIN_NAMES[1:]]
    ik = names.index("kappa")
    ig = [names.index("grain_A_uPa"), names.index("grain_B_uPa_per_100s")]
    en = na.extra_names(model_cfg)
    ie = [names.index(n) for n in en]
    ao = aorig if model in ("background", "logspread") else None

    def fields_of(xv, want, mode=0):
        a = xv[ia] * ascale
        m = np.concatenate([[1.0], xv[im]])
        extra = na.extra_from_values(model_cfg, xv[ie], ascale) if en else {}
        return nm.host_fields(model, Y, a, Bft, cft, A_s, bins, m, s_w, None, want, extra=extra, aorig=ao, mode=mode)

    f0 = fields_of(x, True)
    A = f0["A"]
    a_chk = float(np.max(np.abs(A - A_fit) / np.maximum(np.abs(A_fit), 1e-300))) if DIAG_X == "final" else None
    say("model %s, %d parameters, point %s; A vs nonneg_facets.npz max rel %s" % (model, len(x), DIAG_X, a_chk))
    # columns of dA / dx (optimizer coordinates), in the layout of names
    cols_list, col_of = [], {}
    for j, k in enumerate(ia):
        col_of[k] = len(cols_list)
        cols_list.append(f0["D"] * Y[:, j] * ascale)
    for b, k in enumerate(im):
        col_of[k] = len(cols_list)
        cols_list.append(np.where(bins == b + 1, f0["P"], 0.0))
    if model == "sos":
        r_ = int(cfg["sos_rank"])
        Dm = {jj: fields_of(x, True, mode=jj)["D"] for jj in range(1, r_)}
        for n_, k in zip(en, ie):
            i, jj = int(n_[5:7]), int(n_[8:])
            col_of[k] = len(cols_list)
            cols_list.append(Dm[jj] * Y[:, i] * ascale)
    elif model == "background":
        col_of[ie[0]] = len(cols_list)
        cols_list.append(fields_of(x, True, mode="raw")["D"])
    X = np.column_stack([A] + cols_list)
    del cols_list
    say("activity columns: %d (%.1f GB)" % (X.shape[1], X.nbytes / 1e9))
    fc_, header, cols, rows = NM.base_table(BASE)
    helio = fc_.header_value(header, "HeliocentricModel")
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    assert np.array_equal(ipoint, point_ids)
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    bm, cNG, cRG, _ = NC.gauge_constants(fc_, header, TRG)
    NC._G.update({"ACO2": el["ACO2"], "T": T0, "kappa": 0.0, "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio,
                  "AU": float(fc_.AU), "COS75": float(fc_.COS75), "KT": float(fc_.KT),
                  "fp32_record_fields": device_store.fp32_record_fields})
    _W.update({"ones": np.ones(nf), "zeros": np.zeros(nf)})
    SN, SR, CN, CR, ET = record_pass(X, rows, cols, ipoint, procs)
    say("record pass 1 done (%.0f s)" % (time.time() - t0))
    lay = pp.LateLayout(pp.read_point_table(cfg["point_table"]), pp.read_run_points(BASE), ipoint)
    W = lay.W.tocsr()

    def late(v):
        v = v.copy()
        v[lay.late_rows] = W @ v[lay.grid_rows]
        return v
    t_end = float(json.load(open(fit / "grain_term.json"))["t_end_et"])
    s_ = ET - t_end
    on = s_ >= 0.0
    NNG, NRG = int(use["NG"].sum()), int(use["RG"].sum())
    wN, wR = 1.0 / math.sqrt(NNG), 1.0 / math.sqrt(NRG)

    def pressures(sn_col, sr_col, cn, cr, kappa, gA, gB):
        g = np.maximum(np.where(on, gA * 1e-6 - gB * 1e-8 * np.maximum(s_, 0.0), 0.0), 0.0)
        pN = np.maximum(late(sn_col + kappa * cn) + g, 1e-300)
        pR = np.maximum(late(sr_col + kappa * cr), 1e-300)
        return pN, pR, g

    def resid(pN, pR):
        return np.concatenate([wN * (np.log(pN[use["NG"]]) - np.log(ref["NG"][use["NG"]])),
                               wR * (np.log(pR[use["RG"]]) - np.log(ref["RG"][use["RG"]]))])
    kap, gA, gB = x[ik], x[ig[0]], x[ig[1]]
    pN, pR, g = pressures(SN[:, 0], SR[:, 0], CN, CR, kap, gA, gB)
    rd = resid(pN, pR)
    data_cost = 0.5 * float(rd @ rd)
    # data Jacobian in the layout of names
    J = np.zeros((len(rd), len(x)))
    for k, c in col_of.items():
        dN, dR = late(SN[:, c + 1]), late(SR[:, c + 1])
        J[:, k] = np.concatenate([wN * dN[use["NG"]] / pN[use["NG"]], wR * dR[use["RG"]] / pR[use["RG"]]])
    dN, dR = late(CN), late(CR)
    J[:, ik] = np.concatenate([wN * dN[use["NG"]] / pN[use["NG"]], wR * dR[use["RG"]] / pR[use["RG"]]])
    gpos = on & (gA * 1e-6 - gB * 1e-8 * np.maximum(s_, 0.0) >= 0.0)
    dgA = np.where(gpos, 1e-6, 0.0)
    dgB = np.where(on & (g > 0.0), -1e-8 * np.maximum(s_, 0.0), 0.0)
    J[:, ig[0]] = np.concatenate([wN * dgA[use["NG"]] / pN[use["NG"]], np.zeros(NRG)])
    J[:, ig[1]] = np.concatenate([wN * dgB[use["NG"]] / pN[use["NG"]], np.zeros(NRG)])
    sv = np.linalg.svd(J, compute_uv=False)
    sv_m = np.asarray(opt.get("data_jacobian_singular_values") or [], float)
    sv_rel = float(np.max(np.abs(sv[:len(sv_m)] / sv_m - 1.0))) if len(sv_m) == len(sv) else None
    # penalty rows (option D)
    if model == "logspread":
        S, C, ybar, Wseen = na.spread_matrix(Y, area, seen)
        f_, cref = float(cfg["spread_f"]), float(cfg["spread_cref"])
        Sx = (ascale / A_s) * S
        Jp = np.zeros((25, len(x)))
        Jp[:, ia] = math.sqrt(2.0 * f_ * cref) * Sx
    else:
        Jp = np.zeros((0, len(x)))
    rp = Jp @ x if len(Jp) else np.zeros(0)
    pen = 0.5 * float(rp @ rp)
    F = data_cost + pen
    out = {"fit": str(fit), "model": model, "config": {k: cfg.get(k) for k in ("bg_min", "sos_rank", "spread_f", "spread_cref")},
           "names": names, "x": x.tolist(),
           "checks": {"A_vs_nonneg_facets_max_rel": a_chk, "data_cost_cpu": data_cost,
                      "data_cost_manifest": float(opt["data_cost"]),
                      "data_cost_rel_diff": abs(data_cost / float(opt["data_cost"]) - 1.0),
                      "singular_values_rel_max_vs_manifest": sv_rel,
                      "penalty_cpu": pen, "penalty_manifest": ((rec["model_check"].get("alternative") or {}).get("spread")
                                                               or {}).get("final", {}).get("penalty")},
           "optimizer": {k: opt.get(k) for k in ("cost", "nfev", "message", "status", "optimality")}}
    say("data cost CPU %.12e vs manifest %.12e (rel %.1e); singular values max rel %s; penalty %.6e" % (
        data_cost, float(opt["data_cost"]), out["checks"]["data_cost_rel_diff"], sv_rel, pen))
    # ---- gradient, scaling, bounds
    Jf = np.vstack([J, Jp])
    rf = np.concatenate([rd, rp])
    gd, gp = J.T @ rd, (Jp.T @ rp if len(Jp) else np.zeros(len(x)))
    gr = gd + gp
    cn_ = np.linalg.norm(Jf, axis=0)
    cn_[cn_ == 0] = 1.0
    atlo = np.isfinite(lo) & ((x - lo) <= 1e-8 * np.maximum(1.0, np.abs(np.where(np.isfinite(lo), lo, 0.0))))
    active = atlo & (gr > 0)
    gs = gr / cn_
    free = ~active
    order = np.argsort(-np.abs(gs * free))
    out["point"] = DIAG_X
    out["colnorm"] = cn_.tolist()
    out["gradient"] = {"g_inf_free_scaled": float(np.max(np.abs(gs[free]))) if free.any() else 0.0,
                       "g_inf_raw": float(np.max(np.abs(gr))),
                       "g_scaled_norm_free": float(np.linalg.norm(gs[free])),
                       "relative_to_sqrt_2F": float(np.linalg.norm(gs[free]) / math.sqrt(2 * F)),
                       "bound_active": [names[k] for k in np.flatnonzero(active)],
                       "top_free_scaled": [(names[k], float(gs[k]), float(gd[k] / cn_[k]), float(gp[k] / cn_[k]))
                                           for k in order[:8]],
                       "note": "scaled g_k = (J^T r)_k / ||J_k|| (x_scale='jac'); entries: (name, total, data part, "
                               "penalty part); bound_active: at the lower bound with g > 0 (KKT)"}
    # ---- bound-constrained Gauss-Newton step
    from scipy.optimize import lsq_linear
    Js = Jf / cn_
    lb_s = (lo - x) * cn_
    ub_s = (hi - x) * cn_
    lb_s = np.minimum(lb_s, 0.0)
    ub_s = np.maximum(ub_s, 0.0)
    sol = lsq_linear(Js, -rf, bounds=(lb_s, ub_s), method="bvls", tol=1e-14, max_iter=2000)
    dlt = sol.x / cn_
    pred = F - 0.5 * float(np.sum((Jf @ dlt + rf) ** 2))
    out["gauss_newton"] = {"predicted_decrease": pred, "predicted_relative": pred / F,
                           "step_scaled_norm": float(np.linalg.norm(sol.x)), "lsq_linear_status": int(sol.status),
                           "largest_step_components": [(names[k], float(dlt[k]), float(x[k])) for k in
                                                       np.argsort(-np.abs(sol.x))[:8]]}
    say("F %.12e (data %.6e, penalty %.6e); scaled |g|_inf (free) %.3e; GN predicted decrease %.3e (relative %.2e)" % (
        F, data_cost, pen, out["gradient"]["g_inf_free_scaled"], pred, pred / F))
    # ---- trial points: along the GN step and along the projected scaled steepest descent
    dsd = -(gr * free) / cn_ ** 2
    Jd_ = Jf @ dsd
    t_c = float(-(gr @ dsd) / max(Jd_ @ Jd_, 1e-300)) if np.any(dsd) else 0.0
    tg = [float(v) for v in os.environ.get("TRIALS_GN", "1:0.5:0.25:0.1").split(":") if v != "none"]
    tsd = [float(v) for v in os.environ.get("TRIALS_SD", "0.5:1:2").split(":") if v != "none"]
    if not tg and not tsd:
        out["trials"] = []
        json.dump(out, open(OUTJ, "x"), indent=1, default=float)
        say("wrote %s without trials (%.0f s)" % (OUTJ, time.time() - t0))
        return 0
    trials = [("GN t=%g" % t, np.maximum(x + t * dlt, lo)) for t in tg]
    trials += [("SD %gx Cauchy" % t, np.maximum(x + t * t_c * dsd, lo)) for t in tsd]
    Xt = np.column_stack([fields_of(xt, False)["A"] for _, xt in trials])
    SN2, SR2, CN2, CR2, _ = record_pass(Xt, rows, cols, ipoint, procs)
    res = []
    for i, (lab, xt) in enumerate(trials):
        pN_, pR_, _ = pressures(SN2[:, i], SR2[:, i], CN2, CR2, xt[ik], xt[ig[0]], xt[ig[1]])
        r_ = resid(pN_, pR_)
        dc = 0.5 * float(r_ @ r_)
        pc = 0.5 * float(np.sum((Jp @ xt) ** 2)) if len(Jp) else 0.0
        lin = F - 0.5 * float(np.sum((Jf @ (xt - x) + rf) ** 2))
        res.append({"trial": lab, "F": dc + pc, "data": dc, "penalty": pc, "actual_decrease": F - dc - pc,
                    "relative_actual_decrease": (F - dc - pc) / F, "predicted_linear_model": lin})
        say("  %-16s F %.12e (data %.6e, penalty %.6e): actual decrease %+.3e (relative %+.2e), linear-model %.3e" % (
            lab, dc + pc, dc, pc, F - dc - pc, (F - dc - pc) / F, lin))
    out["trials"] = res
    out["cauchy_step"] = t_c
    out["wall_seconds"] = time.time() - t0
    json.dump(out, open(OUTJ, "x"), indent=1, default=float)
    say("wrote %s (%.0f s)" % (OUTJ, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
