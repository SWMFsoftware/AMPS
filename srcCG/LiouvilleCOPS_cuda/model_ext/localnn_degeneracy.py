#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_degeneracy.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_degeneracy.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: fdb6cf306e9208678f68a8f673af5296ffd9517f6c7e64ee8dd708532427161d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: analysis (option C, lead 2026-11-14 16:20 CST): m_b vs c_j degeneracy of a final state: effective activity per
#   slope bin, signal share per bin, node dominance, compensation rho of each multiplier by the c_j (full Jacobian and
#   data rows), the global level direction (all c down, multipliers up), conditioning; P_m rows included when the fit had them.
# RUN: tests/localnn_20261114/degeneracy_run.sh (f10, f0.1); jobs/localnn_degeneracy.sbatch (P_m check, model_ext_frozen/localnn_mp_20261115T020233Z).
# ==== end of header ====
"""Option C (local non-negative basis): is a fit degenerate between the slope-bin multipliers m_b and the local
coefficients c_j? (lead 2026-11-14 16:20 CST, after the f = 10 anchored fit ended with m_80-90 = 3617, m_70-80 = 1608)

Both scale the same facets: A_e = m_b(e) (sum_j c_j phi_j(e) + TS_e). Diagnostics at the fit's final state x:
1. effective activity per bin (seen facets of the bin, area weighted): m_b, <u>_b (u = sum_j c_j phi_j), <u + TS>_b,
   <A>_b = m_b <u + TS>_b, the median of A, relative to bin 0-40; the same <A>_b for D2 and the softplus c9 state
2. the bin's share of the modelled surface NG / RG (design matrices: m_b (G^(b) c + T^(b)) / total; the 95 late
   observations with W10x5, and every used row)
3. node dominance: per node j the share of its seen support (sum over seen facets of area phi_j) in each bin; the nodes
   whose seen support is >= 50 % / >= 80 % in one bin
4. compensation of a multiplier by the local coefficients: for the column j_b = dr/d(log m_b) of the full residual vector
   (data + crossing + smoothness rows, as the fit), rho_b = |(I - P_c) j_b| / |j_b|, P_c = the projector onto the span of
   the c columns (least squares); rho_b << 1: a change of m_b is almost absorbed by the c_j (the data cannot tell them
   apart). Also with the data rows only (no smoothness), where the 2745 c columns exceed the 2084 data rows.
5. conditioning: singular values of the full Jacobian with columns scaled as scipy TRF x_scale = 'jac' (column norms),
   condition number, the 12 smallest singular values and the weight of the six multipliers, kappa and grain in their
   singular vectors; the local Gauss-Newton relative sigma of each multiplier (as localnn_report.gn_sigma).

CLI: localnn_degeneracy.py --runroot RR --design DESIGN --grid GRID --base BASE --d2 D2 --point-table PT --c9-ft C9FT
     -> RR/analysis_degeneracy/{degeneracy.json, degeneracy.txt}
"""

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402

NB = lp.NBINS


def wmean(v, w):
    return float(np.sum(v * w) / np.sum(w)) if np.sum(w) > 0 else float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("runroot", "design", "grid", "base", "d2", "point-table", "c9-ft"):
        ap.add_argument("--" + k, required=True)
    args = ap.parse_args(argv)
    t0 = time.time()
    RR = Path(args.runroot).resolve()
    seg = (RR / "fit" / "final").resolve()
    out = RR / "analysis_degeneracy"
    if out.exists():
        raise SystemExit("%s exists" % out)
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        x = np.asarray(z["x"], float)
        meta = json.loads(str(z["meta_json"]))
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        A, u, TS, seen, area = (np.asarray(z[k]) for k in ("A", "u", "TS", "seen", "area"))
        bins = np.asarray(z["bin"]).astype(np.intp)
    seen = seen.astype(bool)
    design = lp.load_design(args.design)
    grid = lp.load_grid(args.grid)
    D = lp.Data(args.base, args.d2, args.point_table, design)
    P = lp.Problem(D, design, grid, meta["templates"] == "on", meta["slope_mult_min"], design["meta"]["A_scale"],
                   meta["lambda_C"], mult_penalty_f=float((meta.get("mult_penalty") or {}).get("f_m", 0.0)))
    c, m, kappa, gA, gB = P.unpack(x)
    # 1. effective activity per bin
    with np.load(Path(args.d2) / "fit_restart.npz", allow_pickle=False) as z:
        A2 = np.maximum(np.asarray(z["activity_linear"], np.float64), 0.0)
    with np.load(Path(args.c9_ft) / "finetune_state.npz", allow_pickle=False) as z:
        A9 = np.asarray(z["A_final"], np.float64)
    eff = []
    for b in range(NB):
        sel = seen & (bins == b)
        w = area[sel]
        eff.append({"bin": lp.BIN_NAMES[b], "m": float(m[b]), "seen_area_m2": float(w.sum()), "seen_facets": int(sel.sum()),
                    "mean_u": wmean(u[sel], w), "mean_u_plus_TS": wmean(u[sel] + TS[sel], w), "mean_A": wmean(A[sel], w),
                    "median_A": float(np.median(A[sel])), "mean_A_D2": wmean(A2[sel], w), "mean_A_c9": wmean(A9[sel], w),
                    "seen_production_share_A_area": float(np.sum(A[sel] * w) / np.sum(A[seen] * area[seen]))})
    for e in eff:
        e["mean_A_rel_bin0"] = e["mean_A"] / eff[0]["mean_A"]
        e["mean_u_rel_bin0"] = e["mean_u"] / eff[0]["mean_u"]
    del A2, A9
    # 2. signal share per bin (design matrices), late window and all used rows
    shares = {}
    for g, G, T, C in (("NG", design["G_NG"], design["T_NG"], design["C_NG"]), ("RG", design["G_RG"], design["T_RG"],
                                                                                 design["C_RG"])):
        per = np.stack([m[b] * (G[:, b, :] @ c + (T[:, b] if P.templates else 0.0)) for b in range(NB)], axis=1)
        tot = per.sum(axis=1) + kappa * C
        perW = np.stack([D.E @ per[:, b] for b in range(NB)], axis=1)
        totW = D.E @ tot
        used = D.rows_ng if g == "NG" else D.rows_rg
        shares[g] = {"late_window_mean": (perW[D.late_rows] / totW[D.late_rows, None]).mean(axis=0).tolist(),
                     "all_used_mean": (perW[used] / totW[used, None]).mean(axis=0).tolist()}
    # 3. node dominance
    fn, fw = grid["facet_node"], grid["facet_w"]
    sw = np.zeros((grid["n"], NB))
    for k in range(3):
        for b in range(NB):
            sel = seen & (bins == b)
            sw[:, b] += np.bincount(fn[sel, k], weights=area[sel] * fw[sel, k], minlength=grid["n"])
    tot_sw = sw.sum(axis=1)
    has = tot_sw > 0
    frac = np.zeros_like(sw)
    frac[has] = sw[has] / tot_sw[has, None]
    dom = {lp.BIN_NAMES[b]: {"nodes_ge_50pct": int(np.sum(frac[:, b] >= 0.5)), "nodes_ge_80pct": int(np.sum(frac[:, b] >= 0.8)),
                             "max_share": float(frac[:, b].max()),
                             "c_at_0_among_ge_50pct": int(np.sum((frac[:, b] >= 0.5) & (x[P.ic] <= 1e-8)))}
           for b in range(NB)}
    # 4. compensation
    J = P.jac(x)
    nd = D.NNG + D.NRG
    Jc = J[:, P.ic]
    comp = {}
    for label, rows in (("full (data + crossing + smoothness)", slice(None)), ("data rows only", slice(0, nd))):
        Jcr = Jc[rows]
        Uc, sc, _ = np.linalg.svd(Jcr, full_matrices=False)        # SVD projector (rank-revealing)
        keep = sc > 1e-12 * sc.max()
        Qk = Uc[:, keep]
        rr = {}
        for k in range(NB - 1):
            jb = J[rows, P.im.start + k] * m[k + 1]
            res = jb - Qk @ (Qk.T @ jb)
            rr[lp.BIN_NAMES[k + 1]] = float(np.linalg.norm(res) / max(np.linalg.norm(jb), 1e-300))
        comp[label] = {"rho": rr, "c_column_rank": int(keep.sum()), "rows": int(Jcr.shape[0])}
    # 4b. the global level direction: all c_j down by a common factor, the six free multipliers up by it (A unchanged
    # on bins 40-180 except the template part; bin 0-40, m held at 1, loses activity); |J e| relative to the same
    # multiplier step alone, and the share of the change that comes from each bin
    e = np.zeros(len(x))
    e[P.ic] = -x[P.ic]
    e[P.im] = x[P.im]
    em = np.zeros(len(x))
    em[P.im] = x[P.im]
    glob = {"rel_norm_full": float(np.linalg.norm(J @ e) / np.linalg.norm(J @ em)),
            "rel_norm_data": float(np.linalg.norm(J[:nd] @ e) / np.linalg.norm(J[:nd] @ em)),
            "abs_norm_full_per_unit_log": float(np.linalg.norm(J @ e)),
            "definition": "d log c_j = -1 (all j), d log m_b = +1 (b = 40-50 .. 90-180); norms of J e over J (multipliers alone)"}
    # 5. conditioning (TRF x_scale = 'jac')
    cn = np.linalg.norm(J, axis=0)
    cn[cn == 0] = 1.0
    Js = J / cn
    U, s, Vt = np.linalg.svd(Js, full_matrices=False)
    idx = np.argsort(s)[:12]
    small = []
    for i in idx:
        v = Vt[i]
        small.append({"sigma": float(s[i]), "weight_multipliers": float(np.sum(v[P.im] ** 2)),
                      "weight_kappa_grain": float(np.sum(v[P.ik:] ** 2)), "weight_c": float(np.sum(v[P.ic] ** 2)),
                      "largest_multiplier": lp.BIN_NAMES[1 + int(np.argmax(np.abs(v[P.im])))]})
    mult_in_weak = {lp.BIN_NAMES[1 + k]: float(np.sum((Vt[idx][:, P.im.start + k]) ** 2)) for k in range(NB - 1)}
    # GN sigma (as localnn_report.gn_sigma)
    Jd = P.data_jacobian(x)
    rd = P.data_residuals(x)
    S = np.zeros((P.R.shape[0], len(x)))
    S[:, P.ic] = np.sqrt(2.0 * P.lam) * P.R
    H = Jd.T @ Jd + S.T @ S
    Hi = np.linalg.pinv(H, rcond=1e-14)
    peff = float(np.trace(Hi @ (Jd.T @ Jd)))
    v = float(rd @ rd) / max(len(rd) - peff, 1.0)
    sig = np.sqrt(np.maximum(np.diag(v * Hi), 0.0))
    gn = {lp.BIN_NAMES[1 + k]: {"m": float(m[k + 1]), "sigma": float(sig[P.im][k]), "rel": float(sig[P.im][k] / m[k + 1])}
          for k in range(NB - 1)}
    rep = {"created_utc": lp.utc_now(), "host": socket.gethostname(), "command": sys.argv, "runroot": str(RR),
           "segment": str(seg), "lambda_C": meta["lambda_C"], "costs": meta["costs"], "effective_activity": eff,
           "signal_share_per_bin": shares, "node_dominance": dom, "compensation_rho": comp, "global_level_direction": glob,
           "conditioning": {"columns": int(J.shape[1]), "rows": int(J.shape[0]), "sigma_max": float(s.max()),
                            "sigma_min": float(s.min()), "condition_number": float(s.max() / s.min()),
                            "smallest": small, "multiplier_weight_in_the_12_weakest": mult_in_weak},
           "gn_sigma_multipliers": gn, "p_eff": peff, "code_sha256": lp.sha256_file(__file__),
           "wall_seconds": time.time() - t0}
    stage = RR / ("analysis_degeneracy.incomplete_%s" % lp.utc_now().replace(":", "").replace("-", "")[:15])
    stage.mkdir()
    with open(stage / "degeneracy.json", "x") as f:
        json.dump(rep, f, indent=1, default=float)
    L = ["degeneracy diagnostics %s (lambda_C %.6e, data cost %.10e)" % (RR.name, meta["lambda_C"], meta["costs"]["data_cost"]),
         "#bin m <u>_seen <u+TS>_seen <A>_seen <A>/<A>_0-40 median_A <A>_D2 <A>_c9 share_seen_prod late_NG late_RG "
         "rho_full rho_data nodes>=50% gn_rel_sigma"]
    for b, e in enumerate(eff):
        nm = e["bin"]
        L.append("%s %.4g %.4e %.4e %.4e %.4g %.4e %.4e %.4e %.4f %.4f %.4f %s %s %d %s" % (
            nm, e["m"], e["mean_u"], e["mean_u_plus_TS"], e["mean_A"], e["mean_A_rel_bin0"], e["median_A"], e["mean_A_D2"],
            e["mean_A_c9"], e["seen_production_share_A_area"], shares["NG"]["late_window_mean"][b],
            shares["RG"]["late_window_mean"][b],
            "%.3e" % comp["full (data + crossing + smoothness)"]["rho"][nm] if b else "-",
            "%.3e" % comp["data rows only"]["rho"][nm] if b else "-", dom[nm]["nodes_ge_50pct"],
            "%.3g" % gn[nm]["rel"] if b else "fixed"))
    L.append("condition number (x_scale jac) %.3e; 12 smallest singular values %s" % (
        rep["conditioning"]["condition_number"], ["%.2e" % d["sigma"] for d in small]))
    L.append("multiplier weight in the 12 weakest singular vectors %s" % json.dumps({k: round(v, 4) for k, v in mult_in_weak.items()}))
    L.append("global level direction (all c down, six multipliers up by the same factor): |J e| / |J e_m| = %.3e (full), %.3e "
             "(data rows)" % (glob["rel_norm_full"], glob["rel_norm_data"]))
    L.append("note: gn_rel_sigma = local Gauss-Newton sigma / m (pinv of a badly conditioned H; indicative only)")
    with open(stage / "degeneracy.txt", "x") as f:
        f.write("\n".join(L) + "\n")
    stage.rename(out)
    print("\n".join(L), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
