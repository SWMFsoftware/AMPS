#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_slopetest_state_check.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_slopetest_state_check.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3f39699f36767800baa2ea3ed015f4cfe95472b2d15d629a5149b9bc7df273b3
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: independent check of a control fit's state in the reference's problem: old code, direct sums, paths (3c)
# RUN: jobs/q002_slopetest_state_check.sbatch
# ==== end of header ====
"""Independent check that a control fit's saved state beats the reference in the REFERENCE's own problem (2026-11-21;
the user: "where the error is found?"). Evaluations only, no optimizer.

Uses only the reference's frozen code (localnn_q002_20261118T162327Z: six multipliers free, >= 1), so the --slope-test
change is not involved:
  1 both states are inside that problem's bounds; their objective and its terms by that code's residual vector
  2 the candidate's data cost by direct facet sums over all kernel records (a route without the design matrices)
  3 the objective along two paths from the reference to the candidate (is there a barrier between the two solutions?)
  4 what the candidate looks like: switched-off share (decision 1b), Q, its seen share, the share above the energy limit

    q002_slopetest_state_check.py --checkpoint <fit/checkpoint.npz> --out <new folder> [--procs 16]
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import shutil
import sys
from pathlib import Path

import numpy as np

R = Path("/data001/highresolu_vt/cops_output")
OLD = R / "model_ext_frozen/localnn_q002_20261118T162327Z"
sys.path.insert(0, str(OLD))
import localnn_design as ld        # noqa: E402
import localnn_problem as lp       # noqa: E402

sys.path.append(str(Path(__file__).resolve().parent))
import q002_slope_activity as Q    # noqa: E402

N = R / "runs/shap8_012m_lategrid_dz1.5"
L = N / "nonneg/localnn"
REF = L / "stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit/final"
DESIGN = L / "grid_150m/design_h2oqmin0.02_d2q002"
BASE = N / "base_q002"
AA = R / "figures_cluster/12M/q002_decisions_20261121/active_area.json"
F_E = 1.0936266398665254e21                          # q002_plausibility/metrics.json, 10:33:10 UTC
LAM, FM = 8.249495752351283e-05, 0.01
_S = {}


def work(task):
    r, ip = task
    return r, ld.direct_row(_S["CS"], BASE, r, ip, _S["A"], _S["kappa"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir()
    ck = out / "checkpoint_copy.npz"
    shutil.copy2(a.checkpoint, ck)                   # the running fit replaces its checkpoint
    with np.load(ck, allow_pickle=False) as z:
        xb, saved, nfev = np.asarray(z["x"], float), float(z["objective"]), int(z["nfev_total"])
    with np.load(REF / "localnn_state.npz", allow_pickle=False) as z:
        xr = np.asarray(z["x"], float)
    design = lp.load_design(DESIGN)
    D = lp.Data(BASE, N / "clip_refit/D2_q002_warm/seg01", N / "late_grid_points.txt", design)
    grid = lp.load_grid(L / "grid_150m")
    P = lp.Problem(D, design, grid, False, 1.0, design["meta"]["A_scale"], LAM, mult_penalty_f=FM)

    def objective(x):
        r = P.fun(x)
        return 0.5 * float(r @ r)

    rec = {"code": str(OLD), "checkpoint": a.checkpoint, "checkpoint_nfev": nfev, "checkpoint_objective": saved,
           "checkpoint_sha256": hashlib.sha256(ck.read_bytes()).hexdigest()}
    with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
        seen, bins = np.asarray(z["seen"], bool), np.asarray(z["bin"]).astype(int)
    nf, cen, area = Q.elements()
    del cen
    T = Q.map_times()
    label = [k for k, v in T.items() if v["utc"].endswith("10:33:10")][0]
    q = Q.qfac(label, T[label]["rh_m"], nf)[0]
    thr = json.loads(AA.read_text())["T"]
    states = {}
    for name, x in (("reference", xr), ("candidate", xb)):
        c, m, kappa, gA, gB = P.unpack(x)
        cost, _ = P.evaluate(x)
        u = lp.facet_field(grid, c)
        A = m[bins] * np.maximum(u, 0.0)
        Fa = A * q * area
        off = seen & (A <= thr)
        states[name] = dict(A=A, kappa=kappa, gA=gA, gB=gB)
        rec[name] = {"inside_reference_bounds": bool(np.all(x >= P.lo) and np.all(x <= P.hi)),
                     "objective": objective(x), "data_cost": cost["data_cost"], "crossing": cost["crossing_cost"],
                     "smoothness": cost["penalty"], "P_m": P.mult_penalty(x), "multipliers": m.tolist(),
                     "late_rms_NG": cost["late_rms_ln_NG"], "late_rms_RG": cost["late_rms_ln_RG"],
                     "rms_NG": cost["rms_ln_NG"], "rms_RG": cost["rms_ln_RG"],
                     "off_seen_1b": float(area[off].sum() / area[seen].sum()),
                     "Q_103310": float(Fa.sum()), "Q_seen_share": float(Fa[seen].sum() / Fa.sum()),
                     "Q_share_over_F_E": float(Fa[A * q > F_E].sum() / Fa.sum()),
                     "F_max_over_F_E": float((A * q).max() / F_E),
                     "seen_production_share_by_bin": [float(Fa[seen & (bins == b)].sum() / Fa[seen].sum()) for b in range(7)]}
        print(name, json.dumps(rec[name]), flush=True)

    s = np.linspace(0.0, 1.0, 101)
    im = P.im

    def straight(t):
        return (1.0 - t) * xr + t * xb

    def geometric(t):
        x = straight(t)
        x[im] = np.exp((1.0 - t) * np.log(xr[im]) + t * np.log(xb[im]))
        return x

    o0 = objective(xr)
    rec["paths"] = {}
    for name, path in (("straight", straight), ("multipliers_geometric", geometric)):
        v = np.array([objective(path(t)) for t in s])
        k = int(v.argmax())
        rec["paths"][name] = {"s": s.tolist(), "objective": v.tolist(), "max_over_reference": float(v[k] / o0),
                              "s_at_max": float(s[k]), "change_at_s_1e-4": objective(path(1e-4)) - o0}
        print("path %s: max %.4f x reference at s %.2f; change at s=1e-4 %+.3e" % (
            name, v[k] / o0, s[k], rec["paths"][name]["change_at_s_1e-4"]), flush=True)

    b = states["candidate"]
    CS = ld.contrib_setup(BASE, R / "model_ext_frozen/nonneg_contrib_20261118T072729Z", lp.USER_T0)
    _S.update(CS=CS, A=b["A"], kappa=b["kappa"])
    sN, sR = np.zeros(D.nrow), np.zeros(D.nrow)
    with mp.get_context("fork").Pool(a.procs) as pool:
        for r, (n_, r_) in pool.imap_unordered(work, [(r, int(D.ipoint[r])) for r in range(D.nrow)], chunksize=1):
            sN[r], sR[r] = n_, r_
    rows = D.model_rows(sN, sR, b["gA"], b["gB"])
    direct = D.costs(rows["pNG"], rows["pRG"])["data_cost"]
    c = rec["candidate"]
    rec["direct"] = {"data_cost": direct, "rel_diff": abs(direct / c["data_cost"] - 1.0)}
    rec["confirmed"] = bool(c["inside_reference_bounds"] and rec["reference"]["inside_reference_bounds"]
                            and c["objective"] < rec["reference"]["objective"] and rec["direct"]["rel_diff"] < 1e-10
                            and abs(c["objective"] / saved - 1.0) < 1e-10)
    (out / "state_check.json").write_text(json.dumps(rec, indent=1, default=float) + "\n")
    print("direct", rec["direct"], "objective vs saved %.1e" % abs(c["objective"] / saved - 1.0), "CONFIRMED", rec["confirmed"])
    return 0 if rec["confirmed"] else 1


if __name__ == "__main__":
    sys.exit(main())
