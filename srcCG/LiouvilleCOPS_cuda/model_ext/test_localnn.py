#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_localnn.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_localnn.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 962ac4213c8056517e53d30688805695fab8074ff270d9296c8e02335b320d05
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tests (option C, CPU): partition of unity, stiffness, design vs direct sums (all rows; independent 4 rows),
#   fast vs all-rows path, Jacobian vs finite differences, bounds, NNLS projection, D2 / c9 reproduction.
# RUN: jobs/localnn_build.sbatch STEP tests (frozen model_ext_frozen/localnn_20261114T081912Z).
# ==== end of header ====
"""Tests of option C (local non-negative basis; localnn_mesh.py, localnn_design.py, localnn_problem.py,
localnn_fit.py). CPU only (a compute node: the independent direct sums read four kernel files of the base).

    test_localnn.py --grid GRID --design DESIGN --base BASE --d2 D2 --point-table PT --c9-ft C9FT --c9-fit C9FIT
                    --contrib-code FROZEN_CONTRIB [--start-file START_DIR] [--out report.json]

 1 partition of unity: sum_j phi_j(e) = 1 (<= 4 eps) and phi >= 0 at every facet; node indices valid; every node has
   support (mass row sum > 0); 1^T M 1 = total facet area
 2 stiffness: Lg 1 = 0, Lg symmetric PSD, Dg^T Dg = Lg; for a linear function u = a . x at the nodes the energy equals
   sum_T area_T |a - (a . n_T) n_T|^2 (the projected gradient)
 3 design route vs direct facet sums: (a) the two test activities of the design pass at every row (<= 1e-12);
   (b) an INDEPENDENT activity (new seed, with templates) at four rows (early observation, last pre-gap observation,
   a late observation, a grid row) by the frozen record_contributions, kernel files read here (<= 1e-12)
 4 fast path (effective used rows) vs all-rows path (design -> W10x5 -> used rows): data residuals (<= 1e-13)
 5 Jacobian vs central finite differences (24 random node columns, every multiplier, kappa, grain A, B; step 1e-4
   max(|x|, 1) for c and m, 1e-4 max(|x|, 1e-3) for kappa and grain) at the start (95 crossing rows active) and at a
   state with multipliers 1.2 .. 4 and grain A 12, B 0.4 (away from the grain's one-sided kink), lambda_C > 0
   (max rel <= 1e-5)
 6 bounds: lower bounds (c 0, multipliers 1, kappa 0, grain 0), a 3-evaluation TRF run stays inside them
 7 start: the NNLS projection of a function in the span (c >= 0) recovers c (rel <= 1e-8)
 8 reproduction (localnn_fit.verification): D2's and the c9 final state's data costs (rel < 1e-9) and the printed pNG,
   pRG of their tables
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402
import localnn_fit as lf          # noqa: E402
import localnn_design as ld       # noqa: E402
import localnn_mesh as lm         # noqa: E402

RESULTS = []


def check(name, ok, **info):
    RESULTS.append({"test": name, "passed": bool(ok), **info})
    print("%s %s %s" % ("PASS" if ok else "FAIL", name, json.dumps(info, default=float)), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    for k in ("grid", "design", "base", "d2", "point-table", "c9-ft", "c9-fit", "contrib-code"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--start-file")
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    t0 = time.time()
    grid = lp.load_grid(args.grid)
    design = lp.load_design(args.design)
    el = lm.read_elements(Path(args.base) / "elements.bin")
    area = el["area"]
    n = grid["n"]
    fn, fw = grid["facet_node"], grid["facet_w"]
    # 1
    s = fw.sum(axis=1)
    rowsum = np.asarray(grid["M"].sum(axis=1)).ravel()
    check("1 partition of unity", np.max(np.abs(s - 1)) <= 4 * np.finfo(float).eps and fw.min() >= 0 and fn.min() >= 0
          and fn.max() < n and rowsum.min() > 0 and abs(rowsum.sum() / area.sum() - 1) < 1e-12,
          max_abs_sum_minus_1=float(np.max(np.abs(s - 1))), min_phi=float(fw.min()), min_support=float(rowsum.min()),
          mass_total_rel=float(rowsum.sum() / area.sum() - 1))
    # 2
    Lg, Dg = grid["Lg"], grid["Dg"]
    Ld = Lg.toarray()
    ev = np.linalg.eigvalsh(0.5 * (Ld + Ld.T))
    pos, T = grid["node_xyz"], grid["tri"]
    a = np.array([0.3, -0.7, 0.2])
    u = pos @ a
    nrm = np.cross(pos[T[:, 1]] - pos[T[:, 0]], pos[T[:, 2]] - pos[T[:, 0]])
    nh = nrm / np.linalg.norm(nrm, axis=1)[:, None]
    proj = a[None, :] - (nh @ a)[:, None] * nh
    e_ref = float(np.sum(grid["tri_area"] * np.sum(proj ** 2, axis=1)))
    e = float(u @ (Lg @ u))
    check("2 stiffness", np.max(np.abs(Lg @ np.ones(n))) < 1e-12 * np.max(np.abs(Ld)) and
          np.max(np.abs(Ld - Ld.T)) == 0 and ev.min() > -1e-10 * ev.max() and
          abs(Dg.T @ Dg - Lg).max() < 1e-12 * np.max(np.abs(Ld)) and abs(e / e_ref - 1) < 1e-10,
          Lg1=float(np.max(np.abs(Lg @ np.ones(n)))), min_eig=float(ev.min()), max_eig=float(ev.max()),
          linear_energy_rel=float(e / e_ref - 1))
    # 3a
    tc = design["meta"]["design_vs_direct_test"]
    check("3a design vs direct (test activities, all rows)", all(v <= 1e-12 for v in tc.values()), **tc)
    # 3b independent activity at four rows
    D = lp.Data(args.base, args.d2, args.point_table, design)
    with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
        TS = np.asarray(z["TS"], np.float64)
        bins = np.asarray(z["bin"]).astype(np.intp)
    rng = np.random.default_rng(7)
    A_s = design["meta"]["A_scale"]
    c = A_s * rng.uniform(0.0, 3.0, n) * (rng.uniform(size=n) > 0.3)
    m = np.concatenate([[1.0], rng.uniform(1.0, 8.0, 6)])
    kappa = 0.37
    A = m[bins] * (lp.facet_field(grid, c) + TS)
    late_pre = int(D.rows_ng[np.flatnonzero(D.et[D.rows_ng] <= D.et_before)[-1]])
    rows = [0, late_pre, int(D.late_rows[50]), int(D.grid_rows[200])]
    S = ld.contrib_setup(args.base, args.contrib_code, lp.USER_T0)
    sN, sR = lp.all_rows_surface(design, c, m, kappa, True)
    worst = 0.0
    det = {}
    for r in rows:
        dn, dr = ld.direct_row(S, args.base, r, int(D.ipoint[r]), A, kappa)
        en, er = abs(sN[r] / dn - 1), abs(sR[r] / dr - 1)
        det["row%d" % r] = [en, er]
        worst = max(worst, en, er)
    check("3b design vs independent direct facet sum (4 rows)", worst <= 1e-12, max_rel=worst, rows=det)
    # 4 fast path vs all rows; 5 Jacobian; 6 bounds
    P = lp.Problem(D, design, grid, True, 1.0, A_s, smooth_lambda=1e-6)
    x = np.concatenate([c / A_s, m[1:], [kappa, 6.9, 1e-3]])
    r_fast = P.data_residuals(x)
    rows_all = D.model_rows(sN, sR, 6.9, 1e-3)
    r_all = D.costs(rows_all["pNG"], rows_all["pRG"])["r"]
    check("4 fast path vs all-rows path", np.max(np.abs(r_fast - r_all)) <= 1e-13,
          max_abs=float(np.max(np.abs(r_fast - r_all))))
    if args.start_file:
        with np.load(Path(args.start_file) / "start.npz", allow_pickle=False) as z:
            xs = np.asarray(z["x0"], float)
            tmpl = str(z["templates"]) == "on"
        P2 = lp.Problem(D, design, grid, tmpl, 1.0, A_s, smooth_lambda=1e-6)
        xs = np.minimum(np.maximum(xs, P2.lo + 1e-9), P2.hi)
    else:
        P2, xs = P, x
    xt = xs.copy()
    xt[P2.im] = np.array([1.5, 2.0, 3.0, 1.2, 4.0, 2.5])
    xt[P2.iA], xt[P2.iB] = 12.0, 0.4
    for label, xx in (("start", xs), ("multipliers 1.2..4, grain A 12, B 0.4", xt)):
        P2.deadline = None
        rr = P2.data_residuals(xx)
        active = int(np.sum(D.crossing(rr)[1] > 0))
        chk = lf.jacobian_check(P2, xx)
        check("5 Jacobian vs finite differences at %s" % label, chk["passed"], max_rel=chk["max_rel_error"],
              columns=chk["columns"], crossing_rows_active=active)
    lo_ok = (np.all(P.lo[P.ic] == 0) and np.all(P.lo[P.im] == 1.0) and P.lo[P.ik] == 0 and P.lo[P.iA] == 0 and
             P.lo[P.iB] == 0 and np.all(np.isinf(P.hi)))
    from scipy.optimize import least_squares
    P.nfev, P.best, P.deadline = 0, None, None
    x0 = np.minimum(np.maximum(x, P.lo), P.hi)
    x0[P.im] = 1.0
    res = least_squares(P.fun, x0, jac=P.jac, bounds=(P.lo, P.hi), x_scale="jac", method="trf", tr_solver="exact",
                        max_nfev=3)
    check("6 bounds", lo_ok and np.all(res.x >= P.lo) and np.all(res.x <= P.hi), lower_bounds_ok=lo_ok,
          min_c=float(res.x[P.ic].min()), min_m=float(res.x[P.im].min()), nfev=int(res.nfev))
    # 7 NNLS projection recovers a function of the span
    c_true = A_s * rng.uniform(0.0, 1.0, n) * (rng.uniform(size=n) > 0.5)
    c_rec, rep = lp.nnls_projection(grid, area, lp.facet_field(grid, c_true))
    rel = float(np.max(np.abs(c_rec - c_true)) / np.max(c_true))
    check("7 NNLS projection of a function in the span", rel <= 1e-8, max_rel=rel, nodes_at_zero=rep["nodes_at_zero"])
    # 8 reproduction
    ns = argparse.Namespace(design=args.design, base=args.base, d2=args.d2, point_table=args.point_table,
                            c9_ft=args.c9_ft, c9_fit=args.c9_fit)
    v = lf.verification(ns)
    for name in ("D2", "c9"):
        vv = v[name]
        check("8 reproduction %s" % name, vv["passed"], data_cost=vv["data_cost"], reported=vv["reported_data_cost"],
              rel=vv["data_cost_rel_diff"], tables={k: (vv[k]["max_rel"], vv[k]["printed_differences"]) for k in vv
                                                    if isinstance(vv[k], dict) and "max_rel" in vv[k]})
    ok = all(r["passed"] for r in RESULTS)
    print("%d tests, %d failed (%.0f s)" % (len(RESULTS), sum(not r["passed"] for r in RESULTS), time.time() - t0))
    if args.out:
        with open(args.out, "x") as f:
            json.dump({"results": RESULTS, "passed": ok, "argv": sys.argv, "code": lf.code_record()}, f, indent=1,
                      default=float)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
