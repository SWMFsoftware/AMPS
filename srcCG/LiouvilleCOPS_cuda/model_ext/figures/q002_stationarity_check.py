#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_stationarity_check.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_stationarity_check.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b84145474677b03bc8be0cf826737acb0ecce348b54ecd7aa08c0f82bf08daa6
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: is an option C end state a minimum or an early ftol stop; evaluations only
# RUN: jobs/q002_stationarity_check.sbatch
# CHANGED: 2026-11-22 states with classes held at a value (cliff test)
# ==== end of header ====
"""Is an option C end state a minimum of its problem, or did ftol stop the optimizer early? (2026-11-22, multi-start
refits: the first one stopped after 6 evaluations.) Evaluations only, no optimizer.

Problem: the reference's (C f10, P_m 0.01, templates off, six multipliers >= 1), its frozen code. At each state x:
  * gradient g = J^T r, per column norm (cosine between the residual and that column): for a parameter at its bound only
    a negative g allows descent
  * the bounded Gauss-Newton subproblem min 0.5 |J d + r|^2, d >= lo - x (columns scaled to norm 1; lsq_linear): the
    decrease the linearized model promises
  * the true objective along that step (alpha 1 ... 1e-4): the decrease that is really available in that direction
  * with --solver nnls also SECOND order: the exact Hessian on the free set (the parameters the step leaves off their
    bounds), checked against differences of the gradient. No descent + smallest eigenvalue > 0 = a local minimum; the
    Gauss-Newton test alone cannot tell (ln m has negative curvature, m x c is bilinear).
With --own-problem a control fit is tested in its own problem (its held multipliers and lower bound).
A state with an available decrease far above ftol x objective is not a converged minimum, whatever the fit printed.

    q002_stationarity_check.py --out <new folder> name=<localnn_state.npz> [...]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear, nnls

R = Path("/data001/highresolu_vt/cops_output")
sys.path.insert(0, str(R / "model_ext_frozen/localnn_q002_20261118T162327Z"))
import localnn_problem as lp       # noqa: E402

N = R / "runs/shap8_012m_lategrid_dz1.5"
G = N / "nonneg/localnn/grid_150m"
ALPHAS = (1.0, 0.5, 0.25, 0.1, 0.03, 0.01, 3e-3, 1e-3, 1e-4)


def second_order(P, x, r, J, free):
    """Exact Hessian of 0.5 r.r on the free set, columns scaled to norm 1; (smallest eigenvalues, check against FD)."""
    D = P.D
    pN, pR, _ = P.pressures_used(x)
    nN, nd = len(pN), len(pN) + len(pR)
    w = np.concatenate([np.broadcast_to(D.wNG, pN.shape), np.broadcast_to(D.wRG, pR.shape)]).astype(float)
    Jf = J[:, free]
    weight = np.ones(len(r))
    weight[:nd] -= r[:nd] / w                                        # from d2(ln p) = -(dp/p)^2
    H = Jf.T @ (weight[:, None] * Jf)
    pos = np.full(len(x), -1)
    pos[free] = np.arange(int(free.sum()))
    sN = r[:nN] * np.broadcast_to(D.wNG, pN.shape) / pN              # from d2 p: only c x m
    sR = r[nN:nd] * np.broadcast_to(D.wRG, pR.shape) / pR
    ic = pos[P.ic]
    for b in range(1, lp.NBINS):
        k = P.im.start + b - 1
        if free[k]:
            cross = P.A_s * (P.GN[b].T @ sN + P.GR[b].T @ sR)
            H[ic[ic >= 0], pos[k]] += cross[ic >= 0]
            H[pos[k], ic[ic >= 0]] += cross[ic >= 0]
            if P.fm > 0:                                             # from d2(ln m) = -1/m^2
                H[pos[k], pos[k]] += P.wm * np.log(x[k]) * (-P.wm / x[k] ** 2)
    n = np.linalg.norm(Jf, axis=0)
    Hs = H / np.outer(n, n)
    v = np.random.default_rng(1).standard_normal(len(n))
    v /= np.linalg.norm(v)
    d = np.zeros(len(x))
    d[free] = v / n
    d[free & (np.abs(d) * 1e-5 > x - P.lo)] = 0.0                   # both difference points inside the bounds
    best = np.inf
    for h in (1e-5, 1e-6, 1e-7):
        xp, xm = np.maximum(x + h * d, P.lo), np.maximum(x - h * d, P.lo)
        num = ((P.jac(xp).T @ P.fun(xp)) - (P.jac(xm).T @ P.fun(xm)))[free] / (2 * h)
        an = H @ d[free]
        best = min(best, float(np.linalg.norm((num - an) / n) / np.linalg.norm(an / n)))
    return np.linalg.eigvalsh(Hs)[:3].tolist(), best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-iter", type=int, default=40)
    ap.add_argument("--solver", choices=("trf", "nnls"), default="trf",
                    help="nnls: the subproblem exactly (active set; every bound is a finite lower bound), for states "
                         "where the interior trf reaches max_iter")
    ap.add_argument("--own-problem", action="store_true",
                    help="a control fit (meta slope_test) is tested with its held multipliers and its lower bound")
    ap.add_argument("states", nargs="+")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir()
    design = lp.load_design(G / "design_h2oqmin0.02_d2q002")
    D = lp.Data(N / "base_q002", N / "clip_refit/D2_q002_warm/seg01", N / "late_grid_points.txt", design)
    grid = lp.load_grid(G)
    P = lp.Problem(D, design, grid, False, 1.0, design["meta"]["A_scale"], 8.249495752351283e-05, mult_penalty_f=0.01)
    rec = {}
    lo0, fm0, wm0 = P.lo.copy(), P.fm, P.wm
    for item in a.states:
        name, path = item.split("=", 1)
        t0 = time.time()
        with np.load(path, allow_pickle=False) as z:
            x = np.asarray(z["x"], float)
            test = json.loads(str(z["meta_json"])).get("slope_test") if a.own_problem else None
        P.lo, held = lo0.copy(), []
        P.fm, P.wm = fm0, wm0
        if test:
            P.lo[P.im] = test["slope_mult_min"]
            at_value = test.get("held_values", {})           # cliff test: held at a value below 1, fitted with f_m 0
            held = [P.im.start + list(lp.BIN_NAMES[1:]).index(b) for b in list(test["held_at_1"]) + list(at_value)]
            if at_value:
                P.fm, P.wm = 0.0, 0.0
                P.lo[held] = x[held]
        x = np.maximum(x, P.lo)
        r, J = P.fun(x), P.jac(x)
        f0 = 0.5 * float(r @ r)
        nr = np.linalg.norm(J, axis=0)
        use = nr > 0
        use[held] = False
        cos = np.zeros(len(x))
        cos[use] = (J[:, use].T @ r) / (nr[use] * np.linalg.norm(r))
        S, low = J[:, use] / nr[use], (P.lo[use] - x[use]) * nr[use]
        if a.solver == "nnls":                       # y = d - low >= 0
            y, rnorm = nnls(S, -(r + S @ low), maxiter=20 * S.shape[1])
            sol = argparse.Namespace(x=y + low, cost=0.5 * rnorm * rnorm, status=1, nit=-1)
        else:
            sol = lsq_linear(S, -r, bounds=(low, np.inf), method="trf", lsq_solver="exact", tol=1e-12,
                             max_iter=a.max_iter)
        d = np.zeros(len(x))
        d[use] = sol.x / nr[use]
        pred = f0 - float(sol.cost)
        trials = []
        for al in ALPHAS:
            rt = P.fun(np.maximum(x + al * d, P.lo))
            trials.append([al, f0 - 0.5 * float(rt @ rt)])
        best = max(trials, key=lambda t: t[1])
        m = x[P.im]
        rec[name] = dict(state=path, solver=a.solver, conclusive=bool(sol.status > 0), objective=f0, predicted_decrease_rel=pred / f0, available_decrease_rel=best[1] / f0,
                         alpha=best[0], trials=trials, lsq_status=int(sol.status), lsq_iterations=int(sol.nit),
                         multipliers=m.tolist(), multiplier_gradient_cos=cos[P.im].tolist(),
                         multiplier_gn_step=d[P.im].tolist(), seconds=time.time() - t0)
        if a.solver == "nnls":
            free = np.zeros(len(x), bool)
            free[np.where(use)[0][y > 0]] = True
            eig, fd = second_order(P, x, r, J, free)
            rec[name].update(problem="own" if test else "all free >= 1", free=int(free.sum()),
                             hessian_smallest_eigenvalues=eig, hessian_vs_fd=fd,
                             local_minimum=bool(best[1] / f0 < 1e-8 and pred / f0 < 1e-8 and eig[0] > 0))
            print("%s: SECOND ORDER free %d, smallest eigenvalues %s, Hessian vs FD %.1e, local minimum %s" % (
                name, free.sum(), " ".join("%.3e" % e for e in eig), fd, rec[name]["local_minimum"]), flush=True)
        print("%s: objective %.10e | Gauss-Newton promises %.3e, available on that step %.3e (alpha %g), relative; "
              "lsq status %d, %d iterations, %.0f s" % (name, f0, pred / f0, best[1] / f0, best[0], sol.status, sol.nit,
                                                        time.time() - t0), flush=True)
        print("   m    " + " ".join("%11.4g" % v for v in m))
        print("   cos  " + " ".join("%11.3e" % v for v in cos[P.im]))
        print("   step " + " ".join("%11.3e" % v for v in d[P.im]), flush=True)
    (out / "stationarity.json").write_text(json.dumps(rec, indent=1) + "\n")


if __name__ == "__main__":
    main()
