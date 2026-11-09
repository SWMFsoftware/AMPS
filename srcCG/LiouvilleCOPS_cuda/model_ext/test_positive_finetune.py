#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_positive_finetune.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_positive_finetune.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 00e384b0dbe66f1be7a14595b1bee66b5e4df506a082699dab99aa7835d33959
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. CPU tests of positive_finetune.py (W10x5 adjoint, grain/crossing vs grain_crossing, bounds, adjoint and
#   FD on a CPU stand-in, zero-sensitivity optimizer path, wall-time stop, writers, positive_table.py). user decisions 2026-11-07/08.
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
"""CPU tests of model_ext/positive_finetune.py (per-element fine tuning of the positive-activity fits; 2026-11-07).

    .venv/bin/python -W ignore test_positive_finetune.py [--scratch DIR]      (login node: < 1 min, < 1 GB)

The GPU part (forward = the unchanged evaluator, J v = eval_fused_m_c<k>, J^T u = kernels/finetune.cu ft_jtu_m) is
validated on the real 12M data by positive_finetune.py --validate (checks 1-6, jobs/positive_finetune.sbatch).
Here a CPU stand-in with the same record semantics (exp activity, per-record H2O floor F = max(q A, F_min), records at
the floor without activity derivative) replaces the GPU engine, so everything on the host side is tested:
 1. Rows.window_T is the adjoint of positive_problem.WindowOperator.apply (W10x5 operator of positive_problem);
 2. Rows.grain == grain_crossing's GrainProblem.grain (the module installed on a stand-in fitter);
 3. Rows.penalty and the crossing Jacobian rows == grain_crossing's own fun_aug / jac_aug (its least_squares wrapper
    driven by a stub), including active (wrong-side) pairs;
 4. lower_bounds: min(0, ln(r A_orig / A)), none for A_orig = 0, the start-below-bound mask;
 5. Objective/Jacobian with the stand-in engine, lambda 0 and 0.3: adjoint <J v, u> == <v, J^T u>, J v and J^T u
    against central differences (random directions, single facets, grain), cost decomposition (data, Tikhonov =
    lambda/2 sum area (s - s0)^2 / sum area, crossing);
 6. optimize(): least_squares TRF/LSMR with bounds decreases the cost and respects the bounds; a wall-time stop
    returns the best evaluated state;
 7. writers: COPS_model.dat (only pNG/pRG tokens change, note line inserted), late_window_model.dat,
    pressure_comparison.txt, finetune_state.npz (held_next = held_prev + delta, s_final = s_fit + held_next), and
    figures_cluster/scripts/positive_table.py runs on the output;
 8. the worker request handler: ft_* requests are added to cops_cuda.transport.WorkerState.handle, others go to the
    original handler; importing the module installs nothing;
 9. FACETS late_only (2026-11-08): the parameter set equals a brute-force first-record classification over the
    objective points, and a random change of its parameters leaves every row before the window start bitwise
    unchanged (strict support) while rows inside change;
10. --facets is required (no silent default), late_only requires --window-min;
11. zero sensitivity (2026-11-08, the lambda-0 cycle-1 case): every parameter facet at the H2O floor -> J v of the
    facet block and J^T u on it exactly 0, zero_sensitivity() counts all of them, the optimizer (TRF/LSMR, zero
    columns, x_scale 1, delta at its lower bound 0 with zero gradient) runs, leaves delta at the start (lambda 0)
    and refits the grain; with lambda > 0 the Tikhonov term alone moves delta.
12. CO2 flux floor (2026-11-08): check_environment (in subprocesses with positive_activity installed) accepts an
    earlier fit (no COPS_POS_FLOOR_CO2 recorded) only without the variable and a ratio-1e-5 fit (recorded "1") only with
    it, the configurations' co2_floor must agree; the job script's environment snippet exports the variable only when
    the fit recorded it; the H2O code that kernels/finetune.cu ft_jtu_m transposes (static_factors_m, the H2O floor and
    the activity/template column code of fused_body_m) is the same text in the 2026-11-07 pressure.cu (frozen c1
    snapshot) and the current one, and PRESSURE_COMPATIBLE_SHA256 lists the current pressure.cu.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
for p in (str(HERE), str(APP / "tools"), str(APP)):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("COPS_OUTPUT_ROOT", "/data001/highresolu_vt/cops_output")
POSITIVE_TABLE = Path("/data001/highresolu_vt/cops_output/figures_cluster/scripts/positive_table.py")
RESULTS = []


def report(name, ok, detail):
    RESULTS.append((name, bool(ok), detail))
    print("%-62s %s  %s" % (name, "PASS" if ok else "FAIL", detail), flush=True)


# ------------------------------------------------------------------------------------------------ stand-ins
class CpuEngine:
    """CPU stand-in of positive_finetune.GpuEngine with the GPU kernels' record semantics."""

    def __init__(self, n, nf, rng, A_orig, floor_ratio, kmax):
        p, f = [], []
        for i in range(n):
            k = int(rng.integers(40, 300))
            p.append(np.full(k, i))
            f.append(rng.choice(nf, size=k, replace=False))
        self.p, self.f = np.concatenate(p), np.concatenate(f)
        R = len(self.p)
        self.n, self.nf = n, nf
        self.aN = rng.uniform(0, 1, R) * (rng.uniform(size=R) > 0.15)
        self.aR = rng.uniform(0, 1, R) * (rng.uniform(size=R) > 0.15)
        self.qmax = rng.uniform(1.0, 2.0, n)
        self.q = rng.uniform(0.02, 1.0, R) * self.qmax[self.p]
        self.fmin = self.qmax[self.p] * floor_ratio * A_orig[self.f]
        self.cN = rng.uniform(0.5, 1.5, n) * 2e-8
        self.cR = rng.uniform(0.5, 1.5, n) * 2e-8
        self.bN = rng.uniform(0.5, 1.5, n) * 1e-7
        self.bR = rng.uniform(0.5, 1.5, n) * 1e-7
        self.kmax = kmax
        self.idx = None
        self.counter = 0
        self.reuploads = 0
        self.worker_seconds = {"ft_eval": [], "ft_jv": [], "ft_jtu": []}

    def forward(self, s):
        A = np.exp(s)
        F = np.maximum(self.q * A[self.f], self.fmin)
        pNG = self.cN * np.bincount(self.p, self.aN * F, minlength=self.n) + self.bN
        pRG = self.cR * np.bincount(self.p, self.aR * F, minlength=self.n) + self.bR
        self.counter += 1
        return self.counter, pNG, pRG

    def set_parameters(self, idx):
        self.idx = np.asarray(idx, np.int64)
        self.v = np.zeros((self.kmax, len(self.idx)))          # parameter space, as GpuEngine
        return [{"param_faces": len(self.idx)}]

    def direction_buffer(self, k):
        return self.v[:k]

    def _facets(self, vp):
        vf = np.zeros(self.nf)
        vf[self.idx] = vp
        return vf

    def _on(self, s):
        qa = self.q * np.exp(s)[self.f]
        return qa, ~(qa < self.fmin)                       # pressure.cu: at_floor = q A < F_min

    def jv(self, st, k):
        qa, on = self._on(st.s)
        JNG, JRG = np.empty((k, self.n)), np.empty((k, self.n))
        for j in range(k):
            d = np.where(on, qa * self._facets(self.v[j])[self.f], 0.0)
            JNG[j] = self.cN * np.bincount(self.p, self.aN * d, minlength=self.n)
            JRG[j] = self.cR * np.bincount(self.p, self.aR * d, minlength=self.n)
        return JNG, JRG

    def jtu(self, st, uN, uR):
        qa, on = self._on(st.s)
        w = np.where(on, (uN[self.p] * self.cN[self.p] * self.aN + uR[self.p] * self.cR[self.p] * self.aR) * qa, 0.0)
        return np.bincount(self.f, w, minlength=self.nf)[self.idx]

    def seen(self, mask, structural=False):
        hit = ((self.aN != 0) | (self.aR != 0)) if structural else ((self.aN > 0) | (self.aR > 0))
        m = np.asarray(mask, bool)[self.p] & hit
        out = np.zeros(self.nf, bool)
        out[self.f[m]] = True
        return out


def synthetic(seed=3, nf=3000, lam=0.0, floor_ratio=0.02):
    """Rows (60 observations with an NG gap, late window from the first post-gap NG point, 152 grid rows), a stand-in
    engine, and the per-facet fields."""
    import grain_crossing as gc
    import positive_finetune as pf
    import positive_problem as pp
    rng = np.random.default_rng(seed)
    n_obs = 60
    et_obs = 528503000.0 + 10.0 * np.arange(n_obs)
    grid_et = 528503000.0 + 440.0 + np.arange(152.0)            # 1-s grid covering the late windows
    et = np.concatenate([et_obs, grid_et])
    n = len(et)
    s_true = rng.normal(0.0, 0.7, nf)
    A_orig = np.exp(s_true + rng.normal(0, 0.5, nf)) * (rng.uniform(size=nf) > 0.05)
    s_start = s_true + rng.normal(0, 0.3, nf)
    s_start[:40] = np.log(0.5 * floor_ratio * np.maximum(A_orig[:40], 1e-3))     # some facets start below the floor
    area = rng.uniform(0.5, 1.5, nf)
    s0 = s_start + rng.normal(0, 0.2, nf)
    engine = CpuEngine(n, nf, rng, A_orig, floor_ratio, pf.KMAX)
    _, pNG_t, pRG_t = engine.forward(s_true)
    refNG = pNG_t * np.exp(rng.normal(0, 0.1, n))
    refRG = pRG_t * np.exp(rng.normal(0, 0.1, n))
    refNG[n_obs:] = refRG[n_obs:] = 0.0
    use_ng = np.ones(n, int)
    use_rg = np.ones(n, int)
    use_ng[41:46] = 0                                            # the NG data gap (et 410 .. 450)
    use_ng[n_obs:] = use_rg[n_obs:] = 0
    late = np.flatnonzero(et_obs >= et_obs[46])
    grid = np.arange(n_obs, n)
    op = pp.WindowOperator.__new__(pp.WindowOperator)
    op.late_pos, op.grid_pos = late, grid
    op.W = pp.w10x5_operator(et[grid], et[late] + 1.0)
    op.nrows, op.calls = n, 0
    # make several crossing pairs wrong-sided in the start model (NG above RG before the gap)
    rows = pf.Rows(np.arange(n), et, refNG, refRG, use_ng, use_rg, op, 10.0, gc.GAP_SEARCH_S, (gc.A_UNIT, gc.B_UNIT))
    seen_mask = engine.seen(rows.objective_points()[0])
    return types.SimpleNamespace(rows=rows, engine=engine, s_start=s_start, s0=s0, area=area, A_orig=A_orig,
                                 seen=seen_mask, lam=lam, floor_ratio=floor_ratio, rng=rng, et=et, op=op,
                                 n_obs=n_obs)


def make_objective(sy, lam=None):
    import positive_finetune as pf
    return pf.Objective(sy.rows, sy.engine, sy.s_start, sy.s0, sy.area, sy.A_orig, sy.seen,
                        sy.lam if lam is None else lam, sy.floor_ratio)


# ------------------------------------------------------------------------------------------------ 1 window
def test_window():
    sy = synthetic()
    rng = np.random.default_rng(11)
    worst = 0.0
    for _ in range(5):
        x, u = rng.standard_normal(sy.rows.n), rng.standard_normal(sy.rows.n)
        a, b = float(sy.rows.window(x) @ u), float(x @ sy.rows.window_T(u))
        worst = max(worst, abs(a - b) / max(abs(a), abs(b)))
    # the window replaces exactly the late rows
    x = rng.standard_normal(sy.rows.n)
    y = sy.rows.window(x)
    keep = np.setdiff1d(np.arange(sy.rows.n), sy.op.late_pos)
    report("1 window_T adjoint of WindowOperator.apply", worst < 1e-13 and np.array_equal(x[keep], y[keep]),
           "rel %.1e" % worst)


# ------------------------------------------------------------------------------------------------ 2/3 grain_crossing
def install_grain_crossing(P, weight=10.0):
    import grain_crossing as gc
    captured = {}

    class Base:
        pass

    def stub(fun, x0, jac=None, **kw):
        captured["fun"], captured["jac"] = fun, jac
        r = fun(x0)
        return types.SimpleNamespace(fun=r, jac=None, cost=0.5 * float(r @ r))

    fake = types.SimpleNamespace(Problem=Base, write_fit=lambda *a, **k: None, least_squares=stub, FLOOR=1e-300, P=P)
    record = {}
    gc.install(fake, weight, record)
    return fake, captured


def test_grain_and_crossing():
    sy = synthetic()
    rows = sy.rows
    P = types.SimpleNamespace(M=types.SimpleNamespace(rows=np.arange(rows.n)), useNG=rows.useNG, useRG=rows.useRG,
                              NNG=rows.NNG, NRG=rows.NRG, wNG=rows.wNG, wRG=rows.wRG, et_fit=rows.et,
                              et_end=rows.et_end, et_before=rows.et_before, lnNG=rows.lnNG, lnRG=rows.lnRG)
    fake, cap = install_grain_crossing(P)
    gp = object.__new__(fake.Problem)
    gp.et_end = rows.et_end
    worst = 0.0
    for x in ([3.0, 0.0], [3.0, 2.5], [0.0, 1.0], [5.0, 40.0]):
        g1 = gp.grain(np.array([9.9] + x), rows.et)
        g2 = rows.grain(x[0], x[1])
        worst = max(worst, max(float(np.max(np.abs(a - b))) for a, b in zip(g1, g2)))
    report("2 Rows.grain == GrainProblem.grain", worst == 0.0, "max abs diff %.1e" % worst)
    # crossing: a linear data model r(x) = r0 + Jd x driven through grain_crossing's least_squares wrapper
    rng = np.random.default_rng(5)
    obj = make_objective(sy)
    st = obj.evaluate(obj.start_vector(3.0, 0.0))
    r0 = st.rd
    Jd = rng.standard_normal((rows.nd, 4)) * 0.05

    def fun(x):
        return r0 + Jd @ x

    def jac(x):
        return Jd

    fake.least_squares(fun, np.zeros(4), jac=jac)
    worst_q, worst_j, active = 0.0, 0.0, 0
    for k in range(4):
        x = rng.standard_normal(4) * 0.3
        fa = cap["fun"](x)
        ja = cap["jac"](x)
        r = fun(x)
        q, d = rows.penalty(r)
        coef = rows.c_w * rows.c_sign * (d > 0.0)
        mine = coef[:, None] * (Jd[rows.c_iNG] / rows.wNG - Jd[rows.c_iRG] / rows.wRG)
        worst_q = max(worst_q, float(np.max(np.abs(fa[rows.nd:] - q))))
        worst_j = max(worst_j, float(np.max(np.abs(ja[rows.nd:] - mine))))
        active += int((d > 0).sum())
    report("3 crossing rows/Jacobian == grain_crossing's", worst_q == 0.0 and worst_j == 0.0 and active > 0,
           "pairs %d, active %d over 4 states, max diff %.1e / %.1e" % (rows.npen, active, worst_q, worst_j))


# ------------------------------------------------------------------------------------------------ 4 bounds
def test_bounds():
    import positive_finetune as pf
    s = np.array([0.0, 1.0, -10.0, 2.0])
    A = np.array([1.0, 0.0, 1.0, 100.0])
    lb, below = pf.lower_bounds(s, A, 0.002)
    exp = np.array([math.log(0.002), -np.inf, 0.0, min(0.0, math.log(0.2) - 2.0)])
    ok = np.allclose(lb[np.isfinite(exp)], exp[np.isfinite(exp)], rtol=0, atol=1e-15) and np.isneginf(lb[1]) \
        and list(below) == [False, False, True, False]
    report("4 lower_bounds rule", ok, "lb %s below %s" % (np.round(lb, 4).tolist(), below.tolist()))


# ------------------------------------------------------------------------------------------------ 5 objective
def test_objective():
    for lam in (0.0, 0.3):
        sy = synthetic(lam=lam)
        obj = make_objective(sy)
        rows = obj.rows
        rng = np.random.default_rng(7)
        z = obj.start_vector(3.0, 0.5)
        st = obj.evaluate(z)
        J = obj.jac(z)
        # cost decomposition
        D2 = float(np.sum(sy.area * (st.s - sy.s0) ** 2) / np.sum(sy.area))
        okc = (abs(st.costs["regularization"] - 0.5 * lam * D2) <= 1e-12 * max(1.0, D2) and
               abs(st.cost - st.costs["total"]) <= 1e-14 * max(1.0, st.cost))
        # adjoint
        worst = 0.0
        for _ in range(4):
            v, u = rng.standard_normal(obj.n), rng.standard_normal(obj.m)
            a, b = float(J.matvec(v) @ u), float(v @ J.rmatvec(u))
            worst = max(worst, abs(a - b) / max(abs(a), abs(b)))
        # J v vs FD
        fd_worst = 0.0
        for name in ("facets", "grain_A", "grain_B", "all"):
            v = np.zeros(obj.n)
            if name in ("facets", "all"):
                v[:obj.npar] = rng.standard_normal(obj.npar)
            if name in ("grain_A", "all"):
                v[-2] = 1.0
            if name in ("grain_B", "all"):
                v[-1] = 1.0
            h = 1e-6
            fd = (obj.evaluate(z + h * v).r - obj.evaluate(z - h * v).r) / (2 * h)
            Jv = J.matvec(v)
            fd_worst = max(fd_worst, float(np.max(np.abs(fd - Jv)) / np.max(np.abs(fd))))
        # single facets: (J^T r)_k vs d cost / d delta_k
        g = J.rmatvec(st.r)
        top = np.argsort(-np.abs(g[:obj.npar]))[:4]
        f_worst = 0.0
        for k in list(top) + [obj.n - 2, obj.n - 1]:
            e = np.zeros(obj.n)
            e[k] = 1.0
            h = 1e-6
            dc = (obj.evaluate(z + h * e).cost - obj.evaluate(z - h * e).cost) / (2 * h)
            f_worst = max(f_worst, abs(dc - g[k]) / max(abs(dc), abs(g[k])))
        # matmat == column-wise matvec
        V = rng.standard_normal((obj.n, 3))
        mm = float(np.max(np.abs(J.matmat(V) - np.column_stack([J.matvec(V[:, j]) for j in range(3)]))))
        ok = okc and worst < 1e-12 and fd_worst < 1e-6 and f_worst < 1e-5 and mm == 0.0 and obj.m == (
            rows.nd + (obj.nf if lam > 0 else 0) + rows.npen)
        report("5 objective lambda %g: costs, adjoint, J v / J^T u vs FD" % lam, ok,
               "adjoint %.1e, Jv-FD %.1e, facet grad-FD %.1e, matmat %.0e, pars %d (of %d facets), rows %d, "
               "active pairs %d" % (worst, fd_worst, f_worst, mm, obj.npar, obj.nf, obj.m, int((st.dpen > 0).sum())))


# ------------------------------------------------------------------------------------------------ 6 optimizer
def test_optimize():
    import positive_finetune as pf
    args = types.SimpleNamespace(x_scale_grain=1.0, lsmr_maxiter=20, max_nfev=25, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    for lam in (0.0, 0.3):
        sy = synthetic(lam=lam)
        obj = make_objective(sy)
        z0 = obj.start_vector(3.0, 0.5)
        c0 = obj.evaluate(z0).cost
        import io
        import contextlib
        with contextlib.redirect_stdout(io.StringIO()):
            st, rec = pf.optimize(obj, z0, args)
        ok = st.cost < 0.9 * c0 and np.all(st.z >= obj.lb) and rec["nfev"] <= 25
        report("6 optimize lambda %g (TRF/LSMR, bounds)" % lam, ok,
               "cost %.4e -> %.4e (%s, nfev %d, J v %d, J^T u %d), at lb %d" % (
                   c0, st.cost, rec["termination"], rec["nfev"], rec["counters"]["jv_columns"],
                   rec["counters"]["jtu"], obj.summary(st)["facets_at_lower_bound"]))

    # wall-time stop: the best evaluated state comes back
    class Stopper(pf.Tracker):
        def __init__(self, n):
            super().__init__()
            self.left = n

        def check(self):
            self.left -= 1
            if self.left < 0:
                raise pf.StopFit("test stop")

    sy = synthetic()
    obj = pf.Objective(sy.rows, sy.engine, sy.s_start, sy.s0, sy.area, sy.A_orig, sy.seen, 0.0, sy.floor_ratio,
                       Stopper(60))
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        st, rec = pf.optimize(obj, obj.start_vector(3.0, 0.5), args)
    best = min(h["cost"] for h in obj.tracker.history)
    report("6 optimize wall-time stop returns the best state", rec["termination"] == "wall_time" and st.cost == best,
           "stopped after %d evaluations, cost %.4e (best evaluated %.4e)" % (rec["nfev"], st.cost, best))


# ------------------------------------------------------------------------------------------------ 7 writers
def test_writers(scratch):
    import positive_finetune as pf
    sy = synthetic()
    obj = make_objective(sy)
    rows = obj.rows
    z_fit = obj.start_vector(3.0, 0.5)
    st_fit = obj.evaluate(z_fit)
    z = z_fit.copy()
    z[:obj.npar] = np.random.default_rng(2).normal(0, 0.1, obj.npar)
    z[:obj.npar] = np.maximum(z[:obj.npar], obj.lb_delta)
    st = obj.evaluate(z)
    root = Path(tempfile.mkdtemp(prefix="ptest_", dir=scratch))
    fitdir, out = root / "fit", root / "out"
    fitdir.mkdir()
    out.mkdir()
    utc = ["2016-09-30T10:%02d:%06.3f" % divmod(t - 528502800.0, 60) for t in rows.et]
    names = "iPoint UTC et x y z altitude inside rh_AU VdotB_RG nVisible Q_H2O Q_CO2 n_H2O n_CO2 PhiRG_H2O PhiRG_CO2 " \
            "pNG_H2O pNG_CO2 pNG pRG_H2O pRG_CO2 pRG refNG refRG NGBelowDeckFraction T_RG useNG useRG"
    lines = ["# LiouvilleCOPS FIT PREDICTION (synthetic)\n", "# ProductKind=fit_prediction\n",
             "# PositiveFitGridRows (synthetic): the grid rows are model-only points\n", "#" + names + "\n"]
    for i in range(rows.n):
        lines.append("%d %s %.6f 1 2 3 4 0 3.8 0.3 100 nan nan nan nan nan nan nan nan %.6e nan nan %.6e %.6e %.6e "
                     "nan 250.0 %d %d\n" % (rows.ipoint[i], utc[i], rows.et[i], st_fit.pNG[i], st_fit.pRG[i],
                                            rows.refNG[i], rows.refRG[i], int(rows.useNG[i]), int(rows.useRG[i])))
    (fitdir / "COPS_model.dat").write_text("".join(lines))
    layout = types.SimpleNamespace(late_rows=sy.op.late_pos, grid_rows=sy.op.grid_pos,
                                   stop_et=rows.et[sy.op.late_pos] + 1.0, summary={"late_start_utc": utc[46]})
    ctx = types.SimpleNamespace(files={"COPS_model.dat": fitdir / "COPS_model.dat"}, layout=layout, fit_dir=fitdir,
                                utc=utc, held_prev=np.zeros(obj.nf), s_param=sy.s_start.copy())
    pf.write_tables(out, ctx, obj, st_fit, st, "2026-11-07T00:00:00+00:00")
    old = (fitdir / "COPS_model.dat").read_text().splitlines()
    new = (out / "COPS_model.dat").read_text().splitlines()
    hdr = [line for line in new if line.startswith("#")]
    c = {n: i for i, n in enumerate(names.split())}
    ok = len(new) == len(old) + 1 and new[0] == old[0] and hdr[1].startswith("# PositiveFineTuning")
    worst = 0.0
    k = 0
    for a, b in zip([line for line in old if not line.startswith("#")], [line for line in new if not line.startswith("#")]):
        va, vb = a.split(), b.split()
        same = [j for j in range(len(va)) if j not in (c["pNG"], c["pRG"])]
        ok &= all(va[j] == vb[j] for j in same)
        worst = max(worst, abs(float(vb[c["pNG"]]) / st.pNG[k] - 1), abs(float(vb[c["pRG"]]) / st.pRG[k] - 1))
        k += 1
    lw = [line.split() for line in (out / "late_window_model.dat").read_text().splitlines() if not line.startswith("#")]
    kinds = {v[3] for v in lw}
    ok &= len(lw) == len(sy.op.late_pos) + len(sy.op.grid_pos) and kinds == {"grid", "late_obs"}
    pc = [line for line in (out / "pressure_comparison.txt").read_text().splitlines() if not line.startswith("#")]
    ok &= len(pc) == rows.n
    # finetune_state.npz
    path, chk = pf.write_state(out, ctx, obj, st, {"fit_restart_sha256": "x"})
    with np.load(path) as zz:
        keys = set(zz.files)
        need = {"delta", "held_prev", "held_next", "s_final", "A_final", "lower_bound", "seen", "grain_A_uPa",
                "grain_B_uPa_per_100s", "meta_json"}
        ok &= need <= keys and np.array_equal(zz["held_next"], zz["held_prev"] + zz["delta"])
        ok &= bool(np.all(zz["delta"][~zz["seen"]] == 0.0)) and bool(np.all(np.isnan(zz["lower_bound"][~zz["seen"]])))
        ok &= float(np.max(np.abs(zz["s_final"] - (ctx.s_param + zz["held_next"])))) < 1e-12
    # the figure helper
    pt = subprocess.run([sys.executable, str(POSITIVE_TABLE), str(out), str(root / "obs_table.dat")],
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    ok &= pt.returncode == 0
    report("7 writers + positive_table.py on the output", ok and worst < 1e-6,
           "table rel %.1e, positive_table: %s" % (worst, pt.stdout.strip().splitlines()[-1] if pt.stdout else "-"))


# ------------------------------------------------------------------------------------------------ 8 handlers
def test_handlers():
    import positive_activity
    import positive_finetune as pf
    from cops_cuda import transport
    before = transport.WorkerState.handle
    not_installed = not getattr(before, "_positive_finetune", False) and not positive_activity._STATE["installed"]
    pf.install_worker_handlers()
    again = pf.install_worker_handlers()
    w = transport.WorkerState(0)
    r1 = w.handle(("ft_eval", 1, None))
    r2 = w.handle(("close",))
    ok = (not_installed and getattr(transport.WorkerState.handle, "_positive_finetune", False) and again is False
          and r1[0] == "error" and "ft_attach" in r1[1] and r2 == ("closed",))
    transport.WorkerState.handle = before
    report("8 worker handlers (ft_* added, others delegated)", ok, "%s / %s" % (r1[0], r2))


# ------------------------------------------------------------------------------------------------ 9/10 late_only
def test_late_only():
    import positive_finetune as pf
    sy = synthetic()
    rows, eng = sy.rows, sy.engine
    own, late_used = rows.objective_points()
    utc = ["2016-09-30T10:%02d:%06.3f" % divmod(t - 528502800.0, 60) for t in rows.et]
    ctx = types.SimpleNamespace(area=sy.area, A_orig=sy.A_orig, utc=utc)
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        late, rec, support = pf.late_only_facets(eng, rows, ctx, 5.0, own, late_used)
    # brute force: earliest structural record time per facet over the objective points
    used = rows.useNG | rows.useRG
    t_end = rows.et[used].max()
    t0 = t_end - 300.0
    pts = (used & ~rows.is_grid) | (own & rows.is_grid)
    struct = (eng.aN != 0) | (eng.aR != 0)
    m = struct & pts[eng.p]
    inside = m & (rows.et[eng.p] >= t0) & (rows.et[eng.p] <= t_end)
    early = m & (rows.et[eng.p] < t0)
    ref = np.zeros(eng.nf, bool)
    ref[eng.f[inside]] = True
    ref[eng.f[early]] = False
    ok = np.array_equal(late, ref) and late.sum() > 0
    sy2 = synthetic()
    sy2.seen = late
    obj = make_objective(sy2)
    z = obj.start_vector(3.0, 0.5)
    a = obj.evaluate(z)
    z2 = z.copy()
    z2[:obj.npar] += np.random.default_rng(4).uniform(-1, 1, obj.npar)
    b = obj.evaluate(z2)
    before = support["before_rows"]
    same = np.array_equal(a.pNG_s[before], b.pNG_s[before]) and np.array_equal(a.pRG_s[before], b.pRG_s[before])
    changed = int(np.sum(support["inside_rows"] & (a.pNG_s != b.pNG_s)))
    report("9 FACETS late_only: set == brute force, strict support", ok and same and changed > 0,
           "%d facets (window %s .. %s), %d rows before unchanged %s, %d rows inside changed" % (
               late.sum(), rec["t_start_utc"], rec["t_end_utc"], before.sum(), same, changed))


def test_arguments():
    import positive_finetune as pf
    import contextlib
    import io
    base = ["--fit", "x", "--out", "y", "--floor-ratio", "0.002", "--devices", "0", "--max-nfev", "5"]
    res = []
    for extra in ([], ["--facets", "late_only"], ["--facets", "all_seen", "--window-min", "5"]):
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                pf.parse_args(base + extra)
            res.append("accepted")
        except SystemExit:
            res.append("refused")
    a = pf.parse_args(base + ["--facets", "late_only", "--window-min", "5"])
    ok = res == ["refused"] * 3 and a.facets == "late_only" and a.window_min == 5.0
    report("10 --facets required, late_only needs --window-min", ok, "%s; late_only 5 accepted" % res)


def test_zero_sensitivity():
    import positive_finetune as pf
    import contextlib
    import io
    args = types.SimpleNamespace(x_scale_grain=1.0, lsmr_maxiter=20, max_nfev=25, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    out = []
    for lam in (0.0, 0.3):
        sy = synthetic(lam=lam)
        pos = sy.A_orig > 0
        sy.s_start = np.where(pos, np.log(0.1 * sy.floor_ratio * np.where(pos, sy.A_orig, 1.0)), sy.s_start)
        sy.seen = pos & sy.seen
        obj = make_objective(sy)
        z0 = obj.start_vector(3.0, 0.5)
        st = obj.evaluate(z0)
        J = obj.jac(z0)
        v = np.zeros(obj.n)
        v[:obj.npar] = np.random.default_rng(1).standard_normal(obj.npar)
        jv = J.matvec(v)[:obj.rows.nd + obj.nreg + obj.rows.npen]
        data_zero = not np.any(jv[:obj.rows.nd]) and not np.any(jv[obj.rows.nd + obj.nreg:])
        zs = pf.zero_sensitivity(obj, st)
        with contextlib.redirect_stdout(io.StringIO()):
            stf, rec = pf.optimize(obj, z0, args)
        moved = float(np.max(np.abs(stf.z[:obj.npar] - z0[:obj.npar])))
        ok = (data_zero and zs["zero_sensitivity"] == obj.npar and stf.cost < st.cost and stf.z[-2] != z0[-2]
              and np.all(stf.z >= obj.lb) and ((moved < 1e-9) if lam == 0 else (moved > 1e-3)))
        out.append(ok)
        report("11 zero sensitivity lambda %g (all facets at the floor)" % lam, ok,
               "J v data block zero %s, zero-sensitivity %d of %d, cost %.4e -> %.4e (%s), max |delta change| %.1e, "
               "grain A %.3f -> %.3f" % (data_zero, zs["zero_sensitivity"], obj.npar, st.cost, stf.cost,
                                         rec["termination"], moved, z0[-2], stf.z[-2]))


# ------------------------------------------------------------------------------------------------ 12 CO2 floor
OLD_C1 = Path("/data001/highresolu_vt/cops_output/model_ext_frozen/c1_20261108T000327Z")
SBATCH = Path("/data001/highresolu_vt/cops_output/jobs/positive_finetune.sbatch")
CHECK_ENV_CHILD = r"""
import json, os, sys, types
for _p in (sys.argv[1], sys.argv[3] + "/tools", sys.argv[3]):          # as positive_finetune.py's own sys.path
    sys.path.insert(0, _p)
import positive_activity
positive_activity.install()
import positive_finetune as pf
case = json.loads(sys.argv[2])
ctx = types.SimpleNamespace(environment=case["fit_env"], environment_in_process=case["fit_env_in_process"],
                            fit_positive_activity={"config": case["fit_config"], "kernel_sha256": {}})
try:
    out = pf.check_environment(ctx, float(case["fit_env"]["COPS_POS_FLOOR_RATIO"]))
    print("RESULT accepted co2_floor=%s compat=%s" % (out["co2_floor"], out["pressure_sha256_compatibility"][:20]))
except SystemExit as exc:
    print("RESULT refused: %s" % str(exc)[:200])
"""


def _fit_record(ratio, co2):
    import positive_activity as pa
    env = {"COPS_POS_ENABLE": "1", "COPS_POS_NIGHT_RATIO": ratio, "COPS_POS_ACTIVITY": "exp",
           "COPS_POS_FLOOR_RATIO": ratio}
    if co2:
        env["COPS_POS_FLOOR_CO2"] = "1"
    cfg = pa._read_environment(env)
    config = {k: cfg[k] for k in ("enabled", "modified", "night_ratio", "qmin_ratio", "activity", "floor_ratio")}
    if co2:
        config["co2_floor"] = True                     # an earlier fit's record has no co2_floor key
    in_process = dict(env)
    if not co2:
        in_process.pop("COPS_POS_FLOOR_CO2", None)
    return env, in_process, config


def test_co2_floor_environment(scratch):
    import positive_finetune as pf
    ok_keys = "COPS_POS_FLOOR_CO2" in pf.ENV_KEYS and pf.PRESSURE_SOURCE_SHA256 in pf.PRESSURE_COMPATIBLE_SHA256
    cur = pf.sha256(HERE / "kernels" / "pressure.cu")
    ok_keys &= cur in pf.PRESSURE_COMPATIBLE_SHA256
    cases = {}
    for name, (ratio, co2_fit, co2_proc) in {"0.002 fit, unset": ("0.002", False, False),
                                             "0.002 fit, FLOOR_CO2=1": ("0.002", False, True),
                                             "1e-5 fit, FLOOR_CO2=1": ("1e-5", True, True),
                                             "1e-5 fit, unset": ("1e-5", True, False)}.items():
        env, in_process, config = _fit_record(ratio, co2_fit)
        proc = {k: v for k, v in os.environ.items() if not k.startswith("COPS_POS_")}
        proc.update({k: v for k, v in env.items() if k != "COPS_POS_FLOOR_CO2"})
        if co2_proc:
            proc["COPS_POS_FLOOR_CO2"] = "1"
        case = {"fit_env": env, "fit_env_in_process": in_process, "fit_config": config}
        r = subprocess.run([sys.executable, "-W", "ignore", "-c", CHECK_ENV_CHILD, str(HERE), json.dumps(case), str(APP)],
                           env=proc, capture_output=True, text=True, timeout=300)
        lines = [l for l in r.stdout.splitlines() if l.startswith("RESULT")]
        cases[name] = lines[0][7:] if lines else "no result (exit %d): %s" % (r.returncode, r.stderr[-400:])
    expect = {"0.002 fit, unset": "accepted co2_floor=False", "0.002 fit, FLOOR_CO2=1": "refused",
              "1e-5 fit, FLOOR_CO2=1": "accepted co2_floor=True", "1e-5 fit, unset": "refused"}
    ok = ok_keys and all(cases[k].startswith(v) for k, v in expect.items())
    report("12 CO2 floor: check_environment (fit's COPS_POS_FLOOR_CO2 and co2_floor)", ok,
           "ENV_KEYS/compatible SHA %s; %s" % (ok_keys, cases))
    # the job script's environment snippet (positive_finetune.sbatch ENVLINES)
    text = SBATCH.read_text()
    i = text.index("ENVLINES=$($PY - \"$FIT\" <<'PY'\n") + len("ENVLINES=$($PY - \"$FIT\" <<'PY'\n")
    snippet = text[i:text.index("\nPY\n", i)]
    root = Path(scratch) / "ft_envlines"
    outs = {}
    for name, co2 in (("old", False), ("r1e-5", True)):
        d = root / name
        d.mkdir(parents=True, exist_ok=True)
        env, _, _ = _fit_record("0.002" if not co2 else "1e-5", co2)
        env = dict(env, ratio_decision="text")
        (d / "fit_manifest.json").write_text(json.dumps({"model_extension": {"positive": {"environment": env}}}))
        r = subprocess.run([sys.executable, "-c", snippet, str(d)], capture_output=True, text=True, timeout=60)
        outs[name] = r.stdout.split()
    ok2 = ("COPS_POS_FLOOR_CO2=1" not in outs["old"] and len(outs["old"]) == 8 and
           "COPS_POS_FLOOR_CO2=1" in outs["r1e-5"] and len(outs["r1e-5"]) == 10)
    report("12 CO2 floor: job script exports COPS_POS_FLOOR_CO2 only when the fit recorded it", ok2, str(outs))


def _block(text, start, end, after="__device__ __forceinline__ void fused_body_m("):
    """The text from ``start`` to the next ``end``, searched after ``after`` (the mixed-precision body; the FP64
    fused_body earlier in the file has a block of the same shape)."""
    i = text.index(start, text.index(after) if after in text and start not in ("__device__ __forceinline__ void "
                                                                               "static_factors_m(",) else 0)
    return text[i:text.index(end, i) + len(end)]


def _code(text):
    """C++ text without // comments and blank lines, whitespace-normalized per line."""
    out = []
    for line in text.splitlines():
        line = line.split("//")[0].rstrip()
        if line.strip():
            out.append(" ".join(line.split()))
    return "\n".join(out)


def test_co2_floor_kernel_text():
    old = (OLD_C1 / "kernels" / "pressure.cu").read_text()
    new = (HERE / "kernels" / "pressure.cu").read_text()
    ft = (HERE / "kernels" / "finetune.cu").read_text()
    sf = ("__device__ __forceinline__ void static_factors_m(", "\n}\n")
    h2o = ("        // record_factors with the recomputed xlit and lit", "        const double dlx = dl * x;\n")
    cols = ("        if (G > 0) {\n", "                acc[base + G + g] = acc[base + G + g] + rg * mb;\n")
    same = {"static_factors_m": _block(old, *sf) == _block(new, *sf),
            "h2o_floor_and_temperature_law": _block(old, *h2o) == _block(new, *h2o),
            "activity_template_columns": _block(old, *cols) == _block(new, *cols),
            "finetune_copy_of_static_factors_m": _code(_block(ft, *sf)) == _code(_block(new, *sf))}
    import positive_finetune as pf
    same["finetune_cu_unchanged"] = pf.sha256(HERE / "kernels" / "finetune.cu") == \
        "b4aa9009400d01b0bc0939c4df0dbd293e64e2c2939a04584be76c917ac0bb15"
    report("12 CO2 floor: the H2O code ft_jtu_m transposes is unchanged", all(same.values()), str(same))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=None)
    args = ap.parse_args()
    scratch = args.scratch or tempfile.gettempdir()
    t0 = time.time()
    test_window()
    test_grain_and_crossing()
    test_bounds()
    test_objective()
    test_optimize()
    test_writers(scratch)
    test_handlers()
    test_late_only()
    test_arguments()
    test_zero_sensitivity()
    test_co2_floor_environment(scratch)
    test_co2_floor_kernel_text()
    bad = [r for r in RESULTS if not r[1]]
    print("%d tests, %d failed (%.1f s)" % (len(RESULTS), len(bad), time.time() - t0), flush=True)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
