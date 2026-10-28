"""Dust-grain term in the nude gauge and an NG = RG crossing constraint (user decisions 2026-10-27).

Local, untracked extension of the unchanged fitter (``fit_cops``) and the CUDA driver (``tools/gpu_fit.py``):
nothing in the git tree is edited. ``grain_fit.py`` installs it before ``gpu_fit.main()`` runs.

Model (user 2026-10-27):
* RG = surface sources (the fitter's model, unchanged).
* NG = surface sources + grain(t). The NG jump across the data gap (no NG data ~10:24-10:33 UTC) is a dust grain
  that entered the nude gauge and sublimates:  grain(t) = max(0, A - B (t - t_end)) for t >= t_end, else 0;
  t_end = the first fitted NG point after the gap (10:33:10 UTC on 12M). A >= 0 [Pa], B >= 0 [Pa/s].
  Fitted as two extra parameters ``grain_A_uPa`` (A in 1e-6 Pa) and ``grain_B_uPa_per_100s`` (B in 1e-6 Pa per
  100 s), appended after kappa; start A = observed jump (first post-gap NG datum minus the surface model at the
  fitter's start there, clipped at 0), B = 0; lower bounds 0 (user 2026-10-27: "A = jump, B = 0"; "the first
  reading after the gap is attributed to the surface effect and the dust grain").
* Crossing (user 2026-10-27: "crossing of model vs data", weight 10): the data have NG < RG before the gap and
  NG > RG after it, so their crossing lies in the gap. At every point where both gauges are fitted, the model
  (full NG) must have the same sign: before the gap  q = MU w max(0, ln pNG - ln pRG),  from t_end on
  q = MU w max(0, ln pRG - ln pNG),  w = sqrt(w_NG w_RG) (the fitter's gauge weights), exact Jacobian.
  Like the monotonic penalty, it is removed from ``fun``/``jac``/``cost`` before the result goes back to the
  fitter, so every reported cost is the fitter's data objective; the penalty at the solution is recorded.
* The monotonic penalty (when used) sees the full NG, because it reads the fitter's residuals (user 2026-10-27).
"""

import datetime as dt
from pathlib import Path

import numpy as np

POLICY = ("user 2026-10-27: NG = surface + grain, grain(t) = max(0, A - B (t - t_end)) after the NG data gap, "
          "A, B >= 0, start A = observed jump, B = 0; RG = surface only; crossing constraint: model sign of "
          "ln(pNG/pRG) equal to the data's before the gap (NG <= RG) and after it (NG >= RG), weight MU * "
          "sqrt(w_NG w_RG), exact Jacobian, reported costs exclude it")
A_UNIT = 1e-6            # Pa per unit of grain_A_uPa
B_UNIT = 1e-6 / 100.0    # Pa/s per unit of grain_B_uPa_per_100s
NAMES = ("grain_A_uPa", "grain_B_uPa_per_100s")
GAP_SEARCH_S = 3600.0    # the gap is the largest NG data gap in the last hour before the last NG point


def et_by_point(run):
    """{iPoint: ephemeris time [s]} from the base run's COPS_model.dat."""
    names, out = None, {}
    with open(Path(run) / "COPS_model.dat", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#iPoint"):
                names = line[1:].split()
                continue
            if line.startswith("#") or not line.strip() or names is None:
                continue
            v = line.split()
            out[int(v[names.index("iPoint")])] = float(v[names.index("et")])
    return out


def utc_of_et(et, et_ref, utc_ref):
    return (utc_ref + dt.timedelta(seconds=float(et - et_ref))).isoformat(timespec="seconds")


def install(fit_cops, crossing_weight, record):
    """Replace ``fit_cops.Problem`` (grain) and wrap ``fit_cops.least_squares`` (crossing) and ``write_fit``
    (the written NG prediction includes the grain). ``record`` (a dict) receives the provenance."""
    base_problem = fit_cops.Problem
    record.update(policy=POLICY, crossing_weight=crossing_weight)

    class GrainProblem(base_problem):
        def __init__(self, M, el, coef, args):
            super().__init__(M, el, coef, args)
            ets = et_by_point(M.run)
            self.et_all = np.array([ets[int(i)] for i in M.ipoint])           # every row of the run
            et = self.et_all[M.rows]                                          # the fit rows
            et_ng = np.sort(et[self.useNG])
            idx = np.nonzero(et_ng >= et_ng[-1] - GAP_SEARCH_S)[0]       # NG points of the last hour
            k = idx[int(np.argmax(np.diff(et_ng[idx])))]                    # the largest gap among them
            self.et_before, self.et_end = float(et_ng[k]), float(et_ng[k + 1])
            self.et_fit = et
            # start: A = observed jump at t_end (first post-gap NG datum minus the surface model at the start)
            a, b, T, kappa = self.unpack(self.x0)
            pNG, _ = M.evaluate(a, T, kappa, b=b if self.nb else None,
                                enhancement_basis=self.enhancement_basis, enhancement_scale=self.enhancement_scale)
            i_end = int(np.nonzero(self.useNG & (et == self.et_end))[0][0])
            jump = max(float(M.refNG[M.rows][i_end] - pNG[i_end]), 0.0)
            self.ngrain = 2
            self.names = list(self.names) + list(NAMES)
            self.x0 = np.concatenate([self.x0, [jump / A_UNIT, 0.0]])
            self.lo = np.concatenate([self.lo, [0.0, 0.0]])
            self.hi = np.concatenate([self.hi, [np.inf, np.inf]])
            record.update(t_end_et=self.et_end, last_ng_before_gap_et=self.et_before,
                          gap_seconds=self.et_end - self.et_before,
                          start={"grain_A_uPa": jump / A_UNIT, "grain_B_uPa_per_100s": 0.0,
                                 "observed_ng_at_t_end_Pa": float(M.refNG[M.rows][i_end]),
                                 "surface_model_ng_at_t_end_Pa": float(pNG[i_end])})

        # grain(t) at given ephemeris times, and its derivatives with respect to the two parameters
        def grain(self, x, et):
            A, B = x[-2] * A_UNIT, x[-1] * B_UNIT
            s = et - self.et_end
            on = s >= 0.0
            g = np.where(on, A - B * np.maximum(s, 0.0), 0.0)
            # one-sided derivatives at the kink g = 0: raising A grows the grain (dA active for g >= 0, so the
            # grain can start from A = 0), raising B cannot lower a grain that is already 0 (dB only for g > 0)
            dA = np.where(on & (g >= 0.0), A_UNIT, 0.0)
            dB = np.where(on & (g > 0.0), -B_UNIT * np.maximum(s, 0.0), 0.0)
            g = np.maximum(g, 0.0)
            return g, dA, dB

        def data_residuals(self, x):
            a, b, T, kappa = self.unpack(x)
            pNG, pRG = self.M.evaluate(a, T, kappa, b=b if self.nb else None,
                                       enhancement_basis=self.enhancement_basis,
                                       enhancement_scale=self.enhancement_scale)
            pNG = pNG + self.grain(x, self.et_fit)[0]
            return np.concatenate([self.wNG * (np.log(np.maximum(pNG[self.useNG], fit_cops.FLOOR)) - self.lnNG),
                                   self.wRG * (np.log(np.maximum(pRG[self.useRG], fit_cops.FLOOR)) - self.lnRG)])

        def data_jacobian(self, x):
            a, b, T, kappa = self.unpack(x)
            pNG, pRG, J = self.M.evaluate(a, T, kappa, b=b if self.nb else None,
                                          enhancement_basis=self.enhancement_basis,
                                          enhancement_scale=self.enhancement_scale, jac=True,
                                          basis_index=self.fit_j, enhancement_index=list(range(self.nb)))
            g, dA, dB = self.grain(x, self.et_fit)
            pNG, pRG = np.maximum(pNG + g, fit_cops.FLOOR), np.maximum(pRG, fit_cops.FLOOR)
            uN, uR = self.useNG, self.useRG
            col = lambda dn, dp: np.concatenate([self.wNG * dn[uN] / pNG[uN], self.wRG * dp[uR] / pRG[uR]])  # noqa: E731
            cols = [col(*J[("a", j)]) * self.ascale for j in self.fit_j]
            cols += [col(*J[("b", j)]) for j in range(self.nb)]
            if not self.args.fix_temperature:
                JT = [col(*J[("T", k)]) for k in range(3)]
                if self.ordered:
                    JT = [JT[0] + JT[1] + JT[2], JT[0] + JT[1], JT[0]]
                cols += JT
            if not self.args.fix_kappa:
                cols.append(col(*J[("kappa", 0)]))
            zero = np.zeros(len(pRG))
            cols += [col(dA, zero), col(dB, zero)]
            return np.column_stack(cols)

        def activity_linear_and_derivative(self, x, parameter_index):
            if parameter_index >= len(self.x0) - self.ngrain:
                a, b, _, _ = self.unpack(x)
                total = self.M.Y @ a
                if self.nb:
                    total += self.enhancement_scale * (self.enhancement_basis @ b)
                return total, None
            return super().activity_linear_and_derivative(x, parameter_index)

    fit_cops.Problem = GrainProblem

    # ---- written prediction: NG of every row includes the grain of the result being written
    original_write = fit_cops.write_fit

    def write_fit(out, r, *rest):
        MF, P = fit_cops.MF, fit_cops.P
        original_evaluate = MF.evaluate
        g_all = P.grain(r["x"], P.et_all)[0]

        def evaluate(*a, **k):
            res = original_evaluate(*a, **k)
            return (res[0] + g_all,) + tuple(res[1:])
        MF.evaluate = evaluate
        try:
            diag = original_write(out, r, *rest)
        finally:
            del MF.evaluate
        write_grain_products(out, P, r, MF, g_all, record)
        return diag
    fit_cops.write_fit = write_fit

    # ---- crossing penalty (wraps least_squares like cops_cuda/monotonic.py)
    state = {}

    def pairs():
        P = fit_cops.P
        rows = np.asarray(P.M.rows)
        both = P.useNG & P.useRG
        iNG = np.cumsum(P.useNG) - 1                       # residual index of each fit row's NG datum
        iRG = np.cumsum(P.useRG) - 1 + P.NNG
        k = np.nonzero(both)[0]
        et = P.et_fit[k]
        sign = np.where(et >= P.et_end, -1.0, np.where(et <= P.et_before, 1.0, 0.0))   # +1: want NG <= RG
        keep = sign != 0.0
        k = k[keep]
        return {"iNG": iNG[k], "iRG": iRG[k], "sign": sign[keep], "rows": rows[k],
                "w": crossing_weight * np.sqrt(P.wNG * P.wRG),
                "lnNG": P.lnNG[iNG[k]], "lnRG": P.lnRG[iRG[k] - P.NNG]}

    def penalty(r, c):
        P = fit_cops.P
        lnNG = r[c["iNG"]] / P.wNG + c["lnNG"]
        lnRG = r[c["iRG"]] / P.wRG + c["lnRG"]
        d = c["sign"] * (lnNG - lnRG)                       # > 0: wrong side
        return c["w"] * np.maximum(d, 0.0), d

    def penalty_jacobian(r, J, c):
        P = fit_cops.P
        _, d = penalty(r, c)
        rows = c["w"] * c["sign"][:, None] * (J[c["iNG"]] / P.wNG - J[c["iRG"]] / P.wRG)
        rows[d <= 0.0] = 0.0
        return rows

    def least_squares_with(original):
        def least_squares(fun, x0, jac="2-point", **kwargs):
            if not callable(jac):
                raise ValueError("the crossing penalty needs the fitter's analytic Jacobian")
            c = pairs()
            n_pen = len(c["sign"])

            def fun_aug(x):
                r = np.asarray(fun(x), dtype=float)
                state["cache"] = (np.array(x, dtype=float, copy=True), r)
                return np.concatenate([r, penalty(r, c)[0]])

            def jac_aug(x):
                cache = state.get("cache")
                r = cache[1] if cache is not None and np.array_equal(cache[0], x) else np.asarray(fun(x), float)
                J = np.asarray(jac(x), dtype=float)
                return np.vstack([J, penalty_jacobian(r, J, c)])

            result = original(fun_aug, x0, jac=jac_aug, **kwargs)
            n = len(result.fun) - n_pen
            r = np.asarray(result.fun[:n], dtype=float)
            q, d = penalty(r, c)
            before, after = c["sign"] > 0, c["sign"] < 0
            record["final"] = {
                "penalty_cost": 0.5 * float(q @ q), "objective_cost_with_penalty": float(result.cost),
                "pairs_before_gap": int(before.sum()), "pairs_after_gap": int(after.sum()),
                "wrong_side_before_gap": int(np.sum(d[before] > 0)),
                "wrong_side_after_gap": int(np.sum(d[after] > 0)),
                "largest_wrong_side_factor": float(np.exp(d.max())) if len(d) and d.max() > 0 else 1.0}
            result.fun = r
            if result.jac is not None:
                result.jac = result.jac[:n]
            result.cost = 0.5 * float(r @ r)
            return result
        return least_squares

    fit_cops.least_squares = least_squares_with(fit_cops.least_squares)
    return GrainProblem


def write_grain_products(out, P, r, MF, g_all, record):
    """<out>/grain_term.txt (parameters, crossing summary) and <out>/grain_components.dat (per row)."""
    import json
    x = r["x"]
    k = len(P.x0) - 2
    sig = r.get("sig")
    A, B = x[k] * A_UNIT, x[k + 1] * B_UNIT
    t0 = None
    try:
        t0 = A / B if B > 0 else float("inf")
    except ZeroDivisionError:
        pass
    rec = dict(record)
    rec["result"] = {"A_Pa": A, "B_Pa_per_s": B,
                     "sigma_A_Pa": float(sig[k]) * A_UNIT if sig is not None else None,
                     "sigma_B_Pa_per_s": float(sig[k + 1]) * B_UNIT if sig is not None else None,
                     "grain_lifetime_s": t0}
    with open(Path(out) / "grain_term.json", "w") as f:
        json.dump(rec, f, indent=1, default=float)
    with open(Path(out) / "grain_term.txt", "w") as f:
        f.write("# dust grain in the nude gauge (model_ext/grain_crossing.py)\n# %s\n" % POLICY)
        f.write("t_end_et %.3f\nlast_ng_before_gap_et %.3f\n" % (P.et_end, P.et_before))
        f.write("A_Pa %.6e\nB_Pa_per_s %.6e\ngrain_reaches_zero_after_s %s\n" % (A, B, "%.1f" % t0 if t0 != float("inf") else "inf"))
        if sig is not None:
            f.write("sigma_A_Pa %.3e\nsigma_B_Pa_per_s %.3e\n" % (sig[k] * A_UNIT, sig[k + 1] * B_UNIT))
        fin = record.get("final") or {}
        for key, v in fin.items():
            f.write("crossing_%s %s\n" % (key, v))
    with open(Path(out) / "grain_components.dat", "w") as f:
        f.write("# per row of the run: ephemeris time, grain contribution to NG [Pa] (pNG in COPS_model.dat "
                "includes it; surface-only NG = pNG - grain)\n#iPoint et grain_NG_Pa\n")
        for ip, et, g in zip(MF.ipoint, P.et_all, g_all):
            f.write("%d %.3f %.6e\n" % (ip, et, g))
