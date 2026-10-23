"""Soft monotonicity of the modelled gauge pressures along the descent (user decision 2026-10-22).

The fitter (``fit_cops``) stays unchanged. ``tools/gpu_fit.py --monotonic-weight MU`` replaces
``fit_cops.least_squares``, the name ``fit_one`` calls, by ``Monotonic.least_squares_with(original)``,
which appends one penalty residual per pair of consecutive fitted points of each gauge (time order;
the spacecraft altitude decreases strictly along the end-of-mission trajectory):

    q_k = MU * w_g * max(0, ln p_k - ln p_{k+1})          (g = NG or RG, w_g the fitter's gauge weight)

so a modelled pressure that decreases along the descent costs like a log misfit of MU times the
decrease. ``ln p`` comes from the fitter's own data residuals, ``r_i = w_g (ln p_i - ln p_meas_i)``,
and the Jacobian rows are exact: ``dq_k/dx = MU (J_k - J_{k+1})`` where the decrease is positive, else 0.
Before the result goes back to ``fit_one``, the penalty is removed from ``fun``/``jac`` and ``cost`` is
recomputed from the fitter's residuals alone (data plus any regularization), so every reported cost
is the fitter's objective; the penalty at the solution is recorded separately (``record()``).
"""

import numpy as np

POLICY = ("user 2026-10-22: soft monotonicity of the modelled NG and RG pressures along the descent; "
          "penalty MU * w_gauge * max(0, ln p_k - ln p_k+1) over consecutive fitted points (time order), "
          "exact Jacobian; reported costs exclude the penalty")


class Monotonic:
    def __init__(self, weight, problem):
        """``problem``: the fitter's ``Problem`` or a callable returning it at fit time (``fit_cops.P``)."""
        if not weight > 0.0:
            raise ValueError("the monotonic weight must be positive")
        self.weight = float(weight)
        self._problem = problem
        self.final = None
        self._cache = None

    # ------------------------------------------------------------------ the pairs
    def _pairs(self):
        """[(gauge, first residual index, ln-offset per residual, weight, i, j)] for consecutive points."""
        P = self._problem() if callable(self._problem) else self._problem
        rows = np.asarray(P.M.rows)
        out = []
        start = 0
        for gauge, use, w, lnm in (("NG", P.useNG, P.wNG, P.lnNG), ("RG", P.useRG, P.wRG, P.lnRG)):
            index = np.arange(int(np.sum(use))) + start
            order = np.argsort(rows[use], kind="stable")          # time order of the fitted points
            out.append((gauge, float(w), np.asarray(lnm, float), start, index[order[:-1]], index[order[1:]]))
            start += int(np.sum(use))
        return out

    def penalty(self, r):
        """Penalty residuals and their (gauge, decrease) details from the fitter's residual vector."""
        q, details = [], []
        for gauge, w, lnm, start, i, j in self.pairs:
            lnp_i = r[i] / w + lnm[i - start]
            lnp_j = r[j] / w + lnm[j - start]
            dec = np.maximum(lnp_i - lnp_j, 0.0)
            q.append(self.weight * w * dec)
            details.append((gauge, dec))
        return np.concatenate(q), details

    def penalty_jacobian(self, r, J):
        rows = []
        for gauge, w, lnm, start, i, j in self.pairs:
            active = (r[i] / w + lnm[i - start]) > (r[j] / w + lnm[j - start])
            block = self.weight * (J[i] - J[j])
            block[~active] = 0.0
            rows.append(block)
        return np.vstack(rows)

    # ------------------------------------------------------------------ the wrapper
    def least_squares_with(self, original):
        """A drop-in ``least_squares`` minimizing the fitter's objective plus the penalty."""
        def least_squares(fun, x0, jac="2-point", **kwargs):
            if not callable(jac):
                raise ValueError("the monotonic penalty needs the fitter's analytic Jacobian")
            self.pairs = self._pairs()

            def fun_aug(x):
                r = np.asarray(fun(x), dtype=float)
                self._cache = (np.array(x, dtype=float, copy=True), r)
                return np.concatenate([r, self.penalty(r)[0]])

            def jac_aug(x):
                if self._cache is not None and np.array_equal(self._cache[0], x):
                    r = self._cache[1]
                else:
                    r = np.asarray(fun(x), dtype=float)
                J = np.asarray(jac(x), dtype=float)
                return np.vstack([J, self.penalty_jacobian(r, J)])

            result = original(fun_aug, x0, jac=jac_aug, **kwargs)
            n = len(result.fun) - sum(len(p[4]) for p in self.pairs)
            r = np.asarray(result.fun[:n], dtype=float)
            q, details = self.penalty(r)
            self.final = {"penalty_cost": 0.5 * float(q @ q), "objective_cost": float(result.cost),
                          "gauges": {g: {"pairs": int(len(dec)), "decreasing_pairs": int(np.sum(dec > 0)),
                                         "largest_decrease_factor": float(np.exp(dec.max())) if len(dec) else 1.0,
                                         "total_log_decrease": float(dec.sum())}
                                     for g, dec in details}}
            result.fun = r
            if result.jac is not None:
                result.jac = result.jac[:n]
            result.cost = 0.5 * float(r @ r)
            return result
        return least_squares

    def record(self):
        return {"weight": self.weight, "policy": POLICY, "final": self.final}
