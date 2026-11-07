# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/late_weights.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/late_weights.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: d52cc3d5fdd075e59a06c27c1cde24b983a4b8e51d2edd0f2f050b073de1bd08
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Multiplies the cost terms of the last N used points of the selected gauges by W (residual and
#   Jacobian rows by sqrt(W)); wraps fit_cops.least_squares inside the crossing wrapper; returns unweighted results.
# USED BY: grain_fit_weighted.py only. Not run directly.
# CALL: late_weights.install(fit_cops, {"last": 6, "factor": 20, "profile": "linear", "start": 1, "gauges": ["NG","RG"]},
#   rec) BEFORE grain_crossing.install; optional "values": [W1, ..., WN] (time order) overrides factor/profile/start.
# HISTORY: "start" and "values" added 2026-11-06 22:13 UTC; with start = 1 and no values the weights are those of the
#   first version (job 122236), checked.
# ==== end of header ====
"""Larger weights for the last N data points of the fit (user request 2026-11-06: "the weight for the last 5 points
will be increased from 1 to 10, where 1 corresponds to the current weight assigned to each point ... make the
selection of the number of points from the end of the dataset and how their importance is increased parameters
defined with CLI"; answers: both gauges, the weight multiplies the point's cost term, constant profile).

Definition. fit_cops' data residual of point i of gauge g is r_i = w_g (ln p_model - ln p_obs), w_g = 1/sqrt(N_g),
and its cost term is r_i^2 / 2 (weight 1). With weight W_i the cost term becomes W_i r_i^2 / 2: the residual and its
Jacobian row are multiplied by sqrt(W_i). W_i = 1 except for the last N used points (time order) of each selected
gauge:
    profile constant: W = F for all N;
    profile linear:   W = S + (F - S) k/(N - 1), k = 0..N-1: S at the N-th point from the end, F at the last
                      (S = --weight-start, default 1; added 2026-11-06, user: "starting from 1 for point 6 from the
                      end to 20 for point 1 from the end ... linearly"; with S = 1 the profile is the earlier one);
    explicit values:  W_1..W_N given in time order (the last value for the last point), any shape.
The NG = RG crossing penalty (grain_crossing.py) is computed from the unweighted residuals, unchanged.

Implementation: a wrapper of fit_cops.least_squares installed BEFORE grain_crossing.install (so the crossing
wrapper sits outside it and sees unweighted residuals). The optimizer minimizes the weighted objective; on return
result.fun/jac are unweighted again and result.cost is the unweighted cost, so every reported cost, rms, sigma and
the restart file are on the usual (unweighted) scale and comparable with the other fits. The weighted objective is
recorded in fit_manifest.json["model_extension"]["late_weights"]["final"].
"""

import numpy as np


def weights(n, factor, profile, start=1.0, values=None):
    """The N weights in time order (the last one for the last point)."""
    if values is not None:
        w = np.asarray(values, float)
        if len(w) != n or np.any(w <= 0):
            raise SystemExit("late_weights: --weight-values needs N = %d positive values, got %s" % (n, list(w)))
        return w
    if n <= 0:
        return np.zeros(0)
    if profile == "constant" or n == 1:
        return np.full(n, float(factor))
    if profile == "linear":
        return float(start) + (float(factor) - float(start)) * np.arange(n) / (n - 1.0)
    raise ValueError("late_weights: profile must be constant or linear")


def install(fit_cops, cfg, record):
    """cfg: {"last": N, "factor": F, "profile": "constant"|"linear", "gauges": ["NG", "RG"],
    optional "start": S (linear profile, default 1), "values": [W_1, ..., W_N] (overrides factor/profile/start)}."""
    n, factor, profile, gauges = int(cfg["last"]), float(cfg["factor"]), cfg["profile"], list(cfg["gauges"])
    start, values = float(cfg.get("start", 1.0)), cfg.get("values")
    if n < 1 or factor <= 0 or not set(gauges) <= {"NG", "RG"} or not gauges:
        raise SystemExit("late_weights: need --weight-last >= 1, --weight-factor > 0, --weight-gauges from ng,rg")
    record.update(last=n, factor=factor, profile=profile, start=start, values=values, gauges=gauges,
                  definition="cost term of each of the last N used points of each gauge multiplied by W "
                             "(residual and Jacobian row by sqrt(W)); W = %s" % [float(v) for v in weights(n, factor, profile, start, values)])
    original = fit_cops.least_squares

    def scale_vector(P, nrows):
        nd = int(P.NNG + P.NRG)
        s = np.ones(nrows)
        rows = {}
        for g, use, off, N in (("NG", P.useNG, 0, int(P.NNG)), ("RG", P.useRG, int(P.NNG), int(P.NRG))):
            if g not in gauges:
                continue
            et = P.et_fit[use] if hasattr(P, "et_fit") else None
            if et is not None and not np.all(np.diff(et) > 0):
                raise SystemExit("late_weights: the %s residual rows are not in time order" % g)
            k = min(n, N)
            w = weights(n, factor, profile, start, values)[n - k:]
            s[off + N - k:off + N] = np.sqrt(w)
            idx = np.nonzero(use)[0][N - k:]
            rows[g] = {"fit_rows": [int(i) for i in idx], "W": [float(v) for v in w],
                       "point_ids": [int(np.asarray(P.M.rows)[i]) for i in idx]}
        if nrows < nd:
            raise SystemExit("late_weights: %d residual rows < %d data rows" % (nrows, nd))
        return s, rows

    def least_squares(fun, x0, jac="2-point", **kwargs):
        if not callable(jac):
            raise ValueError("late_weights needs the fitter's analytic Jacobian")
        P = fit_cops.P
        r0 = np.asarray(fun(x0), float)
        S, rows = scale_vector(P, len(r0))
        record["weighted_rows"] = rows

        def fun_w(x):
            return np.asarray(fun(x), float) * S

        def jac_w(x):
            return np.asarray(jac(x), float) * S[:, None]

        record["start"] = {"unweighted_cost": 0.5 * float(r0 @ r0),
                           "weighted_cost": 0.5 * float((r0 * S) @ (r0 * S))}
        result = original(fun_w, x0, jac=jac_w, **kwargs)
        weighted_cost = float(result.cost)
        result.fun = np.asarray(result.fun, float) / S
        if result.jac is not None:
            result.jac = np.asarray(result.jac, float) / S[:, None]
        result.cost = 0.5 * float(result.fun @ result.fun)
        nd = int(P.NNG + P.NRG)
        rd = result.fun[:nd]
        record["final"] = {"weighted_cost": weighted_cost, "unweighted_cost": float(result.cost),
                           "unweighted_data_cost": 0.5 * float(rd @ rd)}
        lr = {}
        for g, off, w in (("NG", 0, P.wNG), ("RG", int(P.NNG), P.wRG)):
            if g in rows:
                N = int(P.NNG) if g == "NG" else int(P.NRG)
                k = len(rows[g]["W"])
                lr[g] = [float(v) for v in result.fun[off + N - k:off + N] / w]
        record["final"]["ln_model_over_obs_weighted_points"] = lr
        print("late_weights: weighted objective %.6e (unweighted %.6e); ln(model/obs) at the weighted points %s"
              % (weighted_cost, result.cost, {g: ["%.4f" % v for v in lr[g]] for g in lr}), flush=True)
        return result

    fit_cops.least_squares = least_squares
