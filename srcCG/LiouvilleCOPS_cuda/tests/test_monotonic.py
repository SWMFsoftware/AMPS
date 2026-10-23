"""Soft monotonicity penalty around least_squares (cops_cuda/monotonic.py; user decision 2026-10-22)."""

from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares

from cops_cuda import monotonic


def toy_problem(seed=1):
    """Two gauges with the fitter's residual layout: r = w_g (ln p - ln p_meas), NG block then RG block.

    ln p(t) = x0 + x1 t + x2 sin(6 t): the sine lets the model decrease along the descent; the
    measurements carry a dip, so the unpenalized fit follows it."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 1.0, 60)
    rows = np.arange(len(t))
    useNG = rng.random(len(t)) > 0.2
    useRG = np.ones(len(t), bool)
    truth = 0.5 + 2.0 * t + 0.8 * np.sin(6 * t)
    lnNG, lnRG = truth[useNG], truth[useRG] - 0.3
    P = SimpleNamespace(M=SimpleNamespace(rows=rows), useNG=useNG, useRG=useRG,
                        wNG=1 / np.sqrt(useNG.sum()), wRG=1 / np.sqrt(useRG.sum()), lnNG=lnNG, lnRG=lnRG)

    def lnp(x):
        return x[0] + x[1] * t + x[2] * np.sin(6 * t)

    def fun(x):
        m = lnp(x)
        return np.concatenate([P.wNG * (m[useNG] - lnNG), P.wRG * ((m[useRG] - 0.3) - lnRG)])   # RG: model - 0.3

    def jac(x):
        cols = np.column_stack([np.ones_like(t), t, np.sin(6 * t)])
        return np.vstack([P.wNG * cols[useNG], P.wRG * cols[useRG]])

    return P, fun, jac, lnp


def test_penalty_jacobian_matches_finite_differences():
    P, fun, jac, _ = toy_problem()
    mono = monotonic.Monotonic(5.0, P)
    seen = {}

    def capture(f, x0, jac=None, **kw):
        seen["f"], seen["j"] = f, jac
        return least_squares(f, x0, jac=jac, max_nfev=1)
    mono.least_squares_with(capture)(fun, np.array([0.1, 1.0, 0.8]), jac=jac)
    x = np.array([0.2, 0.7, 0.9])                  # a state with decreasing stretches
    f0, J = seen["f"](x), seen["j"](x)
    assert len(f0) > len(fun(x)) and np.any(f0[len(fun(x)):] > 0)
    for k in range(3):
        h = np.zeros(3)
        h[k] = 1e-6
        fd = (seen["f"](x + h) - seen["f"](x - h)) / 2e-6
        assert np.allclose(fd, J[:, k], atol=1e-6)


def test_reported_cost_is_the_fitters_and_the_fit_becomes_monotonic():
    P, fun, jac, lnp = toy_problem()
    plain = least_squares(fun, np.zeros(3), jac=jac)
    assert np.any(np.diff(lnp(plain.x)) < -1e-3)                   # follows the dip
    mono = monotonic.Monotonic(50.0, P)
    res = mono.least_squares_with(least_squares)(fun, np.zeros(3), jac=jac)
    assert len(res.fun) == len(fun(res.x)) and res.jac.shape[0] == len(res.fun)
    assert np.isclose(res.cost, 0.5 * np.sum(fun(res.x) ** 2))
    assert res.cost >= plain.cost - 1e-12                           # the data fit pays for monotonicity
    total = lambda x: np.sum(np.maximum(-np.diff(lnp(x)), 0.0))     # summed log decrease along the descent
    assert total(res.x) < 0.2 * total(plain.x)
    rec = mono.record()["final"]
    assert set(rec["gauges"]) == {"NG", "RG"} and rec["penalty_cost"] >= 0.0
    assert rec["gauges"]["NG"]["pairs"] == int(P.useNG.sum()) - 1


def test_pairs_follow_time_order_not_residual_order():
    P, fun, jac, _ = toy_problem()
    P.M.rows = P.M.rows[::-1].copy()                               # rows listed backwards in time
    mono = monotonic.Monotonic(1.0, P)
    pairs = mono._pairs()
    ng = pairs[0]
    assert np.all(ng[4] > ng[5])                                   # earlier point has the larger index
