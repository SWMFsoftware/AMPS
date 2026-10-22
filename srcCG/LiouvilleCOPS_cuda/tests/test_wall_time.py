"""Wall-time stop of a fit (cops_cuda/wall_time.py, tools/gpu_fit.py; user decision 2026-10-22). CPU only.

A fake clock makes the stop deterministic:

* without a deadline nothing is installed, and the clock with no deadline leaves the optimizer bitwise
  unchanged;
* with a deadline the optimizer stops before an evaluation would end inside the reserve (time to save
  plus one residual and one Jacobian), at the best evaluated state (lowest cost), with status -2 and the
  residual evaluations used;
* a second segment from that state with the remaining evaluations continues the fit; together the two
  use at most the budget. (It is a new least_squares run, so not identical to one continuous run.)
* the first evaluation without time to save raises ``WallTimeExhausted``;
* the real ``fit_one`` (synthetic CPU problem) reports the stopped result, and gpu_fit's segment record
  carries index, status, nfev and cost;
* deadlines from HH:MM:SS, Slurm's end time (environment or squeue) and an absolute value.
"""

from pathlib import Path
import subprocess
import sys
import types

import numpy as np
import pytest
from scipy.optimize import least_squares

import test_backend_contract as contract
import gpu_fit
from cops_cuda import wall_time

cpu_model = contract.cpu_model
TOOL = Path(__file__).resolve().parents[1] / "tools" / "gpu_fit.py"


class FakeClock(object):
    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t


class Rosenbrock(object):
    """Residuals and Jacobian that advance the fake clock by fixed costs and log every residual."""

    def __init__(self, clock, fun_seconds=1.0, jac_seconds=2.0):
        self.clock, self.fun_seconds, self.jac_seconds = clock, fun_seconds, jac_seconds
        self.evaluated = []
        self.ends = []

    def fun(self, x):
        self.clock.t += self.fun_seconds
        self.ends.append(self.clock.t)
        f = np.array([10.0 * (x[1] - x[0] ** 2), 1.0 - x[0]])
        self.evaluated.append((0.5 * float(f @ f), np.array(x)))
        return f

    def jac(self, x):
        self.clock.t += self.jac_seconds
        self.ends.append(self.clock.t)
        return np.array([[-20.0 * x[0], 10.0], [-1.0, 0.0]])


X0 = np.array([-1.2, 1.0])
BOUNDS = ([-2.0, -2.0], [2.0, 2.0])


def solve(problem, clock_state, x0=X0, max_nfev=200):
    ls = clock_state.least_squares_with(least_squares) if clock_state is not None else least_squares
    return ls(problem.fun, x0, jac=problem.jac, bounds=BOUNDS, method="trf", x_scale="jac",
              max_nfev=max_nfev, ftol=1e-12, xtol=1e-12, gtol=1e-12)


def test_no_deadline_changes_nothing():
    namespace = types.SimpleNamespace(least_squares=least_squares)
    assert gpu_fit.install_wall_clock(namespace, None, 600.0) is None
    assert namespace.least_squares is least_squares
    plain = solve(Rosenbrock(FakeClock()), None)
    clock = FakeClock()
    timed = solve(Rosenbrock(clock), wall_time.WallClock(None, 600.0, clock))
    assert plain.x.tobytes() == timed.x.tobytes() and plain.nfev == timed.nfev and plain.status == timed.status


def test_stop_before_the_deadline_at_the_best_state():
    clock = FakeClock()
    problem = Rosenbrock(clock)
    deadline, save = 30.0, 3.0
    state = wall_time.WallClock(deadline, save, clock)
    res = solve(problem, state)
    assert res.status == wall_time.STOP_STATUS and not res.success
    assert "wall-time deadline" in res.message
    assert res.nfev == len(problem.evaluated) and state.stopped is not None
    reserve = save + 1.0 + 2.0
    assert max(problem.ends) <= deadline - reserve          # nothing ran into the reserve
    assert deadline - reserve - max(problem.ends) < 2.0 + 1.0     # and no evaluation was skipped needlessly
    best_cost, best_x = min(problem.evaluated, key=lambda e: e[0])
    assert res.cost == best_cost and res.x.tobytes() == best_x.tobytes()
    assert res.optimality >= 0.0 or res.optimality == -1.0
    record = state.record()
    assert record["nfev"] == res.nfev and record["best_cost"] == best_cost
    assert record["evaluation_seconds"]["fun"]["max"] == 1.0 and record["evaluation_seconds"]["jac"]["max"] == 2.0


def test_resume_with_the_remaining_evaluations():
    budget = 200
    clock = FakeClock()
    first = solve(Rosenbrock(clock), wall_time.WallClock(30.0, 3.0, clock), max_nfev=budget)
    assert first.status == wall_time.STOP_STATUS
    clock2 = FakeClock()
    second = solve(Rosenbrock(clock2), wall_time.WallClock(None, 3.0, clock2), x0=first.x,
                   max_nfev=budget - first.nfev)
    assert second.nfev <= budget - first.nfev
    assert second.cost < first.cost and second.cost < 1e-20 and second.status > 0
    continuous = solve(Rosenbrock(FakeClock()), None, max_nfev=budget)
    assert continuous.cost < 1e-20                     # same minimum; the paths differ (not asserted equal)


def test_no_time_for_the_first_evaluation():
    clock = FakeClock()
    with pytest.raises(wall_time.WallTimeExhausted):
        solve(Rosenbrock(clock), wall_time.WallClock(2.0, 3.0, clock))


def test_fit_one_reports_the_stopped_state_and_the_segment(cpu_model):
    fit_cops = cpu_model["fit_cops"]
    P = contract.make_problem(cpu_model, contract.CONFIGS[0])
    P.args.max_nfev = 50
    ticks = types.SimpleNamespace(t=0.0)

    def clock():                                       # every clock read advances one second
        ticks.t += 1.0
        return ticks.t
    saved = getattr(fit_cops, "P", None), fit_cops.least_squares
    fit_cops.P = P
    try:
        state = gpu_fit.install_wall_clock(fit_cops, 40.0, 5.0, clock)
        assert fit_cops.least_squares is not saved[1]
        started = 0.0
        result = fit_cops.fit_one(0.0)
    finally:
        fit_cops.P, fit_cops.least_squares = saved
    assert state.stopped is not None
    assert result["status"] == wall_time.STOP_STATUS and not result["success"]
    assert result["nfev"] == state.nfev and 1 <= result["nfev"] < 50
    assert result["cost"] == pytest.approx(state.best[0], rel=0, abs=0)
    assert np.array_equal(result["x"], state.best[1])
    record = gpu_fit.segment_record(2, result, 50, state, started)
    assert record["index"] == 2 and record["status"] == "interrupted"
    assert record["nfev"] == result["nfev"] and record["max_nfev"] == 50 and record["cost"] == result["cost"]
    done = gpu_fit.segment_record(3, dict(result, status=1, message="converged"), 30, None, started)
    assert done["status"] == "complete" and done["deadline_utc"] is None


def test_deadlines():
    assert wall_time.parse_duration("01:30:00") == 5400.0
    assert wall_time.parse_duration("2-00:00:01") == 2 * 86400.0 + 1.0
    assert wall_time.parse_duration("45:00") == 2700.0 and wall_time.parse_duration("90") == 90.0
    with pytest.raises(ValueError):
        wall_time.parse_duration("1:2:3:4")
    assert wall_time.resolve_deadline(None) is None and wall_time.resolve_deadline("none") is None
    assert wall_time.resolve_deadline("00:10:00", margin=60.0, start=1000.0) == 1000.0 + 600.0 - 60.0
    assert wall_time.resolve_deadline(deadline=5000.0, margin=300.0) == 4700.0
    assert wall_time.resolve_deadline("slurm", margin=0.0, environ={"SLURM_JOB_END_TIME": "123456"}) == 123456.0

    def squeue(cmd, **kw):
        assert cmd[:4] == ["squeue", "-h", "-j", "77"]
        return types.SimpleNamespace(stdout="2026-10-22T14:00:00\n")
    end = wall_time.slurm_end_time({"SLURM_JOB_ID": "77"}, run=squeue)
    import time
    assert end == time.mktime((2026, 10, 22, 14, 0, 0, 0, 0, -1))
    with pytest.raises(ValueError):
        wall_time.slurm_end_time({})


def test_gpu_fit_refuses_two_deadlines_and_bad_wall_times():
    def run(*args):
        return subprocess.run([sys.executable, str(TOOL)] + list(args), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True)
    r = run("--wall-time", "01:00:00", "--deadline", "5", "--", "run")
    assert r.returncode != 0 and "not both" in r.stderr
    r = run("--wall-time", "soon", "--", "run")
    assert r.returncode != 0 and "wall time must be" in r.stderr
