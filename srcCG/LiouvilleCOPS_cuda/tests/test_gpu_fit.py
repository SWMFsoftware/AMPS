"""tools/gpu_fit.py: start override and run-setting refusals (CPU only).

* ``install_template_start`` changes only the optimizer start of the template
  amplitudes (to the user's c = 1). The prior ``P.x0`` and everything derived
  from it stay the fitter's, and the change is recorded.
* The driver refuses conflicting run settings before importing the fitter.
* ``--tbounds 30 250`` is the default unless the fitter arguments give bounds (user 2026-10-22).
"""

from pathlib import Path
import subprocess
import sys

import numpy as np

import test_backend_contract as contract
import gpu_fit

cpu_model = contract.cpu_model
TOOL = Path(__file__).resolve().parents[1] / "tools" / "gpu_fit.py"


def template_config():
    return next(c for c in contract.CONFIGS if c["templates"])


def run_fit_one(env, start):
    fit_cops = env["fit_cops"]
    P = contract.make_problem(env, template_config())
    P.args.max_nfev = 5
    x0 = P.x0.copy()
    record = {}
    saved_P, saved_fit_one = getattr(fit_cops, "P", None), fit_cops.fit_one
    fit_cops.P = P
    seen = {}
    original_ls = fit_cops.least_squares

    def spy(fun, x_start, **kw):
        seen["x_start"] = np.array(x_start, copy=True)
        return original_ls(fun, x_start, **kw)
    fit_cops.least_squares = spy
    try:
        gpu_fit.install_template_start(fit_cops, start, record)
        result = fit_cops.fit_one(0.0)
    finally:
        fit_cops.P, fit_cops.fit_one, fit_cops.least_squares = saved_P, saved_fit_one, original_ls
    return P, x0, record, seen["x_start"], result


def test_template_start_moves_only_the_optimizer_start(cpu_model):
    P, x0, record, x_start, result = run_fit_one(cpu_model, 1.0)
    first = P.na + P.nh
    t = slice(first, first + P.nt)
    assert P.nt == 2 and np.all(x_start[t] == 1.0)
    other = np.ones(len(x0), bool)
    other[t] = False
    assert np.array_equal(x_start[other], x0[other])
    assert np.array_equal(P.x0, x0) and np.all(P.x0[t] == 0.0)        # prior unchanged
    assert record["applied"] and record["template_parameters"] == P.names[t]
    assert record["x_start_scaled"] == [float(v) for v in x_start]
    # A different start gives a different (short) trajectory than the CPU default.
    _, _, record0, x_start0, result0 = run_fit_one(cpu_model, 0.0)
    assert not record0 and np.all(x_start0[t] == 0.0)
    assert not np.array_equal(result["x"], result0["x"])


def run_tool(*args):
    return subprocess.run([sys.executable, str(TOOL)] + list(args),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def test_conflicting_run_settings_are_refused():
    r = run_tool("--", "run", "--max-nfev", "10")
    assert r.returncode != 0 and "give --max-nfev to gpu_fit.py" in r.stderr
    r = run_tool("--cpu-comparison", "--max-nfev", "2000", "--", "run")
    assert r.returncode != 0 and "--cpu-comparison fixes" in r.stderr
    r = run_tool("--template-start", "-1", "--", "run")
    assert r.returncode != 0 and "template-start" in r.stderr
    r = run_tool("--", "run", "--lambda", "0", "1e-3")
    assert r.returncode != 0 and "one --lambda" in r.stderr


def test_default_temperature_bounds():
    args, added = gpu_fit.with_default_tbounds(["run", "--lambda", "0"])
    assert added and args[-3:] == ["--tbounds", "30", "250"]
    args, added = gpu_fit.with_default_tbounds(["run", "--tbounds", "30", "400"])
    assert not added and args == ["run", "--tbounds", "30", "400"]
    args, added = gpu_fit.with_default_tbounds(["run"], cpu_comparison=True)
    assert not added and args == ["run"]                   # the CPU fitter's own default
