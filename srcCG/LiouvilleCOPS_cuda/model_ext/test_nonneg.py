#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_nonneg.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_nonneg.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 860ab10ba97bdf246150e39a929fa47e89472a4a4530dff3cb952edb02bf5d14
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. Tests of the softplus / squared models (unit: softplus value and derivative, host fields, threads bitwise,
#   one-job requirements A and B; workers; e2e on the synthetic grid base: D2 -> fits with full Jacobian checks,
#   softplus s = 0 == clipped model bitwise, squared start, fine tunings with checks 1-8, cycle-2 start, main() twice in
#   one process) and the GPU helper --gpu-independent. 2026-11-11.
# USED BY: jobs/nonneg_cpu.sbatch (CPU job 124128: 20 tests, 0 failed), jobs/nonneg_test.sbatch (124130-124133).
# 2026-11-13: six multmin tests (option absent = previous snapshot bitwise; bounds [1, inf); TRF start on the bound;
#   records; marker refusals); CPU job 124658: 26 tests, 0 failed.
# ==== end of header ====
"""Tests of the non-negative activity models (nonneg_mult.py, nonneg_problem.py, nonneg_fit.py, nonneg_start.py,
nonneg_finetune.py; user decisions 2026-11-11).

CPU (a Slurm CPU job; a few minutes):
    .venv/bin/python -W ignore test_nonneg.py --scratch DIR [--only unit,workers,e2e]
 unit:    softplus value / derivative (stable form, s = 0 hard clip, large |L/s|) vs an independent expression and
          central FD; host fields of both models; chunked host steps (2, 8 threads) == 1 thread bitwise; set_L_threads
          idempotent (requirement B); the stage held term (set_stage_held) changes the fields and never the installed
          configuration (requirement A); the CPU stand-in's columns of every parameter type vs FD (incl. facets with
          s = 0).
 workers: the fine tuning's worker requests through the real MultiGpuBackend wrappers / pipe / WorkerState.handle.
 e2e:     on the synthetic grid base (with A_orig = 0 on a facet subset): D2 (clip_refit), softplus fit from D2 with the
          full Jacobian check of every parameter, the softplus -> clip check (s = 0, m = 1) vs an independent plain
          clipped evaluation (bitwise), the start cost of both models vs an independent evaluation, the squared start
          (nonneg_start.py) and fit, fine tuning with validation 1-8 for both, threads bitwise (fit and fine tuning),
          cycle-2 start cost = fine-tuning final, nonneg_fit.main twice in one process == separate processes.
GPU (inside a GPU job; jobs/nonneg_test.sbatch): --gpu-independent OUT_DIR: an UNPATCHED process evaluates the model's
          activity A (recomputed independently on the host) as the only template (basis A, scale 1, a = 0: act_lin = A)
          and compares the cost with the driver's start check.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
REF = APP.parent / "LiouvilleCOPS"
for p in (str(HERE), str(APP / "tools"), str(APP)):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("COPS_OUTPUT_ROOT", "/data001/highresolu_vt/cops_output")
OUT_ROOT = Path("/data001/highresolu_vt/cops_output")
RESULTS = []
FIXT = ["--fix", "Tmax=182.09922", "--fix", "T75=136.11149", "--fix", "Tmin=133.0"]


def report(name, ok, detail):
    RESULTS.append((name, ok, detail))
    print("%-72s %s  %s" % (name, "PASS" if ok else "FAIL", detail), flush=True)


# ------------------------------------------------------------------------------------------------ unit
def softplus_independent(L, s):
    """s logaddexp(0, L / s) (s > 0), max(0, L) (s = 0): an independent implementation."""
    L, s = np.asarray(L, float), np.asarray(s, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(s > 0, s * np.logaddexp(0.0, L / np.where(s > 0, s, 1.0)), np.maximum(L, 0.0))
    return v


def test_unit():
    import nonneg_mult as nm
    rng = np.random.default_rng(5)
    n = 20011
    s = rng.uniform(0.5, 2.0, n) * 1e12
    s[rng.random(n) < 0.05] = 0.0
    L = rng.normal(0.0, 3.0, n) * np.where(s > 0, s, 1e12)
    L[:5] = [0.0, 1e5 * 1e12, -1e5 * 1e12, 700e12, -700e12]
    g, d = nm.softplus(L, s)
    gi = softplus_independent(L, s)
    rel = np.abs(g - gi) / np.maximum(np.abs(gi), 1e-300)
    ok_v = bool(np.all(np.isfinite(g)) and np.all(g >= 0) and np.max(rel[gi > 0]) < 1e-13 and np.all(g[gi == 0] == 0))
    h = 1e-6 * np.where(s > 0, s, 1.0)
    sel = (s > 0) & (np.abs(L / np.where(s > 0, s, 1)) < 30)
    fd = (nm.softplus(L + h, s)[0] - nm.softplus(L - h, s)[0]) / (2 * h)
    ok_d = bool(np.max(np.abs(fd[sel] - d[sel])) < 1e-6 and np.all(d[s == 0] == (L[s == 0] > 0)) and
                np.all((d >= 0) & (d <= 1)))
    report("unit: softplus (stable form) vs independent logaddexp and FD; s = 0 hard clip", ok_v and ok_d,
           "max rel value %.1e over %d facets; max |FD - sigmoid| %.1e over %d; s = 0: %d facets" % (
               float(np.max(rel[gi > 0])), int(np.sum(gi > 0)), float(np.max(np.abs(fd[sel] - d[sel]))), int(sel.sum()),
               int(np.sum(s == 0))))
    # host fields, chunked == 1 thread
    k = 25
    Y = np.ascontiguousarray(rng.normal(size=(n, k)))
    Bft = np.ascontiguousarray(np.where(rng.random((n, 3)) < 0.1, rng.uniform(0, 0.3, (n, 3)), 0.0))
    a = rng.normal(size=k) * 1e17
    c = np.ones(3)
    scale = 2.2e18
    bins = rng.integers(0, 7, n).astype(np.int8)
    m = np.array([1.0, 0.3, 2.0, 0.0, 1.5, 0.7, 4.0])
    gh = np.exp(rng.normal(0, 0.3, n))
    same = []
    for model in nm.MODELS:
        nm.set_L_threads(None)
        ref = nm.host_fields(model, Y, a, Bft, c, scale, bins, m, s if model == "softplus" else None, gh, True)
        for th in (2, 8):
            nm.set_L_threads(th)
            got = nm.host_fields(model, Y, a, Bft, c, scale, bins, m, s if model == "softplus" else None, gh, True)
            same.append(all(np.array_equal(ref[x], got[x]) for x in ("A", "D", "E", "P", "u", "TS")))
        nm.set_L_threads(None)
        # the definitions
        u, TS = Y @ a, scale * (Bft @ c)
        mrow = m[bins.astype(np.intp)]
        if model == "softplus":
            want = mrow * (softplus_independent(u + TS, s) * gh)
        else:
            want = mrow * (((u * u) / scale + TS) * gh)
        same.append(bool(np.max(np.abs(ref["A"] - want) / np.maximum(np.abs(want), 1e-300)) < 1e-12))
    report("unit: host fields of both models; chunked (2, 8 threads) == 1 thread (bitwise)", all(same),
           "checks %s" % same)
    # requirement B: set_L_threads idempotent
    r1 = nm.set_L_threads(4)
    p1 = nm._STATE["pool"]
    r2 = nm.set_L_threads(4)
    p2 = nm._STATE["pool"]
    r3 = nm.set_L_threads(3)
    p3 = nm._STATE["pool"]
    ok_b = (p1 is p2 and r2["reused"] and p3 is not p1 and p1._shutdown and not r3["reused"])
    nm.set_L_threads(None)
    report("unit: set_L_threads idempotent (same pool for the same n; the old pool shut down otherwise)", ok_b,
           "%s %s %s" % (r1, r2, r3))


def test_cpu_columns(scratch):
    """The CPU stand-in's columns of every parameter type (harmonics, multipliers, templates, T, kappa) vs central FD,
    for both models, on the synthetic base with A_orig = 0 on a facet subset (softplus hard clip)."""
    import test_clip_refit as tcr
    import test_slope_mult as tsm
    import nonneg_mult as nm
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    root = Path(scratch) / "cols"
    run_g, run_n, bundle, ft, points = tcr.make_bases(root)
    zero = zero_ah2o(run_g, 0.1, 3)
    header, names, rows = fit_cops.read_model_file(str(run_g / "COPS_model.dat"))
    cols = {nm_: i for i, nm_ in enumerate(names)}
    el = fit_cops.read_elements(str(run_g / "elements.bin"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    at = reference.load("activity_templates")
    bins = nm.sm.bin_index(np.random.default_rng(31).uniform(0, 180, len(el["area"])))
    Bft = at.select_templates(at.load_bundle(str(ft), el["centroid"], el["area"]), ["ft_c1"])
    basis = np.ascontiguousarray(np.column_stack([nm.sm.indicator_basis(bins), Bft]))
    meta = {"slope_multiplier_bins": {"names": list(nm.BIN_NAMES), "edges_deg": list(nm.BIN_EDGES_DEG),
                                      "nfaces": len(bins), "counts": [int(v) for v in np.bincount(bins, minlength=7)],
                                      "bins_sha256": nm.sha256_bytes(bins),
                                      "ah2o_sha256": nm.sha256_bytes(np.ascontiguousarray(el["AH2O"])),
                                      "mesh_identity": "test"}}
    worst = {}
    for model in nm.MODELS:
        cfg = nm.process_config(model, meta, [{"name": "ft_c1", "bundle": {"sha256": "x"}}])
        nm._STATE["config"] = cfg
        Cls = nm.cpu_model_class(fit_cops, cfg)
        M = Cls(str(run_g), header, cols, rows, el, coef)
        nm.set_stage_held(None)
        ascale = float(np.max(np.abs(coef[0])))
        rng = np.random.default_rng(7)
        a = coef[0][:25] * (1.0 if model == "softplus" else 0.0) + rng.normal(0, 0.3, 25) * ascale
        if model == "squared":
            a = rng.normal(0, 0.3, 25) * 2e18 * 0.3
        b = np.concatenate([[1.0, 0.5, 2.0, 1.3, 0.8, 1.7, 0.2], [1.0]])
        T = np.array([182.09922, 136.11149, 133.0])
        scale = 2e18
        pNG, pRG, J = M.evaluate(a, T, 1.2, b=b, enhancement_basis=basis, enhancement_scale=scale, jac=True,
                                 basis_index=list(range(25)), enhancement_index=list(range(8)))
        w = 0.0
        for key in [("a", j) for j in (0, 3, 7, 12, 24)] + [("b", k) for k in range(8)] + [("T", 0), ("kappa", 0)]:
            def ev(dx):
                a2, b2, T2, kap = a.copy(), b.copy(), T.copy(), 1.2
                if key[0] == "a":
                    a2[key[1]] += dx
                elif key[0] == "b":
                    b2[key[1]] += dx
                elif key[0] == "T":
                    T2[key[1]] += dx
                else:
                    kap += dx
                return M.evaluate(a2, T2, kap, b=b2, enhancement_basis=basis, enhancement_scale=scale)
            base_v = {"a": a, "b": b, "T": T, "kappa": [1.2]}[key[0]][key[1]]
            hh = 1e-5 * max(abs(float(base_v)), 1.0 if key[0] != "a" else ascale * 1e-3)
            pp_, pm_ = ev(hh), ev(-hh)
            fd = (pp_[0] - pm_[0]) / (2 * hh)
            col = J[key][0]
            sc = max(float(np.max(np.abs(fd))), 1e-300)
            w = max(w, float(np.max(np.abs(fd - col)) / sc))
            worst["%s %s%s" % (model, key[0], key[1])] = float(np.max(np.abs(fd - col)) / sc)
        nm._STATE["config"] = None
    bad = {k: v for k, v in worst.items() if v > 1e-5}
    report("unit: CPU stand-in columns of every parameter type vs FD (both models; A_orig = 0 on %d facets)" %
           len(zero), not bad, "worst %.1e; above 1e-5: %s" % (max(worst.values()), bad))


def zero_ah2o(run, frac, seed):
    """A_orig = 0 on a random subset of the synthetic base's facets (elements.bin column 32), in place."""
    path = Path(run) / "elements.bin"
    raw = bytearray(path.read_bytes())
    import struct
    version, n, nb = struct.unpack("<iii", bytes(raw[:12]))
    d = np.frombuffer(bytes(raw[12:]), dtype="<f8").reshape(n, 34).copy()
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(rng.random(n) < frac)
    d[idx, 32] = 0.0
    path.write_bytes(bytes(raw[:12]) + d.astype("<f8").tobytes())
    return idx


def test_held_requirement_a():
    import nonneg_mult as nm
    import tempfile
    rng = np.random.default_rng(9)
    n = 5000
    d = Path(tempfile.mkdtemp(prefix="nonneg_held_"))
    h1, h2 = rng.normal(0, 0.2, n), rng.normal(0, 0.2, n)
    np.savez(str(d / "h1.npz"), held_next=h1)
    np.savez(str(d / "h2.npz"), held_next=h2)
    cfg0 = {"model": "squared", "nfaces": n}
    nm._STATE["config"] = cfg0
    text0 = nm.canonical(cfg0)
    r1 = nm.set_stage_held(str(d / "h1.npz"), nfaces=n)
    g1 = nm.held_g().copy()
    r2 = nm.set_stage_held(str(d / "h2.npz"), nfaces=n)
    g2 = nm.held_g().copy()
    r0 = nm.set_stage_held(None)
    ok = (np.array_equal(g1, np.exp(h1)) and np.array_equal(g2, np.exp(h2)) and r0 is None and nm.held_g() is None and
          nm.canonical(nm._STATE["config"]) == text0 and r1["sha256"] == nm.sha256_bytes(h1) and
          r2["sha256"] != r1["sha256"])
    try:
        nm.set_stage_held(str(d / "h1.npz"), nfaces=n, expect_sha256=r2["sha256"])
        refused = False
    except ValueError:
        refused = True
    nm._STATE["config"] = None
    nm.set_stage_held(None)
    shutil.rmtree(d)
    report("unit: requirement A: the stage held term (host) changes per call, never the installed configuration",
           ok and refused, "held sha %s -> %s -> none; digest mismatch refused %s" % (r1["sha256"][:12],
                                                                                   r2["sha256"][:12], refused))


# ------------------------------------------------------------------------------------------------ worker requests
WORKER_CHILD = r"""
import json, pickle, sys, traceback
import numpy as np
moddir, app = sys.argv[1], sys.argv[2]
for p in (app, app + "/tools", moddir):
    sys.path.insert(0, p)
import nonneg_mult as nm
cfg = {"version": 1, "model": "squared", "s_ratio": None, "test_only": False, "bin_names": list(nm.BIN_NAMES),
       "nfaces": 12, "templates": [], "n_templates": 0}
nm.install_gpu(cfg)
from cops_cuda import backend as bk, evaluator, multi_gpu, transport
from cops_cuda.errors import BackendError
B, MG, WS = evaluator.CudaBackend, multi_gpu.MultiGpuBackend, transport.WorkerState
NF, NPT = 12, 6
out = {"installed_B": {k: callable(getattr(B, m, None)) for k, m in nm.WORKER_REQUESTS.items()},
       "installed_MG": {m: callable(getattr(MG, m, None)) for m in nm.WORKER_REQUESTS.values()}, "checks": {},
       "errors": {}}
calls = []


def stub_backend(wid, positions):
    be = object.__new__(B)
    def seen_rows(row_mask, structural=False):
        row_mask = np.asarray(row_mask, bool)
        calls.append(("seen", wid, row_mask.tolist(), bool(structural)))
        m = np.zeros(NF, bool)
        m[wid::3] = bool(row_mask.any())
        return m
    def sensitivities(state, facets):
        calls.append(("sens", wid, [int(f) for f in facets], float(state.kappa)))
        ng = np.outer(np.asarray(positions, float) + 1.0, np.arange(1.0, len(facets) + 1.0)) * 10.0 ** wid
        return ng, -2.0 * ng
    def detach():
        calls.append(("detach", wid))
        return True
    be.nonneg_seen_rows, be.nonneg_sensitivities, be.nonneg_detach = seen_rows, sensitivities, detach
    return be


class Channel:
    def __init__(self, ws, label):
        self.ws, self.label, self.exitcode, self.reply = ws, label, None, None
    def send(self, message):
        self.reply = pickle.loads(pickle.dumps(WS.handle(self.ws, pickle.loads(pickle.dumps(message)))))
    def poll(self, timeout=None):
        return self.reply is not None
    def recv(self):
        r, self.reply = self.reply, None
        return r
    def alive(self):
        return True


mg = object.__new__(MG)
mg.timeout_s, mg.broken, mg.rows = 10.0, False, np.arange(NPT)
pos = [np.arange(0, 4), np.arange(4, NPT)]
mg.workers = []
for wid in (0, 1):
    ws = object.__new__(WS)
    ws.device_id, ws.activity, ws.backend = wid, None, stub_backend(wid, pos[wid])
    mg.workers.append({"device": wid, "channel": Channel(ws, "device %d" % wid), "positions": pos[wid]})


def check(name, fn):
    try:
        out["checks"][name] = bool(fn())
    except Exception:
        out["checks"][name] = False
        out["errors"][name] = traceback.format_exc()[-1500:]


row_mask = np.array([1, 0, 0, 1, 0, 1], bool)
def seen():
    calls.clear()
    got = mg.nonneg_seen_rows(row_mask, True)
    want = np.zeros(NF, bool)
    for wid in (0, 1):
        want[wid::3] |= bool(row_mask[pos[wid]].any())
    return np.array_equal(got, want) and calls == [("seen", w, row_mask[pos[w]].tolist(), True) for w in (0, 1)]
check("seen", seen)
def sens():
    calls.clear()
    state = bk.ParameterState(np.zeros(25), np.array([182.0, 136.0, 133.0]), 1.5, np.ones(7))
    ng, rg = mg.nonneg_sensitivities(state, np.array([3, 7, 9], np.int64))
    ok = ng.shape == (NPT, 3)
    for wid in (0, 1):
        w = np.outer(pos[wid] + 1.0, np.arange(1.0, 4.0)) * 10.0 ** wid
        ok = ok and np.array_equal(ng[pos[wid]], w) and np.array_equal(rg[pos[wid]], -2.0 * w)
    return ok and calls == [("sens", w, [3, 7, 9], 1.5) for w in (0, 1)]
check("sensitivities", sens)
def detach():
    calls.clear()
    return mg.nonneg_detach() == [True, True] and calls == [("detach", 0), ("detach", 1)]
check("detach", detach)
def error_path():
    def boom(row_mask, structural=False):
        raise ValueError("nonneg worker test failure")
    mg.workers[1]["channel"].ws.backend.nonneg_seen_rows = boom
    try:
        mg.nonneg_seen_rows(row_mask, False)
    except BackendError as exc:
        return "device 1 failed during nonneg_seen" in str(exc) and "nonneg worker test failure" in str(exc)
    finally:
        mg.broken = False
    return False
check("error_reported", error_path)
def fall_through():
    ch = mg.workers[0]["channel"]
    ch.send(("no_such_request",))
    r = ch.recv()
    return r[0] == "error" and "unknown request" in r[1]
check("other_requests_to_the_original_handler", fall_through)
out["env"] = {"frozen": bool(__import__("os").environ.get(nm.ENV_FROZEN)), "shm": bool(__import__("os").environ.get(nm.ENV_SHM))}
print("RESULT " + json.dumps(out), flush=True)
sys.exit(0 if all(out["checks"].values()) and all(out["installed_B"].values()) and all(out["installed_MG"].values())
         else 1)
"""


def test_workers(scratch):
    env = {k: v for k, v in os.environ.items() if not k.startswith("COPS_POS_") and not k.startswith("COPS_NONNEG")
           and k != "COPS_SLOPEMUL_FROZEN_CONFIG"}
    r = subprocess.run([sys.executable, "-W", "ignore", "-c", WORKER_CHILD, str(HERE), str(APP)], env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (Path(scratch) / "workers.log").write_text(r.stdout)
    res = [json.loads(l[7:]) for l in r.stdout.splitlines() if l.startswith("RESULT ")]
    res = res[0] if res else None
    report("workers: nonneg_* requests: MG wrapper -> pipe -> WorkerState.handle -> backend", r.returncode == 0 and
           res is not None, "rc %s; %s" % (r.returncode, json.dumps(res)[:900] if res else r.stdout[-1500:]))


# ------------------------------------------------------------------------------------------------ e2e harness
def cpu_drive(args, argv, fit_args_, devices, max_nfev, template_start, comm, deadline=None):
    import test_clip_refit as tcr
    import nonneg_mult as nm
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    cfg = nm.installed_config()
    if cfg is not None:
        fit_cops.Model = nm.cpu_model_class(fit_cops, cfg)
    return tcr.cpu_drive(args, argv, fit_args_, devices, max_nfev, template_start, comm, deadline)


def install_xlog(path):
    """Test instrumentation (innermost, below every hook): record the multipliers of every point at which scipy's
    least_squares evaluates fun / jac, and the bounds it was given."""
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    original = fit_cops.least_squares

    def least_squares(fun, x0, jac="2-point", **kwargs):
        P = fit_cops.P
        idx = [k for k, n in enumerate(P.names) if n.startswith("c_slope_")]
        bounds = kwargs.get("bounds")
        rec = {"names": [P.names[k] for k in idx], "x0": [float(x0[k]) for k in idx],
               "lower": [float(bounds[0][k]) for k in idx] if bounds is not None else None,
               "upper": [float(bounds[1][k]) for k in idx] if bounds is not None else None, "fun": [], "jac": []}

        def f(x, *a, **k):
            rec["fun"].append([float(x[i]) for i in idx])
            return fun(x, *a, **k)

        def j(x, *a, **k):
            rec["jac"].append([float(x[i]) for i in idx])
            return jac(x, *a, **k)
        try:
            return original(f, x0, jac=j if callable(jac) else jac, **kwargs)
        finally:
            Path(path).write_text(json.dumps(rec) + "\n")
    fit_cops.least_squares = least_squares


def child(spec):
    import gpu_fit
    gpu_fit.drive = cpu_drive
    import nonneg_mult as nm
    nm.install_gpu = lambda cfg: nm._STATE.update(installed=True, config=dict(cfg), role="cpu-test") or True
    if os.environ.get("NONNEG_TEST_XLOG"):
        install_xlog(os.environ["NONNEG_TEST_XLOG"])
    import nonneg_fit
    rcs = []
    for argv in spec["argvs"]:
        try:
            rc = nonneg_fit.main(argv)
        except SystemExit as exc:
            rc = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
            if not isinstance(exc.code, int) and exc.code is not None:
                print("CHILD_SYSTEMEXIT %s" % exc.code, flush=True)
        rcs.append(rc)
        print("CHILD_RC %s" % rc, flush=True)
    return 0 if all(r == 0 for r in rcs) else 1


def run_child(root, name, argvs, script=None, env_extra=None):
    log = Path(root) / (name + ".log")
    spec = json.dumps({"argvs": [[str(a) for a in argv] for argv in argvs]})
    env = dict(os.environ, **(env_extra or {}))
    rc = subprocess.call([sys.executable, "-W", "ignore", str(script or __file__), "--child", spec], stdout=log.open("w"),
                         stderr=subprocess.STDOUT, env=env)
    return rc, log.read_text()


V_PREVIOUS = OUT_ROOT / "model_ext_frozen/nonneg_20261111T044920Z"     # the snapshot before --slope-mult-min


def same_products(a, b):
    """Fit products of two runs: evaluations, final x, COPS_model.dat rows, nonneg_facets.npz arrays (bitwise)."""
    out = {"evaluations": evals(a) == evals(b)}
    za, zb = np.load(str(Path(a) / "fit_restart.npz")), np.load(str(Path(b) / "fit_restart.npz"))
    out["x"] = bool(np.array_equal(za["x"], zb["x"]))
    rows = lambda d: [l for l in open(Path(d) / "COPS_model.dat") if l.strip() and not l.startswith("#")]  # noqa: E731
    out["COPS_model.dat"] = rows(a) == rows(b)
    fa, fb = np.load(str(Path(a) / "nonneg_facets.npz")), np.load(str(Path(b) / "nonneg_facets.npz"))
    out["facets"] = all(np.array_equal(fa[k], fb[k]) for k in fa.files if k != "meta_json")
    tab = lambda d: [l for l in open(Path(d) / "nonneg_bins.txt") if not l.startswith("#")]  # noqa: E731
    out["nonneg_bins.txt rows"] = tab(a) == tab(b)
    return out


def test_multmin(root, run_g, points, bins_npz, D2, head_sp, SP):
    """--slope-mult-min 1 (user 2026-11-12, softplus): option absent == the previous snapshot (bitwise); option present:
    the six free multipliers have lower bound 1 in the optimizer, every optimizer point >= 1, the start check at m = 1
    equals the unbounded fit's, scipy TRF's shift of a start on the bound, the records; the restart marker refuses a
    mismatch (resume and chain, both directions)."""
    # (a) absent == previous snapshot
    OLD = root / "sp_c1_previous_snapshot"
    argv = head_sp + fit_tail(run_g, OLD, 8)
    rc, t = run_child(root, "sp_c1_previous_snapshot", [argv], script=V_PREVIOUS / "test_nonneg.py")
    NEW = root / "sp_c1_absent"
    rc2, t2 = run_child(root, "sp_c1_absent", [head_sp + fit_tail(run_g, NEW, 8)])
    sp_ = same_products(OLD, NEW) if rc == 0 and rc2 == 0 else {"exit": (rc, rc2)}
    report("multmin: option absent == the previous snapshot nonneg_20261111T044920Z (bitwise products)",
           rc == 0 and rc2 == 0 and all(sp_.values()), "%s" % sp_)
    # (b) present: bounds, every optimizer point >= 1, start, records
    SB = root / "sp_c1_multmin"
    xlog = root / "sp_c1_multmin_xlog.json"
    rc, t = run_child(root, "sp_c1_multmin", [head_sp + ["--slope-mult-min", "1"] + fit_tail(run_g, SB, 12)],
                      env_extra={"NONNEG_TEST_XLOG": str(xlog)})
    if rc:
        report("multmin: bounded fit", False, t[-3000:])
        return
    xl = json.loads(xlog.read_text())
    pts = xl["fun"] + xl["jac"]
    mins = [min(v) for v in pts]
    m, rec = manifest_record(SB)
    mU, recU = manifest_record(NEW)
    names = m["parameters"]["names"] if "parameters" in m and "names" in m["parameters"] else None
    lb = dict(zip(rec["model_check"]["layout"], (m.get("optimizer_settings") or {}).get("lower_bounds") or
                  m.get("parameters", {}).get("lower_bounds") or []))
    zB = np.load(str(SB / "fit_restart.npz"))
    zU = np.load(str(NEW / "fit_restart.npz"))
    xb = dict(zip(rec["model_check"]["layout"], zB["x"]))
    xu = dict(zip(recU["model_check"]["layout"], zU["x"]))
    below_unb = [n for n in xl["names"] if xu[n] < 1.0]
    meta = json.loads(str(zB["meta_json"]))
    ok = (xl["lower"] == [1.0] * 6 and all(np.isinf(xl["upper"])) and xl["x0"] == [1.0] * 6 and min(mins) >= 1.0 and
          all(xb[n] >= 1.0 for n in xl["names"]) and
          rec["clip_core"]["start_check"]["start"]["cost"] == recU["clip_core"]["start_check"]["start"]["cost"] and
          meta["model_extension"].get("slope_mult_min") == 1.0 and
          rec["model_check"]["multiplier_lower_bound"]["value"] == 1.0)
    report("multmin: --slope-mult-min 1: optimizer bounds [1, inf), every optimizer point >= 1, start check = unbounded",
           ok, "bounds %s..%s; %d optimizer points, min multiplier %.17g (first fun point %s); final %s; unbounded fit "
           "below 1: %s; start cost %.15e == %.15e; marker %s" % (
               xl["lower"], xl["upper"], len(pts), min(mins), xl["fun"][0], {n: float(xb[n]) for n in xl["names"]},
               below_unb, rec["clip_core"]["start_check"]["start"]["cost"],
               recU["clip_core"]["start_check"]["start"]["cost"], meta["model_extension"].get("slope_mult_min")))
    first = xl["fun"][0]
    report("multmin: scipy TRF moves the start on the bound to lb + 1e-10 before its first evaluation",
           all(v == 1.0 + 1e-10 for v in first), "first optimizer point %s (start %s)" % (first, xl["x0"]))
    txt = (SB / "NONNEG_MODEL.txt").read_text()
    bins_txt = (SB / "nonneg_bins.txt").read_text()
    rep_ok = ("MULTIPLIER LOWER BOUND" in txt and "multiplier lower bound 1.0" in bins_txt and
              "slope_mult_min" in json.dumps(rec["config"]) and lb.get("c_slope_40_50", 1.0) == 1.0)
    report("multmin: records (manifest config, model_check, marker, NONNEG_MODEL.txt, nonneg_bins.txt)", rep_ok,
           "bins note: %s" % [l for l in bins_txt.splitlines() if "lower bound" in l][:1])
    # (c) marker refuses a mismatch
    ftgrain = None
    h_res = ["--model", "softplus"] + d2_head(points, bins_npz, D2) + ["--start", "d2", "--restart-from",
                                                                       SB / "fit_restart.npz"]
    R1 = root / "sp_resume_without"
    rc1, t1 = run_child(root, "sp_resume_without", [h_res + fit_tail(run_g, R1, 3)])
    h_res2 = ["--model", "softplus"] + d2_head(points, bins_npz, D2) + ["--start", "d2", "--slope-mult-min", "1",
                                                                        "--restart-from", NEW / "fit_restart.npz"]
    R2 = root / "sp_resume_with_from_unbounded"
    rc2, t2 = run_child(root, "sp_resume_with_from_unbounded", [h_res2 + fit_tail(run_g, R2, 3)])
    R3 = root / "sp_resume_with"
    rc3, t3 = run_child(root, "sp_resume_with", [["--model", "softplus"] + d2_head(points, bins_npz, D2) + [
        "--start", "d2", "--slope-mult-min", "1", "--restart-from", SB / "fit_restart.npz"] + fit_tail(run_g, R3, 3)])
    ok = (rc1 != 0 and "multiplier lower bound" in t1 and rc2 != 0 and "multiplier lower bound" in t2 and rc3 == 0)
    report("multmin: the restart marker refuses a resume with another --slope-mult-min (both directions); same value ok",
           ok, "without on a bounded restart: exit %s; with on an unbounded restart: exit %s; same value: exit %s" % (
               rc1, rc2, rc3))
    return SB


def run_child_unused():
    return None


def manifest_record(out):
    m = json.loads((Path(out) / "fit_manifest.json").read_text())
    return m, (m.get("model_extension") or {}).get("nonneg") or {}


def evals(seg):
    return [(l.split()[1], l.split()[3]) for l in open(Path(seg) / "evaluations.txt") if not l.startswith("#")]


def independent_cost(run, point_table, npz_path, x_by_name, t_end_et, A):
    """The UNPATCHED reference Model with the activity A as its only template (basis A, scale 1, a = 0: act_lin = A
    exactly), W10x5 (positive_problem.install_window), grain and data as grain_crossing / fit_cops."""
    import positive_problem as pp
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    header, names, rows = fit_cops.read_model_file(str(Path(run) / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    el = fit_cops.read_elements(str(Path(run) / "elements.bin"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    M = fit_cops.Model(str(run), header, cols, rows, el, coef)
    with np.load(str(npz_path), allow_pickle=False) as z:
        d = {k: np.array(z[k]) for k in z.files}
    table = pp.read_point_table(point_table)
    layout = pp.LateLayout(table, pp.read_run_points(M.run), M.ipoint)
    pp.install_window(M, layout)
    T = np.array([182.09922, 136.11149, 133.0])
    pNG_s, pRG_s = M.evaluate(np.zeros(25), T, x_by_name["kappa"], b=np.array([1.0]),
                              enhancement_basis=np.ascontiguousarray(np.asarray(A, float)[:, None]),
                              enhancement_scale=1.0)
    Ag, Bg = x_by_name["grain_A_uPa"] * 1e-6, x_by_name["grain_B_uPa_per_100s"] * (1e-6 / 100.0)
    s = d["et"] - t_end_et
    on = s >= 0.0
    g = np.maximum(np.where(on, Ag - Bg * np.maximum(s, 0.0), 0.0), 0.0)
    useNG, useRG = d["useNG"].astype(bool), d["useRG"].astype(bool)
    sel = np.asarray(M.rows)
    lnNG, lnRG = np.log(M.refNG[sel][useNG]), np.log(M.refRG[sel][useRG])
    wNG, wRG = 1.0 / np.sqrt(useNG.sum()), 1.0 / np.sqrt(useRG.sum())
    r = np.concatenate([wNG * (np.log(np.maximum((pNG_s + g)[useNG], fit_cops.FLOOR)) - lnNG),
                        wRG * (np.log(np.maximum(pRG_s[useRG], fit_cops.FLOOR)) - lnRG)])
    info = {"surface_bitwise": bool(np.array_equal(pNG_s, d["surface_pNG"]) and np.array_equal(pRG_s, d["surface_pRG"])),
            "surface_max_rel": float(max(np.max(np.abs(pNG_s / d["surface_pNG"] - 1)),
                                         np.max(np.abs(pRG_s / d["surface_pRG"] - 1))))}
    return 0.5 * float(r @ r), info


def independent_A(model, run, x_by_name, ascale, A_s, ft_bundles, s_ratio=1e-5, held=None, m=None):
    """The model's activity from its definition, independently of nonneg_mult (unchunked numpy; softplus via
    logaddexp)."""
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    at = reference.load("activity_templates")
    el = fit_cops.read_elements(str(Path(run) / "elements.bin"))
    a = np.array([x_by_name["a%02d" % j] for j in range(25)]) * ascale
    T = np.column_stack([at.select_templates(at.load_bundle(str(p), el["centroid"], el["area"]), [nm_])
                         for p, nm_ in ft_bundles])
    u = el["Y"] @ a
    TS = A_s * (T @ np.ones(T.shape[1]))
    if model == "softplus":
        q = softplus_independent(u + TS, s_ratio * np.asarray(el["AH2O"], float))
    else:
        q = (u * u) / A_s + TS
    if held is not None:
        q = q * np.exp(held)
    return q if m is None else m * q


def d2_head(points, bins_npz, D2):
    return ["--point-table", points, "--slope-bundle", bins_npz, "--ft-templates-from", D2 / "fit_manifest.json"]


def fit_tail(run, out, nfev):
    import test_slope_mult as tsm
    return FIXT + ["--crossing-weight", "10", "--check-jacobian", "--devices", "0", "--max-nfev", str(nfev), "--"] + \
        tsm.slope_fit_args(run, out)


class CpuEngine:
    """nonneg_finetune's engine interface on the CPU stand-in (tests)."""

    def __init__(self, ctx):
        import nonneg_mult as nm
        from cops_cuda import reference
        reference.bind()
        fit_cops = reference.load("fit_cops")
        nm._STATE.update(installed=True, config=dict(ctx.model_cfg), role="cpu-test")
        held = ctx.held_record
        nm.set_stage_held(held["path"] if held else None, (held or {}).get("key", "held_next"), ctx.nf,
                          expect_sha256=(held or {}).get("sha256"))
        nm.set_override(None)
        Cls = nm.cpu_model_class(fit_cops, ctx.model_cfg)
        self.M = Cls(str(ctx.run), ctx.header, ctx.cols, ctx.table_rows, ctx.el, ctx.coef)
        self.ctx = ctx
        self.nm = nm
        self.calls = {"forward": 0}

    def forward(self):
        c = self.ctx
        self.calls["forward"] += 1
        return self.M.evaluate(c.a, c.T, c.kappa, b=c.b, enhancement_basis=c.basis, enhancement_scale=c.escale)

    def seen(self, mask, structural=False):
        M = self.M
        sel = np.asarray(mask, bool)[M.pt]
        on = ((M.Kn != 0) | (M.KPhi != 0) | (M.KV != 0)) if structural else ((M.Kn > 0) | (M.KPhi > 0) | (M.KV > 0))
        out = np.zeros(M.nfaces, bool)
        out[M.face[sel & on]] = True
        return out

    def sensitivities(self, facets):
        return self.M.sensitivities(self.ctx.T, facets)

    def forward_override(self, facets, g):
        self.nm.set_override(facets, g)
        try:
            return self.forward()
        finally:
            self.nm.set_override(None)

    def close(self):
        self.nm.set_override(None)


# ------------------------------------------------------------------------------------------------ e2e
def test_e2e(scratch):
    import test_clip_refit as tcr
    import test_slope_mult as tsm
    import nonneg_mult as nm
    import nonneg_finetune as nft
    root = Path(scratch) / "e2e"
    run_g, run_n, bundle, ft, points = tcr.make_bases(root)
    late_f = [3, 7, 11, 13, 17, 19, 23, 29]
    t_end, changed = tsm.make_late_only(run_g, late_f, 399)
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    el = fit_cops.read_elements(str(run_g / "elements.bin"))
    bins = nm.sm.bin_index(np.random.default_rng(31).uniform(0.0, 180.0, len(el["area"])))
    bins_npz = root / "slope_bins.npz"
    tsm.write_bins_bundle(bins_npz, el, bins)
    gff_head = FIXT + ["--fix", "c_ft_c1=1", "--crossing-weight", "10", "--devices", "0", "--template-start", "0"]
    S = root / "S0"
    rcs, text = tsm.run_child(root, "S0", "gff", gff_head + ["--max-nfev", "8", "--"] + tcr.fit_args(run_n, S, bundle, ft))
    D2 = root / "D2"
    rcd, text = tsm.run_child(root, "D2", "clip", ["--variant", "D2", "--point-table", points, "--start-restart",
                                                   S / "fit_restart.npz"] + gff_head + ["--max-nfev", "6", "--"]
                              + tcr.fit_args(run_g, D2, bundle, ft))
    if rcs or rcd:
        report("e2e: sources (S0, D2)", False, "S0 %s D2 %s\n%s" % (rcs, rcd, text[-3000:]))
        return
    zero = np.flatnonzero(np.asarray(el["AH2O"]) == 0.0)      # the base's A_orig = max(0, srcCG field) (checked)
    report("e2e: synthetic sources (A_orig = 0 on %d facets: softplus hard clip there; late-only facets, D2)" %
           len(zero), True,
           "%d records moved off %d facets" % (changed, len(late_f)))
    allp = ",".join(["a%02d" % j for j in range(25)] + ["c_" + n for n in nm.BIN_NAMES[1:]] + ["kappa", "grain_A_uPa",
                                                                                              "grain_B_uPa_per_100s"])
    ftb = [(ft, "ft_c1")]
    out = {}
    # ---- softplus fit from D2: full Jacobian check of every parameter
    SP = root / "sp_c1"
    head_sp = ["--model", "softplus"] + d2_head(points, bins_npz, D2) + ["--start", "d2", "--start-restart",
                                                                         D2 / "fit_restart.npz"]
    rc, t = run_child(root, "sp_c1", [head_sp + ["--check-full-jacobian", "--check-parameters", allp] +
                                      fit_tail(run_g, SP, 8)])
    if rc:
        report("e2e: softplus fit (d2 start)", False, t[-4000:])
        return
    m, rec = manifest_record(SP)
    jc = rec["clip_core"]["jacobian_check"]["columns"]
    worst = max(v["data"] for v in jc.values() if isinstance(v, dict))
    report("e2e: softplus fit from D2: full Jacobian check (25 a, 6 m, kappa, grain), products, marker", worst < 1e-5 and
           len(jc) == 34 and (SP / "nonneg_facets.npz").is_file() and rec["start_restart"]["start_choice"] == "d2",
           "worst data %.2e over %d columns; start cost %.9e" % (worst, len(jc),
                                                                 rec["clip_core"]["start_check"]["start"]["cost"]))
    # start cost vs an independent evaluation (softplus: logaddexp; rel 1e-12)
    xb = dict(zip(rec["model_check"]["layout"], rec["start_restart"]["x_start"]))
    ascale = float(json.loads(str(np.load(str(D2 / "fit_restart.npz"))["meta_json"]))["ascale"])
    A_s = rec["model_check"]["A_scale"]
    A_ind = independent_A("softplus", run_g, xb, ascale, A_s, ftb)
    cost, info = independent_cost(run_g, points, SP / "start_check_model.npz", xb, rec["grain_window"]["t_end_et"],
                                  A_ind)
    c_drv = rec["clip_core"]["start_check"]["start"]["data_cost"]
    report("e2e: softplus start cost == independent evaluation (logaddexp form, rel <= 1e-12)",
           abs(cost / c_drv - 1) <= 1e-12, "driver %.15e independent %.15e (rel %.1e); %s" % (
               c_drv, cost, abs(cost / c_drv - 1), info))
    # softplus -> clip (s = 0) with m = 1: the plain clipped model with the ft templates, bitwise
    SC = root / "sp_clip"
    rc, t = run_child(root, "sp_clip", [head_sp + ["--s-ratio", "0", "--test-clip", "--stop-after-checks"] +
                                        fit_tail(run_g, SC, 4)])
    if rc:
        report("e2e: softplus s = 0 (clip) start", False, t[-3000:])
    else:
        m2, rec2 = manifest_record(SC)
        xb2 = dict(zip(rec2["model_check"]["layout"], rec2["start_restart"]["x_start"]))
        a = np.array([xb2["a%02d" % j] for j in range(25)]) * ascale
        at = reference.load("activity_templates")
        Tft = at.select_templates(at.load_bundle(str(ft), el["centroid"], el["area"]), ["ft_c1"])
        lin = el["Y"] @ a + A_s * (np.ascontiguousarray(Tft) @ np.ones(1))      # the reference's act_lin expression
        A_clip = np.maximum(lin, 0.0)
        cost2, info2 = independent_cost(run_g, points, SC / "start_check_model.npz", xb2,
                                        rec2["grain_window"]["t_end_et"], A_clip)
        c2 = rec2["clip_core"]["start_check"]["start"]["data_cost"]
        report("e2e: softplus -> clip (s = 0, m = 1) == plain clipped model Y a + A_s T (independent, bitwise)",
               cost2 == c2 and info2["surface_bitwise"], "driver %.15e independent %.15e; %s" % (c2, cost2, info2))
    # ---- squared: start (nonneg_start.py) and fit
    import nonneg_start
    SQS = root / "squared_start.npz"
    rc_s = nonneg_start.main(["--run", str(run_g), "--d2", str(D2 / "fit_restart.npz"), "--ft-templates-from",
                              str(D2), "--out", str(SQS), "--threads", "2", "--max-iter", "200"])
    smeta = json.loads(str(np.load(str(SQS))["meta_json"]))
    q = smeta["quality"]
    report("e2e: squared start (area-weighted LSQ to D2's activity)", rc_s == 0 and q["cost_final"] <= q["cost_start"],
           "%s" % {k: q[k] for k in ("relative_rms_residual", "total_activity_ratio", "corr_area_weighted",
                                     "iterations", "converged", "stop", "cost_start", "cost_final")})
    SQ = root / "sq_c1"
    head_sq = ["--model", "squared"] + d2_head(points, bins_npz, D2) + ["--start", "lsq", "--start-restart",
                                                                        D2 / "fit_restart.npz", "--start-file", SQS]
    rc, t = run_child(root, "sq_c1", [head_sq + ["--check-full-jacobian", "--check-parameters", allp] +
                                      fit_tail(run_g, SQ, 8)])
    if rc:
        report("e2e: squared fit (lsq start)", False, t[-4000:])
        return
    m3, rec3 = manifest_record(SQ)
    jc3 = rec3["clip_core"]["jacobian_check"]["columns"]
    worst3 = max(v["data"] for v in jc3.values() if isinstance(v, dict))
    xb3 = dict(zip(rec3["model_check"]["layout"], rec3["start_restart"]["x_start"]))
    A3 = independent_A("squared", run_g, xb3, ascale, A_s, ftb)
    cost3, info3 = independent_cost(run_g, points, SQ / "start_check_model.npz", xb3,
                                    rec3["grain_window"]["t_end_et"], A3)
    c3 = rec3["clip_core"]["start_check"]["start"]["data_cost"]
    report("e2e: squared fit from the LSQ start: full Jacobian check, start cost == independent (bitwise)",
           worst3 < 1e-5 and cost3 == c3, "worst data %.2e over %d columns; driver %.15e independent %.15e; %s" % (
               worst3, len(jc3), c3, cost3, info3))
    # ---- threads: fit with --l-threads 2 == without (bitwise evaluations)
    SP2 = root / "sp_c1_l2"
    rc, t = run_child(root, "sp_c1_l2", [head_sp + ["--l-threads", "2", "--check-full-jacobian", "--check-parameters",
                                                    allp] + fit_tail(run_g, SP2, 8)])
    report("e2e: softplus fit with --l-threads 2 == 1 thread (every evaluation bitwise)", rc == 0 and
           evals(SP) == evals(SP2), "%d evaluations" % len(evals(SP)))
    # ---- fine tuning of both fits (CPU engine), validation 1-8; threads bitwise
    for model, FIT in (("softplus", SP), ("squared", SQ)):
        FT = root / ("%s_ft" % model)
        args = nft.parse_args(["--fit", str(FIT), "--out", str(FT), "--facets", "late_only", "--window-min", "5",
                               "--max-nfev", "20", "--validate", "--l-threads", "2"])
        try:
            rc = nft.run(args, engine_factory=lambda ctx: CpuEngine(ctx), command=["test"])
        except SystemExit as exc:
            rc = str(exc)
        v = json.loads((FT / "validation.json").read_text()) if (FT / "validation.json").is_file() else {}
        st_ = {}
        if (FT / "finetune_state.npz").is_file():
            with np.load(str(FT / "finetune_state.npz"), allow_pickle=False) as z:
                st_ = {k: np.array(z[k]) for k in z.files if k != "meta_json"}
                fmeta = json.loads(str(z["meta_json"]))
        ok = (rc == 0 and v.get("passed") and np.all(np.isneginf(st_["lower_bound"][st_["seen"]])) and
              np.array_equal(st_["held_next"], st_["held_prev"] + st_["delta"]) and
              set(np.flatnonzero(st_["delta"]).tolist()) <= set(late_f) and
              fmeta["costs_final"]["cost_total"] <= fmeta["costs_initial"]["cost_total"])
        report("e2e: %s fine tuning: validation 1-8, no lower bound, held_next = held_prev + delta" % model, bool(ok),
               "exit %s; checks %s; parameters %s; cost %s -> %s" % (
                   rc, {k: (vv.get("passed") if isinstance(vv, dict) else vv) for k, vv in v.items()},
                   int(st_["seen"].sum()) if st_ else None,
                   fmeta["costs_initial"]["cost_total"] if st_ else None, fmeta["costs_final"]["cost_total"] if st_ else
                   None))
        out[model] = (FT, st_, fmeta if st_ else None)
        FT1 = root / ("%s_ft_1thread" % model)
        args1 = nft.parse_args(["--fit", str(FIT), "--out", str(FT1), "--facets", "late_only", "--window-min", "5",
                                "--max-nfev", "20", "--validate"])
        nft.run(args1, engine_factory=lambda ctx: CpuEngine(ctx), command=["test"])
        with np.load(str(FT1 / "finetune_state.npz"), allow_pickle=False) as z:
            st1 = {k: np.array(z[k]) for k in z.files if k != "meta_json"}
        report("e2e: %s fine tuning with 2 threads == 1 thread (bitwise)" % model,
               all(np.array_equal(st_[k], st1[k]) for k in ("delta", "held_next", "A_final", "grain_A_uPa",
                                                            "grain_B_uPa_per_100s")), "delta, held_next, A_final, grain")
    # ---- cycle 2 (chain, held) for both models: start cost = fine tuning final
    for model, FIT, head in (("softplus", SP, head_sp), ("squared", SQ, head_sq)):
        FT, st_, fmeta = out[model]
        C2 = root / ("%s_c2" % model)
        h2 = ["--model", model] + d2_head(points, bins_npz, D2) + [
            "--start", "chain", "--start-restart", FIT / "fit_restart.npz", "--held", FT / "finetune_state.npz",
            "--grain-from", FT / "finetune_state.npz"]
        rc, t = run_child(root, "%s_c2" % model, [h2 + fit_tail(run_g, C2, 3)])
        if rc:
            report("e2e: %s cycle-2 fit (chain, held)" % model, False, t[-3000:])
            continue
        m4, rec4 = manifest_record(C2)
        c2s = rec4["clip_core"]["start_check"]["start"]["data_cost"]
        ftf = fmeta["costs_final"]["cost_data"]
        with np.load(str(C2 / "nonneg_facets.npz"), allow_pickle=False) as z:
            has_held = "held" in z.files and np.array_equal(np.asarray(z["held"]), st_["held_next"])
        report("e2e: %s cycle-2 start cost = the fine tuning's final data cost (1e-8); held in the products" % model,
               abs(c2s / ftf - 1) <= 1e-8 and has_held, "c2 start %.15e, fine tuning final %.15e (rel %.1e); held %s"
               % (c2s, ftf, abs(c2s / ftf - 1), has_held))
        out[model + "_c2"] = (C2, h2)
    # ---- nonneg_fit.main twice in one process (c1 then c2 of softplus, then c1 again) == separate processes
    if "softplus_c2" in out:
        C2s, h2 = out["softplus_c2"]
        A1, A2, A3_ = root / "inproc_c1", root / "inproc_c2", root / "inproc_c1b"
        rc, t = run_child(root, "inproc", [head_sp + ["--check-full-jacobian", "--check-parameters", allp] +
                                           fit_tail(run_g, A1, 8), h2 + fit_tail(run_g, A2, 3),
                                           head_sp + ["--check-full-jacobian", "--check-parameters", allp] +
                                           fit_tail(run_g, A3_, 8)])
        ok = rc == 0 and evals(A1) == evals(SP) and evals(A2) == evals(C2s) and evals(A3_) == evals(SP)
        report("e2e: nonneg_fit.main 3x in one process (c1, c2 held, c1) == separate processes (bitwise evaluations)",
               ok, "exit %s; c1 %s, c2 %s, c1 again %s" % (rc, evals(A1) == evals(SP) if rc == 0 else None,
                                                         evals(A2) == evals(C2s) if rc == 0 else None,
                                                         evals(A3_) == evals(SP) if rc == 0 else None))
    test_multmin_chain(root, run_g, points, bins_npz, D2, head_sp, SP, out)


def test_multmin_chain(root, run_g, points, bins_npz, D2, head_sp, SP, out):
    SB = test_multmin(root, run_g, points, bins_npz, D2, head_sp, SP)
    if SB is None or "softplus" not in out:
        return
    FT = out["softplus"][0]
    # chain start: a bounded cycle-k fit from an unbounded previous fit is refused, and vice versa
    h2 = ["--model", "softplus"] + d2_head(points, bins_npz, D2) + [
        "--start", "chain", "--start-restart", SP / "fit_restart.npz", "--held", FT / "finetune_state.npz",
        "--grain-from", FT / "finetune_state.npz", "--slope-mult-min", "1"]
    rc, t = run_child(root, "sp_c2_bounded_from_unbounded", [h2 + fit_tail(run_g, root / "sp_c2_bfu", 3)])
    report("multmin: a chain start with --slope-mult-min from a fit without it is refused", rc != 0 and
           "multiplier lower bound" in t, "exit %s; %s" % (rc, [l for l in t.splitlines() if "lower bound" in l][-1:]))


# ------------------------------------------------------------------------------------------------ GPU independent check
def gpu_independent(out_dir, devices):
    """UNPATCHED process: the model's activity A recomputed independently on the host (independent_A) as the only
    template (basis A, scale 1, a = 0) on the 8 local GPUs (mixed precision, host activity, W10x5, grain, data as the
    fitter); compares the cost with the driver's start check."""
    import positive_problem as pp
    from cops_cuda import gpu_model, precision, reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    out_dir = Path(out_dir)
    m, rec = manifest_record(out_dir)
    cfg = rec["config"]
    run = Path(cfg["run"])
    xb = dict(zip(rec["model_check"]["layout"], rec["start_restart"]["x_start"]))
    with np.load(str(cfg["start_restart"]), allow_pickle=False) as z:
        ascale = float(json.loads(str(z["meta_json"]))["ascale"])
    A_s = rec["model_check"]["A_scale"]
    ftb = [(t["bundle"]["path"], t["name"]) for t in cfg["ft_templates"]]
    A = independent_A(cfg["model"], run, xb, ascale, A_s, ftb, cfg["s_ratio"] or 0.0)
    os.environ[precision.RECORD_PRECISION_ENVIRONMENT] = "fp32"
    gpu_model.REGISTRY.configure(devices[0], "tree", 4096, devices=devices, transport="local", activity="host")
    fit_cops.Model = gpu_model.gpu_model_class()
    cost, info = independent_cost(run, cfg["point_table"], out_dir / "start_check_model.npz", xb,
                                  rec["grain_window"]["t_end_et"], A)
    c_drv = rec["clip_core"]["start_check"]["start"]["data_cost"]
    rel = abs(cost / c_drv - 1.0)
    tol = 0.0 if cfg["model"] == "squared" else 1e-12
    report("GPU: %s start cost == UNPATCHED evaluation of the independently computed A (%s)" % (
        cfg["model"], "bitwise" if tol == 0 else "rel <= 1e-12"), rel <= tol,
        "driver %.15e independent %.15e (rel %.1e); %s" % (c_drv, cost, rel, info))
    gpu_model.REGISTRY.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scratch")
    p.add_argument("--child")
    p.add_argument("--only", default=None, help="comma list of: unit,workers,e2e")
    p.add_argument("--gpu-independent")
    p.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    args = p.parse_args()
    if args.child:
        return child(json.loads(args.child))
    if args.gpu_independent:
        gpu_independent(args.gpu_independent, [int(d) for d in args.devices.split(",")])
    else:
        scratch = Path(args.scratch)
        if scratch.exists():
            raise SystemExit("scratch exists: %s" % scratch)
        scratch.mkdir(parents=True)
        only = set(args.only.split(",")) if args.only else None
        for name, fn in (("unit", lambda: (test_unit(), test_held_requirement_a(), test_cpu_columns(scratch))),
                         ("workers", lambda: test_workers(scratch)), ("e2e", lambda: test_e2e(scratch))):
            if only is None or name in only:
                try:
                    fn()
                except Exception:
                    import traceback
                    report("%s: exception" % name, False, traceback.format_exc()[-3000:])
    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    print("SUMMARY: %d tests, %d failed" % (len(RESULTS), n_fail), flush=True)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
