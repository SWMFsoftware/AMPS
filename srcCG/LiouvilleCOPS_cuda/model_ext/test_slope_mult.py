#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_slope_mult.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_slope_mult.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b1dabddb8759437224f9cd00a8529e44f5f1f46f790a826c9dcf928ca60e0f3e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. CPU tests of the slope-multiplier and bump models (bins, floor, CPU stand-in vs FD, absent option,
#   arguments, end to end through the real driver chain, the fine tuning's worker requests, checks 1-8 with wrong-G /
#   wrong-J cases) and the GPU-job helpers --gpu-independent (unpatched independent evaluation), --compare-runs,
#   --blas-bench. 2026-11-10.
# USED BY: jobs/slope_mult_cpu.sbatch STEP=tests (this version, v8: CPU job 124062, 32 tests, 0 failed),
#   jobs/slope_bumps_test.sbatch, jobs/slope_mult_threads_test.sbatch.
# ==== end of header ====
"""Tests of the slope-multiplier model (model_ext/slope_mult.py, slope_mult_problem.py, slope_mult_fit.py; 2026-11-10).

CPU (a Slurm CPU job; ~3-5 min, < 4 GB):
    .venv/bin/python -W ignore test_slope_mult.py --scratch DIR
 1. bins: edges [0,40) ... [80,90), [90,180] (both sides of every edge), refusals (NaN, < 0, > 180); indicator basis
    round trip and refusal of a non-one-hot basis;
 2. floor: A = max(F, m L), active strictly above the floor (equality = at the floor), F = 0 where A_orig = 0;
 3. the CPU stand-in (slope_mult.cpu_model_class) on a synthetic base (test_positive_fit.make_base, 400 facets):
    every a, multiplier, T and kappa column vs central finite differences (facets that cross the floor within the
    step are counted and reported: the kink); b = None == the reference Model bitwise; floor 0 + all m = 1 == the
    reference clipped model (pressures and a columns) bitwise;
 4. absent option: importing slope_mult / slope_mult_problem changes no evaluator, transport or fit_cops attribute;
    the existing drivers' modules are byte-identical to their recorded copies (RUN_COPIES of the frozen clip snapshot);
 5. arguments: inserted options, refusals;
 6. end to end with the REAL driver chain slope_mult_fit.main -> grain_fit_fixed -> grain_fit -> gpu_fit.main (only
    gpu_fit.drive replaced by a CPU stand-in with the CPU model), each fit in a subprocess, on the synthetic grid base:
    source D2 (clip_refit.py --variant D2) -> slope fit from D2 (start check, full Jacobian check incl. every
    multiplier, products, positive_table.py), floor-0 check vs an independent evaluation (reference Model b=None +
    W10x5 + grain + data), resume (segment 2 reproduces segment 1's cost).
 7. workers (2026-11-10, after job 123985): the fine tuning's worker requests (slopemul_seen / _override / _sens) through
    the real MultiGpuBackend wrappers, _send/_receive, a pickling pipe stand-in and the patched WorkerState.handle, in a
    fresh process with the patch installed; only the three CudaBackend methods are stubs (their installed names are
    checked); the same check on the v6 snapshot must reproduce job 123985's AttributeError.

GPU (inside a GPU job; jobs/slope_mult_test.sbatch):
    --gpu-independent OUT_DIR   independent evaluation of the plain clipped harmonic model (UNPATCHED process: the
                                original GPU model with b = None, 8 local GPUs, W10x5 by positive_problem.install_window,
                                grain and data as the fitter) at the state of OUT_DIR/start_check_model.npz; compares the
                                per-row surface pressures and the cost with the driver's start check (fit_manifest.json).
"""

import argparse
import json
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
os.environ.setdefault("COPS_DATA_ROOT", "/data001/highresolu_vt/cops_cluster_transfer_20261026/"
                                        "rosetta_data_end_of_mission_cops")
OUT_ROOT = Path("/data001/highresolu_vt/cops_output")
RESULTS = []
FIXT = ["--fix", "Tmax=182.09922", "--fix", "T75=136.11149", "--fix", "Tmin=133.0"]


def report(name, ok, detail):
    RESULTS.append((name, ok, detail))
    print("%-66s %s  %s" % (name, "PASS" if ok else "FAIL", detail), flush=True)


# ------------------------------------------------------------------------------------------------ 1-2 bins, floor
def test_bins():
    import slope_mult as sm
    s = np.array([0.0, 39.999999, 40.0, 49.9, 50.0, 59.99, 60.0, 69.9, 70.0, 79.9, 80.0, 89.999, 90.0, 120.0, 180.0])
    want = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6], dtype=np.int8)
    got = sm.bin_index(s)
    ok = np.array_equal(got, want) and got.dtype == np.int8
    bad = []
    for v in ([np.nan], [-1e-9], [180.0000001]):
        try:
            sm.bin_index(np.array(v))
            bad.append(v)
        except ValueError:
            pass
    B = sm.indicator_basis(want)
    rt = np.array_equal(sm.bins_from_basis(B), want)
    nonhot = []
    for M in (np.column_stack([B[:, :6], B[:, 6] * 0.5]), np.column_stack([B[:, :6], B[:, 6] + B[:, 0]]),
              B[:, :6]):
        try:
            sm.bins_from_basis(M)
            nonhot.append(M.shape)
        except ValueError:
            pass
    report("bins: edges, refusals, indicator round trip", ok and not bad and rt and not nonhot,
           "got %s; not refused %s; round trip %s; non-one-hot accepted %s" % (got.tolist(), bad, rt, nonhot))


def test_floor():
    import slope_mult as sm
    L = np.array([2.0, -1.0, 1e-6, 1e-5, 3.0, -2.0, 0.0])
    bins = np.array([0, 1, 2, 2, 3, 4, 5], dtype=np.int8)
    F = 1e-5 * np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    m = np.array([1.0, 2.0, 1.0, 0.0, 1.0, 3.0, 1.0])          # per BIN
    s, A, act = sm.activity(L, m, bins, F)
    wantA = np.array([2.0, 1e-5, 1e-5, 1e-5, 0.0, 0.0, 0.0])
    wantact = np.array([True, False, False, False, False, False, False])
    ok = np.array_equal(A, wantA) and np.array_equal(act, wantact)
    report("floor: A = max(F, m L), active = s > F (equality at the floor)", ok,
           "A %s active %s (s = F = 1e-5 at index 3: at the floor; m = 0 (index 4) or F = 0 with L <= 0: A = 0)"
           % (A.tolist(), act.astype(int).tolist()))


# ------------------------------------------------------------------------------------------------ chunked L
def test_chunked_L():
    """User 2026-11-10 'Bitwise chunked version': slope_mult.chunked_matvec / _matvec_L (mode chunked) == the 1-thread
    Y @ a bitwise, for random a (signed, the fitter's magnitudes) and N in {2, 8, 16, 32, 64}, with n mod 4 tails 0..3,
    odd chunk ends and n smaller than one chunk; chunk bounds cover [0, n) with starts multiples of 4096; OpenBLAS stays
    at 1 thread; absent option = the plain call."""
    import slope_mult as sm
    from concurrent.futures import ThreadPoolExecutor
    rng = np.random.default_rng(23)
    info0 = sm.blas_info()
    sizes = (4096 * 128 + 2, 4096 * 37 + 3, 100001, 4096 * 16, 5000, 4097, 8193, 4096 * 33 + 2, 4096 * 8 + 1, 3, 1)
    bad, nchk, bounds_bad = [], 0, []
    for n in sizes:
        Y = np.ascontiguousarray(rng.normal(0.0, 0.3, (n, 25)))
        Y[:, 0] = 1.0
        for trial in range(3):
            a = rng.normal(0.0, 1.0, 25) * 10.0 ** rng.uniform(16, 21, 25) * rng.choice([-1.0, 1.0], 25)
            ref = Y @ a
            for N in (2, 8, 16, 32, 64):
                b = sm.chunk_bounds(n, N)
                if (b[0][0] != 0 or b[-1][1] != n or any(c0 % sm.L_ALIGN_ROWS for c0, _ in b) or len(b) > N or
                        any(b[k][1] != b[k + 1][0] for k in range(len(b) - 1)) or
                        (len(b) > 1 and min(c1 - c0 for c0, c1 in b) < sm.L_ALIGN_ROWS)):
                    bounds_bad.append((n, N, b[:3]))
                with ThreadPoolExecutor(N) as pool:
                    got = sm.chunked_matvec(Y, a, pool, N)
                nchk += 1
                if not np.array_equal(got, ref) or got.dtype != ref.dtype:
                    bad.append((n, N, trial, int(np.sum(got != ref))))
    # the driver's path: set_L_threads + _matvec_L, then back to the plain call
    n = 4096 * 64 + 1
    Y = np.ascontiguousarray(rng.normal(0.0, 0.3, (n, 25)))
    a = rng.normal(0.0, 1e18, 25)
    plain = sm._matvec_L(Y, a)                        # absent option: the plain Y @ a
    rec = sm.set_L_threads(16, "chunked")
    via = sm._matvec_L(Y, a)
    info1 = sm.blas_info()
    sm._STATE.update(L_mode=None, L_threads=None)
    sm._STATE["L_pool"].shutdown()
    ok = (not bad and not bounds_bad and np.array_equal(plain, Y @ a) and np.array_equal(via, plain) and
          info0["threads"] == 1 and info1["threads"] == 1 and rec["bitwise_by_construction"])
    report("chunked L == 1-thread L (bitwise), bounds, OpenBLAS 1 thread", ok,
           "%d comparisons (n %s, N 2..64, 3 random a each); differing %s; bad bounds %s; driver path equal %s; OpenBLAS "
           "threads before/after %s/%s" % (nchk, list(sizes), bad[:4], bounds_bad[:2], np.array_equal(via, plain),
                                           info0["threads"], info1["threads"]))


# ------------------------------------------------------------------------------------------------ 3 CPU model
def synthetic_model(root, floor_ratio, seed=5):
    import test_positive_fit as tpf
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    run, bundle, points = tpf.make_base(Path(root), seed=seed)
    header, names, rows = fit_cops.read_model_file(str(run / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    el = fit_cops.read_elements(str(run / "elements.bin"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    import slope_mult as sm
    Cls = sm.cpu_model_class(fit_cops, floor_ratio)
    M = Cls(str(run), header, cols, rows, el, coef)
    Mref = fit_cops.Model(str(run), header, cols, rows, el, coef)
    return fit_cops, M, Mref, el, coef, run, bundle, points


def test_cpu_model(scratch):
    import slope_mult as sm
    fit_cops, M, Mref, el, coef, run, _, _ = synthetic_model(Path(scratch) / "cpu_model", sm.FLOOR_RATIO)
    rng = np.random.default_rng(3)
    n = len(el["area"])
    bins = rng.integers(0, 7, n).astype(np.int8)
    B = sm.indicator_basis(bins)
    a = coef[0] * rng.uniform(0.5, 1.5, 25)
    m = rng.uniform(0.5, 2.0, 7)
    T = np.array([182.09922, 136.11149, 133.0])
    kappa = 0.7
    F = sm.FLOOR_RATIO * el["AH2O"]
    pNG, pRG, J = M.evaluate(a, T, kappa, b=m, enhancement_basis=B, enhancement_scale=2.0, jac=True,
                             basis_index=list(range(25)), enhancement_index=list(range(7)))

    def P(a_, m_, T_, k_):
        return np.concatenate(M.evaluate(a_, T_, k_, b=m_, enhancement_basis=B, enhancement_scale=2.0))
    worst = {"a": 0.0, "b": 0.0, "T": 0.0, "kappa": 0.0}
    kink = {}
    s0 = m[bins] * (el["Y"] @ a)
    for kind, count in (("a", 25), ("b", 7), ("T", 3), ("kappa", 1)):
        for j in range(count):
            if kind == "a":
                h = 1e-4 * max(abs(a[j]), 1.0)          # piecewise linear: a larger step, crossings counted
                ap, am = a.copy(), a.copy()
                ap[j] += h
                am[j] -= h
                fd = (P(ap, m, T, kappa) - P(am, m, T, kappa)) / (2 * h)
                sp, sm_ = m[bins] * (el["Y"] @ ap), m[bins] * (el["Y"] @ am)
            elif kind == "b":
                h = 1e-4
                mp, mm = m.copy(), m.copy()
                mp[j] += h
                mm[j] -= h
                fd = (P(a, mp, T, kappa) - P(a, mm, T, kappa)) / (2 * h)
                sp, sm_ = mp[bins] * (el["Y"] @ a), mm[bins] * (el["Y"] @ a)
            elif kind == "T":
                h = 1e-4
                Tp, Tm = T.copy(), T.copy()
                Tp[j] += h
                Tm[j] -= h
                fd = (P(a, m, Tp, kappa) - P(a, m, Tm, kappa)) / (2 * h)
                sp = sm_ = s0
            else:
                h = 1e-6
                fd = (P(a, m, T, kappa + h) - P(a, m, T, kappa - h)) / (2 * h)
                sp = sm_ = s0
            crossing = int(np.sum((sp > F) != (sm_ > F)))
            col = np.concatenate(J[(kind, j)])
            err = float(np.max(np.abs(col - fd)) / max(float(np.max(np.abs(fd))), 1e-300))
            if crossing:
                kink["%s%d" % (kind, j)] = (crossing, err)
            else:
                worst[kind] = max(worst[kind], err)
    ok = all(v < 1e-6 for v in worst.values())
    report("CPU model: a, m, T, kappa columns vs central FD", ok,
           "worst rel. error %s; columns whose FD step crosses the floor (excluded; kink): %s" % (
               {k: "%.1e" % v for k, v in worst.items()}, {k: "%d facets, err %.1e" % v for k, v in kink.items()}))
    # b None == reference, bitwise; floor 0 and m = 1 == reference clipped model, bitwise
    r1 = M.evaluate(a, T, kappa, jac=True, basis_index=list(range(25)))
    r0 = Mref.evaluate(a, T, kappa, jac=True, basis_index=list(range(25)))
    same = (np.array_equal(r1[0], r0[0]) and np.array_equal(r1[1], r0[1]) and
            all(np.array_equal(r1[2][k][g], r0[2][k][g]) for k in r0[2] for g in (0, 1)))
    report("CPU model: b = None == the reference Model (bitwise)", same, "pressures and every column")
    M0 = sm.cpu_model_class(fit_cops, 0.0)(str(run), *M_args(fit_cops, run, el, coef))
    q1 = M0.evaluate(a, T, kappa, b=np.ones(7), enhancement_basis=B, jac=True, basis_index=list(range(25)),
                     enhancement_index=list(range(7)))
    same0 = (np.array_equal(q1[0], r0[0]) and np.array_equal(q1[1], r0[1]) and
             all(np.array_equal(q1[2][("a", j)][g], r0[2][("a", j)][g]) for j in range(25) for g in (0, 1)) and
             all(np.array_equal(q1[2][("T", k)][g], r0[2][("T", k)][g]) for k in range(3) for g in (0, 1)))
    report("CPU model: floor 0, m = 1 == the plain clipped model (bitwise)", same0,
           "pressures, a, T columns equal to fit_cops.Model with b = None")


def M_args(fit_cops, run, el, coef):
    header, names, rows = fit_cops.read_model_file(str(Path(run) / "COPS_model.dat"))
    return header, {n: i for i, n in enumerate(names)}, rows, el, coef


# ------------------------------------------------------------------------------------------------ 4 absent option
def test_absent_option():
    from cops_cuda import evaluator, multi_gpu, transport, reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    B, MG = evaluator.CudaBackend, multi_gpu.MultiGpuBackend
    before = {n: getattr(B, n) for n in ("__init__", "activity_linear", "_upload_activity", "_upload_columns",
                                         "identity", "release", "close")}
    before.update({"MG." + n: getattr(MG, n) for n in ("activity_linear", "_evaluate", "identity")})
    before["pipe_worker_main"] = transport.pipe_worker_main
    before["fit_cops"] = {n: getattr(fit_cops, n) for n in ("Model", "Problem", "least_squares", "write_fit",
                                                             "activity_diagnostics", "save_active_manifest")}
    import slope_mult  # noqa: F401
    import slope_mult_problem  # noqa: F401
    after = {n: getattr(B, n) for n in ("__init__", "activity_linear", "_upload_activity", "_upload_columns",
                                        "identity", "release", "close")}
    after.update({"MG." + n: getattr(MG, n) for n in ("activity_linear", "_evaluate", "identity")})
    after["pipe_worker_main"] = transport.pipe_worker_main
    after["fit_cops"] = {n: getattr(fit_cops, n) for n in before["fit_cops"]}
    same = all(before[k] is after[k] for k in before if k != "fit_cops") and all(
        before["fit_cops"][k] is after["fit_cops"][k] for k in before["fit_cops"])
    report("absent option: importing the modules patches nothing", same, "evaluator/multi_gpu/transport/fit_cops "
           "attributes identical objects")
    # the existing drivers' modules are those that ran D2 (frozen clip snapshot) or unchanged since
    import hashlib
    snap = OUT_ROOT / "model_ext_frozen/clip_20261108T234013Z"
    diffs = {}
    for line in (snap / "SHA256SUMS").read_text().splitlines():
        h, name = line.split()
        if name == "positive_problem.py" or name.startswith("test_"):
            continue                          # edited 2026-11-09 for the positive fits (functions used: identical)
        path = HERE / name if (HERE / name).is_file() else OUT_ROOT / "model_ext" / name   # a snapshot has a subset
        now = hashlib.sha256(path.read_bytes()).hexdigest()
        if now != h:
            diffs[name] = now
    report("absent option: clip_refit chain modules = D2's frozen copies", not diffs,
           "differ: %s" % diffs if diffs else "clip_refit.py, clip_problem.py, grain_*.py, fixed_parameters.py, "
                                               "late_response_eval.py identical (positive_problem.py: see W10x5 check)")


# ------------------------------------------------------------------------------------------------ 5 arguments
def test_arguments(scratch):
    import slope_mult_problem as smp
    pt = Path(scratch) / "pt.txt"
    pt.write_text("# index UTC kind original_index L2_start L2_stop\n")
    bundle = Path(scratch) / "bins.npz"
    bundle.write_bytes(b"x")
    rs = Path(scratch) / "r.npz"
    rs.write_bytes(b"x")
    tail = ["--", "/run", "--lambda", "0", "--degree", "4", "--enhancement-degree", "4"]
    head = ["--point-table", str(pt), "--slope-bundle", str(bundle), "--start", "srccg", "--start-restart", str(rs)] + FIXT
    cfg, argv, notes = smp.rewrite_arguments(head + ["--devices", "0"] + tail)
    h = argv[:argv.index("--")]
    t = argv[argv.index("--") + 1:]
    ok = ("--template-start" in h and h[h.index("--template-start") + 1] == "0" and "c_slope_00_40=1.0" in h and
          t.count("--activity-template") == 7 and t[-1] == "slope_90_180" and cfg["floor_ratio"] == 1e-5 and
          cfg["transport"] == "local" and cfg["start_choice"] == "srccg")
    report("arguments: inserted --template-start 0, --fix c_slope_00_40, 7 templates", ok, "head %s" % h)
    cfg2, argv2, _ = smp.rewrite_arguments(head + ["--transport", "mpi"] + tail)
    ok2 = "--activity" in argv2 and argv2[argv2.index("--activity") + 1] == "host"
    report("arguments: --transport mpi -> --activity host inserted", ok2, str(argv2[:argv2.index("--")][-6:]))
    bad = [(head + tail + ["--activity-template", "x", "y"], "activity-template"),
           (head + ["--floor-ratio", "1e-4"] + tail, "not a decision"),
           (head + ["--floor-ratio", "0"] + tail, "not a decision"),
           (head + ["--test-floor-zero"] + tail, "needs --floor-ratio 0"),
           (head + ["--fix", "c_slope_00_40=2"] + tail, "held at"),
           (head + ["--fix", "c_slope_40_50=1"] + tail, "only Tmax"),
           (head + ["--activity", "device"] + tail, "activity device"),
           (["--point-table", str(pt), "--slope-bundle", str(bundle)] + FIXT + tail, "give --start srccg|d2"),
           (head[:8] + ["--fix", "T75=133"] + head[10:] + tail, "temperatures must be held"),
           (head[:4] + head[6:] + tail, "there is no default"),
           (head[:5] + ["exp"] + head[6:] + tail, "choose one of"),
           (head + ["--template-start", "1"] + tail, "template-start"),
           (head + tail + ["--fix-kappa"], "fix-kappa"),
           (head + ["--", "/run", "--lambda", "1e-3", "--degree", "4", "--enhancement-degree", "4"], "--lambda 0")]
    fails = []
    for argv_, why in bad:
        try:
            smp.rewrite_arguments(argv_)
            fails.append(why)
        except SystemExit as exc:
            if why not in str(exc):
                fails.append("%s: %s" % (why, exc))
    report("arguments: refusals", not fails, "not refused / wrong message: %s" % fails if fails else
           "%d refused" % len(bad))


# ------------------------------------------------------------------------------------------------ 6 end to end
def write_bins_bundle(path, el, bins, gravity_meta=None):
    """A slope_bins.py-style bundle for a synthetic mesh (same metadata keys the driver reads)."""
    import slope_mult as sm
    from cops_cuda import reference
    reference.bind()
    at = reference.load("activity_templates")
    gc, gt = reference.load("gravity_cache"), reference.load("gravity_topology")
    ah2o = np.ascontiguousarray(np.asarray(el["AH2O"], np.float64))
    meta = {"templates": {n: {"definition": "synthetic slope bin %s (test; effective_gravitational_slope)" % n}
                          for n in sm.BIN_NAMES},
            "gravity": gravity_meta or {"cache": {"definition": {"algorithm": gc.ALGORITHM,
                                                                 "orientation_policy": gt.ORIENTATION_POLICY}},
                                        "integrity_check": "DISABLE after topology audit"},
            "slope_multiplier_bins": {"names": list(sm.BIN_NAMES), "edges_deg": list(sm.BIN_EDGES_DEG),
                                      "nfaces": int(len(bins)), "counts": [int(v) for v in np.bincount(bins, minlength=7)],
                                      "bins_sha256": sm.sha256_bytes(bins), "ah2o_sha256": sm.sha256_bytes(ah2o),
                                      "mesh_identity": at.mesh_identity(el["centroid"], el["area"])}}
    at.write_bundle(str(path), list(sm.BIN_NAMES), sm.indicator_basis(bins), el["centroid"], el["area"], meta)


def cpu_drive(args, argv, fit_args_, devices, max_nfev, template_start, comm, deadline=None):
    """gpu_fit.drive without a GPU (test_clip_refit.cpu_drive) with fit_cops.Model = the slope model's CPU stand-in when
    the slope patch configuration is installed (slope_mult._STATE), else the reference Model."""
    import test_clip_refit as tcr
    import slope_mult as sm
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    cfg = sm.installed_config()
    if cfg is not None:
        fit_cops.Model = sm.cpu_model_class(fit_cops, cfg["floor_ratio"], nbumps=int(cfg.get("n_bumps", 0)),
                                            held_g=sm.held_factor(cfg)[0])
    return tcr.cpu_drive(args, argv, fit_args_, devices, max_nfev, template_start, comm, deadline)


def child(spec):
    import gpu_fit
    gpu_fit.drive = cpu_drive
    if spec["mode"] == "slope":
        import slope_mult as sm
        sm.install_gpu = lambda cfg: sm._STATE.update(installed=True, config=dict(cfg), role="cpu-test") or True
        import slope_mult_fit
        sys.argv = [str(HERE / "slope_mult_fit.py")] + spec["argv"]
        return slope_mult_fit.main()
    if spec["mode"] == "clip":
        import clip_refit
        sys.argv = [str(HERE / "clip_refit.py")] + spec["argv"]
        return clip_refit.main()
    import grain_fit_fixed
    sys.argv = [str(HERE / "grain_fit_fixed.py")] + spec["argv"]
    return grain_fit_fixed.main()


def run_child(root, name, mode, argv):
    log = Path(root) / (name + ".log")
    spec = json.dumps({"mode": mode, "argv": [str(a) for a in argv]})
    rc = subprocess.call([sys.executable, "-W", "ignore", __file__, "--child", spec], stdout=log.open("w"),
                         stderr=subprocess.STDOUT)
    return rc, log.read_text()


def slope_fit_args(run, out):
    return [str(run), "--allow-legacy-base", "--degree", "4", "--enhancement-degree", "4", "--lambda", "0", "--jobs",
            "1", "--tbounds", "30", "400", "--ftol", "1e-10", "--xtol", "1e-10", "--gtol", "1e-10",
            "--check-tolerance", "inf", "--out", str(out)]


def independent_cost(run, point_table, npz_path, x_by_name, t_end_et, model_class, configure=None, bumps=None):
    """The plain clipped harmonic model (model_class with b = None: unpatched), W10x5 by positive_problem.install_window
    on the all-rows model, grain and data exactly as grain_crossing / fit_cops. Returns (cost, pNG_s, pRG_s, info)."""
    import positive_problem as pp
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    header, names, rows = fit_cops.read_model_file(str(Path(run) / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    el = fit_cops.read_elements(str(Path(run) / "elements.bin"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    if configure is not None:
        configure()
    M = model_class(str(run), header, cols, rows, el, coef)
    with np.load(str(npz_path), allow_pickle=False) as z:
        d = {k: np.array(z[k]) for k in z.files}
    if not np.array_equal(np.asarray(M.rows), d["rows"]):
        raise SystemExit("independent: the start-check rows are not all rows of the base")
    table = pp.read_point_table(point_table)
    layout = pp.LateLayout(table, pp.read_run_points(M.run), M.ipoint)
    pp.install_window(M, layout)
    a = coef[0].copy()
    ascale = float(np.max(np.abs(coef[0])))
    for j in range(25):
        a[j] = x_by_name["a%02d" % j] * ascale
    T = np.array([182.09922, 136.11149, 133.0])
    kappa = x_by_name["kappa"]
    if bumps is None:
        pNG_s, pRG_s = M.evaluate(a, T, kappa)
    else:      # the plain clipped model with the bump columns as its only templates, a = 0: max(0, A_scale Phi w)
        Phi, w, scale = bumps
        pNG_s, pRG_s = M.evaluate(np.zeros_like(a), T, kappa, b=np.asarray(w, float),
                                  enhancement_basis=np.ascontiguousarray(Phi), enhancement_scale=scale)
    A, B = x_by_name["grain_A_uPa"] * 1e-6, x_by_name["grain_B_uPa_per_100s"] * (1e-6 / 100.0)
    et = d["et"]
    s = et - t_end_et
    on = s >= 0.0
    g = np.maximum(np.where(on, A - B * np.maximum(s, 0.0), 0.0), 0.0)
    useNG, useRG = d["useNG"].astype(bool), d["useRG"].astype(bool)
    sel = np.asarray(M.rows)
    lnNG, lnRG = np.log(M.refNG[sel][useNG]), np.log(M.refRG[sel][useRG])
    wNG, wRG = 1.0 / np.sqrt(useNG.sum()), 1.0 / np.sqrt(useRG.sum())
    pNG = pNG_s + g
    r = np.concatenate([wNG * (np.log(np.maximum(pNG[useNG], fit_cops.FLOOR)) - lnNG),
                        wRG * (np.log(np.maximum(pRG_s[useRG], fit_cops.FLOOR)) - lnRG)])
    cost = 0.5 * float(r @ r)
    rel = lambda u, v: float(np.max(np.abs(u / v - 1.0)))         # noqa: E731
    info = {"surface_pNG_max_rel_diff": rel(pNG_s, d["surface_pNG"]), "surface_pRG_max_rel_diff": rel(pRG_s, d["surface_pRG"]),
            "surface_bitwise": bool(np.array_equal(pNG_s, d["surface_pNG"]) and np.array_equal(pRG_s, d["surface_pRG"])),
            "grain_bitwise": bool(np.array_equal(g, d["grain"])), "NNG": int(useNG.sum()), "NRG": int(useRG.sum()),
            "late_observations": int(len(layout.late_rows)), "grid_rows": int(len(layout.grid_rows))}
    return cost, pNG_s, pRG_s, info


def manifest_record(out):
    m = json.loads((Path(out) / "fit_manifest.json").read_text())
    return m, (m.get("model_extension") or {}).get("slope_multiplier") or {}


def x_names(out):
    m, rec = manifest_record(out)
    core = rec.get("clip_core") or {}
    return m, rec, core


def test_e2e(scratch):
    import test_clip_refit as tcr
    import slope_mult as sm
    root = Path(scratch) / "e2e"
    run_g, run_n, bundle, ft, points = tcr.make_bases(root)
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    el = fit_cops.read_elements(str(run_g / "elements.bin"))
    rng = np.random.default_rng(17)
    slope = rng.uniform(0.0, 180.0, len(el["area"]))
    bins = sm.bin_index(slope)
    bins_npz = root / "slope_bins.npz"
    write_bins_bundle(bins_npz, el, bins)
    gff_head = FIXT + ["--fix", "c_ft_c1=1", "--crossing-weight", "10", "--devices", "0", "--template-start", "0"]
    S = root / "S0"
    rc, text = run_child(root, "S0", "gff", gff_head + ["--max-nfev", "8", "--"] + tcr.fit_args(run_n, S, bundle, ft))
    if rc:
        report("e2e source S0", False, text[-3000:])
        return
    D2 = root / "D2"
    rc, text = run_child(root, "D2", "clip", ["--variant", "D2", "--point-table", points, "--start-restart",
                                              S / "fit_restart.npz"] + gff_head + ["--max-nfev", "6", "--"]
                         + tcr.fit_args(run_g, D2, bundle, ft))
    if rc:
        report("e2e source D2 (clip_refit D2 on the grid base)", False, text[-3000:])
        return
    report("e2e source D2 (clip_refit D2 on the grid base)", True, "")
    head = ["--point-table", points, "--slope-bundle", bins_npz, "--start", "d2", "--start-restart",
            D2 / "fit_restart.npz"] + FIXT + ["--crossing-weight", "10", "--check-jacobian", "--devices", "0"]
    # a. the fit from D2: start check, full Jacobian check of every multiplier, products
    A1 = root / "slope_seg01"
    check = ",".join(["c_" + n for n in sm.BIN_NAMES[1:]] + ["a00", "a02", "a07", "kappa", "grain_A_uPa"])
    rc, t1 = run_child(root, "slope_seg01", "slope", head[:8] + ["--check-full-jacobian", "--check-parameters", check]
                       + head[8:] + ["--max-nfev", "5", "--"] + slope_fit_args(run_g, A1))
    if rc:
        report("e2e slope fit from D2", False, "exit %d\n%s" % (rc, t1[-4000:]))
        return
    m1, rec1, core1 = x_names(A1)
    jc = core1["jacobian_check"]["columns"]
    worst = {k: v.get("data") for k, v in jc.items() if isinstance(v, dict)}
    st = core1["start_check"]
    start = rec1["start_restart"]
    products = all((A1 / f).is_file() for f in ("slope_mult_facets.npz", "slope_mult_bins.txt",
                                                "SLOPE_MULTIPLIER_MODEL.txt", "late_window_model.dat", "evaluations.txt"))
    note = [l for l in open(A1 / "COPS_model.dat") if l.startswith("# SlopeMultiplierModel")]
    with np.load(str(A1 / "fit_restart.npz"), allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
    pt_out = A1 / "positive_table_out.dat"
    rc_pt = subprocess.call([sys.executable, str(OUT_ROOT / "figures_cluster/scripts/positive_table.py"), str(A1),
                             str(pt_out)], stdout=subprocess.DEVNULL)
    with np.load(str(A1 / "slope_mult_facets.npz"), allow_pickle=False) as z:
        fac = {k: np.array(z[k]) for k in z.files}
    ok = (all(v is not None and v < 1e-5 for v in worst.values()) and products and len(note) == 1 and rc_pt == 0 and
          meta["model_extension"]["parameterization"] == "slope_multiplier_floor_w10x5_v1" and
          start["multipliers_set_to"] == 1.0 and np.array_equal(fac["bin"], bins) and
          st["residuals_vs_direct_model_max_abs"] == 0.0)
    report("e2e a: fit from D2: Jacobian check (multipliers, a, kappa, grain), products", ok,
           "Jacobian data worst %s; start cost %.9e (D2 final %.9e); products %s; table note %d; positive_table.py "
           "exit %d; marker %s" % ({k: "%.1e" % v for k, v in worst.items()}, st["start"]["cost"],
                                   st["source"]["cost"], products, len(note), rc_pt,
                                   meta["model_extension"]["parameterization"]))
    for line in [l for l in t1.splitlines() if l.startswith("slope_mult:") or l.startswith("clip_refit:")][:14]:
        print("   | " + line[:240])
    bins_txt = (A1 / "slope_mult_bins.txt").read_text().splitlines()
    print("   | " + "\n   | ".join(l[:200] for l in bins_txt if not l.startswith("# ")))
    # b. floor 0 + all m = 1 == the plain clipped model: start check vs an independent evaluation
    Z = root / "slope_floor0"
    rc, t2 = run_child(root, "slope_floor0", "slope", head[:8] + ["--floor-ratio", "0", "--test-floor-zero",
                                                                 "--stop-after-checks"] + head[8:] +
                       ["--max-nfev", "5", "--"] + slope_fit_args(run_g, Z))
    if rc:
        report("e2e b: floor 0 start check", False, "exit %d\n%s" % (rc, t2[-3000:]))
    else:
        m2, rec2, core2 = x_names(Z)
        xs = rec2["start_restart"]["x_start"]
        names = rec2["model_check"]["layout"]
        xb = dict(zip(names, xs))
        cost, _, _, info = independent_cost(run_g, points, Z / "start_check_model.npz", xb,
                                            rec2["grain_window"]["t_end_et"], fit_cops.Model)
        c_drv = core2["start_check"]["start"]["data_cost"]
        ok = cost == c_drv and info["surface_bitwise"] and info["grain_bitwise"]
        report("e2e b: floor 0, m = 1 == plain clipped model (independent)", ok,
               "driver start data cost %.15e, independent %.15e (rel %.1e); %s" % (c_drv, cost, cost / c_drv - 1, info))
    # c. resume: segment 2 from segment 1 reproduces its cost
    A2 = root / "slope_seg02"
    rc, t3 = run_child(root, "slope_seg02", "slope", head[:4] + head[8:] + ["--restart-from", A1 / "fit_restart.npz",
                                                                           "--max-nfev", "3", "--"]
                       + slope_fit_args(run_g, A2))
    if rc:
        report("e2e c: resume", False, "exit %d\n%s" % (rc, t3[-3000:]))
    else:
        m3, rec3, core3 = x_names(A2)
        rr = core3["start_check"].get("resume_reproduction") or {}
        report("e2e c: resume (segment 2 reproduces segment 1's cost)", bool(rr.get("ok")),
               "%s; chain %s" % (rr, list((rec3.get("fit_chain") or {}).keys())))
    # d. refusals at run time: a D2 restart as --restart-from (resume) is refused
    B4 = root / "slope_bad_resume"
    rc, t4 = run_child(root, "slope_bad_resume", "slope", head[:4] + head[8:] + ["--restart-from", D2 / "fit_restart.npz",
                                                                               "--max-nfev", "3", "--"]
                       + slope_fit_args(run_g, B4))
    report("e2e d: a D2 restart as --restart-from is refused", rc != 0 and "not a segment of this fit" in t4,
           "exit %d" % rc)
    # e. --start srccg (production, user 2026-11-10): a = the fitter's prior x0, multipliers 1, kappa 1, D2's grain
    S1 = root / "slope_srccg_seg01"
    hs = [str(v) for v in head]
    hs[hs.index("d2")] = "srccg"
    rc, t5 = run_child(root, "slope_srccg_seg01", "slope", hs[:8] + ["--check-full-jacobian", "--check-parameters", check]
                       + hs[8:] + ["--max-nfev", "4", "--"] + slope_fit_args(run_g, S1))
    if rc:
        report("e2e e: --start srccg", False, "exit %d\n%s" % (rc, t5[-4000:]))
        return
    m5, rec5, core5 = x_names(S1)
    sr = rec5["start_restart"]
    names = rec5["model_check"]["layout"]
    xs = dict(zip(names, sr["x_start"]))
    with np.load(str(D2 / "fit_restart.npz"), allow_pickle=False) as z:
        dmeta = json.loads(str(z["meta_json"]))
        dx = dict(zip(dmeta["names"], [float(v) for v in z["x"]]))
    header, cn, rws = fit_cops.read_model_file(str(run_g / "COPS_model.dat"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    ascale = float(np.max(np.abs(coef[0])))
    prior_a = [coef[0][j] / ascale for j in range(25)]
    rep_ = rec5.get("start_report") or {}
    with np.load(str(S1 / "fit_restart.npz"), allow_pickle=False) as z:
        meta5 = json.loads(str(z["meta_json"]))
    txt = (S1 / "SLOPE_MULTIPLIER_MODEL.txt").read_text()
    jc5 = core5["jacobian_check"]["columns"]
    worst5 = max(v["data"] for v in jc5.values() if isinstance(v, dict))
    ok = (sr["start_choice"] == "srccg" and sr["a_equals_prior_x0"] and
          [xs["a%02d" % j] for j in range(25)] == prior_a and
          all(xs["c_" + n] == 1.0 for n in sm.BIN_NAMES[1:]) and xs["kappa"] == 1.0 and
          xs["grain_A_uPa"] == dx["grain_A_uPa"] and xs["grain_B_uPa_per_100s"] == dx["grain_B_uPa_per_100s"] and
          (sr.get("grain_term_check") or {}).get("A_ok") and rep_.get("start_choice") == "srccg" and
          rep_.get("cost") == core5["start_check"]["start"]["cost"] and 0.0 <= rep_["floor_share"]["area"] <= 1.0 and
          meta5["model_extension"]["start_choice"] == "srccg" and "--start srccg" in txt and worst5 < 1e-5 and
          "start check (start srccg)" in t5)
    report("e2e e: --start srccg: prior a, m = 1, kappa 1, D2's grain; recorded", ok,
           "start cost %.9e; floor share at the start facets %.4f area %.4f seen %.4f; grain A %r (D2 %r); a vs coef0 %.1e; "
           "Jacobian worst %.1e; marker start %s" % (
               rep_.get("cost", float("nan")), rep_["floor_share"]["facets"], rep_["floor_share"]["area"],
               rep_["floor_share"].get("seen_area", float("nan")), xs["grain_A_uPa"], dx["grain_A_uPa"],
               sr["a_vs_run_coefficients_max_rel_diff"], worst5, meta5["model_extension"]["start_choice"]))
    for line in [l for l in t5.splitlines() if l.startswith("slope_mult: start")][:3]:
        print("   | " + line[:260])
    # f. resume of the srccg fit with --start d2 is refused; with --start srccg it reproduces the cost
    B6 = root / "slope_srccg_bad"
    rc, t6 = run_child(root, "slope_srccg_bad", "slope", head[:4] + ["--start", "d2"] + head[8:] +
                       ["--restart-from", S1 / "fit_restart.npz", "--max-nfev", "3", "--"] + slope_fit_args(run_g, B6))
    S2 = root / "slope_srccg_seg02"
    rc2, t7 = run_child(root, "slope_srccg_seg02", "slope", head[:4] + ["--start", "srccg"] + head[8:] +
                        ["--restart-from", S1 / "fit_restart.npz", "--max-nfev", "3", "--"] + slope_fit_args(run_g, S2))
    rr = (x_names(S2)[2]["start_check"].get("resume_reproduction") or {}) if rc2 == 0 else {}
    # g. --l-threads 4 (chunked; plumbing on the CPU stand-in: records, start-check line) == e's fit; --blas-threads refused
    G = root / "slope_srccg_lthreads"
    rc8, t8 = run_child(root, "slope_srccg_lthreads", "slope", hs[:8] + ["--l-threads", "4", "--check-full-jacobian",
                                                                       "--check-parameters", check]
                        + hs[8:] + ["--max-nfev", "4", "--"] + slope_fit_args(run_g, G))
    rc9, t9 = run_child(root, "slope_srccg_blasthreads", "slope", hs[:8] + ["--blas-threads", "4"] + hs[8:] +
                        ["--max-nfev", "4", "--"] + slope_fit_args(run_g, root / "slope_srccg_blasthreads"))
    if rc8 == 0:
        a8, m8 = (lambda z: ({k: np.array(z[k]) for k in z.files if k != "meta_json"}, json.loads(str(z["meta_json"]))))(
            np.load(str(G / "fit_restart.npz"), allow_pickle=False))
        a5, m5_ = (lambda z: ({k: np.array(z[k]) for k in z.files if k != "meta_json"}, json.loads(str(z["meta_json"]))))(
            np.load(str(S1 / "fit_restart.npz"), allow_pickle=False))
        rec8 = x_names(G)[1]
        okg = (np.array_equal(a8["x"], a5["x"]) and m8["optimizer"]["cost"] == m5_["optimizer"]["cost"] and
               (rec8.get("l_step") or {}).get("mode") == "chunked" and rec8["l_step"]["threads"] == 4 and
               "mode chunked, threads 4" in t8 and rc9 != 0 and "--blas-threads" in t9 and "replaced" in t9)
    else:
        okg = False
    report("e2e g: --l-threads 4 (chunked) == without; recorded; --blas-threads refused", okg,
           "exit %d / %d; %s" % (rc8, rc9, [l for l in t8.splitlines() if "L = Y a" in l][:1]))
    report("e2e f: resume keeps the start choice (d2 refused, srccg reproduces)",
           rc != 0 and "the fit started with --start 'srccg'" in t6 and bool(rr.get("ok")),
           "refusal exit %d; resume exit %d %s" % (rc, rc2, rr))



# ------------------------------------------------------------------------------------------------ 7 bumps, fine tuning
HDR = np.dtype([("version", "<i4"), ("point", "<i4"), ("geometry", "<f8", (10,)), ("count", "<i4")])


def make_late_only(run, facets, replacement, window_s=300.0):
    """Synthetic base: the records of `facets` at every point before [t_end - window_s, t_end] are removed, so that
    these facets are seen only inside the late window (t_end = the last used data point)."""
    import test_positive_fit as tpf
    import grain_crossing as gc
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    header, names, rows = fit_cops.read_model_file(str(Path(run) / "COPS_model.dat"))
    c = {n: i for i, n in enumerate(names)}
    ets = gc.et_by_point(run)
    used = [(int(r[c["useNG"]]) > 0 and float(r[c["refNG"]]) > 0) or (int(r[c["useRG"]]) > 0 and float(r[c["refRG"]]) > 0)
            for r in rows]
    t_end = max(ets[int(r[c["iPoint"]])] for r, u in zip(rows, used) if u)
    changed = 0
    for r in rows:
        ip = int(r[c["iPoint"]])
        if ets[ip] >= t_end - window_s:
            continue
        path = Path(run) / "kernels" / ("point=%04d.bin" % ip)
        raw = path.read_bytes()
        h = np.frombuffer(raw[:HDR.itemsize], HDR)
        rec = np.frombuffer(raw[HDR.itemsize:], tpf.RECORD).copy()
        sel = np.isin(rec["face"], facets)
        if sel.any():                         # drop these records (a point lists every face once)
            rec = rec[~sel]
            h2 = h.copy()
            h2["count"] = len(rec)
            changed += int(sel.sum())
            path.write_bytes(h2.tobytes() + rec.tobytes())
    return t_end, changed


class CpuModelEngine:
    """slope_finetune's engine interface on the CPU stand-in (tests)."""

    def __init__(self, ctx):
        import slope_mult as sm
        from cops_cuda import reference
        reference.bind()
        fit_cops = reference.load("fit_cops")
        g = np.exp(ctx.held_prev) if (ctx.model_cfg.get("held") or {}) else None
        Cls = sm.cpu_model_class(fit_cops, ctx.model_cfg["floor_ratio"], nbumps=int(ctx.model_cfg.get("n_bumps", 0)),
                                 held_g=g)
        self.M = Cls(str(ctx.run), ctx.header, ctx.cols, ctx.table_rows, ctx.el, ctx.coef)
        self.ctx = ctx
        self.calls = {"forward": 0}

    def forward(self):
        c = self.ctx
        self.calls["forward"] += 1
        return self.M.evaluate(c.a, c.T, c.kappa, b=c.b, enhancement_basis=c.basis, enhancement_scale=c.escale)

    def seen(self, mask, structural=False):
        M = self.M
        sel = np.asarray(mask, bool)[M.pt]
        if structural:
            on = (M.Kn != 0) | (M.KPhi != 0) | (M.KV != 0)
        else:
            on = (M.Kn > 0) | (M.KPhi > 0) | (M.KV > 0)
        out = np.zeros(M.nfaces, bool)
        out[M.face[sel & on]] = True
        return out

    def sensitivities(self, facets):
        return self.M.sensitivities(self.ctx.T, facets)

    def forward_override(self, facets, g):
        self.M.set_override(facets, g)
        try:
            return self.forward()
        finally:
            self.M.set_override(None, None)


def test_ft_checks_catch_errors(sft, C1, FT, root, late_f):
    """Validation checks 6 and 7 with the pressure-sized steps (user 2026-11-10 "Fix the check step") still catch a
    deliberately wrong derivative: (a) G scaled by 1 + 1e-3 on one facet that check 7 examines -> check 7 fails at
    ~1e-3; (b) the same on a facet check 7 does not examine -> check 7 passes, check 4 (host model vs direct evaluation
    at a random delta) fails, the validation fails; (c) the host J column of one facet that check 6 examines scaled by
    1 + 1e-3 (HostLinearEngine.jv patched in this process only) -> check 6 fails at ~1e-3. The facets come from the
    correct run's validation.json (same seed, same selection)."""
    v0 = json.loads((FT / "validation.json").read_text())
    f7 = [d["facet"] for d in v0["7_G_vs_gpu_fd"]["per_facet"] if "rel" in d]
    f6 = [d["facet"] for d in v0["6_single_facets_host"]["per_facet"] if "rel" in d]
    report("e2e bumps: checks 6/7 with pressure-sized steps pass on the correct model", v0["6_single_facets_host"][
        "passed"] and v0["7_G_vs_gpu_fd"]["passed"], "check 6 worst %.2e (tau %g, facets %s), check 7 worst %.2e (tau %g, "
        "facets %s)" % (v0["6_single_facets_host"]["worst_rel"], v0["6_single_facets_host"]["tau"], f6,
                        v0["7_G_vs_gpu_fd"]["worst_rel"], v0["7_G_vs_gpu_fd"]["tau"], f7))

    class WrongG(CpuModelEngine):
        facet = None

        def sensitivities(self, facets):
            ng, rg = super().sensitivities(facets)
            j = [int(f) for f in np.asarray(facets)].index(int(self.facet))
            ng, rg = np.array(ng, copy=True), np.array(rg, copy=True)
            ng[:, j] *= 1.0 + 1e-3
            rg[:, j] *= 1.0 + 1e-3
            return ng, rg

    def run_wrong(name, facet=None, jv_facet=None):
        out = root / ("c1_ft_wrong_" + name)
        args = sft.parse_args(["--fit", str(C1), "--out", str(out), "--facets", "late_only", "--window-min", "5",
                               "--validate-only"])
        orig_jv = sft.HostLinearEngine.jv

        def bad_jv(self, st, k):
            q = int(np.flatnonzero(self.idx == jv_facet)[0])
            saved = [float(b[q]) for b in self.bufs[:k]]
            for b in self.bufs[:k]:
                b[q] *= 1.0 + 1e-3
            try:
                return orig_jv(self, st, k)
            finally:
                for b, v in zip(self.bufs[:k], saved):
                    b[q] = v

        def factory(ctx):
            e = WrongG(ctx) if facet is not None else CpuModelEngine(ctx)
            e.facet = facet
            return e
        if jv_facet is not None:
            sft.HostLinearEngine.jv = bad_jv
        try:
            rc = sft.run(args, engine_factory=factory, command=["test"])
        except SystemExit as exc:
            rc = str(exc)
        finally:
            sft.HostLinearEngine.jv = orig_jv
        return rc, json.loads((out / "validation.json").read_text())
    if f7:
        rc, v = run_wrong("g_checked", facet=f7[0])
        w7 = v["7_G_vs_gpu_fd"]["worst_rel"]
        report("e2e bumps: check 7 catches G x (1 + 1e-3) on a checked facet", not v["passed"] and
               not v["7_G_vs_gpu_fd"]["passed"] and 5e-4 <= w7 <= 2e-3, "facet %d: check 7 worst %.2e, check 4 %s "
               "(%.1e), validation passed %s (%s)" % (f7[0], w7, v["4_host_vs_gpu_random_delta"]["passed"],
                                                       v["4_host_vs_gpu_random_delta"]["max_rel"], v["passed"], rc))
    other = [f for f in late_f if f not in f7]
    if other:
        rc, v = run_wrong("g_unchecked", facet=other[0])
        report("e2e bumps: a wrong G on a facet check 7 does not examine is caught by check 4", not v["passed"] and
               v["7_G_vs_gpu_fd"]["passed"] and not v["4_host_vs_gpu_random_delta"]["passed"],
               "facet %d: check 7 %s (%.1e), check 4 %s (%.1e), validation passed %s" % (
                   other[0], v["7_G_vs_gpu_fd"]["passed"], v["7_G_vs_gpu_fd"]["worst_rel"],
                   v["4_host_vs_gpu_random_delta"]["passed"], v["4_host_vs_gpu_random_delta"]["max_rel"], v["passed"]))
    if f6:
        rc, v = run_wrong("j_checked", jv_facet=f6[0])
        w6 = v["6_single_facets_host"]["worst_rel"]
        report("e2e bumps: check 6 catches a J column x (1 + 1e-3) on a checked facet", not v["passed"] and
               not v["6_single_facets_host"]["passed"] and 5e-4 <= w6 <= 2e-3,
               "facet %d: check 6 worst %.2e, validation passed %s" % (f6[0], w6, v["passed"]))


def test_e2e_bumps(scratch):
    """Bump + slope-multiplier model and its fine-tuning cycle on a synthetic base (CPU stand-ins): nnls start, full
    Jacobian check of every weight and multiplier, products (bump table), floor-0 check vs an independent evaluation,
    late-facet fine tuning with its validation, cycle-2 fit with the held term (start cost = the fine tuning's final
    data cost), resume of the held fit."""
    import test_clip_refit as tcr
    import slope_mult as sm
    import slope_bumps
    import slope_finetune as sft
    root = Path(scratch) / "bumps"
    run_g, run_n, bundle, ft, points = tcr.make_bases(root)
    late_f = [3, 7, 11, 13, 17, 19, 23, 29]
    t_end, changed = make_late_only(run_g, late_f, 399)
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    el = fit_cops.read_elements(str(run_g / "elements.bin"))
    rng = np.random.default_rng(31)
    bins = sm.bin_index(rng.uniform(0.0, 180.0, len(el["area"])))
    bins_npz = root / "slope_bins.npz"
    write_bins_bundle(bins_npz, el, bins)
    seen_npz = root / "seen.npz"
    np.savez(str(seen_npz), seen=(np.arange(len(el["area"])) % 3 != 0))
    bumps_npz = root / "bumps6.npz"
    rc = slope_bumps.main(["--run", str(run_g), "--seen", str(seen_npz), "--out", str(bumps_npz), "--n-seen", "4",
                           "--n-unseen", "2", "--threads", "2", "--max-iter", "50"])
    gff_head = FIXT + ["--fix", "c_ft_c1=1", "--crossing-weight", "10", "--devices", "0", "--template-start", "0"]
    S = root / "S0"
    rcs, text = run_child(root, "S0", "gff", gff_head + ["--max-nfev", "8", "--"] + tcr.fit_args(run_n, S, bundle, ft))
    D2 = root / "D2"
    rcd, text = run_child(root, "D2", "clip", ["--variant", "D2", "--point-table", points, "--start-restart",
                                               S / "fit_restart.npz"] + gff_head + ["--max-nfev", "6", "--"]
                          + tcr.fit_args(run_g, D2, bundle, ft))
    if rc or rcs or rcd:
        report("e2e bumps: sources", False, "bumps %s S0 %s D2 %s\n%s" % (rc, rcs, rcd, text[-2000:]))
        return
    report("e2e bumps: late-only facets made, bump bundle, sources", True,
           "%d records moved off %d facets before t_end - 5 min" % (changed, len(late_f)))
    head = ["--point-table", points, "--slope-bundle", bins_npz, "--bump-bundle", bumps_npz, "--start", "nnls",
            "--start-restart", D2 / "fit_restart.npz"] + FIXT + ["--crossing-weight", "10", "--check-jacobian",
                                                                  "--devices", "0"]
    C1 = root / "c1_fit"
    rc, t1 = run_child(root, "c1_fit", "slope", head[:10] + ["--check-full-jacobian"] + head[10:] +
                       ["--max-nfev", "4", "--"] + slope_fit_args(run_g, C1))
    if rc:
        report("e2e bumps: c1 fit (nnls start)", False, "exit %d\n%s" % (rc, t1[-4000:]))
        return
    m1, rec1, core1 = x_names(C1)
    jc = core1["jacobian_check"]["columns"]
    worst = max(v["data"] for v in jc.values() if isinstance(v, dict))
    with np.load(str(C1 / "fit_restart.npz"), allow_pickle=False) as z:
        meta1 = json.loads(str(z["meta_json"]))
    names = rec1["model_check"]["layout"]
    ok = (worst < 1e-5 and len([k for k in jc if k.startswith("c_bump_")]) == 6 and
          len([k for k in jc if k.startswith("c_slope_")]) == 6 and not any(n.startswith("a") and n[1:].isdigit()
                                                                             for n in names) and
          (C1 / "slope_bumps_table.txt").is_file() and
          meta1["model_extension"]["parameterization"] == "bump_slope_multiplier_floor_w10x5_v1" and
          rec1["start_restart"]["start_choice"] == "nnls")
    report("e2e bumps: c1 fit: Jacobian check of every weight + multiplier, products, marker", ok,
           "worst data %.1e over %d columns; layout %s ... %s; start %s" % (worst, len(jc), names[:3], names[-3:],
                                                                          rec1["start_restart"]["start_choice"]))
    # floor 0 + m = 1 == the plain clipped model with the bump columns as templates (independent)
    Z = root / "c1_floor0"
    rc, t2 = run_child(root, "c1_floor0", "slope", head[:10] + ["--floor-ratio", "0", "--test-floor-zero",
                                                                "--stop-after-checks"] + head[10:] +
                       ["--max-nfev", "4", "--"] + slope_fit_args(run_g, Z))
    if rc:
        report("e2e bumps: floor 0 == plain model (independent)", False, t2[-3000:])
    else:
        m2, rec2, core2 = x_names(Z)
        xb = dict(zip(rec2["model_check"]["layout"], rec2["start_restart"]["x_start"]))
        for j in range(25):
            xb["a%02d" % j] = 0.0
        from cops_cuda import figures as fg
        bmeta = fg.bundle_metadata(str(bumps_npz))["slope_bumps"]
        at = reference.load("activity_templates")
        Phi = at.select_templates(at.load_bundle(str(bumps_npz), el["centroid"], el["area"]), bmeta["names"])
        w = [xb["c_" + n] for n in bmeta["names"]]
        cost, _, _, info = independent_cost(run_g, points, Z / "start_check_model.npz", xb,
                                            rec2["grain_window"]["t_end_et"], fit_cops.Model,
                                            bumps=(Phi, w, rec2["model_check"]["A_scale"]))
        c_drv = core2["start_check"]["start"]["data_cost"]
        report("e2e bumps: floor 0, m = 1 == plain clipped model with bump templates", cost == c_drv and
               info["surface_bitwise"], "driver %.15e independent %.15e; %s" % (c_drv, cost, info))
    # late-facet fine tuning of c1 (CPU engine), with validation
    FT = root / "c1_ft"
    args = sft.parse_args(["--fit", str(C1), "--out", str(FT), "--facets", "late_only", "--window-min", "5",
                           "--max-nfev", "20", "--validate", "--l-threads", "2"])
    rc = sft.run(args, engine_factory=lambda ctx: CpuModelEngine(ctx), command=["test"])
    ftj = json.loads((FT / "fine_tuning.json").read_text())
    with np.load(str(FT / "finetune_state.npz"), allow_pickle=False) as z:
        st_ = {k: np.array(z[k]) for k in z.files if k != "meta_json"}
        smeta = json.loads(str(z["meta_json"]))
    v = json.loads((FT / "validation.json").read_text())
    nz = np.flatnonzero(st_["delta"])
    ok = (rc == 0 and v["passed"] and set(np.flatnonzero(st_["late"]).tolist()) == set(late_f) and
          np.array_equal(st_["held_next"], st_["held_prev"] + st_["delta"]) and set(nz.tolist()) <= set(late_f) and
          smeta["costs_final"]["cost_total"] <= smeta["costs_initial"]["cost_total"])
    report("e2e bumps: fine tuning (late_only): validation 1-8, held_next, parameters = the late facets", ok,
           "exit %s; late %s; parameters %d (at floor %d); checks %s; cost %.9e -> %.9e; start vs fit rel %.1e; host "
           "vs model at random delta %.1e; strict support %s" % (
               rc, np.flatnonzero(st_["late"]).tolist(), ftj["seen"]["parameters"], ftj["seen"]["late_facets_at_floor"],
               {k: vv for k, vv in ftj["validation"].items()}, smeta["costs_initial"]["cost_total"],
               smeta["costs_final"]["cost_total"], v["1_start_cost"]["rel_total"],
               v["4_host_vs_gpu_random_delta"]["max_rel"], v["5_strict_support"]))
    # bitwise with 1 thread
    FT1 = root / "c1_ft_1thread"
    args1 = sft.parse_args(["--fit", str(C1), "--out", str(FT1), "--facets", "late_only", "--window-min", "5",
                            "--max-nfev", "20", "--validate"])
    sft.run(args1, engine_factory=lambda ctx: CpuModelEngine(ctx), command=["test"])
    with np.load(str(FT1 / "finetune_state.npz"), allow_pickle=False) as z:
        st1 = {k: np.array(z[k]) for k in z.files if k != "meta_json"}
    same = all(np.array_equal(st_[k], st1[k]) for k in ("delta", "held_next", "grain_A_uPa", "grain_B_uPa_per_100s"))
    report("e2e bumps: fine tuning with 2 threads == 1 thread (bitwise)", same, "delta, held_next, grain")
    test_ft_checks_catch_errors(sft, C1, FT, root, late_f)
    # cycle-2 fit: held term, chain start; start cost = the fine tuning's final data cost
    C2 = root / "c2_fit"
    h2 = [str(x) for x in head]
    h2[h2.index("nnls")] = "chain"
    h2[h2.index(str(D2 / "fit_restart.npz"))] = str(C1 / "fit_restart.npz")
    h2 = h2[:10] + ["--held", str(FT / "finetune_state.npz"), "--grain-from", str(FT / "finetune_state.npz")] + h2[10:]
    rc, t3 = run_child(root, "c2_fit", "slope", h2 + ["--max-nfev", "3", "--"] + slope_fit_args(run_g, C2))
    if rc:
        report("e2e bumps: c2 fit (chain, held)", False, "exit %d\n%s" % (rc, t3[-4000:]))
        return
    m3, rec3, core3 = x_names(C2)
    c2_start = core3["start_check"]["start"]["data_cost"]
    ft_final = smeta["costs_final"]["cost_data"]
    rel_ = abs(c2_start / ft_final - 1.0)
    with np.load(str(C2 / "slope_mult_facets.npz"), allow_pickle=False) as z:
        has_held = "held" in z.files and np.array_equal(np.asarray(z["held"]), st_["held_next"])
    report("e2e bumps: c2 fit start cost = the fine tuning's final data cost (1e-8); held term in the products",
           rel_ <= 1e-8 and has_held, "c2 start %.15e, fine tuning final %.15e (rel %.1e); held recorded %s" % (
               c2_start, ft_final, rel_, has_held))
    # resume of the held fit
    C3 = root / "c2_fit_seg02"
    h3 = [x for x in h2]
    i = h3.index("--start-restart")
    del h3[i:i + 2]
    i = h3.index("--grain-from")
    del h3[i:i + 2]
    rc, t4 = run_child(root, "c2_fit_seg02", "slope", h3 + ["--restart-from", str(C2 / "fit_restart.npz"),
                                                            "--max-nfev", "2", "--"] + slope_fit_args(run_g, C3))
    rr = (x_names(C3)[2]["start_check"].get("resume_reproduction") or {}) if rc == 0 else {}
    report("e2e bumps: resume of the held fit reproduces its cost", bool(rr.get("ok")), "exit %d %s %s" % (
        rc, rr, "" if rc == 0 else t4[-1500:]))

# ------------------------------------------------------------------------------------------------ GPU independent
def gpu_independent(out_dir, devices):
    """Unpatched process: the original GPU model (b = None) at the start state of OUT_DIR."""
    import slope_mult as sm
    if sm._STATE["installed"]:
        raise SystemExit("the slope patch must not be installed in the independent evaluation")
    out = Path(out_dir)
    m, rec, core = x_names(out)
    xb = dict(zip(rec["model_check"]["layout"], rec["start_restart"]["x_start"]))
    run = rec["config"]["run"]
    pts = rec["config"]["point_table"]
    from cops_cuda import gpu_model, precision, reference
    reference.bind()
    os.environ[precision.RECORD_PRECISION_ENVIRONMENT] = "fp32"            # gpu_fit --precision mixed

    def configure():
        gpu_model.REGISTRY.configure(devices[0], "tree", 4096, devices=devices, transport="local", activity="host")
    import time
    t0 = time.perf_counter()
    try:
        cost, pNG, pRG, info = _gpu_cost(run, pts, out, xb, rec, configure)
    finally:
        stats = gpu_model.REGISTRY.stats()
        gpu_model.REGISTRY.close()
    c_drv = core["start_check"]["start"]["data_cost"]
    rel = cost / c_drv - 1.0
    result = {"driver_start_data_cost": c_drv, "independent_cost": cost, "rel_diff": rel, "bitwise_cost": cost == c_drv,
              "info": info, "seconds": time.perf_counter() - t0, "gpu_calls": stats["calls"],
              "driver_dir": str(out), "multipliers_at_start": [xb.get("c_" + n, 1.0) for n in sm.BIN_NAMES]}
    path = out.parent / ("independent_%s.json" % out.name)
    path.write_text(json.dumps(result, indent=2) + "\n")
    report("GPU (a): floor 0, m = 1 == plain clipped harmonic model (independent evaluation)",
           abs(rel) <= 1e-12 and info["surface_pNG_max_rel_diff"] <= 1e-12 and info["surface_pRG_max_rel_diff"] <= 1e-12,
           "driver start data cost %.15e, independent %.15e, rel %.2e (bitwise %s); surface bitwise %s, max rel NG %.1e "
           "RG %.1e; grain bitwise %s; %d NG / %d RG; %s" % (
               c_drv, cost, rel, cost == c_drv, info["surface_bitwise"], info["surface_pNG_max_rel_diff"],
               info["surface_pRG_max_rel_diff"], info["grain_bitwise"], info["NNG"], info["NRG"], path))


def _gpu_cost(run, pts, out, xb, rec, configure):
    from cops_cuda import gpu_model, reference
    Cls = gpu_model.gpu_model_class()
    bumps = None
    if rec["model_check"].get("basis") == "bumps":      # the plain model with the bump columns as its only templates
        reference.bind()
        fit_cops, at = reference.load("fit_cops"), reference.load("activity_templates")
        el = fit_cops.read_elements(str(Path(run) / "elements.bin"))
        bpath = rec["config"]["bump_bundle"]
        names = list(rec["config"]["bump_names"])
        Phi = np.ascontiguousarray(at.select_templates(at.load_bundle(bpath, el["centroid"], el["area"]), names))
        del el
        bumps = (Phi, [xb["c_" + n] for n in names], rec["model_check"]["A_scale"])
        for j in range(25):
            xb.setdefault("a%02d" % j, 0.0)
    return independent_cost(run, pts, out / "start_check_model.npz", xb, rec["grain_window"]["t_end_et"], Cls,
                            configure, bumps=bumps)


# ------------------------------------------------------------------------------------------------ BLAS threads
def blas_bench(run, threads, reps=7):
    """L = Y a with numpy's OpenBLAS at each thread count (the fitter's arrays: fit_cops.read_elements' contiguous Y and
    the run's coefficients); median of `reps` after one warm-up; bitwise comparison with the 1-thread result."""
    import time
    import slope_mult as sm
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    header, names, rows = fit_cops.read_model_file(str(Path(run) / "COPS_model.dat"))
    coef = fit_cops.read_activity_coefficients(str(REF), header)
    el = fit_cops.read_elements(str(Path(run) / "elements.bin"))
    Y, a = el["Y"], np.asarray(coef[0], float)
    out = {"run": str(run), "shape": list(Y.shape), "c_contiguous": bool(Y.flags.c_contiguous), "results": []}
    ref = None
    for n in threads:
        sm.set_blas_threads(n)
        L = Y @ a
        ts = []
        for _ in range(reps):
            t0 = time.perf_counter()
            L = Y @ a
            ts.append(time.perf_counter() - t0)
        if ref is None:
            ref = L.copy()
        r = {"threads": n, "threads_read_back": sm.blas_info()["threads"], "median_s": float(np.median(ts)),
             "min_s": float(np.min(ts)), "bitwise_equal_to_first": bool(np.array_equal(L, ref)),
             "GB_per_s": float(Y.nbytes / np.median(ts) / 1e9)}
        out["results"].append(r)
        print("blas_bench: %2d threads: L = Y a median %.4f s (min %.4f), %.1f GB/s, bitwise equal to %d threads: %s"
              % (n, r["median_s"], r["min_s"], r["GB_per_s"], threads[0], r["bitwise_equal_to_first"]), flush=True)
    sm.set_blas_threads(1)
    from concurrent.futures import ThreadPoolExecutor
    out["chunked"] = []
    for n in threads:
        with ThreadPoolExecutor(n) as pool:
            L = sm.chunked_matvec(Y, a, pool, n)
            ts = []
            for _ in range(reps):
                t0 = time.perf_counter()
                L = sm.chunked_matvec(Y, a, pool, n)
                ts.append(time.perf_counter() - t0)
        r = {"threads": n, "median_s": float(np.median(ts)), "min_s": float(np.min(ts)),
             "bitwise_equal_to_1_thread": bool(np.array_equal(L, ref)), "chunks": len(sm.chunk_bounds(len(Y), n)),
             "GB_per_s": float(Y.nbytes / np.median(ts) / 1e9)}
        out["chunked"].append(r)
        print("blas_bench: chunked %2d threads (OpenBLAS 1, %d chunks of 4096-aligned rows): median %.4f s (min %.4f), "
              "%.1f GB/s, bitwise equal to 1 thread: %s" % (n, r["chunks"], r["median_s"], r["min_s"], r["GB_per_s"],
                                                            r["bitwise_equal_to_1_thread"]), flush=True)
    out["blas"] = sm.blas_info()
    return out


def _evals(seg):
    rows = [l.split() for l in open(Path(seg) / "evaluations.txt") if not l.startswith("#")]
    return [(r[1], r[3]) for r in rows], [float(r[2]) for r in rows if r[1] == "fun"], [float(r[2]) for r in rows
                                                                                         if r[1] == "jac"]


def _table(seg):
    return [l for l in open(Path(seg) / "COPS_model.dat") if l.strip() and not l.startswith("#")]


def compare_runs(ref, runs):
    """Bitwise comparison of short fits (same command but --blas-threads) with the reference: every evaluation's cost
    (evaluations.txt), the final x, the restart's predictions, COPS_model.dat rows, the start cost, the Jacobian check;
    timing per evaluation / Jacobian, the L = Y a step (driver backend identity), the processes' BLAS threads, CPU load."""
    out = {"reference": str(ref), "runs": []}
    re_, rf, rj = _evals(ref)
    with np.load(str(Path(ref) / "fit_restart.npz"), allow_pickle=False) as z:
        rx = {k: np.array(z[k]) for k in ("x", "pred_ng", "pred_rg", "activity_linear")}
    rt = _table(ref)
    rm, rrec, rcore = x_names(ref)
    for run in runs:
        e, f, j = _evals(run)
        with np.load(str(Path(run) / "fit_restart.npz"), allow_pickle=False) as z:
            x = {k: np.array(z[k]) for k in ("x", "pred_ng", "pred_rg", "activity_linear")}
        m, rec, core = x_names(run)
        diff_evals = [i for i, (u, v) in enumerate(zip(re_, e)) if u != v]
        reg = m["compute_backend"]["registry"]
        mext = [(b.get("model_extension") or {}).get("slope_multiplier") or {} for b in reg["backends"]]
        lt = mext[0].get("activity_linear_timing") if mext else None
        workers = [((w.get("model_extension") or {}).get("slope_multiplier") or {}).get("backend") or {}
                   for b in reg["backends"] for w in b.get("workers", [])]
        wthreads = sorted({(w.get("process_blas") or {}).get("threads") for w in workers})
        jc_r = {k: v.get("data") for k, v in rcore["jacobian_check"]["columns"].items() if isinstance(v, dict)}
        jc = {k: v.get("data") for k, v in core["jacobian_check"]["columns"].items() if isinstance(v, dict)}
        load = Path(run).parent.glob("provenance-*/process_load.txt")
        pc, nl = [], []
        for lf in load:
            for line in open(lf):
                v = line.split()
                if len(v) >= 7 and v[0].isdigit() and "python" in v[6]:
                    pc.append(float(v[3]))
                    nl.append(int(v[2]))
        ls_ = rec.get("l_step") or {}
        r = {"run": str(run), "blas": "%s %s" % (ls_.get("mode"), ls_.get("threads")) if ls_ else
             (rec.get("blas", {}).get("L_step") or {}).get("threads_L"),
             "evaluations_equal": len(re_) == len(e) and not diff_evals, "evaluations": len(e),
             "first_differing_evaluation": diff_evals[:3],
             "x_bitwise": bool(np.array_equal(rx["x"], x["x"])),
             "pred_bitwise": bool(np.array_equal(rx["pred_ng"], x["pred_ng"]) and np.array_equal(rx["pred_rg"], x["pred_rg"])),
             "x_max_abs_diff": float(np.max(np.abs(rx["x"] - x["x"]))),
             "pred_max_rel_diff": float(max(np.max(np.abs(x["pred_ng"] / rx["pred_ng"] - 1)),
                                            np.max(np.abs(x["pred_rg"] / rx["pred_rg"] - 1)))),
             "activity_linear_bitwise": bool(np.array_equal(rx["activity_linear"], x["activity_linear"])),
             "cops_model_rows_equal": rt == _table(run),
             "start_cost_equal": rcore["start_check"]["start"]["cost"] == core["start_check"]["start"]["cost"],
             "jacobian_check_equal": jc_r == jc,
             "fun_median_s": float(np.median(f)), "jac_median_s": float(np.median(j)),
             "ref_fun_median_s": float(np.median(rf)), "ref_jac_median_s": float(np.median(rj)),
             "L_step": lt, "L_mean_s": (lt["seconds"] / lt["calls"]) if lt and lt["calls"] else None,
             "worker_blas_threads": wthreads, "fitter_pcpu_max": max(pc) if pc else None,
             "python_nlwp_max": max(nl) if nl else None}
        out["runs"].append(r)
        ok = (r["evaluations_equal"] and r["x_bitwise"] and r["pred_bitwise"] and r["cops_model_rows_equal"] and
              r["start_cost_equal"] and r["jacobian_check_equal"])
        report("L step %s threads vs reference: bitwise" % r["blas"], ok,
               "evaluations equal %s (%d; first differing %s), x %s (max |diff| %.1e), predictions %s (max rel %.1e), "
               "activity_linear %s, COPS_model.dat %s, start cost %s, Jacobian check %s; per evaluation %.3f s (ref %.3f), "
               "per Jacobian %.3f s (ref %.3f); L step mean %s s over %s calls; worker BLAS threads %s; fitter CPU max "
               "%s %%, max threads/process %s" % (
                   r["evaluations_equal"], r["evaluations"], r["first_differing_evaluation"], r["x_bitwise"],
                   r["x_max_abs_diff"], r["pred_bitwise"], r["pred_max_rel_diff"], r["activity_linear_bitwise"],
                   r["cops_model_rows_equal"], r["start_cost_equal"], r["jacobian_check_equal"], r["fun_median_s"],
                   r["ref_fun_median_s"], r["jac_median_s"], r["ref_jac_median_s"],
                   "%.4f" % r["L_mean_s"] if r["L_mean_s"] is not None else "?", (lt or {}).get("calls"), wthreads,
                   r["fitter_pcpu_max"], r["python_nlwp_max"]))
    return out


# ------------------------------------------------------------------------------------------------ 7 worker requests
WORKER_CHILD = r"""
import json, pickle, sys, traceback
import numpy as np
moddir, app = sys.argv[1], sys.argv[2]
for p in (app, app + "/tools", moddir):
    sys.path.insert(0, p)
import slope_mult as sm
cfg = {"version": 1, "model": "slope_multiplier_floor", "test_only": True, "bin_names": list(sm.BIN_NAMES),
       "basis": "bumps", "n_bumps": 2, "bumps_sha256": "worker-handler-test", "nfaces": 12, "floor_ratio": 1e-5}
sm.install_gpu(cfg)
from cops_cuda import backend as bk, evaluator, multi_gpu, transport
from cops_cuda.errors import BackendError
B, MG, WS = evaluator.CudaBackend, multi_gpu.MultiGpuBackend, transport.WorkerState
requests = getattr(sm, "WORKER_REQUESTS", None) or {"slopemul_override": "slopemul_set_override",
                                                    "slopemul_seen": "slopemul_seen_rows",
                                                    "slopemul_sens": "slopemul_sensitivities"}
NF, NPT = 12, 6
out = {"module": sm.__file__, "installed_B": {k: callable(getattr(B, m, None)) for k, m in requests.items()},
       "installed_MG": {m: callable(getattr(MG, m, None)) for m in requests.values()}, "checks": {}, "errors": {}}
calls = []


def stub_backend(wid, positions):
    be = object.__new__(B)                      # the PATCHED CudaBackend class, no device: only the three methods
    def seen_rows(row_mask, structural=False):  # are replaced (instance attributes under the installed names)
        row_mask = np.asarray(row_mask, bool)
        calls.append(("seen", wid, row_mask.tolist(), bool(structural)))
        m = np.zeros(NF, bool)
        m[wid::3] = bool(row_mask.any())
        m[NF - 1 - wid] = True
        return m
    def set_override(facets, g_values):
        calls.append(("override", wid, None if facets is None else [int(f) for f in facets],
                      None if g_values is None else [float(v) for v in g_values]))
        return 0 if facets is None else len(facets) + 10 * wid
    def sensitivities(state, facets):
        calls.append(("sens", wid, [int(f) for f in facets], float(state.kappa)))
        ng = np.outer(np.asarray(positions, float) + 1.0, np.arange(1.0, len(facets) + 1.0)) * 10.0 ** wid
        return ng, -2.0 * ng
    be.slopemul_seen_rows, be.slopemul_set_override, be.slopemul_sensitivities = seen_rows, set_override, sensitivities
    return be


class Channel:                                   # the pipe: messages and replies pickled, WorkerState.handle in between
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
    got = mg.slopemul_seen_rows(row_mask, True)
    want = np.zeros(NF, bool)
    for wid in (0, 1):
        m = np.zeros(NF, bool); m[wid::3] = bool(row_mask[pos[wid]].any()); m[NF - 1 - wid] = True
        want |= m
    return (np.array_equal(got, want) and
            calls == [("seen", w, row_mask[pos[w]].tolist(), True) for w in (0, 1)])
check("seen", seen)
def override():
    calls.clear()
    a = mg.slopemul_set_override(np.array([1, 5], np.int64), np.array([2.0, 3.0]))
    b = mg.slopemul_set_override(None, None)
    return a == [2, 12] and b == [0, 0] and \
        calls == [("override", 0, [1, 5], [2.0, 3.0]), ("override", 1, [1, 5], [2.0, 3.0]),
                  ("override", 0, None, None), ("override", 1, None, None)]
check("override", override)
def sens():
    calls.clear()
    state = bk.ParameterState(np.zeros(25), np.array([182.0, 136.0, 133.0]), 1.5, np.ones(9))
    facets = np.array([3, 7, 9], np.int64)
    ng, rg = mg.slopemul_sensitivities(state, facets)
    ok = ng.shape == (NPT, 3) and rg.shape == (NPT, 3)
    for wid in (0, 1):
        w = np.outer(pos[wid] + 1.0, np.arange(1.0, 4.0)) * 10.0 ** wid
        ok = ok and np.array_equal(ng[pos[wid]], w) and np.array_equal(rg[pos[wid]], -2.0 * w)
    return ok and calls == [("sens", w, [3, 7, 9], 1.5) for w in (0, 1)]
check("sensitivities", sens)
def error_path():
    def boom(row_mask, structural=False):
        raise ValueError("worker-handler-test failure")
    mg.workers[1]["channel"].ws.backend.slopemul_seen_rows = boom
    try:
        mg.slopemul_seen_rows(row_mask, False)
    except BackendError as exc:
        return "device 1 failed during slopemul_seen" in str(exc) and "worker-handler-test failure" in str(exc)
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
print("RESULT " + json.dumps(out), flush=True)
sys.exit(0 if all(out["checks"].values()) and all(out["installed_B"].values()) and all(out["installed_MG"].values())
         else 1)
"""


def _worker_child(moddir, log):
    env = {k: v for k, v in os.environ.items() if not k.startswith("COPS_POS_") and k != "COPS_SLOPEMUL_FROZEN_CONFIG"}
    r = subprocess.run([sys.executable, "-W", "ignore", "-c", WORKER_CHILD, str(moddir), str(APP)], env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log.write_text(r.stdout)
    res = [json.loads(l[7:]) for l in r.stdout.splitlines() if l.startswith("RESULT ")]
    return r.returncode, (res[0] if res else None), r.stdout


def test_worker_handlers(scratch):
    """The fine tuning's worker requests through the REAL code path of a local multi-GPU backend, without a GPU: the
    patch installed in a fresh process (install_gpu), MultiGpuBackend's slopemul_* wrappers with the real _send /
    _receive, a pickling pipe stand-in and the patched transport.WorkerState.handle; only the three CudaBackend methods
    the requests end in are stubs on an instance of the patched class (under the names install_gpu installs, which are
    checked to exist on the class). Job 123985 failed here ('CudaBackend' object has no attribute 'seen_rows'); the
    same check is run against the v6 snapshot and must reproduce that error."""
    d = scratch / "worker_handlers"
    d.mkdir()
    rc, res, text = _worker_child(HERE, d / "current.log")
    ok = rc == 0 and res is not None
    report("worker requests: MG wrapper -> pipe -> WorkerState.handle -> backend", ok,
           "rc %s; installed on CudaBackend %s, on MultiGpuBackend %s; checks %s%s" % (
               rc, res and res["installed_B"], res and res["installed_MG"], res and res["checks"],
               "" if ok else "; errors " + (json.dumps(res["errors"])[:1500] if res else text[-1500:])))
    v6 = OUT_ROOT / "model_ext_frozen/slopemul_20261110T201212Z"
    if v6.is_dir():
        rc6, res6, text6 = _worker_child(v6, d / "v6.log")
        err = (res6 or {}).get("errors", {}).get("seen", "")
        rep = rc6 != 0 and res6 is not None and not res6["checks"].get("seen", True) and \
            "has no attribute 'seen_rows'" in err
        report("worker requests: the same check reproduces job 123985 on v6", rep,
               "v6 rc %s, seen %s, error tail: %s" % (rc6, res6 and res6["checks"].get("seen"),
                                                      err.strip().splitlines()[-1] if err.strip() else text6[-300:]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--blas-bench", help="base run: time L = Y a at --threads")
    p.add_argument("--threads", default="1,8,16,32,64")
    p.add_argument("--compare-runs", nargs="+", help="REF_SEG RUN_SEG [RUN_SEG ...]")
    p.add_argument("--json-out")
    p.add_argument("--scratch")
    p.add_argument("--child")
    p.add_argument("--gpu-independent")
    p.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    p.add_argument("--only", default=None, help="comma list of: bins,floor,chunked,cpu_model,absent,arguments,e2e,bumps,workers")
    args = p.parse_args()
    if args.child:
        return child(json.loads(args.child)) or 0
    if args.blas_bench or args.compare_runs:
        res = (blas_bench(args.blas_bench, [int(v) for v in args.threads.split(",")]) if args.blas_bench else
               compare_runs(args.compare_runs[0], args.compare_runs[1:]))
        if args.json_out:
            Path(args.json_out).write_text(json.dumps(res, indent=2) + "\n")
    elif args.gpu_independent:
        gpu_independent(args.gpu_independent, [int(d) for d in args.devices.split(",")])
    else:
        scratch = Path(args.scratch)
        if scratch.exists():
            raise SystemExit("scratch exists: %s" % scratch)
        scratch.mkdir(parents=True)
        only = set(args.only.split(",")) if args.only else None
        for name, fn in (("bins", lambda: test_bins()), ("floor", lambda: test_floor()),
                         ("chunked", lambda: test_chunked_L()),
                         ("cpu_model", lambda: test_cpu_model(scratch)), ("absent", lambda: test_absent_option()),
                         ("arguments", lambda: test_arguments(scratch)), ("e2e", lambda: test_e2e(scratch)),
                         ("bumps", lambda: test_e2e_bumps(scratch)),
                         ("workers", lambda: test_worker_handlers(scratch))):
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
