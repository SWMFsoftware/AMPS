#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/positive_activity_check.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/positive_activity_check.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 679045ded61b584a2f17f2af3cbfd9bc6567504caf3f7c5535a7d9fdebcf4d3a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test. Checks of positive_activity.py on a fitting base: kernels with the options off bitwise identical to
#   the unchanged ones; NumPy reference of night ratio / exp / floor; Jacobian against finite differences. user decisions 2026-11-07.
# USED BY: jobs/positive_activity_check.sbatch (job 122298).
# ==== end of header ====
"""Validation of model_ext/positive_activity.py on the SHAP8 12M fitting base (mixed precision, FP32 record store).

    .venv/bin/python positive_activity_check.py --out <dir> [--devices 0,1] [--points 0,300,...]

No fit is run. The check binds the GPU evaluator on a subset of the base's rows (default 20 points: 5 before the NG gap,
15 in the late window 10:33:10-10:39:26) and evaluates states of the cycle-3 fit
(runs/shap8_012m_flat_earlier/grain/cycles/c3/C30_n4000_grain_cross10_ftc3/fit_restart.npz: a, the template amplitudes
of cliff_ge_60deg and ft_c1..3 with their columns and A_scale, T = 228.903897/30/30, kappa 0.189178).
Each configuration runs in its own process (the extension is configured through the environment, once per process):

phase ref        COPS_POS_* unset, positive_activity not imported (the unpatched evaluator)
phase def        COPS_POS_ENABLE=1, every other setting at its default
phase mod        COPS_POS_ENABLE=1 NIGHT_RATIO=0.02 ACTIVITY=exp FLOOR_RATIO=0.02 (the user's settings)
phase nofloor    COPS_POS_ENABLE=1 NIGHT_RATIO=0.02 ACTIVITY=exp (no floor: Jacobian without floor switches)
phase clipfloor  COPS_POS_ENABLE=1 NIGHT_RATIO=0.02 FLOOR_RATIO=0.02 (clip activity with the floor)

Backends: cuda1_host (CudaBackend, device 0, host activity), cuda1_device (device activity, activity.cu),
multi2_host / multi2_device (MultiGpuBackend, spawned workers on 2 devices, local transport).

1. Identity: ref vs def, every backend, two states (the fit state and a perturbed one with T75 > Tmin): pressures and
   all 33 Jacobian columns (a0..a24, b0..b3, Tmax, T75, Tmin, kappa) compared BITWISE. Each worker reports the CUDA
   source it compiled (identity + stderr log).
2. Physics: GPU pressures of mod / nofloor / clipfloor (and ref, as a check of the NumPy reference itself) against an
   independent NumPy implementation from the base's kernel files (fit_cops.read_kernel, FP64 records) and from the FP32
   store records (the device inputs; isolates the formulas from the FP32 rounding). The exp state: a = least-squares
   projection of ln(AH2O) onto the 25 harmonics (area-weighted, facets with AH2O > 0); templates with A_scale = 1 and
   amplitudes 0.3, 0.5, 0.5, 0.5 divided by each column's max |B| (bounded exponent change; the fit's A_scale ~2e18
   would be in the exponent; stated simplification).
3. Jacobian: analytic columns of mod / nofloor / clipfloor (cuda1_host) against central finite differences of the GPU
   pressures, two step sizes, for a0, a1, a4, a6, a12, a20, b0 (cliff), b1 (ft_c1), Tmax, T75, Tmin, kappa.
Writes <out>/report.json, <out>/report.txt and the per-phase products (phase_*.npz/json/log) in <out>.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
sys.path.insert(0, str(HERE))

R12 = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_flat_earlier")
DEFAULTS = {"base": R12 / "base", "store": R12 / "store",
            "fit": R12 / "grain/cycles/c3/C30_n4000_grain_cross10_ftc3"}
DEFAULT_POINTS = (0, 300, 600, 900, 1030, 1055, 1062, 1070, 1078, 1086, 1094, 1102, 1110, 1118, 1126, 1134, 1140,
                  1144, 1147, 1149)
PHASES = {
    "ref": None,
    "def": {"COPS_POS_ENABLE": "1"},
    "mod": {"COPS_POS_ENABLE": "1", "COPS_POS_NIGHT_RATIO": "0.02", "COPS_POS_ACTIVITY": "exp",
            "COPS_POS_FLOOR_RATIO": "0.02"},
    "nofloor": {"COPS_POS_ENABLE": "1", "COPS_POS_NIGHT_RATIO": "0.02", "COPS_POS_ACTIVITY": "exp"},
    "clipfloor": {"COPS_POS_ENABLE": "1", "COPS_POS_NIGHT_RATIO": "0.02", "COPS_POS_FLOOR_RATIO": "0.02"},
}
# NumPy model settings of each phase: (night ratio or None = 0.2/0.1, activity, floor ratio, state)
MODEL = {"ref": (None, "clip", 0.0, "fit"), "def": (None, "clip", 0.0, "fit"), "mod": (0.02, "exp", 0.02, "exp"),
         "nofloor": (0.02, "exp", 0.0, "exp"), "clipfloor": (0.02, "clip", 0.02, "fit")}
IDENTITY_BACKENDS = ("cuda1_host", "cuda1_device", "multi2_host", "multi2_device", "cuda1_stream", "cuda1_fp64")
PHYSICS_BACKENDS = {"mod": ("cuda1_host", "cuda1_device", "multi2_host", "cuda1_stream"), "nofloor": ("cuda1_host",),
                    "clipfloor": ("cuda1_host", "multi2_host")}
FD_A = (0, 1, 4, 6, 12, 20)
FD_B = (0, 1)
FD_STEPS = (1e-3, 1e-4, 1e-6)    # exponent change (exp) or fraction of rms(act_lin) (clip); T: x 100 K; kappa: x 10
# cuda1_stream: records streamed from pinned host memory (SCALING_PLAN P2), forced with the evaluator's own test caps
STREAM_ENV = {"COPS_CUDA_MAX_RESIDENT_RECORDS": None, "COPS_CUDA_STREAM_CHUNK_RECORDS": str(1 << 22)}
EXP_TEMPLATE_SCALE = 1.0
# exp state: template amplitude k = EXP_TEMPLATE_EXPONENT[k] / max|B_k| (bounded exponent change; the ft_c columns are
# delta_activity / A_scale carried to the 12M facets and reach large values on small facets: an amplitude 0.2 overflowed
# exp at points 1147/1149 in the first smoke test, job 122293)
EXP_TEMPLATE_EXPONENT = (0.3, 0.5, 0.5, 0.5)


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist() if obj.size <= 64 else "<array %s %s>" % (obj.dtype, obj.shape)
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)


# ------------------------------------------------------------------------------------------------- shared inputs
def reference_model(base, point_ids, gpu=False):
    """(fit_cops, elements, model over the chosen rows): the reference parser; GPU-backed model when gpu."""
    from cops_cuda import gpu_model, host_model, reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    header, names, rows = fit_cops.read_model_file(str(base / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    elements = fit_cops.read_elements(str(base / "elements.bin"))
    coef = fit_cops.read_activity_coefficients(str(reference.bind()), header)
    cls = gpu_model.gpu_model_class() if gpu else host_model._header_model_class()
    M = cls(str(base), header, cols, rows, elements, coef)
    index = {int(p): i for i, p in enumerate(M.ipoint)}
    M.load_kernels(np.array([index[int(p)] for p in point_ids], dtype=int))
    return fit_cops, elements, M


def prepare_inputs(args, out):
    """States (fit, perturbed fit, exp), the basis, the rows; saved to <out>/inputs.npz."""
    from cops_cuda import restart as rs
    fit_cops, el, M = reference_model(args.base, args.points)
    r = rs.read(args.fit / "fit_restart.npz")
    basis = np.ascontiguousarray(np.asarray(r.template_columns, dtype=np.float64))
    scale = float(r.enhancement_scale)
    a_fit = np.asarray(r.a, float)
    b_fit = np.asarray(r.template_values, float)
    T_fit = np.asarray(r.T, float)
    k_fit = float(np.asarray(r.kappa).reshape(-1)[0])
    log("fit state: T %s kappa %.6f templates %s amplitudes %s A_scale %.6e" % (T_fit, k_fit, r.template_names,
                                                                                 b_fit, scale))
    # perturbed fit state (T75 > Tmin: both temperature branches)
    k = np.arange(25)
    a_pert = a_fit * (1.0 + 0.05 * np.sin(1.7 * k + 0.3))
    b_pert = b_fit * np.array([0.8, 1.1, 0.9, 1.2])
    T_pert = np.array([205.0, 150.0, 80.0])
    # exp state: area-weighted least-squares projection of ln(AH2O) on the 25 harmonics (facets with AH2O > 0)
    Y, A, area = el["Y"], el["AH2O"], el["area"]
    G, rhs = np.zeros((25, 25)), np.zeros(25)
    for f0 in range(0, len(A), 1 << 20):
        sl = slice(f0, f0 + (1 << 20))
        pos = A[sl] > 0
        Yc, w = Y[sl][pos], area[sl][pos]
        G += Yc.T @ (Yc * w[:, None])
        rhs += Yc.T @ (w * np.log(A[sl][pos]))
    a_exp = np.linalg.solve(G, rhs)
    lin = Y @ a_exp
    pos = A > 0
    lnres = lin[pos] - np.log(A[pos])
    stats = {"facets_used": int(pos.sum()), "facets_AH2O_zero": int((~pos).sum()),
             "rms_ln_exp_over_AH2O_area_weighted": float(np.sqrt(np.sum(area[pos] * lnres ** 2) / np.sum(area[pos]))),
             "exp_map_min": float(np.exp(lin.min())), "exp_map_max": float(np.exp(lin.max())),
             "AH2O_max": float(A.max()), "AH2O_area_mean": float(np.sum(area * A) / np.sum(area)),
             "exp_map_area_mean": float(np.sum(area * np.exp(lin)) / np.sum(area)), "a_exp": a_exp.tolist()}
    log("exp state: %s" % {k: v for k, v in stats.items() if k != "a_exp"})
    bmax = np.max(np.abs(basis), axis=0)
    b_exp = np.array(EXP_TEMPLATE_EXPONENT, float) / bmax
    full = lin + EXP_TEMPLATE_SCALE * (basis @ b_exp)
    stats.update(template_column_absmax=bmax.tolist(), template_column_nonzero=np.count_nonzero(basis, axis=0).tolist(),
                 template_exponent_max=list(EXP_TEMPLATE_EXPONENT), act_exp_min=float(full.min()),
                 act_exp_max=float(full.max()))
    log("exp state templates: column max %s, amplitudes %s, exponent range [%.3f, %.3f]"
        % (bmax, b_exp, full.min(), full.max()))
    np.savez(out / "inputs.npz", point_ids=np.asarray(args.points, int), basis=basis, scale_fit=scale,
             a_fit=a_fit, b_fit=b_fit, T_fit=T_fit, kappa_fit=k_fit, a_pert=a_pert, b_pert=b_pert, T_pert=T_pert,
             a_exp=a_exp, b_exp=b_exp, scale_exp=EXP_TEMPLATE_SCALE,
             template_names=np.array(r.template_names))
    (out / "inputs.json").write_text(json.dumps(jsonable({
        "base": str(args.base), "store": str(args.store), "fit": str(args.fit), "points": list(args.points),
        "rows": [int(x) for x in M.rows], "utc": [str(x) for x in np.asarray(r.utc)[[list(r.point_ids).index(p) for p in args.points]]],
        "records": [int(c) for c in np.diff(M.indptr)], "fit_state": {"a": a_fit, "b": b_fit, "T": T_fit, "kappa": k_fit,
        "A_scale": scale, "templates": list(r.template_names)}, "perturbed_state": {"a": a_pert, "b": b_pert, "T": T_pert},
        "exp_state": {"b": b_exp, "template_scale": EXP_TEMPLATE_SCALE, "T": T_fit, "kappa": k_fit, **stats}}), indent=1))
    return fit_cops, el, M


def states(inp, which):
    from cops_cuda import backend as bk
    if which == "fit":
        return bk.ParameterState(inp["a_fit"], inp["T_fit"], float(inp["kappa_fit"]), inp["b_fit"])
    if which == "pert":
        return bk.ParameterState(inp["a_pert"], inp["T_pert"], float(inp["kappa_fit"]) * 1.3, inp["b_pert"])
    return bk.ParameterState(inp["a_exp"], inp["T_fit"], float(inp["kappa_fit"]), inp["b_exp"])


# ---------------------------------------------------------------------------------------------------- one phase
def run_phase(name, work, devices):
    """Bind the phase's backends, evaluate, save phase_<name>.npz/json (runs in its own process)."""
    info = {"phase": name, "pid": os.getpid(), "env": {k: v for k, v in os.environ.items()
                                                       if k.startswith("COPS_") or k.startswith("CUDA_V")}}
    if PHASES[name] is not None:
        import positive_activity
        info["installed"] = positive_activity.install()
        info["config"] = positive_activity.config()
    else:
        info["installed"] = False
    from cops_cuda import backend as bk, cuda_sources, evaluator
    from cops_cuda.errors import UnsupportedConfiguration
    from cops_cuda.multi_gpu import MultiGpuBackend
    inp = dict(np.load(work / "inputs.npz"))
    fit_cops, el, M = reference_model(Path(os.environ["CHECK_BASE"]), [int(p) for p in inp["point_ids"]], gpu=True)
    full = bk.ColumnRequest(range(25), range(4), temperature=True, kappa=True)
    night, act, floor, which = MODEL[name]
    exp_state = which == "exp"
    scale = float(inp["scale_exp"] if exp_state else inp["scale_fit"])
    basis = inp["basis"]
    out = {}
    if name in ("ref", "def"):
        backends, state_names = IDENTITY_BACKENDS, ("fit", "pert")
    else:
        backends, state_names = PHYSICS_BACKENDS[name], (which,)
    info["backends"] = {}
    d0, d1 = devices[0], devices[1]
    for be_name in backends:
        t0 = time.perf_counter()
        activity = "device" if be_name.endswith("device") else "host"
        if be_name == "cuda1_stream":
            saved = {k: os.environ.get(k) for k in STREAM_ENV}
            os.environ.update({k: (v if v is not None else str(int(M.indptr[-1]) // 2)) for k, v in STREAM_ENV.items()})
            try:
                be = evaluator.CudaBackend(M, M.run, basis, scale, device_id=d0, activity=activity)
            finally:
                for k, v in saved.items():
                    if v is None:
                        os.environ.pop(k, None)
                    else:
                        os.environ[k] = v
        elif be_name == "cuda1_fp64":
            be = evaluator.CudaBackend(M, M.run, basis, scale, device_id=d0, activity=activity, record_precision="fp64")
        elif be_name.startswith("cuda1"):
            be = evaluator.CudaBackend(M, M.run, basis, scale, device_id=d0, activity=activity)
        else:
            be = MultiGpuBackend(M._spec, M, basis, scale, devices=(d0, d1), activity=activity)
        bind_s = time.perf_counter() - t0
        try:
            for st in state_names:
                s = states(inp, st)
                p = be.evaluate(s)
                j = be.evaluate(s, full)
                out["%s/%s/pNG" % (be_name, st)], out["%s/%s/pRG" % (be_name, st)] = p.pNG, p.pRG
                out["%s/%s/jac_pNG" % (be_name, st)], out["%s/%s/jac_pRG" % (be_name, st)] = j.pNG, j.pRG
                for key, (ng, rg) in j.columns.items():
                    tag = "%s%d" % key
                    out["%s/%s/J/%s/NG" % (be_name, st, tag)] = np.asarray(ng)
                    out["%s/%s/J/%s/RG" % (be_name, st, tag)] = np.asarray(rg)
            # warm timing of the full Jacobian (first state)
            s = states(inp, state_names[0])
            times = []
            for _ in range(5):
                t1 = time.perf_counter()
                be.evaluate(s, full)
                times.append(time.perf_counter() - t1)
            ident = be.identity()
            entry = {"bind_s": bind_s, "jacobian_s_median": float(np.median(times)), "summary": be.summary()}
            if be_name.startswith("cuda1"):
                entry["cuda_source"] = ident.get("cuda_source")
                entry["activity_source"] = ident.get("activity", {}).get("cuda_source")
                entry["model_extension"] = ident.get("model_extension")
                entry["dataset_faces"] = ident["dataset"]["face_compaction"]
                entry["streamed_records"] = ident["dataset"]["streamed_records"]
                entry["record_precision"] = ident["dataset"]["record_precision"]
                if be_name == "cuda1_host":            # rebind_basis keeps the records (and the floor block)
                    be.rebind_basis(basis, scale)
                    st = state_names[0]
                    j = be.evaluate(states(inp, st), full)
                    same = (j.pNG.tobytes() == out["%s/%s/jac_pNG" % (be_name, st)].tobytes() and
                            all(np.asarray(v[0]).tobytes() == out["%s/%s/J/%s%d/NG" % (be_name, st, *key)].tobytes() and
                                np.asarray(v[1]).tobytes() == out["%s/%s/J/%s%d/RG" % (be_name, st, *key)].tobytes()
                                for key, v in j.columns.items()))
                    entry["rebind_basis_bitwise"] = bool(same)
            else:
                entry["workers"] = [{"pid": w.get("pid"), "ready": r, "cuda_source": w.get("cuda_source"),
                                     "activity_source": (w.get("activity") or {}).get("cuda_source"),
                                     "compiled": ((w.get("model_extension") or {}).get("positive_activity") or {}).get("compiled"),
                                     "model_extension_config": {k: ((w.get("model_extension") or {}).get("positive_activity") or {}).get(k)
                                                                for k in ("enabled", "modified", "night_ratio", "qmin_ratio", "activity",
                                                                          "floor_ratio", "role")},
                                     "dataset_floor": ((w.get("model_extension") or {}).get("positive_activity") or {}).get("dataset")}
                                    for w, r in zip(ident["workers"], ident["workers_ready"])]
                entry["model_extension"] = ident.get("model_extension")
            info["backends"][be_name] = entry
            # finite differences (cuda1_host, physics phases)
            if be_name == "cuda1_host" and name not in ("ref", "def"):
                s0 = states(inp, which)
                rms_act = 1.0
                if act == "clip":
                    lin = be.activity_linear(s0)
                    rms_act = float(np.sqrt(np.mean(lin ** 2)))
                Y = M.Y
                params = [("a", k) for k in FD_A] + [("b", k) for k in FD_B] + [("T", 0), ("T", 1), ("T", 2), ("kappa", 0)]
                info["fd"] = {"act_unit": rms_act, "steps": {}}
                for kind, k in params:
                    for delta in FD_STEPS:
                        if kind == "a":
                            h = delta * rms_act / float(np.max(np.abs(Y[:, k])))
                        elif kind == "b":
                            h = delta * rms_act / (abs(scale) * float(np.max(np.abs(basis[:, k]))))
                        elif kind == "T":
                            h = delta * 100.0
                        else:
                            h = delta * 10.0 * s0.kappa
                        vals = []
                        for sign in (1.0, -1.0):
                            a_, b_, T_, kap = s0.a.copy(), s0.b.copy(), s0.T.copy(), s0.kappa
                            if kind == "a":
                                a_[k] += sign * h
                            elif kind == "b":
                                b_[k] += sign * h
                            elif kind == "T":
                                T_[k] += sign * h
                            else:
                                kap += sign * h
                            r = be.evaluate(bk.ParameterState(a_, T_, kap, b_))
                            vals.append((r.pNG, r.pRG))
                        tag = "%s%d" % (kind, k)
                        out["fd/%s/%g/NG" % (tag, delta)] = (vals[0][0] - vals[1][0]) / (2 * h)
                        out["fd/%s/%g/RG" % (tag, delta)] = (vals[0][1] - vals[1][1]) / (2 * h)
                        info["fd"]["steps"]["%s/%g" % (tag, delta)] = h
        finally:
            be.close()
    # refusals of the modified model (FP64 records; restart.cu)
    if name not in ("ref", "def"):
        refusals = {}
        try:
            evaluator.CudaBackend(M, M.run, basis, scale, device_id=d0, record_precision="fp64")
            refusals["fp64"] = "NOT REFUSED"
        except UnsupportedConfiguration as exc:
            refusals["fp64"] = "refused: %s" % exc
        try:
            cuda_sources.function("restart", "facet_point_shares_m", d0)
            refusals["restart"] = "NOT REFUSED"
        except UnsupportedConfiguration as exc:
            refusals["restart"] = "refused: %s" % exc
        info["refusals"] = refusals
    elif name == "def":
        cuda_sources.module("restart", d0)               # defaults: the original restart.cu stays usable
        info["restart_module_default"] = "compiled"
    info["compiled_in_driver"] = cuda_sources.loaded_identities()
    np.savez(work / ("phase_%s.npz" % name), **out)
    (work / ("phase_%s.json" % name)).write_text(json.dumps(jsonable(info), indent=1))
    log("phase %s done" % name)


# --------------------------------------------------------------------------------------------- NumPy reference
def numpy_pressures(fit_cops, el, M, base, store, point_ids, act_lin, night_ratio, activity, floor_ratio, T, kappa,
                    source):
    """pNG, pRG of the modified model at the points (formulas of fit_cops Model.evaluate + the extension)."""
    from cops_cuda.record_store import split_fields, store_file_name
    AU, KT, COS75 = fit_cops.AU, fit_cops.KT, fit_cops.COS75
    p_exp = 4.2143229 if M.helio == "shap5" else 4.8
    rH, rC = (0.2, 0.1) if night_ratio is None else (night_ratio, night_ratio)
    Tmax, T75, Tmin = T
    A = np.exp(act_lin) if activity == "exp" else np.maximum(act_lin, 0.0)
    cNG, cRG, bm = M.cNG, M.cRG, M.bm
    pNG, pRG = np.zeros(len(point_ids)), np.zeros(len(point_ids))
    stats = {"records": 0, "at_floor": 0, "ng_h2o_sum": 0.0, "ng_h2o_floor_sum": 0.0, "per_point_floor_fraction": []}
    for j, ip in enumerate(point_ids):
        _, h, rec = fit_cops.read_kernel(str(base / "kernels" / ("point=%04d.bin" % ip)))
        if source == "fp64":
            face, kn, kphi, kv, cos = rec["face"], rec["Kn"], rec["KPhi"], rec["KV"], rec["cos"]
            night = (cos <= 0.0) | (rec["shadow"] != 0)
        else:
            raw = (store / store_file_name(ip)).read_bytes()
            if len(raw) != 20 * len(rec):
                raise SystemExit("store file of point %d has %d bytes for %d records" % (ip, len(raw), len(rec)))
            f = split_fields(raw, len(rec))
            face = f["face"].astype(np.int64)
            kn, kphi, kv, cos = (f[n].astype(np.float64) for n in ("kn", "kphi", "kv", "cos"))
            night = cos <= 0.0
        rhA = h["rh"] / AU
        qmaxH, qmaxC = 600 / rhA ** p_exp, 600 / rhA ** 2
        qminH, qminC = rH * 600 / rhA ** p_exp, rC * 600 / rhA ** 2
        cz = np.minimum(cos, 1.0)
        qH = np.where(night, qminH, qminH + (qmaxH - qminH) * cz)
        qC = np.where(night, qminC, qminC + (qmaxC - qminC) * cz) * el["ACO2"][face]
        c = np.where(night, -1.0, np.clip(cos, -1.0, 1.0))
        lit = c >= COS75
        xlit = np.where(lit, 1.0 / np.where(lit, c, 1.0) - 1.0, 0.0)
        F = qH * A[face]
        if floor_ratio:
            Fmin = qmaxH * (floor_ratio * el["AH2O"][face])
            on = F < Fmin
            F = np.where(on, Fmin, F)
        else:
            on = np.zeros(len(F), bool)
        Tl = Tmax + KT * (Tmax - T75) * xlit
        use_lit = lit & (Tl > Tmin)
        u = 1.0 / np.sqrt(np.where(use_lit, Tl, Tmin))
        VB = h["VdotB"]
        nH = np.sum(F * kn * u)
        pNG[j] = cNG["H2O"] * bm["H2O"] * nH + kappa * cNG["CO2"] * bm["CO2"] * np.sum(qC * kn * u)
        pRG[j] = (cRG["H2O"][j] * (np.sum(F * kphi) + bm["H2O"] * np.sum(F * VB * kv * u))
                  + kappa * cRG["CO2"][j] * (np.sum(qC * kphi) + bm["CO2"] * np.sum(qC * VB * kv * u)))
        stats["records"] += len(F)
        stats["at_floor"] += int(on.sum())
        stats["ng_h2o_sum"] += float(nH)
        nf = float(np.sum((F * kn * u)[on]))
        stats["ng_h2o_floor_sum"] += nf
        stats["per_point_floor_fraction"].append([int(ip), float(on.mean()) if len(on) else 0.0, nf / nH if nH else 0.0])
    return pNG, pRG, stats


def relmax(x, ref):
    if np.size(ref) == 0:
        return float("nan")
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.max(np.abs(x / ref - 1.0)))


def late_mask(utc):
    """Points in the late window 10:33:10-10:39:26 UTC (time of day of the table's UTC strings)."""
    def seconds(s):
        h, m, sec = s.split("T")[1].split(":")
        return int(h) * 3600 + int(m) * 60 + float(sec)
    t = np.array([seconds(u) for u in utc])
    return (t >= 10 * 3600 + 33 * 60 + 10) & (t <= 10 * 3600 + 39 * 60 + 26)


# ---------------------------------------------------------------------------------------------------- analysis
def analyse(work, args, fit_cops, el, M):
    inp = dict(np.load(work / "inputs.npz"))
    res = {p: dict(np.load(work / ("phase_%s.npz" % p))) for p in PHASES}
    meta = {p: json.loads((work / ("phase_%s.json" % p)).read_text()) for p in PHASES}
    report = {"inputs": json.loads((work / "inputs.json").read_text()), "phases": meta}
    late = late_mask(report["inputs"]["utc"])
    lines = []
    # 1. identity
    ident = {"compared_arrays": 0, "bitwise_equal": 0, "mismatches": []}
    for key, ref in res["ref"].items():
        got = res["def"][key]
        ident["compared_arrays"] += 1
        if got.dtype == ref.dtype and got.shape == ref.shape and got.tobytes() == ref.tobytes():
            ident["bitwise_equal"] += 1
        else:
            ident["mismatches"].append({"key": key, "max_abs_rel": relmax(got, ref)})
    ident["keys_only_in_def"] = sorted(set(res["def"]) - set(res["ref"]))
    ident["bitwise_identical"] = (not ident["mismatches"] and not ident["keys_only_in_def"]
                                  and ident["compared_arrays"] > 0)
    # one GPU vs two GPUs inside each phase (point independence: expected bitwise)
    for p in ("ref", "def"):
        for other in ("multi2_host", "cuda1_stream"):
            eq = all(res[p][k].tobytes() == res[p][k.replace("cuda1_host", other, 1)].tobytes()
                     for k in res[p] if k.startswith("cuda1_host/"))
            ident["%s_cuda1_host_vs_%s_bitwise" % (p, other)] = bool(eq)
    workers = {p: {be: [{"pid": (w.get("ready") or {}).get("pid"), "pressure": (w["cuda_source"] or {}).get("path"),
                         "activity": (w["activity_source"] or {}).get("path"),
                         "compiled": [(c["name"], c["path"], c["sha256"], c["device_id"]) for c in (w["compiled"] or [])]}
                        for w in meta[p]["backends"][be].get("workers", [])]
                   for be in meta[p]["backends"] if be.startswith("multi")} for p in ("ref", "def", "mod")}
    ident["worker_sources"] = workers
    report["identity"] = ident
    lines.append("1. IDENTITY (COPS_POS_ENABLE=1, defaults vs unpatched): %d arrays (backends %s x 2 states x "
                 "(2 pressures x 2 + 33 columns x 2)), bitwise equal %d, mismatches %d -> %s"
                 % (ident["compared_arrays"], "/".join(IDENTITY_BACKENDS), ident["bitwise_equal"], len(ident["mismatches"]),
                    "BITWISE IDENTICAL" if ident["bitwise_identical"] else "DIFFERENT"))
    lines.append("   one GPU vs two GPUs bitwise: ref %s, def %s; resident vs streamed records bitwise: ref %s, def %s"
                 % (ident["ref_cuda1_host_vs_multi2_host_bitwise"], ident["def_cuda1_host_vs_multi2_host_bitwise"],
                    ident["ref_cuda1_host_vs_cuda1_stream_bitwise"], ident["def_cuda1_host_vs_cuda1_stream_bitwise"]))
    lines.append("   streamed records: %s; rebind_basis bitwise: %s" % (
        {p: meta[p]["backends"]["cuda1_stream"]["streamed_records"] for p in ("ref", "def", "mod")},
        {p: meta[p]["backends"]["cuda1_host"].get("rebind_basis_bitwise") for p in PHASES}))
    for p in ("ref", "def", "mod"):
        for be, ws in workers[p].items():
            for w in ws:
                lines.append("   %s %s worker pid %s: pressure.cu %s; compiled %s" % (
                    p, be, w["pid"], w["pressure"], [(c[0], c[1], c[2][:12], c[3]) for c in w["compiled"]] or "-"))
    # 2. physics
    phys = {}
    act_cache = {}
    for p in ("ref", "mod", "nofloor", "clipfloor"):
        night, activity, floor, which = MODEL[p]
        s = states(inp, which)
        scale = float(inp["scale_exp"] if which == "exp" else inp["scale_fit"])
        if which not in act_cache:
            act_cache[which] = M.Y @ s.a + scale * (inp["basis"] @ s.b)
        entry = {}
        for source in ("fp64", "fp32store"):
            t0 = time.perf_counter()
            nNG, nRG, st = numpy_pressures(fit_cops, el, M, args.base, args.store, [int(x) for x in inp["point_ids"]],
                                           act_cache[which], night, activity, floor, s.T, s.kappa,
                                           "fp64" if source == "fp64" else "fp32")
            entry[source] = {"stats": st, "seconds": time.perf_counter() - t0, "backends": {}}
            for be in (PHYSICS_BACKENDS.get(p) or ("cuda1_host",)):
                gNG, gRG = res[p]["%s/%s/pNG" % (be, which)], res[p]["%s/%s/pRG" % (be, which)]
                jNG, jRG = res[p]["%s/%s/jac_pNG" % (be, which)], res[p]["%s/%s/jac_pRG" % (be, which)]
                entry[source]["backends"][be] = {
                    "max_rel_NG": relmax(gNG, nNG), "max_rel_RG": relmax(gRG, nRG),
                    "max_rel_NG_late": relmax(gNG[late], nNG[late]), "max_rel_RG_late": relmax(gRG[late], nRG[late]),
                    "jacobian_call_pressures_bitwise": bool(jNG.tobytes() == gNG.tobytes() and jRG.tobytes() == gRG.tobytes())}
            if source == "fp64":
                entry["numpy_pressures"] = {"NG": nNG.tolist(), "RG": nRG.tolist()}
        if p != "ref":
            host = {k[len("cuda1_host/"):]: v for k, v in res[p].items() if k.startswith("cuda1_host/")}
            for be in PHYSICS_BACKENDS[p][1:]:
                other = {k[len(be) + 1:]: v for k, v in res[p].items() if k.startswith(be + "/")}
                entry["%s_vs_cuda1_host" % be] = {
                    "bitwise": all(other[k].tobytes() == host[k].tobytes() for k in host),
                    "max_rel": max(relmax(other[k], host[k]) for k in host if np.all(host[k] != 0)) if host else None,
                    "max_abs_over_colmax": max(float(np.max(np.abs(other[k] - host[k])) / max(np.max(np.abs(host[k])), 1e-300))
                                               for k in host)}
        phys[p] = entry
    report["physics"] = phys
    lines.append("2. PHYSICS: GPU pressures vs independent NumPy (max |GPU/NumPy - 1| over the %d points)"
                 % len(inp["point_ids"]))
    for p, e in phys.items():
        night, activity, floor, which = MODEL[p]
        st = e["fp64"]["stats"]
        lines.append("   %-9s (night %s, %s, floor %g, state %s; records %d, at floor %d = %.3f%%, NG H2O sum at floor "
                     "%.3f%%)" % (p, night, activity, floor, which, st["records"], st["at_floor"],
                                  100.0 * st["at_floor"] / max(st["records"], 1),
                                  100.0 * st["ng_h2o_floor_sum"] / max(st["ng_h2o_sum"], 1e-300)))
        for source in ("fp64", "fp32store"):
            for be, v in e[source]["backends"].items():
                lines.append("      vs NumPy(%-9s records) %-13s NG %.2e RG %.2e (late window NG %.2e RG %.2e)"
                             % (source, be, v["max_rel_NG"], v["max_rel_RG"], v["max_rel_NG_late"], v["max_rel_RG_late"]))
        for k, v in e.items():
            if k.endswith("_vs_cuda1_host"):
                lines.append("      %s: bitwise %s, max rel %.2e, max abs/colmax %.2e" % (k, v["bitwise"], v["max_rel"] or 0.0,
                                                                                         v["max_abs_over_colmax"]))
    # 3. Jacobian vs finite differences
    jac = {}
    lines.append("3. JACOBIAN (cuda1_host analytic vs central finite differences of GPU pressures): max |Jan - Jfd| / "
                 "max|Jan| (column-scaled) and max pointwise rel (|Jan| > 1e-3 max|Jan|), per step")
    for p in ("mod", "nofloor", "clipfloor"):
        night, activity, floor, which = MODEL[p]
        jac[p] = {}
        lines.append("   %s (%s, floor %g):" % (p, activity, floor))
        for kind, ks in (("a", FD_A), ("b", FD_B), ("T", (0, 1, 2)), ("kappa", (0,))):
            for k in ks:
                tag = "%s%d" % (kind, k)
                row = {}
                for delta in FD_STEPS:
                    for g in ("NG", "RG"):
                        an = res[p]["cuda1_host/%s/J/%s/%s" % (which, tag, g)]
                        fd = res[p]["fd/%s/%g/%s" % (tag, delta, g)]
                        m = float(np.max(np.abs(an)))
                        col = float(np.max(np.abs(an - fd))) / m if m > 0 else float(np.max(np.abs(fd)))
                        sig = np.abs(an) > 1e-3 * m
                        pt = float(np.max(np.abs(an[sig] - fd[sig]) / np.abs(an[sig]))) if sig.any() else 0.0
                        row["%g/%s" % (delta, g)] = {"colscaled": col, "pointwise": pt, "colmax": m}
                jac[p][tag] = row
                lines.append("      %-6s " % tag + "  ".join("h%g %s: col %.1e pt %.1e" % (
                    delta, g, row["%g/%s" % (delta, g)]["colscaled"], row["%g/%s" % (delta, g)]["pointwise"])
                    for delta in FD_STEPS for g in ("NG", "RG")))
    report["jacobian"] = jac
    # refusals and timing
    for p in ("mod", "nofloor", "clipfloor"):
        lines.append("   refusals (%s): %s" % (p, {k: v[:90] for k, v in meta[p]["refusals"].items()}))
    lines.append("4. TIMING (median warm full Jacobian, 33 columns, %d points): %s" % (len(inp["point_ids"]), {
        p: {be: round(v["jacobian_s_median"], 4) for be, v in meta[p]["backends"].items()} for p in PHASES}))
    # pass/fail: identity bitwise; GPU vs NumPy <= 1e-6 (FP64 kernel files) and <= 1e-12 (FP32 store records)
    nan_inf = lambda x: float("inf") if x != x else x     # noqa: E731  (a NaN fails)
    worst = {src: max(nan_inf(v[k]) for e in phys.values() for v in e[src]["backends"].values()
                      for k in ("max_rel_NG", "max_rel_RG")) for src in ("fp64", "fp32store")}
    passed = {"identity_bitwise": bool(ident["bitwise_identical"]),
              "physics_fp64_le_1e-6": bool(worst["fp64"] <= 1e-6), "physics_fp32store_le_1e-12": bool(worst["fp32store"] <= 1e-12)}
    report["passed"] = dict(passed, worst_physics=worst)
    lines.append("PASS/FAIL: %s; worst GPU/NumPy - 1: FP64 records %.2e, FP32 store records %.2e"
                 % (passed, worst["fp64"], worst["fp32store"]))
    report["summary_lines"] = lines
    (work / "report.json").write_text(json.dumps(jsonable(report), indent=1))
    (work / "report.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    return all(passed.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    run = sub.add_parser("run")
    run.add_argument("--out", type=Path, required=True)
    run.add_argument("--devices", default="0,1")
    run.add_argument("--points", default=",".join(str(p) for p in DEFAULT_POINTS))
    run.add_argument("--base", type=Path, default=DEFAULTS["base"])
    run.add_argument("--store", type=Path, default=DEFAULTS["store"])
    run.add_argument("--fit", type=Path, default=DEFAULTS["fit"])
    run.add_argument("--phases", default=",".join(PHASES))
    ph = sub.add_parser("phase")
    ph.add_argument("name", choices=list(PHASES))
    ph.add_argument("--work", type=Path, required=True)
    ph.add_argument("--devices", default="0,1")
    a = ap.parse_args()
    if a.cmd == "phase":
        run_phase(a.name, a.work, [int(d) for d in a.devices.split(",")])
        return 0
    if a.cmd != "run":
        ap.error("use 'run' or 'phase'")
    a.points = [int(p) for p in a.points.split(",")]
    a.out.mkdir(parents=True, exist_ok=False)
    log("inputs")
    fit_cops, el, M = prepare_inputs(a, a.out)
    base_env = {k: v for k, v in os.environ.items() if not k.startswith("COPS_POS_")}
    base_env.update(COPS_CUDA_RECORD_PRECISION="fp32", COPS_CUDA_RECORD_STORE=str(a.store), CHECK_BASE=str(a.base))
    for name in a.phases.split(","):
        env = dict(base_env, **(PHASES[name] or {}))
        cmd = [sys.executable, "-W", "ignore", str(Path(__file__).resolve()), "phase", name, "--work", str(a.out),
               "--devices", a.devices]
        log("phase %s: %s with %s" % (name, " ".join(cmd), PHASES[name]))
        with open(a.out / ("phase_%s.log" % name), "w") as logf:
            rc = subprocess.call(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
        log("phase %s: exit %d" % (name, rc))
        if rc != 0:
            print(open(a.out / ("phase_%s.log" % name)).read()[-6000:], flush=True)
            return 3
    ok = analyse(a.out, a, fit_cops, el, M)
    return 0 if ok else 4


if __name__ == "__main__":
    sys.exit(main())
