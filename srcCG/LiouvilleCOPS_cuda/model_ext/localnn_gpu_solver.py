#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_gpu_solver.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_gpu_solver.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: ea91470bc363ae500de5134085f0387ab018e1baca7bfefd504e69b0cae09928
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: GPU solver path of option C (lead / user 2026-11-15 'option 3'): localnn_fit.py fit --device gpu puts the
#   design products (surface_used, data_jacobian) and scipy TRF's svd(J_augmented) (modes qr = QR + SVD(R), gesvd)
#   on one GPU with CuPy FP64; the rest of TRF unchanged on the CPU; CLI trajectory (CPU / GPU cost trajectory) and
#   verify (fixed state, one TRF step, 50-evaluation trajectory, speed).
# 2026-11-15: GPU solver path of option C (lead / user 2026-11-15 'option 3'): localnn_fit.py fit --device gpu puts
#   the design products (surface_used, data_jacobian) and scipy TRF's svd(J_augmented) (modes qr = QR + SVD(R),
#   gesvd) on one GPU with CuPy FP64; the rest of TRF unchanged on the CPU; CLI trajectory (CPU / GPU cost
#   trajectory) and verify (fixed state, one TRF step, 50-evaluation trajectory, speed).
# RUN 2026-11-15: CPU reference 132047 (11.11 s per evaluation, 16 threads); GPU verification 132810 (frozen
#   localnn_s2gpu_20261115T192024Z; 132048 cancelled before it started).
# ==== end of header ====
"""GPU solver path of option C (lead / user 2026-11-15 ~12:40 CST, "option 3"): localnn_fit.py fit --device gpu.

One GPU, CuPy, FP64. Two parts, both installed by install(P, mode):
  1. the design-matrix products of the Problem on the GPU: surface_used (GN[b] @ c + tN[b], the m_b sum, kappa) and
     data_jacobian (sum_b m_b GN[b] A_s, the m columns uN[b], kappa, grain, the held term GN_h K) - the same formulas
     as localnn_problem.Problem, evaluated with CuPy; the small parts (logs, crossing rows, smoothness rows, P_m rows,
     the held term of the residuals) stay on the CPU;
  2. the TRF linear algebra: scipy.optimize._lsq.trf calls svd(J_augmented, full_matrices=False) once per iteration
     (the exact trust-region solver; J_augmented = [J diag(d); diag(diag_h^0.5)], ~8522 x 2754, 9.6 s of the ~11 s per
     evaluation on 16 CPU threads). The module-level name trf.svd is replaced by a CuPy SVD:
       mode "qr"    : Q, R = qr(A) (cuSOLVER geqrf / orgqr), U_R, s, Vt = svd(R) (gesvd), U = Q U_R
       mode "gesvd" : cupy.linalg.svd(A) directly (cuSOLVER gesvd).
     Everything else in TRF (the trust-region solve on (s, V, U^T f), the scaling, the step control) is scipy's code
     on the CPU, unchanged. The default device (cpu) does not import this module: bitwise the previous path.
uninstall() restores trf.svd. stats() returns the device record (GPU name, CuPy / CUDA versions, SVD calls / seconds).

CLI (verification, lead 2026-11-15: "on a fixed state: residuals, Jacobian and one TRF step equal to rounding; a
50-evaluation fit with the same cost trajectory to ~1e-10 relative; the speed per evaluation"):
  localnn_gpu_solver.py trajectory --device cpu|gpu [--mode qr|gesvd] --case CASE.json --nfev 50 --out DIR
      runs least_squares exactly as localnn_fit.py fit does (trf, exact, x_scale jac, TOL, bounds) from the case's
      start for --nfev evaluations; records every evaluation's objective and the final x (DIR/trajectory.npz, .json)
  localnn_gpu_solver.py verify --case CASE.json --cpu-ref DIR --out DIR2 [--modes qr,gesvd]
      on a GPU node: fixed-state residuals and Jacobian CPU vs GPU, one TRF step (SVD + solve_lsq_trust_region with the
      same Delta) CPU vs GPU per mode, the GPU trajectory per mode vs the CPU reference, seconds per evaluation
  CASE.json: {"state": <localnn_state.npz>, "lambda_C", "mult_penalty_f", "held": <finetune_state.npz> | null,
              "late_kernels", "grain_from": <finetune_state.npz> | null, "grid", "design", "base", "d2", "point_table"}
"""

import argparse
import json
import os
import socket
import sys
import time
import types
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402

_STATE = {"installed": False, "orig_svd": None, "mode": None, "svd_calls": 0, "svd_seconds": 0.0, "products": False}


def _cp():
    import cupy as cp
    return cp


def gpu_svd_factory(mode):
    cp = _cp()

    def svd(a, full_matrices=False, **kw):
        if full_matrices:
            raise ValueError("gpu svd: full_matrices=False only")
        t = time.time()
        A = cp.asarray(np.ascontiguousarray(a, dtype=np.float64))
        if mode == "qr" and A.shape[0] > A.shape[1]:
            Q, R = cp.linalg.qr(A, mode="reduced")
            Ur, s, Vt = cp.linalg.svd(R, full_matrices=False)
            U = Q @ Ur
        elif mode in ("qr", "gesvd"):
            U, s, Vt = cp.linalg.svd(A, full_matrices=False)
        else:
            raise ValueError("unknown svd mode %r" % mode)
        out = cp.asnumpy(U), cp.asnumpy(s), cp.asnumpy(Vt)
        _STATE["svd_calls"] += 1
        _STATE["svd_seconds"] += time.time() - t
        return out
    return svd


# ------------------------------------------------------------------------------------------------ products on the GPU
def _gpu_surface_used(self, c, m, kappa):
    cp = _cp()
    g = self._gpu
    cd = cp.asarray(c)
    uN = cp.empty((lp.NBINS, g["GN"].shape[1]))
    uR = cp.empty((lp.NBINS, g["GR"].shape[1]))
    sN = cp.zeros(g["GN"].shape[1])
    sR = cp.zeros(g["GR"].shape[1])
    for b in range(lp.NBINS):
        uN[b] = g["GN"][b] @ cd + g["tN"][b]
        uR[b] = g["GR"][b] @ cd + g["tR"][b]
        sN += float(m[b]) * uN[b]
        sR += float(m[b]) * uR[b]
    sN = sN + kappa * g["CN"]
    sR = sR + kappa * g["CR"]
    return cp.asnumpy(sN), cp.asnumpy(sR), cp.asnumpy(uN), cp.asnumpy(uR)


def _gpu_data_jacobian(self, x):
    cp = _cp()
    g = self._gpu
    c, m, kappa, A, B = self.unpack(x)
    pNG, pRG, (uN, uR, gr, dA, dB) = self.pressures_used(x)
    D = self.D
    fN = cp.asarray(D.wNG / np.maximum(pNG, lp.FLOOR))
    fR = cp.asarray(D.wRG / np.maximum(pRG, lp.FLOOR))
    md = cp.asarray(m)
    HN = cp.tensordot(md, g["GN"], axes=(0, 0)) * self.A_s
    HR = cp.tensordot(md, g["GR"], axes=(0, 0)) * self.A_s
    n_x = len(x)
    JN = cp.empty((len(pNG), n_x))
    JR = cp.empty((len(pRG), n_x))
    JN[:, self.ic] = HN * fN[:, None]
    JR[:, self.ic] = HR * fR[:, None]
    JN[:, self.im] = cp.asarray(uN[1:].T) * fN[:, None]
    JR[:, self.im] = cp.asarray(uR[1:].T) * fR[:, None]
    JN[:, self.ik] = g["CN"] * fN
    JR[:, self.ik] = g["CR"] * fR
    JN[:, self.iA] = cp.asarray(dA) * fN
    JN[:, self.iB] = cp.asarray(dB) * fN
    JR[:, self.iA] = 0.0
    JR[:, self.iB] = 0.0
    if self.held is not None:
        H = self.held
        qH, _ = self._held_terms(c, m)
        act = (qH > 0.0).astype(float)
        K = (sp.diags(m[H["bins"]] * H["gm1"] * act * self.A_s) @ H["Phi"]).toarray()   # d aH / d x_c (as the CPU)
        Kd = cp.asarray(K)
        JN[:, self.ic] += (g["GNh"] @ Kd) * fN[:, None]
        JR[:, self.ic] += (g["GRh"] @ Kd) * fR[:, None]
        w = qH * H["gm1"]
        for b in range(1, lp.NBINS):
            sel = H["bins"] == b
            if sel.any():
                ws = cp.asarray(w[sel])
                JN[:, self.im.start + b - 1] += (g["GNh"][:, cp.asarray(np.flatnonzero(sel))] @ ws) * fN
                JR[:, self.im.start + b - 1] += (g["GRh"][:, cp.asarray(np.flatnonzero(sel))] @ ws) * fR
    return cp.asnumpy(cp.vstack([JN, JR]))


def install(P, mode="qr", products=True):
    """Put P's design products (products=True) and TRF's SVD on the current GPU. Returns the device record."""
    cp = _cp()
    import scipy.optimize._lsq.trf as trf
    if _STATE["installed"]:
        raise RuntimeError("gpu solver already installed")
    if products:
        P._gpu = {"GN": cp.asarray(P.GN), "GR": cp.asarray(P.GR), "tN": cp.asarray(P.tN), "tR": cp.asarray(P.tR),
                  "CN": cp.asarray(P.CN), "CR": cp.asarray(P.CR)}
        if P.held is not None:
            P._gpu["GNh"] = cp.asarray(np.asarray(P.held["GN"]))
            P._gpu["GRh"] = cp.asarray(np.asarray(P.held["GR"]))
        P.surface_used = types.MethodType(_gpu_surface_used, P)
        P.data_jacobian = types.MethodType(_gpu_data_jacobian, P)
    _STATE.update(installed=True, orig_svd=trf.svd, mode=mode, svd_calls=0, svd_seconds=0.0, products=bool(products))
    trf.svd = gpu_svd_factory(mode)
    return stats()


def uninstall():
    import scipy.optimize._lsq.trf as trf
    if _STATE["installed"]:
        trf.svd = _STATE["orig_svd"]
        _STATE["installed"] = False


def stats():
    cp = _cp()
    dev = cp.cuda.Device()
    props = cp.cuda.runtime.getDeviceProperties(dev.id)
    name = props["name"].decode() if isinstance(props["name"], bytes) else str(props["name"])
    return {"device": "gpu", "gpu_name": name, "gpu_index_in_job": int(dev.id),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"), "pci_bus_id": dev.pci_bus_id,
            "cupy": cp.__version__, "cuda_runtime": int(cp.cuda.runtime.runtimeGetVersion()),
            "cuda_driver": int(cp.cuda.runtime.driverGetVersion()), "svd_mode": _STATE["mode"],
            "products_on_gpu": _STATE["products"], "svd_calls": _STATE["svd_calls"],
            "svd_seconds": _STATE["svd_seconds"], "precision": "FP64",
            "definition": "localnn_gpu_solver.py: design products (surface_used, data_jacobian) and scipy TRF's svd of "
                          "J_augmented on one GPU (CuPy FP64); the rest of TRF unchanged on the CPU"}


# ------------------------------------------------------------------------------------------------ verification CLI
def load_case(path):
    case = json.load(open(path))
    grid = lp.load_grid(case["grid"])
    design = lp.load_design(case["design"])
    D = lp.Data(case["base"], case["d2"], case["point_table"], design)
    with np.load(case["state"], allow_pickle=False) as z:
        x = np.asarray(z["x"], float)
        meta = json.loads(str(z["meta_json"]))
    held = None
    if case.get("held"):
        import localnn_late as llate
        nf = len(grid["facet_node"])
        LK = llate.load(case["late_kernels"], nf)
        with np.load(design["meta"]["paths"]["facets_from"], allow_pickle=False) as z:
            TS, bins = np.asarray(z["TS"], np.float64), np.asarray(z["bin"])
        with np.load(case["held"], allow_pickle=False) as z:
            h = np.asarray(z["held_next"], float)
        held = {"h": h, "LK": LK, "bins": bins, "TS": TS}
    if case.get("grain_from"):
        P0 = lp.Problem(D, design, grid, meta["templates"] == "on", 1.0, design["meta"]["A_scale"])
        with np.load(case["grain_from"], allow_pickle=False) as z:
            x[P0.iA], x[P0.iB] = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])

    def make():
        kw = {}
        if float(case["mult_penalty_f"]) > 0:
            kw["mult_penalty_f"] = float(case["mult_penalty_f"])
        if held is not None:
            kw["held"] = held
        return lp.Problem(D, design, grid, meta["templates"] == "on", float(meta.get("slope_mult_min", 1.0)),
                          design["meta"]["A_scale"], float(case["lambda_C"]), **kw)
    return case, make, x


def run_trajectory(make, x0, nfev, device, mode):
    from scipy.optimize import least_squares
    P = make()
    rec = {"device": "cpu"}
    if device == "gpu":
        rec = install(P, mode)
    traj = []
    f0 = P.fun

    def fun(x):
        r = f0(x)
        traj.append(0.5 * float(r @ r))
        return r
    x0 = np.minimum(np.maximum(x0, P.lo), P.hi)
    t = time.time()
    try:
        res = least_squares(fun, x0, jac=P.jac, bounds=(P.lo, P.hi), x_scale="jac", method="trf", tr_solver="exact",
                            ftol=lp.TOL, xtol=lp.TOL, gtol=lp.TOL, max_nfev=nfev, verbose=0)
    finally:
        sec = time.time() - t
        if device == "gpu":
            rec = stats()
            uninstall()
    return {"objective": np.array(traj), "x": res.x, "seconds": sec, "nfev": int(res.nfev), "njev": int(res.njev),
            "message": res.message, "device_record": rec, "seconds_per_evaluation": sec / max(len(traj), 1)}


def cmd_trajectory(args):
    case, make, x0 = load_case(args.case)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    r = run_trajectory(make, x0, args.nfev, args.device, args.mode)
    np.savez(out / "trajectory.npz", objective=r["objective"], x=r["x"])
    js = {k: v for k, v in r.items() if k not in ("objective", "x")}
    js.update(host=socket.gethostname(), slurm_job_id=os.environ.get("SLURM_JOB_ID"), argv=sys.argv, case=case,
              threads=os.environ.get("OPENBLAS_NUM_THREADS"), created_utc=lp.utc_now(),
              objective_first_last=[float(r["objective"][0]), float(r["objective"][-1])])
    json.dump(js, open(out / "trajectory.json", "x"), indent=1, default=float)
    print("trajectory %s: %d evaluations, %.2f s per evaluation, objective %.12e -> %.12e (%s)" % (
        args.device, len(r["objective"]), r["seconds_per_evaluation"], r["objective"][0], r["objective"][-1],
        r["message"]))
    return 0


def cmd_verify(args):
    from scipy.linalg import svd as cpu_svd
    from scipy.optimize._lsq.common import solve_lsq_trust_region
    cp = _cp()
    case, make, x = load_case(args.case)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    res = {"case": case, "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "created_utc": lp.utc_now(), "argv": sys.argv}
    # 1 fixed state: residuals and Jacobian
    Pc, Pg = make(), make()
    from scipy.optimize._lsq.common import make_strictly_feasible
    x = make_strictly_feasible(np.minimum(np.maximum(x, Pc.lo), Pc.hi), Pc.lo, Pc.hi)   # as least_squares does
    install(Pg, "qr")
    t = time.time(); rc = Pc.fun(x); jc = Pc.jac(x); tc = time.time() - t
    t = time.time(); rg = Pg.fun(x); jg = Pg.jac(x); tg = time.time() - t
    uninstall()
    colmax = np.maximum(np.max(np.abs(jc), axis=0), 1e-300)
    res["1_residuals"] = {"max_abs_diff": float(np.max(np.abs(rg - rc))), "max_abs": float(np.max(np.abs(rc))),
                          "max_rel_to_max": float(np.max(np.abs(rg - rc)) / np.max(np.abs(rc))),
                          "objective_rel": float(abs(rg @ rg / (rc @ rc) - 1)), "tolerance_rel_to_max": 1e-13}
    res["1_jacobian"] = {"max_abs_diff": float(np.max(np.abs(jg - jc))), "max_abs": float(np.max(np.abs(jc))),
                         "max_columnwise_rel": float(np.max(np.max(np.abs(jg - jc), axis=0) / colmax)),
                         "tolerance_columnwise_rel": 1e-12, "fun_jac_seconds_cpu": tc, "fun_jac_seconds_gpu": tg}
    res["1_passed"] = bool(res["1_residuals"]["max_rel_to_max"] <= 1e-13 and
                           res["1_jacobian"]["max_columnwise_rel"] <= 1e-12)
    # 2 one TRF step at the fixed state, from scipy's own formulas (trf_bounds, first iteration)
    from scipy.optimize._lsq.common import CL_scaling_vector, compute_jac_scale
    lb, ub = Pc.lo, Pc.hi
    f = rc
    J = jc
    g = J.T.dot(f)
    scale, scale_inv = compute_jac_scale(J)
    v, dv = CL_scaling_vector(x, g, lb, ub)
    v[dv != 0] *= scale_inv[dv != 0]
    Delta = float(np.linalg.norm(x * scale_inv / v ** 0.5)) or 1.0          # trf_bounds' initial radius
    d = v ** 0.5 * scale
    diag_h = g * dv * scale
    m_, n_ = J.shape
    Ja = np.empty((m_ + n_, n_))
    Ja[:m_] = J * d
    Ja[m_:] = np.diag(diag_h ** 0.5)
    fa = np.zeros(m_ + n_)
    fa[:m_] = f
    t = time.time(); U, s, V = cpu_svd(Ja, full_matrices=False); tsc = time.time() - t
    p_c, a_c, it_c = solve_lsq_trust_region(n_, m_, U.T.dot(fa), s, V.T, Delta, initial_alpha=0.0)
    res["2_trf_step"] = {"J_augmented": [int(Ja.shape[0]), int(Ja.shape[1])], "Delta": Delta, "svd_seconds_cpu": tsc,
                         "cond": float(s[0] / s[-1]), "modes": {}}
    ok2 = True
    for mode in args.modes.split(","):
        sv = gpu_svd_factory(mode)
        sv(Ja[:64, :32])                                              # warm-up (handles, kernels)
        t = time.time(); Ug, sg, Vg = sv(Ja, full_matrices=False); tsg = time.time() - t
        p_g, a_g, it_g = solve_lsq_trust_region(n_, m_, Ug.T.dot(fa), sg, Vg.T, Delta, initial_alpha=0.0)
        rel_p = float(np.linalg.norm(p_g - p_c) / np.linalg.norm(p_c))
        # the step's quality: the predicted reduction of the model with the CPU J (equal to rounding if equivalent)
        def pred(p):
            # predicted reduction of 0.5 |f + J_h p|^2 (the data block; the same J_h for both steps)
            Jp = Ja[:m_].dot(p)
            return float(-(fa[:m_] @ Jp + 0.5 * Jp @ Jp))
        rec = {"svd_seconds_gpu": tsg, "s_max_rel": float(np.max(np.abs(sg - s) / s[0])),
               "s_rel_per_value_max": float(np.max(np.abs(sg - s) / np.maximum(s, 1e-300))),
               "step_rel_l2": rel_p, "alpha": [a_c, a_g], "alpha_iterations": [int(it_c), int(it_g)],
               "step_norm": [float(np.linalg.norm(p_c)), float(np.linalg.norm(p_g))],
               "predicted_reduction": [pred(p_c), pred(p_g)], "tolerance_step_rel_l2": 1e-8}
        rec["passed"] = bool(rel_p <= 1e-8)
        ok2 = ok2 and rec["passed"]
        res["2_trf_step"]["modes"][mode] = rec
    res["2_passed"] = ok2
    # 3 the trajectory per mode vs the CPU reference
    with np.load(Path(args.cpu_ref) / "trajectory.npz", allow_pickle=False) as z:
        oc, xc = np.asarray(z["objective"]), np.asarray(z["x"])
    jref = json.load(open(Path(args.cpu_ref) / "trajectory.json"))
    res["3_trajectory"] = {"cpu_ref": str(Path(args.cpu_ref).resolve()), "cpu_seconds_per_evaluation":
                           jref["seconds_per_evaluation"], "cpu_host": jref["host"], "modes": {}}
    ok3 = True
    for mode in args.modes.split(","):
        r = run_trajectory(make, x, len(oc), "gpu", mode)
        og = r["objective"]
        k = min(len(og), len(oc))
        rel = np.abs(og[:k] / oc[:k] - 1)
        rec = {"evaluations": [int(len(oc)), int(len(og))], "max_rel_objective": float(rel.max()),
               "rel_objective_per_evaluation": [float(v) for v in rel], "final": [float(oc[-1]), float(og[-1])],
               "x_final_rel_l2": float(np.linalg.norm(r["x"] - xc) / np.linalg.norm(xc)),
               "gpu_seconds_per_evaluation": r["seconds_per_evaluation"], "device_record": r["device_record"],
               "message": r["message"], "tolerance_max_rel_objective": 1e-10}
        rec["passed"] = bool(len(og) == len(oc) and rel.max() <= 1e-10)
        ok3 = ok3 and rec["passed"]
        res["3_trajectory"]["modes"][mode] = rec
        np.savez(out / ("trajectory_gpu_%s.npz" % mode), objective=og, x=r["x"])
    res["3_passed"] = ok3
    res["passed"] = bool(res["1_passed"] and res["2_passed"] and res["3_passed"])
    json.dump(res, open(out / "verify.json", "x"), indent=1, default=float)
    short = {"1_res_rel": res["1_residuals"]["max_rel_to_max"], "1_jac_colrel": res["1_jacobian"]["max_columnwise_rel"],
             "2": {m: {k: v for k, v in r.items() if k in ("step_rel_l2", "svd_seconds_gpu", "s_max_rel")}
                   for m, r in res["2_trf_step"]["modes"].items()}, "2_svd_cpu_s": tsc,
             "3": {m: {k: v for k, v in r.items() if k in ("max_rel_objective", "gpu_seconds_per_evaluation", "evaluations")}
                   for m, r in res["3_trajectory"]["modes"].items()},
             "3_cpu_s_per_eval": jref["seconds_per_evaluation"], "passed": res["passed"]}
    print("verify: " + json.dumps(short, default=float))
    return 0 if res["passed"] else 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("trajectory")
    t.add_argument("--device", choices=("cpu", "gpu"), required=True)
    t.add_argument("--mode", default="qr", choices=("qr", "gesvd"))
    t.add_argument("--case", required=True)
    t.add_argument("--nfev", type=int, default=50)
    t.add_argument("--out", required=True)
    v = sub.add_parser("verify")
    v.add_argument("--case", required=True)
    v.add_argument("--cpu-ref", required=True)
    v.add_argument("--modes", default="qr,gesvd")
    v.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    return {"trajectory": cmd_trajectory, "verify": cmd_verify}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
