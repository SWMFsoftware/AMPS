#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nightfix_cpu_gpu_check.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nightfix_cpu_gpu_check.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1531050368ceda6ab0aaa26656a798053539f0040803eaa8241cdfaefd0faf82
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: check (2026-11-18, H2O night-side fix). fit_cops.Model / Problem vs the GPU evaluator with the fixed law on 16
#   used rows of the real 12M base at the B29 fit state: HostConstants Q bitwise, q_min/q_max 0.02, sequential mode
#   bitwise, tree mode within the VALIDATION.md budgets (512 eps M, 1024 eps M_J, residual bound); old-to-new change.
# USED BY: jobs/nightfix_cpu_gpu.sbatch (frozen model_ext_frozen/nightfix_check_20261118T073758Z).
# RUN: .venv/bin/python nightfix_cpu_gpu_check.py --base BASE --restart FIT/fit_restart.npz --rows 16 --out OUT.json (GPU)
# ==== end of header ====
"""H2O night-side fix (2026-11-18): CPU / GPU agreement of the fixed law on REAL data. The reference fitter's model
(fit_cops.Model / Problem, ../LiouvilleCOPS/tools/fit_cops.py, now q_min/q_max = 0.02) against the GPU evaluator
(cops_cuda.evaluator, HostConstants.qmin = 0.02) on a subset of rows of a real 12M fitting base, at a real fit state
(a fit_cops fit restart: x, the 25 degree-4 harmonics + T + kappa, no templates, e.g. the 12M B29 fit), with the
approved budgets of VALIDATION.md / validation/tolerances.json (the checks of tests/test_evaluator.py, C06):
  1 HostConstants: q per record = Model.Q bitwise (both species); q_min/q_max of H2O = 0.02 at every point
  2 sequential mode: pressures and every requested Jacobian column (activity, temperature, kappa) BITWISE = the CPU
  3 tree mode (production): |dp| <= 512 eps M (M = per-point sum of |contribution|); Jacobian |dJ| <= 1024 eps M_J;
    data residuals through the unmodified fit_cops.Problem within w 512 eps M / p
Also reported: the change of the pressures at this state from the old law (H2O q_min = 0.2 q_max; CPU Model with its Q
rebuilt from the old factor) to the new one, per row.
usage (GPU node, cwd $APP): nightfix_cpu_gpu_check.py --base BASE --restart FIT/fit_restart.npz --rows N --out OUT.json
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
for p in (str(APP), str(APP / "tools")):
    if p not in sys.path:
        sys.path.insert(0, p)
EPS = float(np.finfo(float).eps)
T0 = time.time()


def say(t):
    print("nightfix_cpu_gpu_check [%6.0f s]: %s" % (time.time() - T0, t), flush=True)


def magnitudes(P, x, KT):
    """Per-point sums of |contribution| (the M of the budget): copy of tests/test_evaluator.py magnitudes."""
    M = P.M
    a, b, T, kappa = P.unpack(x)
    act = M.Y @ a
    if P.nb:
        act = act + P.enhancement_scale * (P.enhancement_basis @ b)
    af = np.maximum(act, 0.0)[M.face]
    Tmax, T75, Tmin = T
    Tl = Tmax + KT * (Tmax - T75) * M.xlit
    u = 1.0 / np.sqrt(np.where(M.lit & (Tl > Tmin), Tl, Tmin))
    bc = lambda w: np.bincount(M.pt, weights=np.abs(w), minlength=M.npt)  # noqa: E731
    c, r, bm = M.cNG, M.cRG, M.bm
    mNG = abs(c["H2O"] * bm["H2O"]) * bc(M.wn["H2O"] * af * u) + abs(kappa * c["CO2"] * bm["CO2"]) * bc(M.wn["CO2"] * u)
    mRG = np.abs(r["H2O"]) * (bc(M.wp["H2O"] * af) + bm["H2O"] * bc(M.wv["H2O"] * af * u)) + \
        np.abs(kappa * r["CO2"]) * (bc(M.Q["CO2"] * M.ACO2[M.face] * M.KPhi) + bm["CO2"] * bc(M.wv["CO2"] * u))
    return mNG, mRG


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True)
    ap.add_argument("--restart", required=True)
    ap.add_argument("--rows", type=int, default=16)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    import argparse as _ap
    import cupy as cp
    from cops_cuda import reference, backend as bk, evaluator
    from c06_jacobian_budget import jacobian_magnitudes
    fit_cops = reference.load("fit_cops")
    run = str(Path(a.base).resolve())
    header, names, rows = fit_cops.read_model_file(run + "/COPS_model.dat")
    cols = {n: i for i, n in enumerate(names)}
    z = np.load(a.restart, allow_pickle=False)
    x, coef = np.asarray(z["x"], float), np.asarray(z["coef_prior"], float)
    if z["template_values"].size:
        raise SystemExit("the restart has templates; this check takes a template-free fit (B29-like)")
    use = [k for k, r in enumerate(rows) if int(r[cols["useNG"]]) or int(r[cols["useRG"]])]
    pick = sorted(set(np.round(np.linspace(0, len(use) - 1, a.rows)).astype(int).tolist()))
    sub = [rows[use[k]] for k in pick]
    say("base %s: %d rows, %d used; subset %d rows (iPoint %s)" % (run, len(rows), len(use), len(sub),
                                                                  [int(r[cols["iPoint"]]) for r in sub]))
    elements = fit_cops.read_elements(run + "/elements.bin")
    model = fit_cops.Model(run, header, cols, sub, elements, coef)
    model.load_kernels(model.sel)
    say("model: %d records" % len(model.face))
    res = {"schema": "nightfix_cpu_gpu_check_v1", "base": run, "restart": str(Path(a.restart).resolve()),
           "iPoint": [int(r[cols["iPoint"]]) for r in sub], "records": int(len(model.face)), "x": x.tolist(),
           "host": os.uname().nodename, "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "checks": {}}
    # 1 HostConstants vs Model.Q
    c = evaluator.HostConstants(model)
    for s in ("H2O", "CO2"):
        qmin, qmax = c.qmin[s][model.pt], c.qmax[s][model.pt]
        exp = np.where(model.night, qmin, qmin + (qmax - qmin) * np.minimum(model.cos, 1.0))
        res["checks"]["host_constants_Q_bitwise_%s" % s] = bool(exp.tobytes() == model.Q[s].tobytes())
    ratio = c.qmin["H2O"] / c.qmax["H2O"]
    res["checks"]["H2O_qmin_over_qmax_range"] = [float(ratio.min()), float(ratio.max())]
    res["checks"]["H2O_qmin_over_qmax_is_0.02"] = bool(np.all(np.abs(ratio - 0.02) <= 4 * EPS))
    args = _ap.Namespace(degree=4, enhancement_degree=4, fix_legacy_activity=False, fix_temperature=False,
                         fix_kappa=False, ordered_temperature=False, tbounds=(30.0, 400.0), enhancement_scale=None,
                         max_nfev=50, ftol=1e-10, xtol=1e-10, gtol=1e-10, jobs=1)
    P = fit_cops.Problem(model, elements, coef, args)
    if P.x0.size != x.size:
        raise SystemExit("problem has %d parameters, restart %d" % (P.x0.size, x.size))
    x = np.clip(x, P.lo, P.hi)
    state = bk.ParameterState.from_problem(P, x)
    request = bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True)
    ref_b = bk.ReferenceBackend.from_problem(P)
    want = ref_b.evaluate(state, request)
    # 2 sequential
    with evaluator.from_problem(P, run, mode="sequential") as gpu:
        got = gpu.evaluate(state, request)
        seq = {"pNG_bitwise": got.pNG.tobytes() == want.pNG.tobytes(), "pRG_bitwise": got.pRG.tobytes() == want.pRG.tobytes()}
        seq["columns_bitwise"] = all(g.tobytes() == w.tobytes() for key in request.keys()
                                     for g, w in zip(got.columns[key], want.columns[key]))
    res["checks"]["sequential"] = seq
    say("sequential: %s" % seq)
    cp.get_default_memory_pool().free_all_blocks()
    # 3 tree
    mNG, mRG = magnitudes(P, x, fit_cops.KT)
    r_cpu = P.data_residuals(x)
    with evaluator.from_problem(P, run, mode="tree") as gpu:
        ref = ref_b.evaluate(state)
        got = gpu.evaluate(state)
        dNG, dRG = np.abs(got.pNG - ref.pNG), np.abs(got.pRG - ref.pRG)
        tree = {"pressure_NG_within_512epsM": bool(np.all(dNG <= 512 * EPS * mNG)),
                "pressure_RG_within_512epsM": bool(np.all(dRG <= 512 * EPS * mRG)),
                "pressure_NG_max_in_epsM": float(np.max(dNG / (EPS * mNG))), "pressure_RG_max_in_epsM":
                    float(np.max(dRG / (EPS * mRG))), "pressure_max_rel": float(max(np.max(dNG / ref.pNG),
                                                                                   np.max(dRG / ref.pRG)))}
        ref_j = ref_b.evaluate(state, request)
        got_j = gpu.evaluate(state, request)
        worst, ok = 0.0, True
        for key, mags in jacobian_magnitudes(P, P.M, x, request, fit_cops).items():
            for g, w, m in zip(got_j.columns[key], ref_j.columns[key], mags):
                ok &= bool(np.all(np.abs(g - w) <= 1024 * EPS * m))
                nz = m > 0
                if np.any(nz):
                    worst = max(worst, float(np.max(np.abs(g - w)[nz] / (EPS * m[nz]))))
        tree["jacobian_within_1024epsMJ"], tree["jacobian_max_in_epsMJ"] = ok, worst
        bk.install(P, gpu)
        r_gpu = P.data_residuals(x)
    floor = 1e-300
    delta = np.concatenate([P.wNG * 512 * EPS * mNG[P.useNG] / np.maximum(ref.pNG[P.useNG], floor),
                            P.wRG * 512 * EPS * mRG[P.useRG] / np.maximum(ref.pRG[P.useRG], floor)])
    tree["residuals_within_bound"] = bool(np.all(np.abs(r_gpu - r_cpu) <= delta))
    res["checks"]["tree"] = tree
    say("tree: %s" % tree)
    # the old law at the same state (CPU): rebuild Model.Q of H2O with q_min = 0.2 q_max
    q_new = model.Q["H2O"].copy()
    qmax_p = c.qmax["H2O"][model.pt]
    qmin_old = 0.2 * 600 / (model.rh / fit_cops.AU) ** (4.2143229 if model.helio == "shap5" else 4.8)
    qmin_old = qmin_old[model.pt]
    cz = np.minimum(model.cos, 1.0)
    model.Q["H2O"] = np.where(model.night, qmin_old, qmin_old + (qmax_p - qmin_old) * cz)
    P_old = fit_cops.Problem(model, elements, coef, args)
    old = bk.ReferenceBackend.from_problem(P_old).evaluate(bk.ParameterState.from_problem(P_old, x))
    model.Q["H2O"] = q_new
    res["old_to_new_at_this_state"] = {"NG_new_over_old": (ref.pNG / old.pNG).tolist(),
                                       "RG_new_over_old": (ref.pRG / old.pRG).tolist(),
                                       "utc": [r[cols["UTC"]] for r in sub]}
    res["passed"] = bool(res["checks"]["host_constants_Q_bitwise_H2O"] and res["checks"]["host_constants_Q_bitwise_CO2"]
                         and res["checks"]["H2O_qmin_over_qmax_is_0.02"] and all(seq.values()) and
                         tree["pressure_NG_within_512epsM"] and tree["pressure_RG_within_512epsM"] and
                         tree["jacobian_within_1024epsMJ"] and tree["residuals_within_bound"])
    res["wall_seconds"] = time.time() - T0
    json.dump(res, open(out, "x"), indent=1, default=float)
    say("passed %s -> %s" % (res["passed"], out))
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
