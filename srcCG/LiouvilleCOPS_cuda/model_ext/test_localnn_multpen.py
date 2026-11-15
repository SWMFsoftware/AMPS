#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_localnn_multpen.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_localnn_multpen.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 28fd14f7872b76fffda25ae65527bfa429599f01bc1f29acebc4b93fe1822fa3
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tests (option C, CPU) of the multiplier penalty: option absent = bitwise vs the frozen previous code (fun, jac,
#   3-evaluation TRF), P_m rows = f_m C_ref sum (ln m)^2 on the 0.5 r.r scale, Jacobian vs finite differences.
# RUN: srun 2026-11-15 ~01:5x UTC on the f10 end state -> tests/localnn_20261114/multpen/tests_multpen.json (4 / 4 passed).
# ==== end of header ====
"""Tests of the option C multiplier penalty (user 2026-11-14 ~17:30 CST): P_m = f_m C_ref sum_b (ln m_b)^2,
localnn_problem.Problem(..., mult_penalty_f=f_m), localnn_fit.py fit --mult-penalty-f. CPU.

    test_localnn_multpen.py --grid GRID --design DESIGN --base BASE --d2 D2 --point-table PT --frozen OLD_SNAPSHOT
                            --state SEG/localnn_state.npz [--out report.json]

 1 option absent = bitwise: the new Problem without the argument and with mult_penalty_f = 0 vs the frozen previous
   localnn_problem.py (OLD_SNAPSHOT): fun and jac bitwise equal at the given state and at a perturbed state; a 3-evaluation
   TRF run (least_squares as localnn_fit) gives bitwise the same x
 2 P_m rows: with f_m = 0.01 fun has 6 more rows; 0.5 |rows|^2 = P_m = f_m C_ref sum (ln m_b)^2 (rel <= 1e-14);
   Problem.evaluate's objective = data + crossing + penalty + P_m
 3 Jacobian vs central finite differences with P_m on (localnn_fit.jacobian_check: 24 node columns, every multiplier,
   kappa, grain) at the state and at a state with all multipliers x 1.7 (max rel <= 1e-5); the P_m rows' derivative
   sqrt(2 f_m C_ref) / m_b exactly
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402
import localnn_fit as lf          # noqa: E402

RES = []


def check(name, ok, **info):
    RES.append({"test": name, "passed": bool(ok), **info})
    print("%s %s %s" % ("PASS" if ok else "FAIL", name, json.dumps(info, default=float)), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    for k in ("grid", "design", "base", "d2", "point-table", "frozen", "state"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    spec = importlib.util.spec_from_file_location("localnn_problem_frozen", str(Path(args.frozen) / "localnn_problem.py"))
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    design = lp.load_design(args.design)
    grid = lp.load_grid(args.grid)
    D = lp.Data(args.base, args.d2, args.point_table, design)
    with np.load(args.state, allow_pickle=False) as z:
        x = np.asarray(z["x"], float)
        meta = json.loads(str(z["meta_json"]))
    tm, lam, A_s = meta["templates"] == "on", float(meta["lambda_C"]), design["meta"]["A_scale"]
    Po = old.Problem(old.Data(args.base, args.d2, args.point_table, design), design, grid, tm, 1.0, A_s, lam)
    Pn = lp.Problem(D, design, grid, tm, 1.0, A_s, lam)
    Pz = lp.Problem(D, design, grid, tm, 1.0, A_s, lam, mult_penalty_f=0.0)
    rng = np.random.default_rng(3)
    xp = x * (1.0 + 0.01 * rng.standard_normal(len(x)))
    xp = np.minimum(np.maximum(xp, Pn.lo + 1e-9), Pn.hi)
    ok = True
    det = {}
    for lab, xx in (("state", x), ("perturbed", xp)):
        fo, jo = Po.fun(xx), Po.jac(xx)
        for nm_, P in (("new", Pn), ("new f_m=0", Pz)):
            e = bool(np.array_equal(fo, P.fun(xx)) and np.array_equal(jo, P.jac(xx)))
            det["%s %s" % (lab, nm_)] = e
            ok &= e
    runs = {}
    for nm_, P in (("frozen", Po), ("new", Pn), ("new f_m=0", Pz)):
        P.nfev, P.best = 0, None
        r = least_squares(P.fun, x, jac=P.jac, bounds=(P.lo, P.hi), x_scale="jac", method="trf", tr_solver="exact",
                          ftol=lp.TOL, xtol=lp.TOL, gtol=lp.TOL, max_nfev=3)
        runs[nm_] = r.x
    det["trf3_new_equals_frozen"] = bool(np.array_equal(runs["frozen"], runs["new"]))
    det["trf3_fm0_equals_frozen"] = bool(np.array_equal(runs["frozen"], runs["new f_m=0"]))
    ok &= det["trf3_new_equals_frozen"] and det["trf3_fm0_equals_frozen"]
    check("1 option absent = bitwise (fun, jac, 3-evaluation TRF)", ok, **det)
    fm = 0.01
    Pm = lp.Problem(D, design, grid, tm, 1.0, A_s, lam, mult_penalty_f=fm)
    f0, f1 = Pn.fun(x), Pm.fun(x)
    rows = f1[len(f0):]
    lm = np.log(x[Pm.im])
    pm_formula = fm * lp.C_REF_C9 * float(lm @ lm)
    cst, _ = Pm.evaluate(x)
    cst0, _ = Pn.evaluate(x)
    rel1 = abs(0.5 * float(rows @ rows) / pm_formula - 1.0)
    rel2 = abs(cst["mult_penalty"] / pm_formula - 1.0)
    obj_ok = cst["objective"] == cst0["objective"] + cst["mult_penalty"]
    check("2 P_m rows and bookkeeping", len(rows) == 6 and np.array_equal(f1[:len(f0)], f0) and rel1 <= 1e-14 and
          rel2 <= 1e-14 and obj_ok, extra_rows=int(len(rows)), P_m=pm_formula, rel_rows=rel1, rel_evaluate=rel2,
          objective_ok=bool(obj_ok), C_ref=lp.C_REF_C9, multipliers=[float(v) for v in x[Pm.im]])
    Jm = Pm.jac(x)
    exact = np.allclose(Jm[-6:, Pm.im], np.diag(np.sqrt(2 * fm * lp.C_REF_C9) / x[Pm.im]), rtol=0, atol=0) and \
        not np.any(np.delete(Jm[-6:], np.arange(Pm.im.start, Pm.im.stop), axis=1))
    for lab, xx in (("state", x), ("multipliers x 1.7", np.where((np.arange(len(x)) >= Pm.im.start) &
                                                                  (np.arange(len(x)) < Pm.im.stop), x * 1.7, x))):
        Pm.deadline = None
        chk = lf.jacobian_check(Pm, xx)
        check("3 Jacobian vs finite differences with P_m (%s)" % lab, chk["passed"] and exact,
              max_rel=chk["max_rel_error"], columns=chk["columns"], exact_pm_rows=bool(exact),
              multiplier_columns={k: v for k, v in chk["per_column"].items() if k.startswith("c_slope")})
    allok = all(r["passed"] for r in RES)
    print("%d tests, %d failed" % (len(RES), sum(not r["passed"] for r in RES)))
    if args.out:
        with open(args.out, "x") as f:
            json.dump({"results": RES, "passed": allok, "argv": sys.argv}, f, indent=1, default=float)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
