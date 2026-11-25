#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_ac_preview.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_ac_preview.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f6dbf663672719c25885f2d8238fbf6218e0714031571d939f4a933d6d67461e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: what the follow-up a and c can achieve: one bounded step >= 0.01 at the five solutions, late-facet leverage
# RUN: jobs/q002_ac_preview.sbatch
# ==== end of header ====
"""What can the planned follow-up fits a and c achieve? (2026-11-22; the user: "what is expected to achieve with a and
c?".) Evaluations only, no optimizer iteration and no fit.

  a  at each of the five solutions: the bounded Gauss-Newton step with the multipliers' lower bound lowered from 1 to
     0.01 (one linearized step, exact NNLS), the decrease really available along it and the multipliers it points to.
     If promised = available, the step is a fair preview of the minimum nearby; it says nothing about distant minima.
  c  the chain's fine tuning can change only the 271 late facets (factor <= 100): their share of the all-held fit's
     surface pressure at the late rows, i.e. the leverage it has on the late NG misfit.

    q002_ac_preview.py --tools <frozen ac folder> --out <new folder>
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import nnls

R = Path("/data001/highresolu_vt/cops_output")
CODE = R / "model_ext_frozen/localnn_slopetest_20261122T033334Z"
ALPHAS = (1.0, 0.75, 0.5, 0.25, 0.1, 0.03, 0.01)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tools", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir()
    sys.path.insert(0, a.tools)
    import localnn_ac as ac
    sys.path.insert(0, str(CODE))
    import localnn_late as llate
    import localnn_problem as lp
    design = lp.load_design(ac.DESIGN)
    D = lp.Data(ac.BASE, ac.D2, ac.PT, design)
    grid = lp.load_grid(ac.GRID)
    P = lp.Problem(D, design, grid, False, ac.L_MIN, design["meta"]["A_scale"], ac.LAM, mult_penalty_f=ac.FM)
    rec = {"a": {}, "c": {}}
    im = np.arange(P.im.start, P.im.stop)
    for cls in ac.CLASSES:
        with np.load(ac.solution(cls), allow_pickle=False) as z:
            x = np.asarray(z["x"], float)
        r, J = P.fun(x), P.jac(x)
        f0, d0 = 0.5 * float(r @ r), P.evaluate(x)[0]["data_cost"]
        nr = np.linalg.norm(J, axis=0)
        use = nr > 0
        S, low = J[:, use] / nr[use], (P.lo[use] - x[use]) * nr[use]
        y, rn = nnls(S, -(r + S @ low), maxiter=20 * S.shape[1])
        d = np.zeros(len(x))
        d[use] = (y + low) / nr[use]
        tr = []
        for al in ALPHAS:
            xt = np.maximum(x + al * d, P.lo)
            rt = P.fun(xt)
            tr.append((0.5 * float(rt @ rt), al, xt))
        f1, al, x1 = min(tr, key=lambda t: t[0])
        if f1 >= f0:
            f1, al, x1 = f0, 0.0, x
        d1 = P.evaluate(x1)[0]["data_cost"]
        v = dict(objective=f0, promised_rel=(f0 - 0.5 * rn * rn) / f0, available_rel=(f0 - f1) / f0, alpha=al,
                 data_cost=d0, data_cost_change_rel=d1 / d0 - 1.0, multipliers=x[im].tolist(),
                 multipliers_full_step=(x + d)[im].tolist(), multipliers_best=x1[im].tolist())
        rec["a"][cls] = v
        print("a %s: promised %.3e, available %.3e (alpha %g) of the objective; data cost %+.3e" % (
            cls, v["promised_rel"], v["available_rel"], al, v["data_cost_change_rel"]))
        print("   m now       " + " ".join("%10.4g" % q for q in x[im]))
        print("   m full step " + " ".join("%10.4g" % q for q in (x + d)[im]), flush=True)

    seg = (ac.CONTROLS / "none/fit/final").resolve()
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        c, m, kappa = np.asarray(z["c"]), np.asarray(z["m"]), float(z["kappa"])
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        A = np.asarray(z["A"], float)
    LK = llate.load(ac.LATE, len(A))
    sN, sR = lp.all_rows_surface(design, c, m, kappa, False)
    li = LK["late_idx"]
    rows = np.asarray(D.layout.late_rows, int)
    for g, G, s in (("NG", LK["G_NG"], sN), ("RG", LK["G_RG"], sR)):
        sh = (G @ A[li])[rows] / s[rows]
        rec["c"][g] = dict(late_rows=len(rows), late_facets=len(li), mean_share=float(sh.mean()),
                           max_share=float(sh.max()), share_last_row=float(sh[np.argmax([D.layout.et[i] for i in rows])]))
        print("c %s: %d late facets give mean %.3e, max %.3e, last observation %.3e of the surface pressure at %d late rows" % (
            g, len(li), sh.mean(), sh.max(), rec["c"][g]["share_last_row"], len(rows)))
    (out / "ac_preview.json").write_text(json.dumps(rec, indent=1) + "\n")


if __name__ == "__main__":
    main()
