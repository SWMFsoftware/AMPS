#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_multistart.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_multistart.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 09e38f73b82e0f789ce46ecacb17e9bb1a9e0d20fe53424e6d8e250bd9a7e696
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: rule and checks of the option C multi-start refits (user 2026-11-22)
# RUN: jobs/localnn_multistart.sbatch
# ==== end of header ====
"""Multi-start refits of the option C reference problem (user 2026-11-22: "multi-start when ready"): the rule and the
checks. The fits themselves are run by the reference's UNCHANGED frozen localnn_fit.py (start --start state, fit).

Approved: rerun the all-free fit (six multipliers >= 1) from each better end state whose multipliers are all >= 1, up
to five fits.

  gate        may this parent be used? complete (final + passed direct check), objective below the reference's, every
              multiplier >= 1. only_40_50 was not below the reference when the five were proposed: it runs only if
              fewer than five of the other candidates qualify, so there are never more than five fits.
              exit 0 run, 10 skip (the reason is printed), 3 not decidable (incomplete parent)
  check-start the start folder holds the parent's end state unchanged
  preflight   evaluations only: the start-up Jacobian check of the frozen fitter at the given states

    localnn_multistart.py gate --controls <slopetest root> --parent only_60_70
    localnn_multistart.py check-start --state <seg/localnn_state.npz> --start <start folder>
    localnn_multistart.py preflight --code <frozen code> STATE.npz [STATE.npz ...]
"""
import argparse
import json
import sys
from pathlib import Path

R = Path("/data001/highresolu_vt/cops_output")
N = R / "runs/shap8_012m_lategrid_dz1.5"
L = N / "nonneg/localnn"
REF_OBJECTIVE = 6.949871848720e-03          # stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit/final
CANDIDATES = ("only_50_60", "only_60_70", "only_70_80", "only_80_90", "free_min")       # the five proposed
EXTRA = "only_40_50"
MAX_FITS = 5
RUN, SKIP, INCOMPLETE = 0, 10, 3


def decide_one(root, parent):
    fit = Path(root) / parent / "fit"
    final = fit / "final"
    if not final.exists():
        return INCOMPLETE, "%s has no final state" % parent
    seg = final.resolve()
    check = fit / ("direct_check_%s" % seg.name) / "direct_check.json"
    if not check.exists() or not json.loads(check.read_text()).get("passed"):
        return INCOMPLETE, "%s has no passed direct check" % parent
    man = json.loads((seg / "fit_manifest.json").read_text())
    obj, m = float(man["costs"]["objective"]), [float(v) for v in man["multipliers"]]
    if not obj < REF_OBJECTIVE:
        return SKIP, "%s: objective %.10e is not below the reference's %.10e" % (parent, obj, REF_OBJECTIVE)
    if min(m) < 1.0:
        return SKIP, "%s: a multiplier is below 1 (smallest %.6g)" % (parent, min(m))
    return RUN, "%s: objective %.10e, smallest multiplier %.6g" % (parent, obj, min(m))


def decide(root, parent):
    if parent not in CANDIDATES + (EXTRA,):
        return SKIP, "%s is not a candidate" % parent
    code, why = decide_one(root, parent)
    if code != RUN or parent != EXTRA:
        return code, why
    others = [decide_one(root, p) for p in CANDIDATES]
    for c, w in others:
        if c == INCOMPLETE:
            return INCOMPLETE, "cannot count the candidates: " + w
    n = sum(c == RUN for c, _ in others)
    if n >= MAX_FITS:
        return SKIP, "%s: qualifies, but %d others do and the approval is up to %d fits" % (parent, n, MAX_FITS)
    return RUN, why


def cmd_gate(a):
    code, why = decide(a.controls, a.parent)
    print({RUN: "RUN", SKIP: "SKIP", INCOMPLETE: "INCOMPLETE"}[code] + " " + why)
    return code


def cmd_check_start(a):
    import numpy as np
    with np.load(a.state, allow_pickle=False) as z:
        x = np.asarray(z["x"], float)
        want = json.loads(str(z["meta_json"]))["costs"]["data_cost"]
    with np.load(Path(a.start) / "start.npz", allow_pickle=False) as z:
        x0 = np.asarray(z["x0"], float)
        got = json.loads(str(z["meta_json"]))["costs_at_start"]["data_cost"]
    ok = np.array_equal(x, x0) and abs(got / want - 1.0) < 1e-12
    print("start = parent's end state: %s (data cost %.12e, parent %.12e)" % (ok, got, want))
    return 0 if ok else 4


def cmd_preflight(a):
    sys.path.insert(0, a.code)
    import numpy as np
    import localnn_fit as lf
    import localnn_problem as lp
    design = lp.load_design(L / "grid_150m/design_h2oqmin0.02_d2q002")
    D = lp.Data(N / "base_q002", N / "clip_refit/D2_q002_warm/seg01", N / "late_grid_points.txt", design)
    grid = lp.load_grid(L / "grid_150m")
    P = lp.Problem(D, design, grid, False, 1.0, design["meta"]["A_scale"], 8.249495752351283e-05, mult_penalty_f=0.01)
    bad = 0
    for s in a.states:
        with np.load(s, allow_pickle=False) as z:
            x = np.asarray(z["x"], float)
        x0 = np.minimum(np.maximum(x, P.lo), P.hi)                      # as cmd_fit
        r = P.fun(x0)
        chk = lf.jacobian_check(P, x0)
        bad += not chk["passed"]
        print("%s: smallest m %.6g, objective at the start %.10e, Jacobian check %s (max rel %.2e)" % (
            s, float(x[P.im].min()), 0.5 * float(r @ r), "passed" if chk["passed"] else "FAILED", chk["max_rel_error"]),
            flush=True)
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate")
    g.add_argument("--controls", required=True)
    g.add_argument("--parent", required=True)
    c = sub.add_parser("check-start")
    c.add_argument("--state", required=True)
    c.add_argument("--start", required=True)
    p = sub.add_parser("preflight")
    p.add_argument("--code", required=True)
    p.add_argument("states", nargs="+")
    a = ap.parse_args()
    return {"gate": cmd_gate, "check-start": cmd_check_start, "preflight": cmd_preflight}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
