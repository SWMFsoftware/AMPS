#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_ac.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_ac.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: a7c475213e569644c79f36cafdc9272a690415af8451f4ffc4349709dde9a2cd
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: helper of the slope-test follow-up a + c (2026-11-22): starts, preflight, pick, summary
# RUN: jobs/localnn_ac_*.sbatch
# ==== end of header ====
"""Helper of the follow-up fits 'a' and 'c' of the slope test (user 2026-11-22: "as recommended: a, c"; FINDINGS A22).
All of option C f10 / P_m 0.01 / templates off, SHAP8 12M, new H2O law, the fitter and code of the control fits.

  a  six multipliers free, >= 0.01, started from each of the five solutions (are cliffs less active than flat terrain?)
  c  the candidate model with all six multipliers held at 1: fitted from other starts, then a five-cycle chain with fine
     tuning from the lowest converged of these fits and the existing one from D2 (does it stay active and stable?)

Starts of c: the activity map c_j of each of the five solutions with the multipliers set to 1, and a uniform map. Each
is scaled to the seen production of the existing all-held fit at 10:33:10 UTC; kappa and grain from that fit.

  make-starts --out DIR            the derived start states of c (evaluations only)
  preflight --code CODE --starts DIR   evaluations only: the fitter's Jacobian check at all eleven starts in their own
                                   problems, and the reduced Jacobian with a non-zero held fine-tuning term
  after-stage STAGE_FIT            exit 0 finished, 10 stopped at max_nfev (one restart follows)
  pick --runroot R                 the fit the chain starts from; writes R/chain_pick.json, prints the fit folder
  states --runroot R               name=state list of every finished fit, for q002_stationarity_check
  summary --runroot R --out DIR    table of all fits and of the chain
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R = Path("/data001/highresolu_vt/cops_output")
N = R / "runs/shap8_012m_lategrid_dz1.5"
L = N / "nonneg/localnn"
GRID = L / "grid_150m"
DESIGN = GRID / "design_h2oqmin0.02_d2q002"
BASE = N / "base_q002"
D2 = N / "clip_refit/D2_q002_warm/seg01"
PT = N / "late_grid_points.txt"
LATE = L / "late_kernels_150m_h2oqmin0.02.npz"
CONTROLS = L / "slopetest_q002_f10_pm0.01_tploff_20261122T044019Z"
MULTI = L / "multistart_q002_f10_pm0.01_tploff_20261122T065041Z"
REF = L / "stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit/final"
HELD_TEST = L / "stage2r_q002_f10_pm0.01_tplon_20261118T162407Z/c5/ft/seg01/finetune_state.npz"     # preflight only
CLASSES = ("40_50", "50_60", "60_70", "70_80", "80_90")
A_FITS = ["a_from_" + v for v in CLASSES]
C_FITS = ["c_from_" + v for v in CLASSES + ("uniform",)]
LAM, FM, L_MIN = 8.249495752351283e-05, 0.01, 0.01


def solution(cls):
    return (REF if cls == "40_50" else MULTI / ("from_only_" + cls) / "fit/final").resolve() / "localnn_state.npz"


def cmd_make_starts(a):
    sys.path.insert(0, a.code)
    import localnn_problem as lp
    sys.path.append(str(HERE))
    import q002_slope_activity as geo
    out = Path(a.out)
    out.mkdir()
    held = (CONTROLS / "none/fit/final").resolve()
    with np.load(held / "nonneg_facets.npz", allow_pickle=False) as z:
        A0, seen, area = np.asarray(z["A"]), np.asarray(z["seen"], bool), np.asarray(z["area"])
    with np.load(held / "localnn_state.npz", allow_pickle=False) as z:
        base = {k: np.asarray(z[k]) for k in z.files}
    meta = json.loads(str(base["meta_json"]))
    t = geo.map_times()
    label = [k for k, v in t.items() if v["utc"].endswith("10:33:10")][0]
    w = geo.qfac(label, t[label]["rh_m"], len(area))[0] * area * seen
    target = float((A0 * w).sum())
    grid = lp.load_grid(GRID)
    As = float(meta["A_scale"])
    nc = len(base["x"]) - (lp.NBINS - 1) - 3
    rec = {}
    for name in CLASSES + ("uniform",):
        x = np.array(base["x"], float)                      # kappa and grain of the all-held fit
        if name == "uniform":
            src, c = None, np.ones(nc)
        else:
            src = str(solution(name))
            with np.load(src, allow_pickle=False) as z:
                c = np.asarray(z["x"], float)[:nc] * As
        c = c * target / float((np.maximum(lp.facet_field(grid, c), 0.0) * w).sum())
        x[:nc] = c / As
        x[nc:nc + lp.NBINS - 1] = 1.0
        check = float((np.maximum(lp.facet_field(grid, c), 0.0) * w).sum())
        info = dict(name=name, map_from=src, multipliers="all 1", scaled_to_seen_production_103310=target,
                    seen_production=check, kappa_and_grain_from=str(held))
        m = dict(meta, segment=None, costs=None, derived_start=info, created_utc=lp.utc_now())
        m.pop("slope_test", None)
        d = out / name
        d.mkdir()
        np.savez(d / "localnn_state.npz", x=x, names=base["names"], c=c, m=np.ones(lp.NBINS),
                 meta_json=json.dumps(m, sort_keys=True, default=float))
        rec[name] = dict(info, sha256=lp.sha256_file(d / "localnn_state.npz"))
        print("start %s: seen production %.4e (target %.4e)" % (name, check, target), flush=True)
    (out / "starts.json").write_text(json.dumps(rec, indent=1) + "\n")
    return 0


def cmd_preflight(a):
    sys.path.insert(0, a.code)
    import localnn_fit as lf
    import localnn_problem as lp
    design = lp.load_design(DESIGN)
    D = lp.Data(BASE, D2, PT, design)
    grid = lp.load_grid(GRID)
    As = design["meta"]["A_scale"]
    args = argparse.Namespace(held=str(HELD_TEST), late_kernels=str(LATE))
    held, _ = lf.held_setup(args, design)
    cases = [("a_from_" + v, solution(v), L_MIN, [], None) for v in CLASSES]
    cases += [("c_from_" + v, Path(a.starts) / v / "localnn_state.npz", 1.0, "all", None) for v in CLASSES + ("uniform",)]
    cases.append(("chain, held h non-zero", (CONTROLS / "none/fit/final").resolve() / "localnn_state.npz", 1.0, "all", held))
    bad = 0
    for name, sf, lo, hold, h in cases:
        kw = {"held": h} if h is not None else {}
        P = lp.Problem(D, design, grid, False, lo, As, LAM, mult_penalty_f=FM, **kw)
        with np.load(sf, allow_pickle=False) as z:
            x = np.asarray(z["x"], float)
        ok = bool(np.all(x >= P.lo))
        idx = list(range(P.im.start, P.im.stop)) if hold else []
        if hold:
            ok = ok and bool(np.all(x[idx] == 1.0))
        fun, jac, y, bounds, expand = lf.reduced_problem(P, x, idx)
        r = fun(y)
        cols = np.ones(len(x), bool)
        cols[idx] = False
        ok = ok and np.array_equal(jac(y), P.jac(x)[:, cols]) and np.array_equal(expand(y), x)
        cost, _ = P.evaluate(x)
        chk = lf.jacobian_check(P, x)
        ok = ok and chk["passed"]
        bad += not ok
        print("%-24s free %d lower bound %g: objective %.10e data %.10e Jacobian %.1e %s" % (
            name, len(y), lo, 0.5 * float(r @ r), cost["data_cost"], chk["max_rel_error"], "PASS" if ok else "FAIL"),
            flush=True)
    print("ALL PASSED" if not bad else "FAILED %d" % bad)
    return 1 if bad else 0


def manifest(fit):
    return json.loads(((Path(fit) / "final").resolve() / "fit_manifest.json").read_text())


def at_limit(man):
    return "maximum number of function evaluations" in man["message"]


def cmd_after_stage(a):
    return 10 if at_limit(manifest(a.stage_fit)) else 0


def last_fit(root, name):
    """The folder of the fit's last stage (stage2 = the one restart after max_nfev), None if it is not finished."""
    d = Path(root) / name
    fit = d / "stage2/fit" if (d / "stage2").exists() else d / "stage1/fit"
    if not (fit / "final").exists():
        return None
    if not (fit / ("direct_check_%s" % (fit / "final").resolve().name) / "direct_check.json").exists():
        return None
    if fit.parent.name == "stage1" and at_limit(manifest(fit)):
        return None                                         # the restart is due
    return fit


def cmd_pick(a):
    root = Path(a.runroot)
    out = root / "chain_pick.json"
    if out.exists():
        print(json.loads(out.read_text())["fit"])
        return 0
    cand = {"c_from_d2": CONTROLS / "none/fit"}
    for name in C_FITS:
        fit = last_fit(root, name)
        if fit is None:
            print("not finished: %s" % name, file=sys.stderr)
            return 3
        cand[name] = fit
    mans = {k: manifest(v) for k, v in cand.items()}
    obj = {k: float(m["costs"]["objective"]) for k, m in mans.items()}
    best = min(obj, key=obj.get)
    if obj["c_from_d2"] <= obj[best] * (1.0 + 1e-6):        # the same solution: keep the existing fit
        best = "c_from_d2"
    if at_limit(mans[best]):
        print("the lowest fit (%s) is not converged: no chain" % best, file=sys.stderr)
        return 4
    out.write_text(json.dumps(dict(fit=str(cand[best].resolve()), name=best, objectives=obj,
                                   rule="lowest objective; within 1e-6 the existing fit from D2"), indent=1) + "\n")
    print(cand[best].resolve())
    return 0


def cmd_states(a):
    out = []
    for name in A_FITS + C_FITS:
        fit = last_fit(a.runroot, name)
        if fit is not None:
            out.append("%s=%s/localnn_state.npz" % (name, (fit / "final").resolve()))
    print(" ".join(out))
    return 0


def cmd_summary(a):
    sys.path.append(str(HERE))
    import q002_slopetest_report as rep
    out = Path(a.out)
    out.mkdir()
    T = json.loads(rep.AA.read_text())["T"]
    rows, waiting = [], []
    for name in A_FITS + C_FITS:
        fit = last_fit(a.runroot, name)
        if fit is None:
            waiting.append(name)
            continue
        r = rep.row(name, (fit / "final").resolve(), T)
        r["stage"] = fit.parent.name
        rows.append(r)
    rows.append(dict(rep.row("c_from_d2", (CONTROLS / "none/fit/final").resolve(), T), stage="existing"))
    chain = []
    for k in range(1, 6):
        d = Path(a.runroot) / "chain" / ("c%d" % k)
        if not (d / "ft/final").exists():
            waiting.append("chain c%d" % k)
            continue
        ft = (d / "ft/final").resolve()
        cost = json.loads((ft / "fine_tuning.json").read_text())["costs_final"]
        seg = (d / "fit/final").resolve()
        with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
            seen, area = np.asarray(z["seen"], bool), np.asarray(z["area"])
        with np.load(ft / "finetune_state.npz", allow_pickle=False) as z:
            A = np.asarray(z["A_final"])
            hmax = float(np.abs(z["held_next"]).max())
        chain.append(dict(cycle=k, fit_objective=manifest(d / "fit")["costs"]["objective"], data_cost=cost["cost_data"],
                          late_rms_NG=cost["late_window_rms_ln_NG"], late_rms_RG=cost["late_window_rms_ln_RG"],
                          max_abs_h=hmax, off_seen_1b=float(area[seen & (A <= T)].sum() / area[seen].sum())))
    (out / "followup.json").write_text(json.dumps(dict(rows=rows, chain=chain, waiting=waiting), indent=1) + "\n")
    with open(out / "followup.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("%-18s %-8s obj %.6e data %.6e late %.4f/%.4f off %.1f %% maxF/F_E %.3f %s | m %s" % (
            r["variant"], r["stage"], r["objective"], r["data_cost"], r["late_rms_NG"], r["late_rms_RG"],
            100 * r["off_seen_1b"], r["F_max_over_F_E"], r["termination"][:28],
            " ".join("%.4g" % r["m_" + b] for b in rep.BINS)))
    for r in chain:
        print("chain c%d data %.6e late %.4f/%.4f off %.1f %% max|h| %.3g" % (
            r["cycle"], r["data_cost"], r["late_rms_NG"], r["late_rms_RG"], 100 * r["off_seen_1b"], r["max_abs_h"]))
    print("waiting:", waiting or "none")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("make-starts")
    p.add_argument("--out", required=True)
    p.add_argument("--code", required=True)
    p = sub.add_parser("preflight")
    p.add_argument("--code", required=True)
    p.add_argument("--starts", required=True)
    p = sub.add_parser("after-stage")
    p.add_argument("stage_fit")
    for n in ("pick", "states"):
        sub.add_parser(n).add_argument("--runroot", required=True)
    p = sub.add_parser("summary")
    p.add_argument("--runroot", required=True)
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    return {"make-starts": cmd_make_starts, "preflight": cmd_preflight, "after-stage": cmd_after_stage,
            "pick": cmd_pick, "states": cmd_states, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
