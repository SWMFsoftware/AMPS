#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_localnn_slopetest.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_localnn_slopetest.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 55d16441677be93143cbaa08e2dc169bb4667d2a80450e9de9a3d14db925a6a0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tests of localnn_fit.py --slope-test; no optimizer runs (decision 3c)
# RUN: jobs/localnn_slopetest_code_test.sbatch
# ==== end of header ====
"""Tests of localnn_fit.py --slope-test (user 2026-11-21, decision 3c). NO optimizer runs: scipy's least_squares is
replaced by a stub that inspects what the fitter hands to it and returns a given state.

 1 default path: new and old frozen code hand bitwise identical fun(x0), jac(x0), x0 and bounds to the optimizer, and
   write the same costs for the reference state (= the reference fit's recorded costs)
 2 each variant (none, six one-free): the free parameters and bounds; jac = the full Jacobian's free columns and agrees
   with finite differences; held multipliers exactly 1 in the products and not listed as active bounds
 3 lowered bound (TEST value 0.5, no decision) at a state with m = 0.7: bounds, derivatives, products
 4 wall-time stop and resume: refused with other settings, completes with the same
 5 the guards refuse what is not decided

    test_localnn_slopetest.py --code NEW --old-code OLD --reference <fit/segNN> --out DIR -- <the fit's common options>
"""
import argparse
import contextlib
import hashlib
import io
import json
import shutil
import subprocess
import sys
import types
from pathlib import Path

import numpy as np

BINS = ["slope_40_50", "slope_50_60", "slope_60_70", "slope_70_80", "slope_80_90", "slope_90_180"]


def sha(a):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(a, float)).tobytes()).hexdigest()


def worker(job):
    sys.path.insert(0, job["code"])
    import localnn_fit as lf
    import localnn_problem as lp

    rec = {}
    box = {}
    make = lf.make_problem

    def keep(*a, **k):
        box["P"], g = make(*a, **k)
        return box["P"], g

    lf.make_problem = keep
    hold = job.get("hold", [])

    def stub(fun, x0, jac, bounds, **kw):
        P = box["P"]
        free = np.ones(len(P.lo), bool)
        free[[P.im.start + BINS.index(b) for b in hold]] = False
        r, J = fun(x0), jac(x0)
        rec.update(n=len(x0), sha={"x0": sha(x0), "lo": sha(bounds[0]), "hi": sha(bounds[1]), "fun": sha(r),
                                   "jac": sha(J)}, kw={k: str(v) for k, v in sorted(kw.items())})
        rec["mult_lo"] = [float(v) for v in np.asarray(bounds[0])[-(6 - len(hold) + 3):-3]]
        if job.get("wall"):
            raise lp.WallTime()
        with np.load(Path(job["reference"]) / "localnn_state.npz", allow_pickle=False) as z:
            full = np.array(z["x"], float)
        full[~free] = 1.0
        if job.get("m") is not None:
            full[P.im] = np.where(free[P.im], job["m"], 1.0)
        y = full[free]
        if job.get("derivatives"):
            J = jac(y)
            rec["jac_is_free_columns"] = bool(np.array_equal(J, P.jac(full)[:, free]))
            idx = np.where(free)[0]
            probe = [int(np.argmax(y[:P.n]))] + [i for i, k in enumerate(idx) if P.im.start <= k < P.im.stop]
            worst = 0.0
            for i in probe:
                h = 1e-5 * max(abs(y[i]), 1e-3)
                e = np.zeros(len(y))
                e[i] = h
                num = (-fun(y + 2 * e) + 8 * fun(y + e) - 8 * fun(y - e) + fun(y - 2 * e)) / (12 * h)
                den = float(np.max(np.abs(J[:, i])))
                worst = max(worst, float(np.max(np.abs(num - J[:, i]))) / den)
            rec["fd_max_rel"] = worst
        return types.SimpleNamespace(x=np.array(y), message="TEST STUB: no optimizer ran", status=0, njev=0,
                                     optimality=0.0, cost=0.5 * float(np.sum(fun(y) ** 2)))

    lf.least_squares = stub
    out = io.StringIO()
    try:
        with contextlib.redirect_stdout(out):
            rec["rc"] = lf.main(job["argv"])
    except SystemExit as e:
        rec["refused"] = str(e)
    final = Path(job["out"]) / "final"
    if final.exists() and not job.get("guard"):
        with np.load(final / "localnn_state.npz", allow_pickle=False) as z:
            meta = json.loads(str(z["meta_json"]))
        rec.update(costs=meta["costs"], m=meta["multipliers"], bounds_active=meta["active_bounds"],
                   slope_test=meta.get("slope_test"))
        rec["header"] = [x for x in (final / "fit_result.txt").read_text().split("\n")
                         if x.startswith(("# A_e", "c_slope"))]
    print("RESULT " + json.dumps(rec, default=float))


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "--worker":
        return worker(json.loads(Path(sys.argv[2]).read_text()))
    i = sys.argv.index("--")
    common = sys.argv[i + 1:]
    ap = argparse.ArgumentParser()
    for k in ("code", "old-code", "reference", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args(sys.argv[1:i])
    O = Path(a.out)
    O.mkdir(parents=True)
    rep, fails = {}, []

    def run(name, code=None, extra=(), out=None, **kw):
        out = out or str(O / name)
        job = dict(kw, code=code or a.code, reference=a.reference, out=out,
                   argv=["fit"] + common + ["--out", out] + list(extra))
        jp = O / ("job_%s.json" % name)
        jp.write_text(json.dumps(job))
        p = subprocess.run([sys.executable, "-W", "ignore", __file__, "--worker", str(jp)], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, universal_newlines=True)
        got = [x[7:] for x in p.stdout.split("\n") if x.startswith("RESULT ")]
        rep[name] = json.loads(got[-1]) if got else {"crash": p.stdout[-3000:]}
        print(name, json.dumps({k: rep[name].get(k) for k in ("rc", "refused", "n", "fd_max_rel", "crash")}), flush=True)
        return rep[name]

    def check(what, ok):
        print("  %s %s" % ("PASS" if ok else "FAIL", what), flush=True)
        if not ok:
            fails.append(what)

    new, old = run("default_new"), run("default_old", code=a.old_code)
    check("default: identical x0, bounds, residuals, Jacobian", bool(new.get("sha")) and new["sha"] == old.get("sha")
          and new["kw"] == old["kw"])
    check("default: identical costs", bool(new.get("costs")) and new["costs"] == old.get("costs"))
    ref = json.loads((Path(a.reference) / "fit_manifest.json").read_text())["costs"]
    check("default: the reference fit's costs", new.get("costs") == ref)
    nfull = new.get("n", 0)

    variants = [("none", list(BINS))] + [("only_" + b[6:], [x for x in BINS if x != b]) for b in BINS]
    for name, hold in variants:
        extra = ["--slope-test", "--hold-mult", ",".join(hold)]
        if name == "none":
            extra.append("--check-jacobian")
        r = run(name, extra=extra, hold=hold, m=2.0, derivatives=True)
        check(name + ": free parameters and bounds", r.get("n") == nfull - len(hold)
              and r.get("mult_lo") == [1.0] * (6 - len(hold)))
        check(name + ": Jacobian", r.get("jac_is_free_columns") is True and r.get("fd_max_rel", 1) < 1e-6)
        m = dict(zip(BINS, (r.get("m") or [0] * 7)[1:]))
        check(name + ": products", r.get("rc") == 0 and all(m[b] == (1.0 if b in hold else 2.0) for b in BINS)
              and r["bounds_active"]["multipliers_at_lower_bound"] == [] and r["slope_test"]["held_at_1"] == hold
              and sum(x.endswith(" held") for x in r["header"]) == len(hold))

    r = run("free_min_TEST", extra=["--slope-test", "--slope-mult-min", "0.5", "--check-jacobian"], m=0.7,
            derivatives=True)
    check("lowered bound: bounds", r.get("n") == nfull and r.get("mult_lo") == [0.5] * 6)
    check("lowered bound: Jacobian", r.get("jac_is_free_columns") is True and r.get("fd_max_rel", 1) < 1e-6)
    check("lowered bound: products", r.get("rc") == 0 and r.get("m") == [1.0] + [0.7] * 6)

    out = str(O / "resume")
    h = ",".join(BINS)
    r = run("resume_wall", extra=["--slope-test", "--hold-mult", h], out=out, hold=BINS, wall=True)
    check("wall-time stop: exit 75", r.get("rc") == 75)
    r = run("resume_other", extra=["--slope-test", "--hold-mult", BINS[0]], out=out, guard=True)
    check("resume with other settings refused", "resume" in r.get("refused", ""))
    r = run("resume_plain", out=out, guard=True)
    check("resume without --slope-test refused", "resume" in r.get("refused", ""))
    r = run("resume_same", extra=["--slope-test", "--hold-mult", h], out=out, hold=BINS)
    check("resume with the same settings completes", r.get("rc") == 0 and r.get("m") == [1.0] * 7)

    for name, extra in (("hold_without_test", ["--hold-mult", BINS[0]]),
                        ("min_without_test", ["--slope-mult-min", "0.5"]),
                        ("min_zero", ["--slope-test", "--slope-mult-min", "0"]),
                        ("min_above_1", ["--slope-test", "--slope-mult-min", "2"]),
                        ("hold_and_min", ["--slope-test", "--hold-mult", BINS[0], "--slope-mult-min", "0.5"]),
                        ("gpu", ["--slope-test", "--device", "gpu"]),
                        ("unknown_bin", ["--slope-test", "--hold-mult", "slope_00_40"])):
        r = run("guard_" + name, extra=extra, guard=True)
        check("guard " + name, bool(r.get("refused")) and "sha" not in r)

    (O / "test_localnn_slopetest.json").write_text(json.dumps({"failures": fails, "cases": rep}, indent=1) + "\n")
    print("ALL PASSED" if not fails else "FAILED: %s" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
