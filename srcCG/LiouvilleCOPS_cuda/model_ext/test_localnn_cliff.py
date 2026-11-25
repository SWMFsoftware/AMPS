#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_localnn_cliff.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_localnn_cliff.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: cbd51273930c28ff838676ff16fd6475fc2daee270124a2396c7f49272fd8f71
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tests of localnn_fit.py --hold-mult BIN=V; no optimizer runs
# RUN: jobs/localnn_cliff_code_test.sbatch
# ==== end of header ====
"""Tests of the held value (--hold-mult BIN=V, cliff test 2026-11-22). No optimizer runs (stub of test_localnn_slopetest).

 1 held at 1 (as follow-up c): new and previous code hand bitwise identical x0, bounds, residuals and Jacobian
 2 guards: V = 0, V > 1, text, --hold-mult without --slope-test, a start with another value

    test_localnn_cliff.py --code NEW --old-code PREVIOUS --reference <fit/segNN> --out DIR -- <the fit's common options>
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BINS = ["slope_40_50", "slope_50_60", "slope_60_70", "slope_70_80", "slope_80_90", "slope_90_180"]


def main():
    i = sys.argv.index("--")
    common = sys.argv[i + 1:]
    ap = argparse.ArgumentParser()
    for k in ("code", "old-code", "reference", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args(sys.argv[1:i])
    O = Path(a.out)
    O.mkdir(parents=True)
    fails = []

    def run(name, code, extra, hold=(), guard=False, com=None):
        job = dict(code=code, reference=a.reference, out=str(O / name), hold=list(hold), guard=guard,
                   argv=["fit"] + (com or common) + ["--out", str(O / name)] + extra)
        jp = O / ("job_%s.json" % name)
        jp.write_text(json.dumps(job))
        p = subprocess.run([sys.executable, "-W", "ignore", str(HERE / "test_localnn_slopetest.py"), "--worker", str(jp)],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        got = [x[7:] for x in p.stdout.split("\n") if x.startswith("RESULT ")]
        return json.loads(got[-1]) if got else {"crash": p.stdout[-2000:]}

    def check(what, ok, r=None):
        print("%s %s %s" % ("PASS" if ok else "FAIL", what, "" if ok else r), flush=True)
        if not ok:
            fails.append(what)

    extra = ["--slope-test", "--hold-mult", ",".join(BINS)]
    new, old = run("held1_new", a.code, extra, BINS), run("held1_previous", a.old_code, extra, BINS)
    check("held at 1: identical to the previous code", bool(new.get("sha")) and new["sha"] == old.get("sha")
          and new.get("costs") == old.get("costs") and new.get("slope_test") == old.get("slope_test"), new)

    for name, ex in (("value_0", ["--slope-test", "--hold-mult", "slope_60_70=0"]),
                     ("value_2", ["--slope-test", "--hold-mult", "slope_60_70=2"]),
                     ("value_text", ["--slope-test", "--hold-mult", "slope_60_70=x"]),
                     ("without_slope_test", ["--hold-mult", "slope_60_70=0.01"])):
        r = run("guard_" + name, a.code, ex, guard=True)
        check("guard " + name, bool(r.get("refused")) and "sha" not in r, r)

    k = common.index("--start-file")
    d = O / "start_other_value"
    d.mkdir()
    with np.load(Path(common[k + 1]) / "start.npz", allow_pickle=False) as z:
        s = {n: np.asarray(z[n]) for n in z.files}
    x = np.array(s["x0"], float)
    x[-9:-3] = 3.0                                   # the six multipliers
    s["x0"] = x
    np.savez(d / "start.npz", **s)
    com = list(common)
    com[k + 1] = str(d)
    r = run("guard_start", a.code, ["--slope-test", "--hold-mult", "slope_60_70=0.01"], guard=True, com=com)
    check("guard start with another value", "neither 1 nor the held value" in r.get("refused", ""), r)
    print("CLIFF TESTS PASSED" if not fails else "FAILED %s" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
