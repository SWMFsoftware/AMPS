# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/12M/q002_top5_models/slope_activity/shutdown_check/test_seen_gauges.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/shutdown_check/test_seen_gauges.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2cfe8b515434887b00544187abe4cc322879445803654c5a5457b03fc6b6310d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test of the visibility pass on four kernel files (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
"""Test of q002_seen_gauges.work on four real kernel files against option A's saved snapshots (same rows)."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, "/data001/highresolu_vt/cops_output/figures_cluster/scripts")
import q002_seen_gauges as S
R = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5")
d = R / "nonneg/sos/stage2r_q002_r3_20261118T193956Z/contrib_phase_c5"
nf = 12364290; use = np.ones(1551, bool); S.G.update(base=R / "base_q002", nf=nf, useN=use, useR=use)
ok = True
for s in ("103310", "103731_282", "103909_282", "103926"):
    z = np.load(d / ("snap_%s.npz" % s)); f = z["face"]; n, r = z["ng"] > 0, z["rg"] > 0
    M, nr = S.work([int(z["ipoint"])])
    U = {k: np.unpackbits(v, count=nf).astype(bool) for k, v in M.items()}
    for k, e in (("ng_obs", n), ("rg_obs", r), ("both_obs", n & r)):
        want = np.zeros(nf, bool); want[f[e]] = True
        same = np.array_equal(U[k], want); ok &= same
        print(s, k, int(U[k].sum()), "expected", int(want.sum()), "SAME" if same else "DIFFERENT")
    sub = not (U["ng_lit"] & ~U["ng_obs"]).any() and not (U["both_lit"] & ~U["both_obs"]).any(); ok &= sub
print("ALL PASS" if ok else "FAIL")
