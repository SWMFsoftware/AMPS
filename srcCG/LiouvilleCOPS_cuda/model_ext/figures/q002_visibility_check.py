# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/12M/q002_top5_models/slope_activity/shutdown_check/visibility_check.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/shutdown_check/visibility_check.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 329467a04bcaa1394002f98214d03206567476e328bfd85ce9aa99966c73156f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: NG vs RG visibility per phase and at four instants (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
"""2026-11-21: NG vs RG visibility on the current base (option A: A > 0 everywhere, so contribution > 0 = kernel > 0)."""
import sys, numpy as np
d = sys.argv[1]
z = np.load(d + "/phase_means.npz")
for ph in ("pre", "late"):
    n, r = z["NG_" + ph] > 0, z["RG_" + ph] > 0
    print("phase %-5s NG %d RG %d both %d NG-only %d RG-only %d" % (ph, n.sum(), r.sum(), (n & r).sum(), (n & ~r).sum(), (r & ~n).sum()))
for s in ("103310", "103731_282", "103909_282", "103926"):
    z = np.load(d + "/snap_%s.npz" % s)
    n, r = z["ng"] > 0, z["rg"] > 0
    print("instant %-11s row %d NG %d RG %d both %d NG-only %d RG-only %d" % (s, int(z["row"]), n.sum(), r.sum(), (n & r).sum(), (n & ~r).sum(), (r & ~n).sum()))
