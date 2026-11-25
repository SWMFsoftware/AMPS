# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/12M/q002_plausibility/spike_check.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/12M/q002_plausibility/spike_check.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 86be6eaa94df344accd08928b1627d48ba3b544c4d9f9b3ec087541a33cc288c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: facet with the largest activity per fit (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
# 2026-11-21: which facet carries the largest H2O activity in each new-law run (spike check for the plausibility analysis)
import json, numpy as np
from pathlib import Path
res = json.load(open("/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_plausibility/metrics.json"))["runs"]
idx = {}
for k, r in res.items():
    with np.load(Path(r["ft_dir"]) / "finetune_state.npz", allow_pickle=False) as z:
        A = np.asarray(z["A_final"], np.float64)
    i = int(np.argmax(A)); s = np.sort(A)[::-1]
    idx[k] = i
    print(f"{k:30s} argmax facet {i:9d}  A_max {A[i]:.3e}  2nd {s[1]:.3e}  10th {s[9]:.3e}  mean {A.mean():.3e}")
common = sorted(set(idx.values()))
print("distinct argmax facets:", common)
for k, r in res.items():
    with np.load(Path(r["ft_dir"]) / "finetune_state.npz", allow_pickle=False) as z:
        A = np.asarray(z["A_final"], np.float64)
    print(f"{k:30s} A at the facets {common}: " + ", ".join(f"{A[c]:.3e}" for c in common))
