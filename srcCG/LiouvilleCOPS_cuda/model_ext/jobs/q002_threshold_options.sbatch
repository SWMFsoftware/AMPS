#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_threshold_options.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_threshold_options.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: b3bdcb843235a0ab5b69d50bee142ddb6755d66549a725a4dd8ceafb5b925fee
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: switched-off share of the ten new-law fits with the project's existing active-area threshold (2026-11-21)
# RUN: sbatch jobs/q002_threshold_options.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_thresh
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121/logs/%x_%j.log
# Switched-off share of the ten new-law fits with the project's existing threshold (A <= 1 % of the srcCG median AH2O of
# the seen facets; localnn_report.py, nonneg_alt_summary.py) and around it; how much of option B is its imposed floor.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
D=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121
{ echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; git -C /data001/highresolu_vt/AMPS rev-parse HEAD; } > $D/provenance_$SLURM_JOB_ID.txt
cd /data001/highresolu_vt/cops_output/figures_cluster/scripts
/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python - <<'P'
import json
from pathlib import Path
import numpy as np
import q002_slope_activity as Q
D = Path("/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121")
met = json.load(open(Q.PL))["runs"]
n, cen, a = Q.elements(); del cen
mm = np.memmap(Q.BASE / "elements.bin", dtype="<f8", mode="r", offset=12, shape=(n, 34))
orig = np.array(mm[:, 32]); del mm
with np.load(Path(met[Q.A_RUN]["fit_dir"]) / "nonneg_facets.npz") as z:
    seen = z["seen"].astype(bool); b = z["bin"].astype(int); special = (z["TS"] != 0) | (z["held"] != 0)
med = float(np.median(orig[seen]))
T = Q.map_times(); t = [k for k, v in T.items() if v["utc"].endswith("10:33:10")][0]
q, qm, _, _ = Q.qfac(t, T[t]["rh_m"], n)
out = {"srcCG_median_seen": med, "threshold": 0.01 * med, "fits": {}}
sa = a[seen].sum()
print("median AH2O seen %.6e  threshold %.6e  subsolar flux at threshold %.3e" % (med, 0.01 * med, 0.01 * med * qm))
for k in met:
    with np.load(Path(met[k]["ft_dir"]) / "finetune_state.npz") as z:
        A = z["A_final"].astype(float)
    r = {str(f): float(a[seen & (A <= f * med)].sum() / sa) for f in (1e-3, 3e-3, 1e-2, 3e-2, 1e-1)}
    r["median_A_over_orig_median"] = Q.wquant(A[seen], a[seen], [0.5])[0] / med
    if k == "B_q002_reg":
        m = np.array([1, 1, 6.686, 1, 1, 8.099, 13.04])
        row = [x for x in open(str(Q.ROOT / "figures_cluster/12M/model_vs_observations_nonneg_options_stage2_q002/summary.csv"))]
        import csv
        s = {x["key"]: x for x in csv.DictReader(row)}["B_q002"]
        m = np.array([1.0] + [float(s["mult_" + c]) for c in ("40_50", "50_60", "60_70", "70_80", "80_90", "90_180")])
        floor = seen & ~special & (A <= 1.01 * m[b] * 0.01 * orig)
        r["floor_only_seen_area"] = float(a[floor].sum() / sa)
        r["floor_only_share_of_seen_production_103310"] = float((A * q * a)[floor].sum() / (A * q * a)[seen].sum())
    out["fits"][k] = r
    print(k, {x: round(y, 4) if y > 1e-4 else float("%.2g" % y) for x, y in r.items()})
(D / "threshold_options.json").write_text(json.dumps(out, indent=1))
P
