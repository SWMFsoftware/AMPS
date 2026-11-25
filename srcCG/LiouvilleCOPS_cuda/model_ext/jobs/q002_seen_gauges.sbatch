#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_seen_gauges.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_seen_gauges.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 63353b6e02b559f4f6a75ba4bfc9726d1386f4c968a028cc8f2f0e0e66193aba
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job: visibility pass (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
#SBATCH --job-name=cops_seen
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/logs/%x_%j.log
# Per-gauge and concurrent visibility per facet: one pass over the 190-GB kernels of base_q002 (2026-11-21).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5
O=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_seen_gauges.py
echo "host=$(hostname) start=$(date -Is)"; { echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; git -C /data001/highresolu_vt/AMPS rev-parse HEAD; $PY -c "import numpy; print('numpy', numpy.__version__)"; sha256sum $S; } > $O/provenance_seen_$SLURM_JOB_ID.txt
bash -c "echo 1000 > /proc/self/oom_score_adj; exec $PY $S --base $R/base_q002 --regions $R/nonneg/sos/stage2r_q002_r3_20261118T193956Z/contrib_phase_c5/regions.npz --out $O/seen_gauges.npz --procs 16"
echo "end=$(date -Is)"
