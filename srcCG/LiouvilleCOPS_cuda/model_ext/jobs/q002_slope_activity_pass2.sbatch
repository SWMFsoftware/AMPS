#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_slope_activity_pass2.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_slope_activity_pass2.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6c3be545b24877ff2ea11e9bc6127f21f5320f0afbfee1e064af112ecfac61bc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job: slope statistics with the same-observation classes (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
#SBATCH --job-name=cops_slope2
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/logs/%x_%j.log
# Slope / activity statistics with the same-observation visibility classes (after q002_seen_gauges.sbatch), 2026-11-21.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
O=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_slope_activity.py
mkdir -p $O/pass2
{ echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; git -C /data001/highresolu_vt/AMPS rev-parse HEAD; $PY -c "import numpy; print('numpy', numpy.__version__)"; sha256sum $S $O/seen_gauges.npz; } > $O/pass2/provenance_$SLURM_JOB_ID.txt
$PY $S --out $O/pass2 --seen-gauges $O/seen_gauges.npz
