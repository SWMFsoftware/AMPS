#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_touchdown_maps.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_touchdown_maps.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 01b49fdd5aa7bcef1358ec5be40d2f183b5f8d6caf41d5886c723ab30ed89f72
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job: touchdown-vicinity maps (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
#SBATCH --job-name=cops_tdmaps
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/logs/%x_%j.log
# Touchdown-vicinity source-rate maps with the corrected centre (q002_slope_activity.py --maps-only), 2026-11-21.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
O=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/maps_fixed
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_slope_activity.py
mkdir -p $O
{ echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; git -C /data001/highresolu_vt/AMPS rev-parse HEAD; sha256sum $S; } > $O/provenance_$SLURM_JOB_ID.txt
$PY $S --out $O --maps-only
