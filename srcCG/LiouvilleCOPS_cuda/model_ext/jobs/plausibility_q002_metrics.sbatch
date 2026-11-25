#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/plausibility_q002_metrics.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/plausibility_q002_metrics.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: a617471ec373f83f6939ad9cab5daa9b93fe1916effa41acdb67d7412525baf0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm CPU job: figures_cluster/scripts/plausibility_q002_metrics.py --out-dir figures_cluster/12M/q002_plausibility
#   (2026-11-21)
# RUN: sbatch jobs/plausibility_q002_metrics.sbatch (the output files are never overwritten; move them aside first)
# ==== end of header ====
#SBATCH --job-name=cops_plaus
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_plausibility/logs/%x_%j.log
# Activity-spread measures of the ten new-law fits (figures_cluster/scripts/plausibility_q002_metrics.py), 2026-11-21.
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
O=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_plausibility
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
echo "host=$(hostname) start=$(date -Is) job=$SLURM_JOB_ID"
free -g | head -2
{ echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; git -C /data001/highresolu_vt/AMPS rev-parse HEAD;
  git -C /data001/highresolu_vt/AMPS status --short | grep -v '^??' || true; $PY --version; $PY -c "import numpy; print('numpy', numpy.__version__)";
  sha256sum /data001/highresolu_vt/cops_output/figures_cluster/scripts/plausibility_q002_metrics.py; } > $O/provenance_$SLURM_JOB_ID.txt
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec '"$PY"' /data001/highresolu_vt/cops_output/figures_cluster/scripts/plausibility_q002_metrics.py --out-dir '"$O"
echo "end=$(date -Is)"
