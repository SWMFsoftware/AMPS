#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_top5_figures.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_top5_figures.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 267b4e3c5ccdb5110de88aa3bb8aac0e8e673d5d4637f48542659c8ffc1aa8a2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm CPU job: figures_cluster/scripts/q002_top5_figures.py --out figures_cluster/12M/q002_top5_models
#   (2026-11-21; ran as jobs 137802 and 137803)
# RUN: sbatch jobs/q002_top5_figures.sbatch (the figures/ folder is never overwritten; move it aside first)
# ==== end of header ====
#SBATCH --job-name=cops_top5fig
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/logs/%x_%j.log
# Figures and numbers of the five most plausible new-law fits (figures_cluster/scripts/q002_top5_figures.py), 2026-11-21.
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
O=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
echo "host=$(hostname) start=$(date -Is) job=$SLURM_JOB_ID"
{ echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; git -C /data001/highresolu_vt/AMPS rev-parse HEAD;
  git -C /data001/highresolu_vt/AMPS status --short | grep -v '^??' || true; $PY --version;
  $PY -c "import numpy, matplotlib; print('numpy', numpy.__version__, 'matplotlib', matplotlib.__version__)";
  sha256sum /data001/highresolu_vt/cops_output/figures_cluster/scripts/{q002_top5_figures,model_obs_figures,options_residual_figure}.py; } > $O/provenance_$SLURM_JOB_ID.txt
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
$PY /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_top5_figures.py --out $O
echo "end=$(date -Is)"
