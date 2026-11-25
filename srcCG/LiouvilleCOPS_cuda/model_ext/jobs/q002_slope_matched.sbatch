#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_slope_matched.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_slope_matched.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2ab24a4a65cbe31b572938de091c391df639028e76321ae5ecf13c68c8dbbc4e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job: matched slope statistics + report figures (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
#SBATCH --job-name=cops_matched
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/logs/%x_%j.log
# Matched slope statistics + regenerated comparison figures of 12M_q002_top5_models (2026-11-21).
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
O=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
{ echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; git -C /data001/highresolu_vt/AMPS rev-parse HEAD; $PY -c "import numpy, matplotlib; print('numpy', numpy.__version__, 'matplotlib', matplotlib.__version__)"; sha256sum $S/q002_slope_matched.py $S/q002_slope_activity.py $S/q002_top5_figures.py $S/q002_top5_sources_figures.py; } > $O/provenance_$SLURM_JOB_ID.txt
cd $S
$PY q002_slope_matched.py --dir $O/slope_activity
$PY q002_top5_figures.py --out $O
$PY q002_top5_sources_figures.py --out $O
