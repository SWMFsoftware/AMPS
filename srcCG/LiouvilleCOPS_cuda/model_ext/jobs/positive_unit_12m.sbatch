#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_unit_12m.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_unit_12m.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 535c117d8f521cf95448b162345be1cfa25e972e5c8f7c48cac1e5ee32969f1d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. test_positive_fit.py on the 12M base; user decisions 2026-11-07.
# RUN: sbatch positive_unit_12m.sbatch
# ==== end of header ====
# CPU check of model_ext/positive_problem.py on the 12M mesh (2026-11-07, positive-activity fit preparation): the start
# projection of ln A_orig on [Y, cliff_ge_60deg], emission-weighted (w = area A_orig; user 2026-11-07, second run) with
# the area-weighted first rule and a linear-space fit as diagnostics (rms ln residual area/emission weighted, area and
# emission fractions with factor > 2, the cliff coefficient) and the
# log Tikhonov operator (Lt^T Lt == G, D == area-weighted rms of s - s0), with the elements.bin of
# runs/shap8_012m_flat_earlier/base (same mesh as runs/shap8_012m_lategrid_dz1.5) and restart_templates.npz.
# No GPU; ~10 GB RAM, a few minutes. Products: model_ext/logs/positive_unit_12m-<job>.{json,log}.
#   sbatch positive_unit_12m.sbatch
#SBATCH --job-name=cops12m_posunit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_posunit-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONUNBUFFERED=1
X=$COPS_OUTPUT_ROOT/model_ext; L=$X/logs; mkdir -p "$L"
cp -p "$0" "$L/positive_unit_12m-$SLURM_JOB_ID.sbatch"
sha256sum "$X/positive_problem.py" "$X/test_positive_fit.py" > "$L/positive_unit_12m-$SLURM_JOB_ID.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$L/positive_unit_12m-$SLURM_JOB_ID.git.txt" 2>&1
free -g
cd "$APP" || exit 1
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "$APP/.venv/bin/python" -W ignore "$X/test_positive_fit.py" \
    --12m --out "$L/positive_unit_12m-$SLURM_JOB_ID.json" > "$L/positive_unit_12m-$SLURM_JOB_ID.log" 2>&1; rc=$?
tail -50 "$L/positive_unit_12m-$SLURM_JOB_ID.log"
exit $rc
