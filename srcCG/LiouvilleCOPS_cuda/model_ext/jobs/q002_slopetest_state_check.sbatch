#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_slopetest_state_check.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_slopetest_state_check.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d0be2fca1192af4f76e6f6fea01527d9ad739954c575bb2f79cd8e9a55dcb3ab
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job of figures/q002_slopetest_state_check.py
# RUN: sbatch --export=ALL,RUNROOT=<run root>,VARIANT=<folder> jobs/q002_slopetest_state_check.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_slope3c_check
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_slope3c_check-%j.out
# Independent check of a control fit's saved state with the reference's frozen code and direct facet sums
# (figures_cluster/scripts/q002_slopetest_state_check.py); evaluations only, no fit.
#   sbatch --export=ALL,RUNROOT=<run root>,VARIANT=<folder> jobs/q002_slopetest_state_check.sbatch
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
RUNROOT=${RUNROOT:?RUNROOT is required}
VARIANT=${VARIANT:?VARIANT is required}
S=$R/figures_cluster/scripts/q002_slopetest_state_check.py
PY=$APP/.venv/bin/python
O=$RUNROOT/state_check_${VARIANT}_$SLURM_JOB_ID
P=$RUNROOT/provenance-state_check-$SLURM_JOB_ID
mkdir "$P"
cp -p "$0" "$P/q002_slopetest_state_check.sbatch.used"
cp -p "$S" $R/figures_cluster/scripts/q002_slope_activity.py "$P/"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
cd "$APP"
exec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy \
  "$PY" -W ignore "$S" --checkpoint "$RUNROOT/$VARIANT/fit/checkpoint.npz" --out "$O" --procs "$SLURM_CPUS_PER_TASK"
