#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_ac_preview.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_ac_preview.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4f71b5c4842c884a1059deed70a1830481d96cf9dfcde3c3830057fa216161fa
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job of figures/q002_ac_preview.py; no fit
# RUN: sbatch jobs/q002_ac_preview.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_ac_preview
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_ac_preview-%j.out
# One-step preview of the planned follow-up a, and the fine tuning's leverage for c (2026-11-22; evaluations only, no fit;
# figures_cluster/scripts/q002_ac_preview.py).
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
T=$R/model_ext_frozen/ac_20261122T183034Z
S=$R/figures_cluster/scripts/q002_ac_preview.py
D=$R/figures_cluster/12M/q002_decisions_20261121
PY=$APP/.venv/bin/python
P=$D/provenance_$SLURM_JOB_ID
mkdir "$P"
cp -p "$0" "$P/q002_ac_preview.sbatch.used"
cp -p "$S" "$P/"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
cd "$APP"
exec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy \
  "$PY" -W ignore "$S" --tools "$T" --out "$D/ac_preview_$SLURM_JOB_ID"
