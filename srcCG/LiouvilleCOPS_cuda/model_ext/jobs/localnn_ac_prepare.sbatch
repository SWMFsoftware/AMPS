#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_ac_prepare.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_ac_prepare.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: be4cb624ca4bc1ca3dc7008ca7339500c53874c8a5c903b24f77b7dc807e2c54
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: derived starts + preflight of the follow-up; no fit
# RUN: sbatch --export=ALL,TOOLS=<new frozen folder> <this file>
# ==== end of header ====
# Preparation of the slope-test follow-up 'a' and 'c' (user 2026-11-22), evaluations only, no fit: the derived starts of
# c (TOOLS/starts) and the fitter's Jacobian check at all eleven starts in their own problems. Env: TOOLS.
#SBATCH --job-name=cops_ac_prepare
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_ac_prepare-%j.out
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
TOOLS=${TOOLS:?TOOLS is required}
CODE=$R/model_ext_frozen/localnn_slopetest_20261122T033334Z
PY=$APP/.venv/bin/python
P=$TOOLS/provenance-prepare-$SLURM_JOB_ID
mkdir "$P"
cp -p "$0" "$P/localnn_ac_prepare.sbatch.used"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
cd "$APP"
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
heavy "$PY" -W ignore "$TOOLS/localnn_ac.py" make-starts --code "$CODE" --out "$TOOLS/starts"
heavy "$PY" -W ignore "$TOOLS/localnn_ac.py" preflight --code "$CODE" --starts "$TOOLS/starts"
