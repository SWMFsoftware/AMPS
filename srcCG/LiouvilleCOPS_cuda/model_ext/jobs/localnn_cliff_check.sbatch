#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_cliff_check.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_cliff_check.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 56340bbc132070fc13456f7ddf9d70503e38c76ffbf44f78bf76af6c42098448
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: minimum test, table and verdict of the cliff test
# RUN: through jobs/localnn_cliff_submit.sh
# ==== end of header ====
# End checks of the cliff test (user 2026-11-22; PLAN_cliff.md): minimum test of every finished fit in its own problem,
# the table and the reading rules written before the evaluations. Evaluations only.
# Env: RUNROOT, TOOLS, CONTROLS_C (or BASELINE, test only).
#SBATCH --job-name=cops_cliff_check
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_cliff_check-%j.out
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
RUNROOT=${RUNROOT:?RUNROOT is required}
TOOLS=${TOOLS:?TOOLS is required}
PY=$APP/.venv/bin/python
(cd "$TOOLS" && sha256sum -c --quiet SHA256SUMS)
P=$RUNROOT/provenance-check-$SLURM_JOB_ID
mkdir "$P"
cp -p "$0" "$P/localnn_cliff_check.sbatch.used"
cp -p "$TOOLS/SHA256SUMS" "$P/TOOLS_SHA256SUMS"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
cd "$APP"
read -r -a S <<< "$("$PY" "$TOOLS/localnn_cliff.py" states --runroot "$RUNROOT")"
[ "${#S[@]}" -gt 0 ] && bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy \
  "$PY" -W ignore "$TOOLS/q002_stationarity_check.py" --solver nnls --own-problem --out "$RUNROOT/stationarity_$SLURM_JOB_ID" "${S[@]}"
if [ -z "${BASELINE:-}" ]; then BASELINE=$("$PY" "$TOOLS/localnn_cliff.py" baseline --controls-c "${CONTROLS_C:?CONTROLS_C is required}"); fi
"$PY" -W ignore "$TOOLS/localnn_cliff.py" summary --runroot "$RUNROOT" --baseline "$BASELINE" --code "$TOOLS/code" --out "$RUNROOT/report_$SLURM_JOB_ID"
