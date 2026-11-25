#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_multistart_preflight.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_multistart_preflight.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: aa83920a7fe05c6979c1292faf979c83cf2fae5737f50b0f7bc1f8e2f3222652
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: start-up Jacobian check at the parents' states; no fit
# RUN: sbatch --export=ALL,CONTROLS=<root>,GATE=<frozen script> <this file>
# ==== end of header ====
#SBATCH --job-name=cops_multistart_pre
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_multistart_pre-%j.out
# Preflight of the multi-start refits (user 2026-11-22): the frozen fitter's start-up Jacobian check at the finished
# parent's end state and at the running parents' saved states. Evaluations only, no fit. Env: CONTROLS, GATE.
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
CONTROLS=${CONTROLS:?CONTROLS is required}
GATE=${GATE:?GATE is required}
O=$R/tests/localnn_multistart_preflight_$SLURM_JOB_ID
mkdir "$O"
cp -p "$0" "$O/localnn_multistart_preflight.sbatch.used"
cp -p "$GATE" "$O/"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$O/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $APP/.venv/bin/python --version; } > "$O/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$O/module_list.txt" 2>&1 || echo "no module command" >> "$O/module_list.txt"
$APP/.venv/bin/python -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$O/pip_list.txt"
S=()
for v in only_50_60 only_60_70 only_70_80 only_80_90 free_min only_40_50; do
  if [ -e "$CONTROLS/$v/fit/final" ]; then f=$CONTROLS/$v/fit/final/localnn_state.npz
  elif [ -e "$CONTROLS/$v/fit/checkpoint.npz" ]; then f=$CONTROLS/$v/fit/checkpoint.npz
  else continue; fi
  cp "$f" "$O/state_$v.npz"                       # the running fits replace their checkpoints
  S+=("$O/state_$v.npz")
done
cd "$APP"
exec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy \
  "$APP/.venv/bin/python" -W ignore "$GATE" preflight --code "$R/model_ext_frozen/localnn_q002_20261118T162327Z" "${S[@]}"
