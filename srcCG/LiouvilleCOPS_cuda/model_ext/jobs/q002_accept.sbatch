#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_base_20261118/jobs/q002_accept.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_base_20261118/jobs/q002_accept.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: dfbf87d6e00468c844ca5e1e3d822788de7b124cd771905f0bd20024d66600b6
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18, Phase 0 step 3). FP64 base consistency on base_q002 and the old base (nightfix_diag --base-run),
#   gpu_fit (fit_cops.main) with --max-nfev 1 on both: base_q002 must pass the consistency check, the old base fail.
# ==== end of header ====
#SBATCH --job-name=q002_accept
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodelist=gpu002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/q002_base_20261118/logs/%x-%j.out
# Phase 0 step 3 of the H2O night-side fix (2026-11-18): does fit_cops.main's base consistency check pass on base_q002
# with the new law and fail on the old base, and does the GPU path accept base_q002?
#   1 FP64 numbers of the check for all 1551 rows on both bases (nightfix_diag.py --base-run; CPU part of the job)
#   2 gpu_fit (= the unchanged fit_cops.main with the GPU backend) on base_q002: base validation of every file (no audit
#     cache for this path yet), the activity check, the consistency check (--check-tolerance 1e-6), then ONE evaluation
#     (--max-nfev 1) of a plain clipped fit (D2's base arguments, no templates / grain); products in a scratch folder
#   3 the same command on the old base (expected: "the fit model does not reproduce the LiouvilleCOPS pressures")
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
POSUNSET=$(env | grep -o '^COPS_POS_[A-Z_]*' | sort | paste -sd' ' -); [ -n "$POSUNSET" ] && unset $POSUNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
W=$COPS_OUTPUT_ROOT/tests/q002_base_20261118
DIAG=$COPS_OUTPUT_ROOT/model_ext_frozen/nightfix_diag_20261118T074014Z
CONTRIB=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261118T072729Z
PY=$APP/.venv/bin/python
P=$W/provenance/job-$SLURM_JOB_ID; mkdir -p "$P" "$W/accept" || exit 1
stamp() { echo "$(date -u +%FT%TZ) $*"; }
cp -p "$0" "$P/q002_accept.sbatch.used"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1; module list > "$P/module_list.txt" 2>&1; "$PY" -m pip list > "$P/pip_list.txt" 2>&1
cd "$APP" || exit 1
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
stamp "FP64 base check, both bases"
for b in base_q002 base; do
  heavy "$PY" -W ignore "$DIAG/nightfix_diag.py" --base-run --base "$R/$b" --contrib-code "$CONTRIB" \
      --out "$W/accept/base_run_consistency_$b.json" --procs 48 > "$W/accept/base_run_consistency_$b.log" 2>&1
  stamp "$b rc=$?: $(grep -E 'consistency max rel' "$W/accept/base_run_consistency_$b.log" | tr '\n' ' ')"
done
gfit() {  # $1 base name, $2 output name
  local O=$W/accept/gpu_fit_$2
  local CMD=("$PY" -W ignore tools/gpu_fit.py --devices 0,1,2,3,4,5,6,7 --precision mixed --max-nfev 1 --no-figures \
             -- "$R/$1" --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 \
             --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O")
  printf '%q ' "${CMD[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
  stamp "gpu_fit on $1"
  heavy "${CMD[@]}" > "$W/accept/gpu_fit_$2.log" 2>&1; local rc=$?
  stamp "gpu_fit on $1 rc=$rc"
  grep -E 'consistency with LiouvilleCOPS|Error|activity factors reproduced|base validation|validated' "$W/accept/gpu_fit_$2.log" | cut -c1-300
}
gfit base_q002 q002
gfit base old
stamp "done"
exit 0
