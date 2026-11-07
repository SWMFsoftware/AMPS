#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/base_consistency_check.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/base_consistency_check.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1f0675c14c7f01d2d6001bc431acbdb83b0ccd6d2d9bd7763d70afa12156a1ce
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Consistency of a new fitting base with the unchanged GPU forward model (tools/gpu_fit.py, no model_ext,
#   --check-tolerance 1e-6, --data-only); refuses if COPS_POS_* are set (independent review 2026-11-07).
# RUN: sbatch --export=ALL,RUN=<run dir> base_consistency_check.sbatch   (job 122304)
# ==== end of header ====
# Consistency of a new fitting base with the GPU forward model, ORIGINAL physics (independent review 2026-11-07: the
# positive-activity fits switch the fitter's start-up check off with --check-tolerance inf, so the new base
# runs/shap8_012m_lategrid_dz1.5/base must be checked once with the unchanged code). The unchanged tools/gpu_fit.py
# (no model_ext, no COPS_POS_*) evaluates the run's own coefficients and temperatures on every point and compares
# with the C++ pressures of the base's COPS_model.dat (fit_cops.py consistency check, --check-tolerance 1e-6, as in
# every production fit), then writes the fit data and stops (--data-only). Mixed precision, no record store (the
# base's kernels are converted at bind).
#   sbatch --export=ALL,RUN=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5 base_consistency_check.sbatch
#SBATCH --job-name=cops_basecheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=02:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_basecheck-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
RUN=${RUN:?set RUN=<run dir>}; PY="$APP/.venv/bin/python -W ignore"
env | grep -q '^COPS_POS_' && { echo "COPS_POS_* set: this check must run the unchanged physics" >&2; exit 2; }
O=$RUN/checks/consistency-$SLURM_JOB_ID; P=$O.provenance; [ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/base_consistency_check.sbatch"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1; free -g > "$P/free_at_start.txt"
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
CMD=($PY tools/gpu_fit.py --devices "$DEV" --precision mixed -- "$RUN/base" --degree 4 --tbounds 30 400 \
     --check-tolerance 1e-6 --data-only --out "$O")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs; base $RUN/base"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/check.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) exit $rc"; grep -E "consistency with LiouvilleCOPS|activity factors reproduced|Error" "$P/check.log"
exit $rc
