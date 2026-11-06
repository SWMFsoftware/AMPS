#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/area_scaling.sbatch (on Voyager: /data001/highresolu_vt/cops_output/model_ext/area_scaling.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8a6e6c36eb84105f601a63f68003ad9dbe29387be86fa59bc29fb06de9ab0aed
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job for area_scaling.py (one GPU node, minutes).
# RUN: sbatch --export=ALL,COPS_RUN=<run>,FIT=<fit dir>,FT=<fine-tuning dir>,OUT=<new dir>[,ARGS="..."] area_scaling.sbatch
# ==== end of header ====
# Scale the fine-tuned activity on the main-contributor areas of chosen points; pressures of every point (user
# request 2026-11-04; model_ext/area_scaling.py). No restart file is written; nothing is refitted.
#   sbatch --export=ALL,COPS_RUN=<run>,FIT=<fit dir>,FT=<fine-tuning dir>,OUT=<new dir>[,ARGS="<more options>"] area_scaling.sbatch
#SBATCH --job-name=cops_area_scaling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops_area_scaling-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=${COPS_RUN:?}; X=$COPS_OUTPUT_ROOT/model_ext; P=${OUT:?}.provenance-$SLURM_JOB_ID
export COPS_CUDA_RECORD_STORE=$R/store
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$X/area_scaling.py" "$P/"; sha256sum "$X/area_scaling.py" > "$P/area_scaling.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
CMD="$APP/.venv/bin/python -W ignore $X/area_scaling.py --fit ${FIT:?} --finetune ${FT:?} --run $R/base --out $OUT --devices $DEV --precision mixed ${ARGS:-}"
echo "$CMD" > "$P/command.txt"; echo "$(date -u +%FT%TZ) start: $CMD"
$CMD > "$P/area_scaling.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) exit $rc"; tail -12 "$P/area_scaling.log"
exit $rc
