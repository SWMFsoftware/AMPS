#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/gpu_solver/gpu_verify.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/gpu_solver/gpu_verify.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: ff8df7796dd960b9b37c5672c4297050555ecd366f552fae1d70c022cbd4c125
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (1 GPU, not gpu002 / gpu003): localnn_gpu_solver.py verify (fixed state, one TRF step,
#   50-evaluation trajectory per SVD mode, speed); MODEL_EXT = the frozen GPU snapshot.
# RUN 2026-11-15: 132810 (MODEL_EXT=localnn_s2gpu_20261115T192024Z); 132048 cancelled (pending, run copy).
# ==== end of header ====
# GPU solver path of option C (lead 2026-11-15 'option 3'): verification on ONE GPU (localnn_gpu_solver.py verify):
# fixed-state residuals / Jacobian CPU vs GPU, one TRF step (SVD modes qr and gesvd) CPU vs GPU, the 50-evaluation
# trajectory per mode vs the CPU reference cpu_ref_50 (job cpu_ref.sh), the seconds per evaluation.
#SBATCH --job-name=lnn_gpu_verify
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/localnn_20261114/gpu_solver/gpu_verify-%j.out
set -uo pipefail
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=/data001/highresolu_vt/cops_output/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/gpu_solver
export PYTHONPYCACHEPREFIX=$W/pycache_$SLURM_JOB_ID
MODEL_EXT=${MODEL_EXT:-/data001/highresolu_vt/cops_output/model_ext}
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; PY=$APP/.venv/bin/python; cd $APP
P=$W/provenance-$SLURM_JOB_ID; mkdir -p $P
cp -p "$0" $P/gpu_verify.sbatch.used; cp -p $MODEL_EXT/localnn_gpu_solver.py $MODEL_EXT/localnn_problem.py $MODEL_EXT/localnn_fit.py $P/
( cd $P && sha256sum *.py > code_sha256.txt )
{ git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C /data001/highresolu_vt/AMPS status --short; } > $P/git.txt 2>&1
{ type module >/dev/null 2>&1 && module list; } > $P/module_list.txt 2>&1 || echo "module command not available (no modules loaded)" >> $P/module_list.txt
nvidia-smi -q > $P/nvidia_smi_q.txt 2>&1
{ nvcc --version; gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > $P/toolchain.txt 2>&1
~/.local/bin/uv pip list --python $PY > $P/pip_list.txt 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|SLURM_JOB_GPUS|SLURM_JOB_NODELIST|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS)=' | sort > $P/environment.txt
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname), GPU $CUDA_VISIBLE_DEVICES: $(nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader -i 0)"
$PY -W ignore $MODEL_EXT/localnn_gpu_solver.py verify --case $W/case_c2_f10_pm0.01.json --cpu-ref $W/cpu_ref_50 --modes qr,gesvd \
    --out $W/verify_${SLURM_JOB_ID} 2>&1 | grep -v "  eval " | tail -5 | cut -c1-3000; rc=${PIPESTATUS[0]}
echo "$(date -u +%FT%TZ) done rc $rc"
cp -p $W/gpu_verify-$SLURM_JOB_ID.out $P/job.log 2>/dev/null
exit $rc
