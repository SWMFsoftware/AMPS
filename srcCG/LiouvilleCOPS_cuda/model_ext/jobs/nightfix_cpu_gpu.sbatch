#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/nightfix_20261118/jobs/nightfix_cpu_gpu.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/nightfix_20261118/jobs/nightfix_cpu_gpu.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: eb06768e7635688c743dd1f61b0d6ce39c68acef10121a0af2245c8c26388da1
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18). nightfix_cpu_gpu_check.py on gpu002 (1 GPU, 200 GB): CPU vs GPU with the fixed law.
# RUN: sbatch tests/nightfix_20261118/jobs/nightfix_cpu_gpu.sbatch
# ==== end of header ====
#SBATCH --job-name=nightfix_cpu_gpu
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodelist=gpu002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=01:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/nightfix_20261118/jobs/%x-%j.out
# H2O night-side fix (2026-11-18): CPU (fit_cops.Model / Problem) vs GPU (cops_cuda.evaluator) with the fixed law on 16
# used rows of the 12M flat_earlier base at the 12M B29 fit state (VALIDATION.md budgets; nightfix_cpu_gpu_check.py).
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
CODE=/data001/highresolu_vt/cops_output/model_ext_frozen/nightfix_check_20261118T073758Z
W=$COPS_OUTPUT_ROOT/tests/nightfix_20261118/cpu_gpu_check
P=$W/provenance-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cp -p "$0" "$P/nightfix_cpu_gpu.sbatch.used"; cp -p "$CODE/nightfix_cpu_gpu_check.py" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1; module list > "$P/module_list.txt" 2>&1
"$APP/.venv/bin/python" -m pip list > "$P/pip_list.txt" 2>&1
cd "$APP" || exit 1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
CMD=(.venv/bin/python -W ignore "$CODE/nightfix_cpu_gpu_check.py" --base $R/base --restart $R/fits/B29_n2000/seg01/fit_restart.npz --rows 16 --out $W/cpu_gpu_check_B29.json)
printf '%q ' "${CMD[@]}" > "$P/commands.sh"; echo >> "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${CMD[@]}"
