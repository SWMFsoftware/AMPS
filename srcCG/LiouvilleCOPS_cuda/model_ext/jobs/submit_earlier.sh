#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/submit_earlier.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/submit_earlier.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 803735050a8d219e8702647d8d44414554b096addc823f019a35105fa6708d5d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the first 12M "earlier set" campaign submission (B29, N30, C30, NC31; mixed; --tbounds 30 400; fine tuning)
#   through the runner. Superseded by submit_campaign.sh. Login node.
# RUN: submit_earlier.sh
# ==== end of header ====
# Submit the 12M "earlier set" campaign (B29, N30, C30, NC31; mixed; --tbounds 30 400; fine tuning) through the runner.
# The runner's scripts (runs/shap8_012m_flat_earlier/slurm/*.sbatch) are used unchanged. This wrapper only
#  (1) exports the site environment: sbatch's default --export=ALL passes it to prepare.sbatch, and the
#      runner's chained sbatch calls (--export=ALL) pass it on to the fit and finetune jobs and to resubmissions;
#  (2) saves the provenance of the submission next to the products;
#  (3) submits the first job (prepare). The later jobs are submitted by the jobs themselves.
set -euo pipefail
export AMPS=/data001/highresolu_vt/AMPS
export APP=$AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
P=$R/provenance/submit-$(date -u +%Y%m%dT%H%M%SZ); mkdir -p "$P"
cp -p "$0" "$P/"; cp -p "$R"/slurm/*.sbatch "$R/pipeline.json" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
git -C "$AMPS" diff -- srcCG/LiouvilleCOPS/fit_base.in > "$P/fit_base.in.diff" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; mpicc --showme:version; sinfo --version; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
sha256sum "$APP/native/liouville_cops_cuda" "$AMPS/srcCG/LiouvilleCOPS/spice_angular_velocity" > "$P/executables.sha256"
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
cp -p /data001/highresolu_vt/cops_output/build/job-120810/nvidia_smi_q.txt "$P/nvidia_smi_q_gpu002_job120810.txt"   # per-stage hardware: runner stamps
"$APP/.venv/bin/python" -W ignore "$APP/tools/cluster/cops_pipeline.py" status --restart-dir "$R" > "$P/status_before.txt" 2>&1 || true
sbatch "$R/slurm/prepare.sbatch" | tee "$P/sbatch.txt"
