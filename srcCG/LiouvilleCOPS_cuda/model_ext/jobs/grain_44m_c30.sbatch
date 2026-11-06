#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/grain_44m_c30.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/grain_44m_c30.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8626b84e8c96b29c0b0d5dc1ce01284b19ead8fad6aa2684d3c2d6f887390c27
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the best 12M grain fit's exact command (C30_n2000_grain_cross10: default start, free T, --tbounds 30 400,
#   cliff >=60, grain + crossing 10, mixed, max_nfev 2000) on the 44M base, then its fine tuning. Job 122078.
# RUN: sbatch grain_44m_c30.sbatch
# ==== end of header ====
# The best 12M grain fit (C30_n2000_grain_cross10: cliff >= 60 deg template, free temperatures, dust grain in NG,
# crossing weight 10) redone on the 44M base with exactly the 12M conditions (default start, no restart, no fixed
# parameters, no sin(slope)), then its fine tuning (user 2026-11-05: "rerun 44M exactly for the same conditions as
# those of the best fit of 12M").
# The exact command recorded in runs/shap8_012m_flat_earlier/grain/fits/C30_n2000_grain_cross10/fit_manifest.json
# (default start), with only the base, the record store, the template bundle (the 44M cliff_ge_60deg) and --out changed.
# Fine tuning: grain_finetune.py, 5 min, 5 intervals, cpu, mixed.
#SBATCH --job-name=cops44m_c30
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops44m_c30-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_044m_flat_earlier; B=$R/bundle/restart_templates.npz
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
FIT=C30_n2000_grain_cross10; O=$R/grain/fits/$FIT; F=$R/grain/finetune/${FIT}_m5_i5_cpu; P=$R/grain/provenance-$SLURM_JOB_ID
[ -e "$O" ] && { echo "$O exists" >&2; exit 2; }
mkdir -p "$P" "$(dirname "$O")" "$(dirname "$F")"; cd "$APP" || exit 1
cp -p "$0" "$P/grain_44m_c30.sbatch"; cp -p "$X/grain_crossing.py" "$X/grain_fit.py" "$X/grain_finetune.py" "$P/"
sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
echo "$(date -u +%FT%TZ) fit $FIT on 44M: start"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/grain_fit.py" --crossing-weight 10 --check-jacobian \
    --devices "$DEV" --precision mixed --record-store "$R/store" --max-nfev 2000 --wall-time slurm -- "$R/base" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 \
    --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O" \
    --activity-template "$B" cliff_ge_60deg > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit: exit $rc"; tail -6 "$P/fit.log" | cut -c1-250; [ $rc -eq 0 ] || exit $rc
export COPS_CUDA_RECORD_STORE=$R/store
echo "$(date -u +%FT%TZ) fine tuning: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/fit_restart.npz" --minutes 5 --intervals 5 --method cpu \
    --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
