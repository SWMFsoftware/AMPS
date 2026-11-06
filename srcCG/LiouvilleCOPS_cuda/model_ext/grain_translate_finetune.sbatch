#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_translate_finetune.sbatch (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_translate_finetune.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: c08ba22acc63286d51f67939940f5158009af96c2ebb0a3da99387e49111394d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: translate_grain.py (a grain fit -> another base), then the grain fine tuning there (NO_FT=1: translation only).
#   Jobs 121974 (12M -> 20M), 121998 (12M -> 44M).
# RUN: sbatch --export=ALL,SRC_FIT=<fit dir>,COPS_RUN=<target runner dir>,NAME=<new name>[,NO_FT=1]
#        grain_translate_finetune.sbatch
# ==== end of header ====
# Translate a grain fit's restart to another base and run the grain fine tuning there (user request 2026-11-04:
# the 12M fit with T75 = Tmin = 133 K fixed -> the 20M surface model; model_ext/translate_grain.py, then
# grain_finetune.py with the settings of grain_finetune.sbatch: 5 min, 5 intervals, cpu, mixed, max_nfev 2000).
#   sbatch --export=ALL,SRC_FIT=<fit dir>,COPS_RUN=<target runner dir>,NAME=<new name>[,NO_FT=1] grain_translate_finetune.sbatch
# NO_FT=1: translation only (no fine tuning; e.g. the start of fits on the target mesh).
# Products: <run>/grain/translated/<NAME>/ (fit_restart.npz with target predictions, grain_term.json copy),
# <run>/grain/translated/finetune/<NAME>_m5_i5_cpu/, provenance <run>/grain/translated/provenance-<jobid>/.
#SBATCH --job-name=cops_grain_xlate_ft
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops_grain_xlate_ft-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
S=${SRC_FIT:?}; R=${COPS_RUN:?}; N=${NAME:?}
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
G=$R/grain/translated; O=$G/$N; F=$G/finetune/${N}_m5_i5_cpu; P=$G/provenance-$SLURM_JOB_ID
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
mkdir -p "$P" "$G/finetune"; cd "$APP" || exit 1
export COPS_CUDA_RECORD_STORE=$R/store
cp -p "$0" "$X/translate_grain.py" "$X/grain_finetune.py" "$X/grain_crossing.py" "$P/"; sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|SRC_FIT|NAME|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
echo "$(date -u +%FT%TZ) translate: start"
$PY "$X/translate_grain.py" --restart "$S/fit_restart.npz" --to-run "$R/base" --target-bundle "$R/bundle/restart_templates.npz" \
    --out "$O" --devices "$DEV" --precision mixed > "$P/translate.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) translate: exit $rc"; tail -2 "$P/translate.log"
[ $rc -eq 0 ] || exit $rc
[ "${NO_FT:-0}" = 1 ] && exit 0
echo "$(date -u +%FT%TZ) fine tuning: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/fit_restart.npz" --minutes 5 --intervals 5 \
    --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"; tail -4 "$P/finetune.log"
exit $rc
