#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_fixT.sbatch (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_fixT.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3e8a27decafeac5c2d6d464e103bcf48106c634f94915580d52d647a3b7e3500
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: 12M cliff >=60 grain fit with T75 = Tmin = 133 K held (grain_fit_fixed_T.py), then its fine tuning. Job 121932.
# RUN: sbatch --export=ALL,COPS_RUN=<runner restart dir>[,FIXT_TEST=1] grain_fixT.sbatch
# ==== end of header ====
# 12M grain fit with T75 and Tmin fixed (user request 2026-11-04: "run the fit with Tmin and T75 fixed at 133 K";
# answers: the cliff >=60 deg free-temperature fit, default start). Everything else is the setting of
# grain/fits/C30_n2000_grain_cross10 (grain.sbatch GRAIN_SET=free): grain A, B (start A = observed jump, B = 0),
# crossing weight 10, mixed precision, --tbounds 30 400 (Tmax free in it), max_nfev 2000, default start, the run's
# base, store and templates. Driver model_ext/grain_fit_fixed_T.py (+ fixed_temperature.py); git tree unchanged.
#   sbatch --export=ALL,COPS_RUN=<runner restart dir>[,FIXT_TEST=1] grain_fixT.sbatch
# FIXT_TEST=1: max_nfev 4, no fine tuning, into <run>/grain/fixed_T/test-<jobid>/.
# Products: <run>/grain/fixed_T/<fit>/, fine tuning <run>/grain/fixed_T/finetune/<fit>_m5_i5_cpu/ (as
# grain_finetune.sbatch), logs and provenance in <run>/grain/fixed_T/provenance-<jobid>/.
#SBATCH --job-name=cops_grain_fixT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops_grain_fixT-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=${COPS_RUN:?set COPS_RUN}; TEST=${FIXT_TEST:-0}
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
G=$R/grain/fixed_T; O=$G; P=$G/provenance-$SLURM_JOB_ID
NAME=C30_n2000_grain_cross10_T75_133K_Tmin_133K; NFEV=2000
[ "$TEST" = 1 ] && { O=$G/test-$SLURM_JOB_ID; NAME=C30_n4_grain_cross10_T75_133K_Tmin_133K; NFEV=4; }
F=$G/finetune/${NAME}_m5_i5_cpu
[ -e "$O/$NAME" ] && { echo "grain_fixT: $O/$NAME exists" >&2; exit 2; }
mkdir -p "$O" "$P" "$G/finetune"
cd "$APP" || exit 1
cp -p "$0" "$P/grain_fixT.sbatch"; cp -p "$X"/grain_crossing.py "$X"/grain_fit.py "$X"/grain_fit_fixed_T.py \
    "$X"/fixed_temperature.py "$X"/grain_finetune.py "$P/"; sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|FIXT_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
for s in base store bundle; do [ -e "$R/$s" ] || { echo "grain_fixT: $R/$s missing" >&2; exit 2; }; done
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
B=$R/bundle/restart_templates.npz
FITARGS="--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
GA="--crossing-weight 10 --check-jacobian --devices $DEV --precision mixed --record-store $R/store --max-nfev $NFEV --wall-time slurm"
echo "$(date -u +%FT%TZ) fit $NAME: start"
$PY "$X/grain_fit_fixed_T.py" --fix-T T75=133 --fix-T Tmin=133 $GA -- "$R/base" $FITARGS --out "$O/$NAME" \
    --activity-template "$B" cliff_ge_60deg > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit $NAME: exit $rc"
[ $rc -eq 0 ] || exit $rc
[ "$TEST" = 1 ] && exit 0
echo "$(date -u +%FT%TZ) fine tuning: start"
export COPS_CUDA_RECORD_STORE=$R/store
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/$NAME/fit_restart.npz" --minutes 5 --intervals 5 \
    --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"
exit $rc
