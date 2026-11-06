#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_ft20.sbatch (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_ft20.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 80afea0085399f392d98301c9b3dea1f11e24c0eac526f971d2e01814fb0485b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: 20M grain fit from the 12M fine-tuned state (option 2): ft_template.py (12M correction -> 20M template
#   ft12m_late), then grain_fit_fixed.py from the translated fit with c_ft12m_late = 1 held, then fine tuning.
#   Job 121995 (poor result, cost 1.581).
# RUN: sbatch --export=ALL,COPS_RUN=<20M run>,SRC12=<12M fit dir>,FT12=<12M fine-tuning dir>,XL=<translated dir>
#        [,FT20_TEST=1] grain_ft20.sbatch
# ==== end of header ====
# 20M grain fit starting from the 12M fine-tuned state (user request 2026-11-04: "run option 2 on 20M").
# 1. ft_template.py: the 12M fine-tuning activity correction carried to the 20M facets (nearest centroid) as a
#    template "ft12m_late", amplitude held at 1 (fixed added term);
# 2. grain_fit_fixed.py (fixed_parameters.py): the 20M fit from the translated 12M fit's restart
#    (grain/translated/<...>_from_12M/fit_restart.npz) with the settings of the 12M fit: cliff >=60 deg template,
#    T75 = Tmin = 133 K held, Tmax free (--tbounds 30 400), grain + crossing weight 10, mixed, max_nfev 2000;
#    --template-start 0 (all other amplitudes come from the restart; the held one is not in the optimizer vector);
# 3. grain fine tuning of the result (5 min, 5 intervals, cpu), as grain_finetune.sbatch.
#   sbatch --export=ALL,COPS_RUN=<20M run>,SRC12=<12M fit dir>,FT12=<12M fine-tuning dir>,XL=<translated dir>[,FT20_TEST=1] grain_ft20.sbatch
# FT20_TEST=1: max_nfev 4, no fine tuning, into <run>/grain/ft12m/test-<jobid>/.
#SBATCH --job-name=cops_grain_ft20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops_grain_ft20-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=${COPS_RUN:?}; S12=${SRC12:?}; F12=${FT12:?}; XLD=${XL:?}; TEST=${FT20_TEST:-0}
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
G=$R/grain/ft12m; O=$G; NFEV=2000; NAME=C30_n2000_grain_cross10_T75_133K_Tmin_133K_ft12m
[ "$TEST" = 1 ] && { O=$G/test-$SLURM_JOB_ID; NFEV=4; NAME=C30_n4_grain_cross10_T75_133K_Tmin_133K_ft12m; }
P=$G/provenance-$SLURM_JOB_ID; TPL=$O/template_ft12m; F=$G/finetune/${NAME}_m5_i5_cpu
[ -e "$O/$NAME" ] && { echo "exists: $O/$NAME" >&2; exit 2; }
mkdir -p "$O" "$P" "$G/finetune"; cd "$APP" || exit 1
export COPS_CUDA_RECORD_STORE=$R/store
cp -p "$0" "$X"/ft_template.py "$X"/fixed_parameters.py "$X"/grain_fit_fixed.py "$X"/grain_fit.py "$X"/grain_crossing.py \
    "$X"/grain_finetune.py "$P/"; sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|SRC12|FT12|XL|FT20_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
echo "$(date -u +%FT%TZ) template: start"
$PY "$X/ft_template.py" --finetune "$F12" --source-restart "$S12/fit_restart.npz" --target-restart "$XLD/fit_restart.npz" \
    --out "$TPL" > "$P/template.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) template: exit $rc"; [ $rc -eq 0 ] || exit $rc
FITARGS="--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
GA="--crossing-weight 10 --check-jacobian --devices $DEV --precision mixed --record-store $R/store --max-nfev $NFEV --wall-time slurm --restart-from $XLD/fit_restart.npz --template-start 0"
echo "$(date -u +%FT%TZ) fit $NAME: start"
$PY "$X/grain_fit_fixed.py" --fix T75=133 --fix Tmin=133 --fix c_ft12m_late=1 $GA -- "$R/base" $FITARGS --out "$O/$NAME" \
    --activity-template "$R/bundle/restart_templates.npz" cliff_ge_60deg --activity-template "$TPL/ft12m_late.templates.npz" ft12m_late \
    > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit $NAME: exit $rc"; tail -6 "$P/fit.log" | cut -c1-250
[ $rc -eq 0 ] || exit $rc
[ "$TEST" = 1 ] && exit 0
echo "$(date -u +%FT%TZ) fine tuning: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/$NAME/fit_restart.npz" --minutes 5 --intervals 5 \
    --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
