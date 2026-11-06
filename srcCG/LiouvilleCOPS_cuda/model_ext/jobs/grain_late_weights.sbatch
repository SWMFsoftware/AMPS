#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/grain_late_weights.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/grain_late_weights.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: ef079fbf10d6828e6f6da23bce21d2bfc7f750fd018f90152a613c0d6c2fbc01
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: 12M fit from cycle 3 with the last N points weighted (grain_fit_weighted.py), then the fine tuning.
#   Job 122236 (last 5 points x10, NG and RG, 3 GPUs of gpu002).
# RUN: sbatch [--gpus-per-node=3 --mem=80G --time=06:00:00] --export=ALL,WLAST=5,WFACTOR=10,WPROFILE=constant
#        grain_late_weights.sbatch   (WGAUGES: export it in the shell; sbatch splits --export values at commas)
# ==== end of header ====
# 12M fit from the last best fit (cycle 3) with larger weights for the last N data points (user request 2026-11-06:
# "on 6 GPUs of gpu002: from the best 12M fit start a new fit where we weight for the last 5 points will be increased
# from 1 to 10, where 1 corresponds to the current weight assigned to each point ... make the selection of the number
# of points from the end of the dataset and how their importance is increased parameters defined with CLI";
# answers: start = cycle 3 (ft_c1..ft_c3 held at 1, grain start = cycle 3's fine-tuned grain, free T 30-400 K);
# both gauges; W multiplies the point's cost term; constant profile).
# model_ext/grain_fit_weighted.py (late_weights.py + fixed_parameters.py around the unchanged grain_fit.py), the C30
# command as in the cycles (grain + crossing 10, mixed, --tbounds 30 400, degree 4/4, lambda 0, cliff_ge_60deg +
# ft_c1..3), --restart-from, --template-start 0, max_nfev 4000; then grain_finetune.py (5 min, 5 intervals, cpu;
# its own objective, unweighted).
#   sbatch --export=ALL,WLAST=5,WFACTOR=10,WPROFILE=constant grain_late_weights.sbatch   (WGAUGES default ng,rg; to set it, export WGAUGES=ng,rg in the shell and use --export=ALL: sbatch splits --export values at commas)
# TEST=1: max_nfev 4, into <run>/grain/late_weights_test-<jobid>/.
#SBATCH --job-name=cops12m_latew
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=6
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --nodelist=gpu002
#SBATCH --time=03:30:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_latew-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
WLAST=${WLAST:-5}; WFACTOR=${WFACTOR:-10}; WPROFILE=${WPROFILE:-constant}; WGAUGES=${WGAUGES:-ng,rg}
NFEV=4000; TEST=${TEST:-0}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
CY=$R/grain/cycles; PF=$CY/c3/C30_n4000_grain_cross10_ftc3; PFT=$CY/c3/finetune/C30_n4000_grain_cross10_ftc3_m5_i5_cpu
LW=$R/grain/late_weights; [ "$TEST" = 1 ] && { LW=$R/grain/late_weights_test-$SLURM_JOB_ID; NFEV=4; }
G=${WGAUGES//,/}; NAME=C30_n${NFEV}_grain_cross10_c3_last${WLAST}_w${WFACTOR}_${WPROFILE}_${G}
D=$LW/$NAME; O=$D/fit; F=$D/finetune/${NAME}_m5_i5_cpu; P=$D/provenance-$SLURM_JOB_ID
for f in "$PF/fit_restart.npz" "$PFT/grain_finetune.json"; do [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
mkdir -p "$P" "$D/finetune"; cd "$APP" || exit 1
cp -p "$0" "$P/grain_late_weights.sbatch"
cp -p "$X"/late_weights.py "$X"/grain_fit_weighted.py "$X"/restart_grain_start.py "$X"/fixed_parameters.py "$X"/grain_fit.py \
    "$X"/grain_crossing.py "$X"/grain_finetune.py "$P/"; sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|W[A-Z]+|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS)=' | sort > "$P/environment.txt"
free -g > "$P/free_at_start.txt"
# devices as CUDA numbers them inside the job (on a shared node nvidia-smi lists all GPUs; job 122226 failed with them)
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
echo "$(date -u +%FT%TZ) $NG GPUs visible ($DEV); CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
RS=$D/restart_start.npz
if [ ! -e "$RS" ]; then
    $PY "$X/restart_grain_start.py" --restart "$PF/fit_restart.npz" --finetune "$PFT" --out "$RS" > "$P/restart_start.log" 2>&1 \
        || { echo "restart_grain_start failed" >&2; cat "$P/restart_start.log"; exit 3; }
    cat "$P/restart_start.log"
fi
FIX=(); TA=()
for j in 1 2 3; do FIX+=(--fix c_ft_c$j=1); TA+=(--activity-template "$CY/c$j/template_ft_c$j/ft_c$j.templates.npz" ft_c$j); done
echo "$(date -u +%FT%TZ) fit $NAME: start (weights: last $WLAST, factor $WFACTOR, $WPROFILE, $WGAUGES; held ${FIX[*]})"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/grain_fit_weighted.py" \
    --weight-last "$WLAST" --weight-factor "$WFACTOR" --weight-profile "$WPROFILE" --weight-gauges "$WGAUGES" \
    "${FIX[@]}" --crossing-weight 10 \
    --check-jacobian --devices "$DEV" --precision mixed --record-store "$R/store" --max-nfev "$NFEV" --wall-time slurm \
    --restart-from "$RS" --template-start 0 -- "$R/base" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 \
    --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O" \
    --activity-template "$R/bundle/restart_templates.npz" cliff_ge_60deg "${TA[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit: exit $rc"; grep -E '^late_weights:' "$P/fit.log" | cut -c1-300; tail -6 "$P/fit.log" | cut -c1-250
[ $rc -eq 0 ] || exit $rc
export COPS_CUDA_RECORD_STORE=$R/store
echo "$(date -u +%FT%TZ) fine tuning: start"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/fit_restart.npz" \
    --minutes 5 --intervals 5 --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
