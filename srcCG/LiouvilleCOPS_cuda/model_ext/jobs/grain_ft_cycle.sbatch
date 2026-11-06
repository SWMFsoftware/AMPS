#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/grain_ft_cycle.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/grain_ft_cycle.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: cee0636bd58d80a622553603a087b0335742c69a1e6458668320a546b613264d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: 12M fit/fine-tuning cycle k = 1, 2, 3 from the best fit C30: previous fine-tuning correction -> template
#   ft_c<k> (ft_template.py), held at 1 with ft_c1..ft_c<k-1>; grain start fine-tuned (restart_grain_start.py);
#   the C30 command with --restart-from, --template-start 0, max_nfev 4000 (grain_fit_fixed.py); then the fine tuning.
#   Jobs 122080-122082 (cycle 3: cost 2.7151e-3, the lowest 12M fit).
# RUN: c1=$(sbatch --parsable --export=ALL,CYCLE=1 grain_ft_cycle.sbatch); then CYCLE=2, 3 with --dependency=afterok
#   (CYCLE_TEST=1: 4-evaluation test).
# ==== end of header ====
# 12M fit / fine-tuning cycles from the best 12M grain fit (C30_n2000_grain_cross10) (user request 2026-11-05:
# "redo the best fit 12M: restart after the fine tuning, so 4000 iteration fitting, then fine tuning. repeat the
# latter fitting and fine tuning cycle 3 times"; answers: the fine-tuning correction as a template held at
# amplitude 1, corrections accumulated; grain start = the fine-tuned grain; 3 cycles in total).
# Cycle k (CYCLE=1,2,3), from the previous fit P and its fine tuning PFT (cycle 1: the C30 fit and its fine tuning):
#  1. ft_template.py: PFT's delta_activity -> template ft_c<k> on the same 12M mesh (each facet maps to itself);
#  2. restart_grain_start.py: P's restart with the grain start = PFT's tuned grain A, B;
#  3. grain_fit_fixed.py: the C30 command (grain + crossing 10, mixed, --tbounds 30 400, degree 4/4, lambda 0,
#     free temperatures, cliff_ge_60deg free) + templates ft_c1..ft_c<k>, each held at 1, --restart-from 2.,
#     --template-start 0 (all free parameters come from the restart), max_nfev 4000;
#  4. grain_finetune.py (5 min, 5 intervals, cpu, mixed), as for C30.
#   sbatch --export=ALL,CYCLE=<k> grain_ft_cycle.sbatch      (cycles chained with --dependency=afterok)
# CYCLE_TEST=1: max_nfev 4, everything into <run>/grain/cycles_test-<jobid>/ (cycle 1 only).
#SBATCH --job-name=cops12m_ftcycle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:30:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_ftcycle-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
K=${CYCLE:?1, 2 or 3}; NFEV=4000; TEST=${CYCLE_TEST:-0}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
CY=$R/grain/cycles; [ "$TEST" = 1 ] && { CY=$R/grain/cycles_test-$SLURM_JOB_ID; NFEV=4; K=1; }; D=$CY/c$K; NAME=C30_n${NFEV}_grain_cross10_ftc$K
if [ "$K" = 1 ]; then
    PF=$R/grain/fits/C30_n2000_grain_cross10; PFT=$R/grain/finetune/C30_n2000_grain_cross10_m5_i5_cpu
else
    J=$((K - 1)); PF=$CY/c$J/C30_n${NFEV}_grain_cross10_ftc$J; PFT=$CY/c$J/finetune/C30_n${NFEV}_grain_cross10_ftc${J}_m5_i5_cpu
fi
O=$D/$NAME; F=$D/finetune/${NAME}_m5_i5_cpu; P=$D/provenance-$SLURM_JOB_ID
for f in "$PF/fit_restart.npz" "$PFT/fine_tuning_state.npz" "$PFT/grain_finetune.json"; do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
mkdir -p "$P" "$D/finetune"; cd "$APP" || exit 1
cp -p "$0" "$P/grain_ft_cycle.sbatch"
cp -p "$X"/ft_template.py "$X"/restart_grain_start.py "$X"/fixed_parameters.py "$X"/grain_fit_fixed.py "$X"/grain_fit.py \
    "$X"/grain_crossing.py "$X"/grain_finetune.py "$P/"; sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|CYCLE|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
TPL=$D/template_ft_c$K
if [ ! -e "$TPL/ft_c$K.templates.npz" ]; then
    echo "$(date -u +%FT%TZ) cycle $K template ft_c$K from $PFT: start"
    $PY "$X/ft_template.py" --finetune "$PFT" --source-restart "$PF/fit_restart.npz" --target-restart "$PF/fit_restart.npz" \
        --out "$TPL" --name ft_c$K > "$P/template.log" 2>&1 || { echo "template failed" >&2; cat "$P/template.log"; exit 3; }
    cat "$P/template.log"
fi
RS=$D/restart_start.npz
if [ ! -e "$RS" ]; then
    $PY "$X/restart_grain_start.py" --restart "$PF/fit_restart.npz" --finetune "$PFT" --out "$RS" > "$P/restart_start.log" 2>&1 \
        || { echo "restart_grain_start failed" >&2; cat "$P/restart_start.log"; exit 3; }
    cat "$P/restart_start.log"
fi
FIX=(); TA=()
for j in $(seq 1 "$K"); do
    FIX+=(--fix c_ft_c$j=1); TA+=(--activity-template "$CY/c$j/template_ft_c$j/ft_c$j.templates.npz" ft_c$j)
done
echo "$(date -u +%FT%TZ) cycle $K fit $NAME: start (restart $RS; held ${FIX[*]})"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/grain_fit_fixed.py" "${FIX[@]}" --crossing-weight 10 \
    --check-jacobian --devices "$DEV" --precision mixed --record-store "$R/store" --max-nfev "$NFEV" --wall-time slurm \
    --restart-from "$RS" --template-start 0 -- "$R/base" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 \
    --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O" \
    --activity-template "$R/bundle/restart_templates.npz" cliff_ge_60deg "${TA[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) cycle $K fit: exit $rc"; tail -6 "$P/fit.log" | cut -c1-250; [ $rc -eq 0 ] || exit $rc
export COPS_CUDA_RECORD_STORE=$R/store
echo "$(date -u +%FT%TZ) cycle $K fine tuning: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/fit_restart.npz" --minutes 5 --intervals 5 --method cpu \
    --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) cycle $K fine tuning: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
