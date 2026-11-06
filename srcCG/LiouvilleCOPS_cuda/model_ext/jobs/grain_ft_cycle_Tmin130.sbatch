#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/grain_ft_cycle_Tmin130.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/grain_ft_cycle_Tmin130.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: f16c059709fe1755981c5c278d2da8cd37ddac879b64f6c26a7497b3a2db5185
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: as grain_ft_cycle.sbatch, cycles 4-6 from cycle 3 with Tmin = 130 K and Tmax = 231.982 K held, T75 free in
#   30-400 K. Jobs 122117, 122140, 122141.
# RUN: sbatch --export=ALL,CYCLE=4 grain_ft_cycle_Tmin130.sbatch; CYCLE=5, 6 chained with --dependency=afterok
# ==== end of header ====
# 12M fit / fine-tuning cycles 4-6 with Tmin held at 130 K and Tmax held at the original best fit's value (user
# request 2026-11-05: "130 should be the minimum temperature, not a temperature everywhere; set the max temperature
# from the max temperature of the original model"; answers: Tmax = 231.9820451565314 K, the Tmax of
# grain/fits/C30_n2000_grain_cross10 (held); Tmin = 130 K held; T75 free within --tbounds 30 400; the 130 K
# all-temperature cycles (grain_ft_cycle_T130.sbatch) cancelled after cycle 4).
# Temperature law T = Tlit if lit and Tlit > Tmin, else Tmin: Tmin is the floor everywhere.
# Same scheme as jobs/grain_ft_cycle.sbatch (fine-tuning correction as a template held at 1, corrections accumulated;
# grain start = the fine-tuned grain), continuing its numbering: cycle 4 from cycles/c3 (the last best fit) and its
# fine tuning, cycle k from cycles_Tmin130/c<k-1>; ft_c1..ft_c<k> held at 1; max_nfev 4000; then the fine tuning
# (5 min, 5 intervals, cpu). Products: runs/shap8_012m_flat_earlier/grain/cycles_Tmin130/c<k>/.
#   sbatch --export=ALL,CYCLE=<4|5|6> grain_ft_cycle_Tmin130.sbatch      (chained with --dependency=afterok)
# CYCLE_TEST=1: cycle 4 with max_nfev 4, into <run>/grain/cycles_Tmin130_test-<jobid>/.
#SBATCH --job-name=cops12m_ftcycTmin130
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:30:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_ftcycTmin130-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
K=${CYCLE:?4, 5 or 6}; NFEV=4000; TEST=${CYCLE_TEST:-0}; TLO=30; TMIN=130; TMAX=231.9820451565314
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
CY0=$R/grain/cycles; CY=$R/grain/cycles_Tmin130
[ "$TEST" = 1 ] && { CY=$R/grain/cycles_Tmin130_test-$SLURM_JOB_ID; NFEV=4; K=4; }
case $K in 4|5|6) ;; *) echo "CYCLE must be 4, 5 or 6" >&2; exit 2;; esac
dir_of() { if [ "$1" -le 3 ]; then echo "$CY0/c$1"; else echo "$CY/c$1"; fi; }
name_of() { if [ "$1" -le 3 ]; then echo "C30_n4000_grain_cross10_ftc$1"; else echo "C30_n4000_Tmin130_grain_cross10_ftc$1"; fi; }
D=$CY/c$K; NAME=C30_n${NFEV}_Tmin130_grain_cross10_ftc$K
J=$((K - 1)); PD=$(dir_of $J); PN=$(name_of $J); PF=$PD/$PN; PFT=$PD/finetune/${PN}_m5_i5_cpu
O=$D/$NAME; F=$D/finetune/${NAME}_m5_i5_cpu; P=$D/provenance-$SLURM_JOB_ID
for f in "$PF/fit_restart.npz" "$PFT/fine_tuning_state.npz" "$PFT/grain_finetune.json"; do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
[ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
mkdir -p "$P" "$D/finetune"; cd "$APP" || exit 1
cp -p "$0" "$P/grain_ft_cycle_Tmin130.sbatch"
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
FIX=(--fix Tmin=$TMIN --fix Tmax=$TMAX); TA=()
for j in $(seq 1 "$K"); do
    FIX+=(--fix c_ft_c$j=1); TA+=(--activity-template "$(dir_of $j)/template_ft_c$j/ft_c$j.templates.npz" ft_c$j)
done
echo "$(date -u +%FT%TZ) cycle $K fit $NAME: start (restart $RS; held ${FIX[*]})"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/grain_fit_fixed.py" "${FIX[@]}" --crossing-weight 10 \
    --check-jacobian --devices "$DEV" --precision mixed --record-store "$R/store" --max-nfev "$NFEV" --wall-time slurm \
    --restart-from "$RS" --template-start 0 -- "$R/base" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds $TLO 400 --ftol 1e-10 \
    --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O" \
    --activity-template "$R/bundle/restart_templates.npz" cliff_ge_60deg "${TA[@]}" > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) cycle $K fit: exit $rc"; tail -6 "$P/fit.log" | cut -c1-250; [ $rc -eq 0 ] || exit $rc
export COPS_CUDA_RECORD_STORE=$R/store
echo "$(date -u +%FT%TZ) cycle $K fine tuning: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/fit_restart.npz" --minutes 5 --intervals 5 --method cpu \
    --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) cycle $K fine tuning: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
