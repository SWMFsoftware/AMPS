#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_activity_check.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_activity_check.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: e5072c79295105b016e1cd5d7eec3e0340bec859ae833f5e3ceb0830addb1e81
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. positive_activity_check.py (kernel identity, NumPy reference, Jacobian); user decisions 2026-11-07.
# RUN: sbatch positive_activity_check.sbatch   (job 122298)
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
# Validation of model_ext/positive_activity.py (user decisions 2026-11-07: night/shadow production ratio 0.02 for H2O
# and CO2, positive activity af = exp(act_lin), per-record H2O flux floor F = max(q A, 0.02 qmax A_orig)) on the
# SHAP8 12M fitting base, mixed precision (FP32 record store), 2 GPUs. No fit is run.
#   1. identity: COPS_POS_ENABLE=1 with defaults vs the unpatched evaluator, bitwise (1 GPU host/device activity,
#      2 GPUs through spawned workers, host/device activity);
#   2. physics: GPU pressures vs an independent NumPy implementation (FP64 kernel files and FP32 store records);
#   3. Jacobian: analytic columns vs central finite differences of the GPU pressures.
#   sbatch positive_activity_check.sbatch     (POINTS=<colon-separated point ids> for a smaller subset;
#                                              products analysis/positive_activity_check/<job>/)
# 2026-11-08 (CO2 flux floor, COPS_POS_FLOOR_CO2=1): also the phases oldmod/r002/oldr002/r1e5/co2floor/co2floor02 of
# positive_activity_check.py (old 2026-11-07 module and kernels vs new, bitwise; the CO2 floor vs NumPy and finite
# differences; PTX of the old vs new kernels). MODEL_EXT=<dir> (default $COPS_OUTPUT_ROOT/model_ext): the model_ext
# directory the check runs from (a frozen snapshot); OLD_MODEL_EXT=<dir> (default the c1 production snapshot).
#SBATCH --job-name=cops_poscheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=2
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=01:30:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_poscheck-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
unset COPS_POS_ENABLE COPS_POS_NIGHT_RATIO COPS_POS_ACTIVITY COPS_POS_FLOOR_RATIO COPS_POS_FLOOR_CO2 COPS_POS_FROZEN_CONFIG
unset CHECK_MODULE_DIR
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
OLD_MODEL_EXT=$(readlink -f "${OLD_MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext_frozen/c1_20261108T000327Z}")
X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
D=$COPS_OUTPUT_ROOT/analysis/positive_activity_check/$SLURM_JOB_ID; P=$D/provenance
[ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$P"; cd "$X" || exit 1
cp -p "$0" "$P/positive_activity_check.sbatch"
cp -p "$X/positive_activity.py" "$X/positive_activity_check.py" "$P/"; mkdir -p "$P/kernels"; cp -p "$X"/kernels/*.cu "$P/kernels/"
sha256sum "$OLD_MODEL_EXT/positive_activity.py" "$OLD_MODEL_EXT"/kernels/*.cu > "$P/old_model_ext.sha256"
sha256sum "$X/positive_activity.py" "$X/positive_activity_check.py" "$X"/kernels/*.cu "$APP"/cops_cuda/kernels/pressure.cu \
    "$APP"/cops_cuda/kernels/activity.cu "$APP"/cops_cuda/{evaluator,transport,multi_gpu,cuda_sources,device_store}.py > "$P/sha256.txt"
for k in pressure activity; do diff "$APP/cops_cuda/kernels/$k.cu" "$X/kernels/$k.cu" > "$P/kernels/$k.cu.diff"; done
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z0-9_]*|POINTS|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
free -g > "$P/free_at_start.txt"
nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv > "$P/gpus_at_start.txt" 2>&1
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "model_ext run: $X (MODEL_EXT); old module/kernels: $OLD_MODEL_EXT" | tee "$P/model_ext.txt"
ARGS=(run --out "$D/check" --devices 0,1 --old-model-ext "$OLD_MODEL_EXT")
[ -n "${POINTS:-}" ] && ARGS+=(--points "${POINTS//:/,}")      # POINTS colon-separated (sbatch --export splits commas)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"; grep '^export ' "$0";
  echo "cd $X"; printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/positive_activity_check.py" "${ARGS[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/positive_activity_check.py" "${ARGS[@]}"; rc=$?
echo "$(date -u +%FT%TZ) check: exit $rc"
dmesg 2>/dev/null | grep -i -E "out of memory|oom-kill" | tail -5 > "$P/dmesg_oom.txt" || true
exit $rc
