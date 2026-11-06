#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_44m_case.sbatch (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_44m_case.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2a7fd9ebcc203a3203b4a471f958b64f2abd6c5736c10dae428ef472d50d5fa2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: 44M grain fit cases from the translated 12M fit (T75 = Tmin = 133 K held, Tmax free), CASE=cliff40 or cliff60
#   (each + sin(slope)), max_nfev 2000, then the grain fine tuning. Jobs 121999, 122000.
# RUN: sbatch --export=ALL,COPS_RUN=<44M run>,XL=<translated dir>,CASE=cliff40|cliff60[,CASE_TEST=1] grain_44m_case.sbatch
# ==== end of header ====
# 44M grain fit cases (user decisions 2026-11-04): start from the 12M fit (cliff >=60 deg, T75 = Tmin = 133 K)
# translated to 44M without the fine-tuning correction; separate cases cliff >= 40 deg + sin(slope) and
# cliff >= 60 deg + sin(slope), sin(slope) added with a free amplitude; the grain + crossing constraint (weight 10);
# T75 = Tmin = 133 K held, Tmax free (--tbounds 30 400); mixed; max_nfev 2000; template start 1 (user rule;
# amplitudes the restart defines, cliff >= 60, come from it); then the grain fine tuning (5 min, 5 intervals, cpu).
#   sbatch --export=ALL,COPS_RUN=<44M run>,XL=<translated dir>,CASE=cliff40|cliff60[,CASE_TEST=1] grain_44m_case.sbatch
# CASE_TEST=1: max_nfev 4, no fine tuning, into <run>/grain/cases/test-<jobid>/.
#SBATCH --job-name=cops44m_case
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops44m_case-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=${COPS_RUN:?}; XLD=${XL:?}; C=${CASE:?cliff40 or cliff60}; TEST=${CASE_TEST:-0}
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
S=$R/c16/bundle/slope_templates.npz
case $C in
    cliff40) B=$R/bundle_cliff40/restart_templates_cliff40_60.npz; TPL=cliff_ge_40deg; CODE=C40S;;
    cliff60) B=$R/bundle/restart_templates.npz; TPL=cliff_ge_60deg; CODE=C60S;;
    *) echo "CASE must be cliff40 or cliff60" >&2; exit 2;;
esac
G=$R/grain/cases; O=$G; NFEV=2000; NAME=${CODE}_n2000_grain_cross10_T75_133K_Tmin_133K_from12M
[ "$TEST" = 1 ] && { O=$G/test-$SLURM_JOB_ID; NFEV=4; NAME=${CODE}_n4_grain_cross10_T75_133K_Tmin_133K_from12M; }
P=$G/provenance-$SLURM_JOB_ID; F=$G/finetune/${NAME}_m5_i5_cpu
[ -e "$O/$NAME" ] && { echo "exists: $O/$NAME" >&2; exit 2; }
for f in "$B" "$S" "$XLD/fit_restart.npz" "$R/store"; do [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
mkdir -p "$O" "$P" "$G/finetune"; cd "$APP" || exit 1
export COPS_CUDA_RECORD_STORE=$R/store
cp -p "$0" "$X"/fixed_parameters.py "$X"/grain_fit_fixed.py "$X"/grain_fit.py "$X"/grain_crossing.py "$X"/grain_finetune.py "$P/"
sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|XL|CASE|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
FITARGS="--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
GA="--crossing-weight 10 --check-jacobian --devices $DEV --precision mixed --record-store $R/store --max-nfev $NFEV --wall-time slurm --restart-from $XLD/fit_restart.npz"
echo "$(date -u +%FT%TZ) fit $NAME: start"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/grain_fit_fixed.py" --fix T75=133 --fix Tmin=133 $GA -- "$R/base" $FITARGS \
    --out "$O/$NAME" --activity-template "$B" "$TPL" --activity-template "$S" slope_sin > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit $NAME: exit $rc"; tail -6 "$P/fit.log" | cut -c1-250
[ $rc -eq 0 ] || exit $rc
[ "$TEST" = 1 ] && exit 0
echo "$(date -u +%FT%TZ) fine tuning: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/$NAME/fit_restart.npz" --minutes 5 --intervals 5 \
    --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
