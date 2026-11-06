#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/grain_rough.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/grain_rough.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: afc0be21d4a5d6ff3f49885dce0ca2d93312785a796cdb8d996e3e3ea86b7b17
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the best 12M grain fit's exact command on the 12M rough-surface base (prepare_12m_rough.sbatch), then its
#   fine tuning. Job 122034 (worse than flat).
# RUN: sbatch grain_rough.sbatch (after prepare_12m_rough.sbatch)
# ==== end of header ====
# The best 12M grain fit (C30_n2000_grain_cross10: cliff >= 60 deg template, free temperatures, dust grain in NG,
# crossing weight 10) redone on the 12M "rough" base (emission uniform over the hemisphere), then its fine tuning
# (user request 2026-11-05: "yes, run the best grain fit on the rough base").
# The exact command recorded in runs/shap8_012m_flat_earlier/grain/fits/C30_n2000_grain_cross10/fit_manifest.json
# (default start), with only the base, the record store and --out changed; the 12M template bundle is used unchanged
# (same mesh; load_bundle checks the mesh identity). Fine tuning: grain_finetune.py, 5 min, 5 intervals, cpu, mixed.
#SBATCH --job-name=cops12m_rough_fit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_rough_fit-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R=$COPS_OUTPUT_ROOT/runs/shap8_012m_rough_earlier
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
FIT=C30_n2000_grain_cross10; O=$R/grain/fits/$FIT; F=$R/grain/finetune/${FIT}_m5_i5_cpu; P=$R/grain/provenance-$SLURM_JOB_ID
grep -q "^SurfaceBC *= *rough" "$R/base/input.snapshot.in" || { echo "the base at $R is not a rough base" >&2; exit 2; }
[ -e "$O" ] && { echo "$O exists" >&2; exit 2; }
mkdir -p "$P" "$(dirname "$O")" "$(dirname "$F")"; cd "$APP" || exit 1
cp -p "$0" "$P/grain_rough.sbatch"; cp -p "$X/grain_crossing.py" "$X/grain_fit.py" "$X/grain_finetune.py" "$P/"
sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
echo "$(date -u +%FT%TZ) fit $FIT on the rough base: start"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/grain_fit.py" --crossing-weight 10 --check-jacobian \
    --devices "$DEV" --precision mixed --record-store "$R/store" --max-nfev 2000 --wall-time slurm -- "$R/base" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 \
    --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$O" \
    --activity-template "$R12/bundle/restart_templates.npz" cliff_ge_60deg > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit: exit $rc"; tail -6 "$P/fit.log" | cut -c1-250; [ $rc -eq 0 ] || exit $rc
export COPS_CUDA_RECORD_STORE=$R/store
echo "$(date -u +%FT%TZ) fine tuning: start"
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/fit_restart.npz" --minutes 5 --intervals 5 --method cpu \
    --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"; tail -4 "$P/finetune.log" | cut -c1-250
exit $rc
