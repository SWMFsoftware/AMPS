#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/prepare_44m.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/prepare_44m.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8c2d327e4be974b562589308b2385ed2d954d5c3a21db4def4e6e45d9b83cff8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: SHAP8 44M (repaired) preparation: runner prepare group (base, store, visibility, gravity with
#   --allow-genus 6, cache, bundles cliff >=40/60), the sin(slope) templates. Jobs 121996 (OOM in gravity), 122020.
# RUN: sbatch [--export=ALL,PREP_GPUS=6] prepare_44m.sbatch
# ==== end of header ====
# SHAP8 44M (repaired) preparation (user decisions 2026-11-04: run 44M; gravity on genus 6 admitted as an
# experiment (--allow-genus 6); cliff templates >=40 and >=60 deg and sin(slope) built; grain fits to follow).
# 1. runner prepare group (cops_pipeline.py run --group prepare; the runner's fits are not run): base (native GPU
#    base build), record store (mixed), visibility, gravity (mixed, all GPUs), gravity cache, template bundle
#    (cliff >= 45, 60 deg: the runner's default bundle); same campaign settings as 12M/20M;
# 2. an extra bundle with cliff >= 40 and >= 60 deg (c05_write_bundle.py --cliff-threshold 40 --cliff-threshold 60);
# 3. the C16 slope bundle (c16_slope_template.py: slope_sin, late_ng_not_rg), as jobs/c16.sbatch builds it.
# Resubmitting resumes: the runner skips finished stages; steps 2 and 3 are skipped when their bundle exists.
#SBATCH --job-name=cops44m_prepare
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops44m_prepare-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_044m_flat_earlier
SHAPE=$COPS_DATA_ROOT/shape_models/SHAP8_MIARD/cg-spg_mspcd-shap8-v0.4cmb-44Mfacets_repaired.bdf
PY="$APP/.venv/bin/python -W ignore"
P=$R/provenance/prepare-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
sha256sum "$APP/native/liouville_cops_cuda" "$AMPS/srcCG/LiouvilleCOPS/spice_angular_velocity" > "$P/executables.sha256"
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
free -g > "$P/free_before.txt"
# the same campaign settings as runs/shap8_020m_flat_earlier (its pipeline.json submit argv), genus 6, 44M shape
OPTS=(--shape "$SHAPE" --surface-bc flat --max-nfev 2000 --allow-genus 6 --fits B29,N30,C30,NC31 --precision mixed
      "--fit-args=--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
      --output-root "$COPS_OUTPUT_ROOT" --fine-tune --ft-minutes 5 --ft-intervals 5 --ft-method cpu --ft-max-nfev 2000)
# job 121996: the 44M gravity (one worker per GPU, about 30 GB host memory each, parent 23 GB) on 8 GPUs exceeded the
# 251 GB node: one worker OOM-killed at 18:37 CST, the parent waited forever. PREP_GPUS (default 6) limits the
# runner's GPU count (all GPU stages of this resume; base/store/visibility were complete).
OPTS+=(--gpus-per-node "${PREP_GPUS:-6}")
echo "$(date -u +%FT%TZ) runner prepare: start (gpus ${PREP_GPUS:-6})"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY tools/cluster/cops_pipeline.py run --restart-dir "$R" \
    --group prepare --wall-time slurm "${OPTS[@]}" > "$P/runner_prepare.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) runner prepare: exit $rc"; tail -15 "$P/runner_prepare.log"
[ $rc -eq 0 ] || exit $rc
B40=$R/bundle_cliff40/restart_templates_cliff40_60.npz
if [ ! -e "$B40" ]; then
    echo "$(date -u +%FT%TZ) cliff >=40 bundle: start"
    $PY tools/c05_write_bundle.py --run "$R/base" --gravity-cache "$R/cache/admitted_gravity_cache.npz" \
        --gravity-shape "$SHAPE" --allow-genus 6 --cliff-threshold 40 --cliff-threshold 60 --out "$B40" --device 0 \
        > "$P/bundle_cliff40.log" 2>&1; rc=$?
    echo "$(date -u +%FT%TZ) cliff >=40 bundle: exit $rc"; tail -3 "$P/bundle_cliff40.log" | cut -c1-300
    [ $rc -eq 0 ] || exit $rc
fi
S=$R/c16/bundle/slope_templates.npz
if [ ! -e "$S" ]; then
    mkdir -p "$(dirname "$S")"
    echo "$(date -u +%FT%TZ) slope bundle: start"
    $PY tools/c16_slope_template.py --run "$R/base" --bundle "$R/bundle/restart_templates.npz" --out "$S" \
        > "$P/slope_template.log" 2>&1; rc=$?
    echo "$(date -u +%FT%TZ) slope bundle: exit $rc"; tail -3 "$P/slope_template.log" | cut -c1-300
    [ $rc -eq 0 ] || exit $rc
fi
free -g > "$P/free_after.txt"
$PY tools/cluster/cops_pipeline.py status --restart-dir "$R" > "$P/status_after.txt" 2>&1
echo "$(date -u +%FT%TZ) 44M preparation complete"
