#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/prepare_12m_rough.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/prepare_12m_rough.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 822f0f5e9c8f9038bd2aa5e63dc8e4f5772dbf4cb3bcfb7bda14ad92f8fc1237
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: runner prepare group (base + store) for 12M with --surface-bc rough (dF/dOmega = F/2pi) into
#   runs/shap8_012m_rough_earlier. Job 122033 (52 min).
# RUN: sbatch prepare_12m_rough.sbatch
# ==== end of header ====
# 12M base with the "rough" surface boundary (emission uniform over the hemisphere, dF/dOmega = F/(2 pi)) instead of
# "flat" (Lambertian) (user request 2026-11-05: "prepare the rough base for 12M on gpu002"; option A of the roughness
# discussion: brackets the effect of sub-facet roughness on the angular emission law).
# Runner prepare group with --fits B29 (base + record store only; the mesh is the 12M one, so the 12M template bundles
# apply unchanged), the 12M campaign settings (runs/shap8_012m_flat_earlier/pipeline.json) except --surface-bc rough.
#SBATCH --job-name=cops12m_rough
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --nodelist=gpu002
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_rough-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_rough_earlier; PY="$APP/.venv/bin/python -W ignore"
SHAPE=$COPS_DATA_ROOT/shape_models/SHAP8_MIARD/cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf
P=$R/provenance/prepare-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
git -C "$AMPS" diff -- srcCG/LiouvilleCOPS/fit_base.in > "$P/fit_base.in.diff" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
sha256sum "$APP/native/liouville_cops_cuda" "$AMPS/srcCG/LiouvilleCOPS/spice_angular_velocity" > "$P/executables.sha256"
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
free -g > "$P/free_before.txt"
OPTS=(--shape "$SHAPE" --surface-bc rough --max-nfev 2000 --allow-genus 4 --fits B29 --precision mixed
      "--fit-args=--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"
      --output-root "$COPS_OUTPUT_ROOT" --fine-tune --ft-minutes 5 --ft-intervals 5 --ft-method cpu --ft-max-nfev 2000)
echo "$(date -u +%FT%TZ) runner prepare (rough base + store): start"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY tools/cluster/cops_pipeline.py run --restart-dir "$R" \
    --group prepare --wall-time slurm "${OPTS[@]}" > "$P/runner_prepare.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) runner prepare: exit $rc"; grep -v "^\[base\] .*point " "$P/runner_prepare.log" | tail -8
grep -n "SurfaceBC" "$R/base/input.snapshot.in" 2>/dev/null
free -g > "$P/free_after.txt"
exit $rc
