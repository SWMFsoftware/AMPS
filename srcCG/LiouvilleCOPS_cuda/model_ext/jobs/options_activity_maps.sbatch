#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/options_activity_maps.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/options_activity_maps.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: cfc91a0ca0f5598d876faf45342dcf7e135d6a9477c1e3960bb7a563080a82c6
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: GPU job of pub_options_activity_maps.py (one GPU, short; provenance in
#   <folder>/provenance/activity_maps-<job>) (2026-11-14)
# RUN: sbatch --export=ALL,OPT_DIR=<comparison folder> options_activity_maps.sbatch (job 125205)
# ==== end of header ====
# Activity maps of the best stage-1 setting of each non-negative option and the bounded softplus c9, side by side
# (figures_cluster/scripts/pub_options_activity_maps.py; 2026-11-14). One GPU (ground track + lon/lat face images; short).
# Environment: OPT_DIR (the comparison folder made by figures_12m_nonneg_options_stage1.sh; writes OPT_DIR/activity_maps).
#SBATCH --job-name=cops_options_maps
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_options_maps-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
# c9 state checks (A_fit reproduced bitwise) need numpy's AVX2 kernels, as nonneg_maps.sbatch (gpu002 is AVX2 anyway)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; OPT_DIR=${OPT_DIR:?}
C9=$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
MCODE=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z/nonneg_maps.py
P=$OPT_DIR/provenance/activity_maps-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cd "$APP" || exit 1
for v in OPT_DIR C9 MCODE; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/options_activity_maps.sbatch.used"
mkdir "$P/code" && cp -p "$S/pub_options_activity_maps.py" "$S/pub_maps_figures.py" "$S/pub_maps_v2_figures.py" "$MCODE" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
free -g > "$P/free_before.txt" 2>&1
cmd=($PY "$S/pub_options_activity_maps.py" --folder "$OPT_DIR" --base "$R/base" --track-base "$R/../shap8_012m_flat_earlier/base"
     --bundle "$R/../shap8_012m_flat_earlier/bundle/restart_templates.npz" --maps-code "$MCODE"
     --c9-ft "$C9/ft/final" --c9-fit "$C9/fit/final" --procs "$SLURM_CPUS_PER_TASK")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/options_activity_maps.log" 2>&1; rc=$?
tail -n 20 "$P/options_activity_maps.log"
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_options_maps-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
