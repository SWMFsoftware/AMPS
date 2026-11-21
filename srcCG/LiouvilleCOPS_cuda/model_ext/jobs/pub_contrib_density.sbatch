#!/bin/bash
# ==== version-controlled copy: header added 2026-11-20, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/pub_contrib_density.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/pub_contrib_density.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2b26a412e17c5395cc94ccfccc613640884ce445db615ec31da9da0ef02da2d7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Job of pub_contrib_density_figures.py with provenance; 2026-11-20 (full sets of the new-law fits): BASE and
#   MCODE overrides (defaults unchanged)
# RUN: sbatch --nodelist=gpu002 --export=ALL (DATA, OUTF [, BASE, MCODE])
# ==== end of header ====
# Contribution-density maps (figures_cluster/scripts/pub_contrib_density_figures.py; 2026-11-14). One GPU (lon/lat face
# images only; short). Environment: DATA (contrib_phase_data output), OUTF (figure folder), MODEL_NAME.
#SBATCH --job-name=cops_contrib_density
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_contrib_density-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; DATA=${DATA:?}; OUTF=${OUTF:?}
BASE=${BASE:-$R/base}   # 2026-11-20: BASE override (default unchanged)
MODEL_NAME=${MODEL_NAME:-the softplus activity model with slope multipliers >= 1, final state after 9 fit and fine-tuning cycles}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; PY="$APP/.venv/bin/python -W ignore"
MCODE=${MCODE:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z/nonneg_maps.py}   # 2026-11-20: MCODE, BASE overrides
P=$OUTF/provenance-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cd "$APP" || exit 1
for v in DATA OUTF MODEL_NAME MCODE BASE; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/pub_contrib_density.sbatch.used"
mkdir "$P/code" && cp -p "$S/pub_contrib_density_figures.py" "$S/pub_maps_figures.py" "$S/pub_maps_v2_figures.py" "$MCODE" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
cmd=($PY "$S/pub_contrib_density_figures.py" --base "$BASE" --track-base "$R/../shap8_012m_flat_earlier/base"
     --bundle "$R/../shap8_012m_flat_earlier/bundle/restart_templates.npz" --data "$DATA" --maps-code "$MCODE" --out-dir "$OUTF"
     --model-name "$MODEL_NAME")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${cmd[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${cmd[@]}" > "$P/pub_contrib_density.log" 2>&1; rc=$?
tail -n 20 "$P/pub_contrib_density.log"
echo "$(date -u +%FT%TZ) done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_contrib_density-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
