#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_degeneracy.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_degeneracy.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: dc919938375af0176f3bdc7bafacfc6d9465b52aa2ecb040b999f957b9eaa5b2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, option C): localnn_degeneracy.py for one run root.
# RUN: 125362 (P_m check, chained).
# 2026-11-15: 2026-11-17: NPY_DISABLE_CPU_FEATURES (as the chain driver; reproducible on AVX-512 nodes).
# RUN 2026-11-15: the chain posts.
# ==== end of header ====
# Option C (local non-negative basis): m_b vs c_j degeneracy diagnostics of ONE finished fit (model_ext/localnn_degeneracy.py:
# effective activity per bin, signal share per bin, node dominance, compensation rho, the global level direction,
# conditioning) -> $RUNROOT/analysis_degeneracy. CPU (normal128cores).
# Settings: RUNROOT (required), MODEL_EXT (the frozen localnn snapshot), GRID_DIR (default grid_150m).
#SBATCH --job-name=cops_localnn_degen
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_degen-%j.out
set -uo pipefail
# numpy without AVX-512 kernels, as the chain driver (2026-11-17: on an AVX-512 node (node028, EPYC 9554) np.exp differs
# in the last bit at 38 of the 271 held facets and the tools' bitwise reproduction of A failed, post 136094)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
GRID_DIR=${GRID_DIR:-$R/nonneg/localnn/grid_150m}
RUNROOT=$(readlink -f "${RUNROOT:?RUNROOT is required}")
PY=$APP/.venv/bin/python
P=$RUNROOT/provenance-$SLURM_JOB_ID-degeneracy; mkdir "$P" || exit 1
cd "$APP" || exit 1
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): RUNROOT=$RUNROOT MODEL_EXT=$MODEL_EXT"
cp -p "$0" "$P/localnn_degeneracy.sbatch.used"; cp -p "$MODEL_EXT"/localnn_*.py "$MODEL_EXT/positive_problem.py" "$P/"
( cd "$P" && sha256sum *.py localnn_degeneracy.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
CMD=("$PY" -W ignore "$MODEL_EXT/localnn_degeneracy.py" --runroot "$RUNROOT" --design "$GRID_DIR/design" --grid "$GRID_DIR" \
     --base "$R/base" --d2 "$R/clip_refit/D2/final" --point-table "$R/late_grid_points.txt" \
     --c9-ft "$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9/ft/final")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands.sh"
"${CMD[@]}" > "$P/degeneracy.log" 2>&1; rc=$?
tail -n 14 "$P/degeneracy.log" | cut -c1-400
echo "$(date -u +%FT%TZ) done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_degen-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
