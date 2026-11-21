#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_gpu_eval.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_gpu_eval.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: bb825d9edc69beaa9e8d8d2c1b8303cefb902b6d1c2059f8ca1d4d47e09bfb4a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (gpu002, 8 GPUs, 1 h, option C): GPU standard evaluation of one final state (SEG) -> <run root>/gpu_eval;
#   grain free (2026-11-14 22:2x UTC fix).
# RUN: failed with held grain: 125052, 125123; fixed: 125183 (f10, passed), 125185 (f0.1), 125181 (f1, chained).
# 2026-11-18: 2026-11-18: BASE_DIR / D2_DIR, DESIGN_DIR default of the current law.
# RUN 2026-11-18: planned: the cross-check of the q002 final states.
# ==== end of header ====
# Option C (local non-negative basis, lead 2026-11-14): GPU cross-check of ONE localnn fit segment's final per-facet
# activity through the standard fitter path, writing the standard products (model_ext/localnn_gpu_eval.py docstring):
#   1. localnn_gpu_eval.py prep: template bundle (localnn_a = A_final / A_s) + the exact held values
#   2. clip_refit.py --variant D2 (D2's frozen code model_ext_frozen/clip_20261108T234013Z, unchanged grain_fit_fixed /
#      grain_fit / gpu_fit / fit_cops) with a00..a24 = 0, c_localnn_a = 1, kappa, T held; grain A, B free (D2's code
#      cannot hold them: 125052, 125123 failed); --max-nfev 1 (TRF returns its start, the fitter's grain start) -> $OUT/fit
#   3. localnn_gpu_eval.py compare: the CPU state evaluated with the GPU's grain vs the GPU (data cost, tables)
#      -> $OUT/gpu_check.json
# Settings: SEG (required: <run root>/fit/final or segNN), OUT (default <run root>/gpu_eval), GRID_DIR (default grid_150m),
#   FT_DIR (stage 2, 2026-11-15: a localnn fine tuning of SEG, e.g. <chain>/c5/ft/final: its final state A_final and its
#   grain are evaluated; OUT default <run root>/gpu_eval_ft); the node: override on the sbatch command line
#   (--nodelist=<free node>; never gpu003: #SBATCH --exclude=gpu003 since 2026-11-15 14:00 CST, before: --nodelist=gpu002),
#   MODEL_EXT (the localnn
#   frozen snapshot), CLIP_EXT (default model_ext_frozen/clip_20261108T234013Z, the code of D2).
#SBATCH --job-name=cops_localnn_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_gpu-%j.out
set -uo pipefail
# 2026-11-18 (H2O night law q_min = 0.02 q_max): the products of the current law by default (DESIGN_DIR, CONTRIB, LATE,
# NONNEG_CODE); the loads refuse a design / late kernels of another law (localnn_problem.law_check)
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; BASE=$(readlink -f "${BASE_DIR:-$R/base}"); PTS=$R/late_grid_points.txt   # 2026-11-18: BASE_DIR / D2_DIR
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
CX=$(readlink -f "${CLIP_EXT:-$COPS_OUTPUT_ROOT/model_ext_frozen/clip_20261108T234013Z}")
SEG=$(readlink -f "${SEG:?SEG is required}")
RUNROOT=$(dirname "$(dirname "$SEG")")
FT_DIR=${FT_DIR:-}; FTO=()
if [ -n "$FT_DIR" ]; then FT_DIR=$(readlink -f "$FT_DIR"); FTO=(--ft-dir "$FT_DIR"); OUT=${OUT:-$RUNROOT/gpu_eval_ft}; fi
OUT=${OUT:-$RUNROOT/gpu_eval}
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
[ -e "$OUT" ] && { stamp "$OUT exists: not overwritten"; exit 2; }
mkdir -p "$OUT" || exit 1
P=$OUT/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): SEG=$SEG FT_DIR=$FT_DIR OUT=$OUT MODEL_EXT=$MODEL_EXT CLIP_EXT=$CX"
cp -p "$0" "$P/localnn_gpu_eval.sbatch.used"
cp -p "$MODEL_EXT/localnn_gpu_eval.py" "$MODEL_EXT/localnn_problem.py" "$MODEL_EXT/localnn_mesh.py" "$MODEL_EXT/localnn_late.py" "$P/" 2>/dev/null
mkdir "$P/clip_ext" && cp -p "$CX"/*.py "$CX/SHA256SUMS" "$P/clip_ext/"
( cd "$P" && sha256sum *.py clip_ext/*.py localnn_gpu_eval.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|SEG|OUT|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
stamp "1. prep"
rec $PY "$MODEL_EXT/localnn_gpu_eval.py" prep --seg "$SEG" --base "$BASE" "${FTO[@]}" --out "$OUT/inputs"
$PY "$MODEL_EXT/localnn_gpu_eval.py" prep --seg "$SEG" --base "$BASE" "${FTO[@]}" --out "$OUT/inputs" > "$P/prep.log" 2>&1 || { tail -5 "$P/prep.log"; exit 1; }
tail -1 "$P/prep.log" | cut -c1-400
FIX=(); while read -r f; do FIX+=(--fix "$f"); done < "$OUT/inputs/fix_args.txt"
B=$OUT/inputs/localnn_activity.templates.npz
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
stamp "2. GPU evaluation on $NG GPUs ($DEV) -> $OUT/fit"
CMD=($PY "$CX/clip_refit.py" --variant D2 --point-table "$PTS" "${FIX[@]}" --crossing-weight 10 --devices "$DEV" \
    --precision mixed --max-nfev 1 --segment 1 --template-start 0 -- "$BASE" \
    --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 \
    --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$OUT/fit" \
    --activity-template "$B" localnn_a)
rec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/fit.log" 2>&1; rc=$?
stamp "GPU evaluation exit $rc"
grep -E '^(clip_refit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: [0-9])' "$P/fit.log" | cut -c1-400
tail -4 "$P/fit.log" | cut -c1-300
if [ $rc -eq 0 ]; then
  stamp "3. compare"
  GD=${GRID_DIR:-$R/nonneg/localnn/grid_150m}
  CC=(compare --seg "$SEG" --gpu "$OUT/fit" --design "${DESIGN_DIR:-$GD/design_h2oqmin0.02}" --base "$BASE" --d2 "${D2_DIR:-$R/clip_refit/D2/final}" --point-table "$PTS" "${FTO[@]}" --out "$OUT/gpu_check.json")
  rec $PY "$MODEL_EXT/localnn_gpu_eval.py" "${CC[@]}"
  $PY "$MODEL_EXT/localnn_gpu_eval.py" "${CC[@]}" > "$P/compare.log" 2>&1; rc=$?
  tail -1 "$P/compare.log" | cut -c1-1500
fi
stamp "done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_gpu-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
