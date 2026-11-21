#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_late.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_late.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: bbf2fe7bb6f8306ee665a4c9f49804376d62da39980685440aeb399be6121142
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU): the production late kernels runs/.../nonneg/localnn/late_kernels_150m.npz
#   (localnn_late.py).
# 2026-11-18: 2026-11-18: CONTRIB, NONNEG_CODE and OUT defaults of the current law
#   (late_kernels_150m_h2oqmin0.02.npz).
# RUN 2026-11-18: no job since the change (no C fit started).
# ==== end of header ====
# Option C stage 2 (lead 2026-11-15): the late kernels of the CPU fine tuning (model_ext/localnn_late.py): the 9 row
# masks the frozen positive_finetune.late_only_facets asks for (structural / positive seen per row), the 271 late facets
# and G_NG / G_RG = dp/dA of every row on them (exact kernel arithmetic of the frozen nonneg_contrib code). CPU.
# Settings: OUT (default runs/.../nonneg/localnn/late_kernels_150m.npz; never overwritten), MODEL_EXT (frozen localnn
#   snapshot), NONNEG_CODE (default model_ext_frozen/nonneg_20261113T035423Z), CONTRIB, PROCS.
#SBATCH --job-name=cops_localnn_late
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_late-%j.out
set -uo pipefail
# 2026-11-18 (H2O night law q_min = 0.02 q_max): the products of the current law by default (DESIGN_DIR, CONTRIB, LATE,
# NONNEG_CODE); the loads refuse a design / late kernels of another law (localnn_problem.law_check)
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
NONNEG_CODE=$(readlink -f "${NONNEG_CODE:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_20261118T072729Z}")
CONTRIB=$(readlink -f "${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261118T072729Z}")
OUT=${OUT:-$R/nonneg/localnn/late_kernels_150m_h2oqmin0.02.npz}
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-16}}
PY=$APP/.venv/bin/python
[ -e "$OUT" ] && { echo "$OUT exists: not overwritten"; exit 2; }
P=$OUT.provenance-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cd "$APP" || exit 1
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): OUT=$OUT MODEL_EXT=$MODEL_EXT NONNEG_CODE=$NONNEG_CODE"
cp -p "$0" "$P/localnn_late.sbatch.used"; cp -p "$MODEL_EXT/localnn_late.py" "$MODEL_EXT/localnn_problem.py" "$P/"
( cd "$P" && sha256sum *.py localnn_late.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
CMD=("$PY" -W ignore "$MODEL_EXT/localnn_late.py" --base "$R/base" --d2 "$R/clip_refit/D2/final"
     --point-table "$R/late_grid_points.txt" --contrib-code "$CONTRIB" --nonneg-code "$NONNEG_CODE" --out "$OUT" --procs "$PROCS")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/late.log" 2>&1; rc=$?
tail -n 3 "$P/late.log" | cut -c1-600
echo "$(date -u +%FT%TZ) done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_late-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
