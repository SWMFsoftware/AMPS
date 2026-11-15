#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_start.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_start.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 142ac802ee9808e01ee9c91d86e52eefb96eb26529309c4a50d25804cba353cc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (CPU). Start of the squared model (model_ext/nonneg_start.py; user decision 2026-11-11): b = the
#   area-weighted least squares to D2's activity max(0, D2 field). OUT, THREADS (32), MODEL_EXT, D2.
# RUN: sbatch --export=ALL,MODEL_EXT=<snapshot> nonneg_start.sbatch   (test 124129; production 124138 ->
#   runs/shap8_012m_lategrid_dz1.5/nonneg/bundle/squared_start.npz)
# 2026-11-14: MODEL=sos RANK=r [REFINE=lm].
# ==== end of header ====
# Start of the squared-harmonics model (user decision 2026-11-11; model_ext/nonneg_start.py): b = the area-weighted
# least squares of A(m = 1) = A_s ((Y b)^2 + sum ft_c*) to D2's activity max(0, D2 field) over all facets.
# MODEL=sos RANK=r (option A, 2026-11-14): the Gram-matrix start of the sos model (nonneg_start.py --model sos --rank r;
#   default OUT runs/.../nonneg/bundle/sos_start_r<r>.npz). REFINE=lm (sos, proposal 2026-11-14 pending the lead's OK):
#   then Levenberg-Marquardt of the rank-r problem from that start (--refine lm; default OUT .../sos_start_r<r>_lm.npz).
# Settings: OUT (default runs/shap8_012m_lategrid_dz1.5/nonneg/bundle/squared_start.npz; a new file), THREADS (32),
#   MODEL_EXT (the model_ext copy to run; default $COPS_OUTPUT_ROOT/model_ext), D2 (D2's final fit dir).
#SBATCH --job-name=cops12m_nonneg_start
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_nonneg_start-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
MODEL=${MODEL:-squared}; RANK=${RANK:-}; REFINE=${REFINE:-}
case $MODEL in squared) DEF=$R/nonneg/bundle/squared_start.npz; MX=()
         [ -n "$REFINE" ] && { echo "REFINE is a sos setting" >&2; exit 2; };;
    sos) [ -n "$RANK" ] || { echo "MODEL=sos needs RANK" >&2; exit 2; }; DEF=$R/nonneg/bundle/sos_start_r$RANK.npz
         MX=(--model sos --rank "$RANK")
         case $REFINE in "") ;; lm) DEF=$R/nonneg/bundle/sos_start_r${RANK}_lm.npz; MX+=(--refine lm);;
             *) echo "REFINE must be lm (or unset)" >&2; exit 2;; esac;;
    *) echo "MODEL must be squared or sos" >&2; exit 2;; esac
OUT=$(readlink -m "${OUT:-$DEF}"); THREADS=${THREADS:-32}
D2=$(readlink -f "${D2:-$R/clip_refit/D2/final}")
[ -e "$OUT" ] && { echo "exists: $OUT" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"; P=$(dirname "$OUT")/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/nonneg_start.sbatch"; cp -p "$MODEL_EXT"/nonneg_*.py "$MODEL_EXT"/slope_mult.py "$MODEL_EXT"/clip_problem.py \
    "$MODEL_EXT"/positive_problem.py "$P/"; [ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
CMD=("$APP/.venv/bin/python" -W ignore "$P/nonneg_start.py" --run "$R/base" --d2 "$D2/fit_restart.npz"
     --ft-templates-from "$D2/fit_manifest.json" --out "$OUT" --threads "$THREADS" "${MX[@]}")
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"
  printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): nonneg_start -> $OUT (MODEL_EXT $MODEL_EXT)"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/nonneg_start.log" 2>&1; rc=$?
grep -E '^nonneg_start|Traceback|Error' "$P/nonneg_start.log" | cut -c1-900
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
