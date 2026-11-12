#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_start.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_start.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 79a23565b900777d5701baf3c340ddc8cf6669e478517ae395e557d917d293bf
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (CPU). Start of the squared model (model_ext/nonneg_start.py; user decision 2026-11-11): b = the
#   area-weighted least squares to D2's activity max(0, D2 field). OUT, THREADS (32), MODEL_EXT, D2.
# RUN: sbatch --export=ALL,MODEL_EXT=<snapshot> nonneg_start.sbatch   (test 124129; production 124138 ->
#   runs/shap8_012m_lategrid_dz1.5/nonneg/bundle/squared_start.npz)
# ==== end of header ====
# Start of the squared-harmonics model (user decision 2026-11-11; model_ext/nonneg_start.py): b = the area-weighted
# least squares of A(m = 1) = A_s ((Y b)^2 + sum ft_c*) to D2's activity max(0, D2 field) over all facets.
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
OUT=$(readlink -m "${OUT:-$R/nonneg/bundle/squared_start.npz}"); THREADS=${THREADS:-32}
D2=$(readlink -f "${D2:-$R/clip_refit/D2/final}")
[ -e "$OUT" ] && { echo "exists: $OUT" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"; P=$(dirname "$OUT")/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/nonneg_start.sbatch"; cp -p "$MODEL_EXT"/nonneg_*.py "$MODEL_EXT"/slope_mult.py "$MODEL_EXT"/clip_problem.py \
    "$MODEL_EXT"/positive_problem.py "$P/"; [ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
CMD=("$APP/.venv/bin/python" -W ignore "$P/nonneg_start.py" --run "$R/base" --d2 "$D2/fit_restart.npz"
     --ft-templates-from "$D2/fit_manifest.json" --out "$OUT" --threads "$THREADS")
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"
  printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): nonneg_start -> $OUT (MODEL_EXT $MODEL_EXT)"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/nonneg_start.log" 2>&1; rc=$?
grep -E '^nonneg_start|Traceback|Error' "$P/nonneg_start.log" | cut -c1-900
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
