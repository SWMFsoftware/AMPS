#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_alt_summary.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_alt_summary.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: ab6fa6b0bc1944937d673e5de5201dac1b002c419ee3475ca6bcff5ad0e4d000
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. CPU summary of one finished fit (nonneg_alt_summary.py), about 30 s with 64 processes.
# RUN: sbatch --export=ALL,FIT=<root>/fit/final,OUT=<root>/summary.txt,MODEL_EXT=<summary snapshot>
#   nonneg_alt_summary.sbatch (test 125030; submitted by nonneg_alt_stage1.sh after each guard).
# ==== end of header ====
# Summary of ONE finished non-negative-model fit (model_ext/nonneg_alt_summary.py; the stage-1 report numbers of options
# A, B, D, lead 2026-11-14): costs (data, crossing, spread penalty), termination, rms ln NG/RG overall / before the gap /
# late window, parameters at bounds, Q_H2O at 10:33:10, the active seen area (A > 1 % of the srcCG median of the seen
# facets), the seen area holding 50 % / 90 % of the seen production, NG/RG shares from the active area.
# Settings: FIT (the fit dir, e.g. <root>/fit/final; required), OUT (default <root>/summary.txt where <root> = the
#   directory above fit/), MODEL_EXT (the model_ext copy; default $COPS_OUTPUT_ROOT/model_ext), GEOMETRY (default the
#   validated maps geometry job_124775), BASE (default the lategrid base). CPU, about 3 min with 64 processes.
#SBATCH --job-name=cops_nonneg_summary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=01:00:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_nonneg_summary-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
FIT=$(readlink -f "${FIT:?set FIT=<fit dir>}")
ROOT=$(dirname "$(dirname "$(readlink -m "${FIT_LINK:-$FIT}")")")
OUT=$(readlink -m "${OUT:-$ROOT/summary.txt}")
GEOMETRY=$(readlink -f "${GEOMETRY:-$R/maps_geometry/job_124775/geometry}")
BASE=$(readlink -f "${BASE:-$R/base}")
[ -e "$OUT" ] && { echo "exists: $OUT" >&2; exit 2; }
P=$(dirname "$OUT")/$(basename "$OUT" .txt).provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/nonneg_alt_summary.sbatch"
for f in nonneg_alt_summary.py nonneg_maps.py nonneg_contrib.py positive_problem.py; do cp -p "$MODEL_EXT/$f" "$P/"; done
[ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
lscpu > "$P/lscpu.txt" 2>&1
CMD=("$APP/.venv/bin/python" -W ignore "$P/nonneg_alt_summary.py" --fit "$FIT" --base "$BASE" --geometry "$GEOMETRY"
     --out "$OUT" --procs "${SLURM_CPUS_PER_TASK:-64}")
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); MODEL_EXT $MODEL_EXT"; grep '^export ' "$0"
  echo "cd $APP"; printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): summary of $FIT -> $OUT (MODEL_EXT $MODEL_EXT)"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/summary.log" 2>&1; rc=$?
tail -5 "$P/summary.log" | cut -c1-600
[ -f "$OUT" ] && cat "$OUT"
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
