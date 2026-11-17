#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_ft_summary.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_ft_summary.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: c57b58d3128999260f623c183e0113d872dec9d1919d450161b098a0adc45860
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. CPU summary of one fine tuning (nonneg_ft_summary.py), optional reference fine tuning.
# RUN: sbatch --export=ALL,FT=<ft seg dir>,REF=<ref ft seg dir>,OUT=<summary.txt>,MODEL_EXT=<snapshot> nonneg_ft_summary.sbatch.
# ==== end of header ====
# Summary of ONE finished late-facet fine tuning (model_ext/nonneg_ft_summary.py; the delta-regularization tests of the
# user's 2026-11-15 decision): costs, late-window rms ln NG/RG, P_delta, rms / max |h + delta|, facets at the cap, the
# seen area holding 50 % / 90 % of the seen production at 10:33:10, the tuned facets' share of the late signal.
# Settings: FT (the ft seg dir; required), REF (optional reference ft seg dir), OUT (required), MODEL_EXT (default
#   $COPS_OUTPUT_ROOT/model_ext), GEOMETRY (default job_124775), BASE (default the lategrid base). CPU, ~1 min.
#SBATCH --job-name=cops_nonneg_ft_summary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=01:00:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_nonneg_ft_summary-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
FT=$(readlink -f "${FT:?set FT=<ft seg dir>}"); REF=${REF:-}
OUT=$(readlink -m "${OUT:?set OUT}")
GEOMETRY=$(readlink -f "${GEOMETRY:-$R/maps_geometry/job_124775/geometry}")
BASE=$(readlink -f "${BASE:-$R/base}")
[ -e "$OUT" ] && { echo "exists: $OUT" >&2; exit 2; }
P=$(dirname "$OUT")/$(basename "$OUT" .txt).provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/nonneg_alt_summary.sbatch"
for f in nonneg_ft_summary.py nonneg_alt_summary.py nonneg_maps.py nonneg_contrib.py positive_problem.py; do cp -p "$MODEL_EXT/$f" "$P/"; done
[ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
lscpu > "$P/lscpu.txt" 2>&1
CMD=("$APP/.venv/bin/python" -W ignore "$P/nonneg_ft_summary.py" --ft "$FT" --base "$BASE" --geometry "$GEOMETRY"
     --out "$OUT" --procs "${SLURM_CPUS_PER_TASK:-64}")
[ -n "$REF" ] && CMD+=(--reference-ft "$(readlink -f "$REF")")
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); MODEL_EXT $MODEL_EXT"; grep '^export ' "$0"
  echo "cd $APP"; printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): ft summary of $FT${REF:+ (reference $REF)} -> $OUT (MODEL_EXT $MODEL_EXT)"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/summary.log" 2>&1; rc=$?
tail -5 "$P/summary.log" | cut -c1-600
[ -f "$OUT" ] && cat "$OUT"
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
