#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_cpu.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_cpu.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 7c2455cb6f3949186a80cfb6bf6b710850f7216652ee3501a521c9f5e2f69ec7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test, CPU). CPU tests of the softplus / squared models (model_ext/test_nonneg.py: unit, workers, e2e) from
#   MODEL_EXT. 2026-11-11.
# RUN: sbatch --export=ALL,MODEL_EXT=<snapshot> nonneg_cpu.sbatch   (job 124128: 20 tests, 0 failed)
# ==== end of header ====
# CPU tests of the non-negative activity models (model_ext/test_nonneg.py: unit, workers, e2e) from MODEL_EXT.
#SBATCH --job-name=cops_nonneg_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_nonneg_cpu-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
D=$COPS_OUTPUT_ROOT/checks/nonneg_cpu_$SLURM_JOB_ID; mkdir -p "$D"; P=$D/provenance; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/"; cp -p "$MODEL_EXT"/*.py "$P/" 2>/dev/null; [ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
CMD=("$APP/.venv/bin/python" -W ignore "$MODEL_EXT/test_nonneg.py" --scratch "$D/scratch")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): MODEL_EXT=$MODEL_EXT"
"${CMD[@]}"; rc=$?
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
