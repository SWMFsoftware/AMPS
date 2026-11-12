#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_start_determinism.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_start_determinism.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 99c1f7317e16d992178ff3c58718448637e95d93d6d866a07903b5ba666a3cb2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test, CPU). Determinism check of nonneg_start.py after the production start 124138 (node032) differed from
#   the test start 124129 (node023) at 8e-15: the frozen snapshot's nonneg_start.py with THREADS into OUT, with the
#   node's CPU and numpy's OpenBLAS core recorded. Result: bitwise within a CPU type, not across. 2026-11-11.
# RUN: sbatch --nodelist=node023 --export=ALL,MODEL_EXT=<snapshot>,OUT=<new>,THREADS=32 nonneg_start_determinism.sbatch
#   (jobs 124149-124153)
# ==== end of header ====
# Determinism check of model_ext/nonneg_start.py (2026-11-11, the lead: the production start 124138 on node032 differs
# from the test start 124129 on node023 at rounding level, 23 of 25 coefficients, max rel 8e-15). Runs the FROZEN
# snapshot's nonneg_start.py (MODEL_EXT) into OUT with THREADS and records the node's CPU and numpy's runtime (OpenBLAS
# core, SIMD dispatch). Settings: MODEL_EXT (required), OUT (required, new), THREADS (default 32).
#SBATCH --job-name=cops12m_nonneg_start_det
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=00:30:00
#SBATCH --partition=normal128cores
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_nonneg_start_det-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT}"); OUT=$(readlink -m "${OUT:?set OUT}"); THREADS=${THREADS:-32}
D2=$R/clip_refit/D2/final
[ -e "$OUT" ] && { echo "exists: $OUT" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"; P=$(dirname "$OUT")/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/"; cp -p "$MODEL_EXT/nonneg_start.py" "$MODEL_EXT/SHA256SUMS" "$P/"
lscpu > "$P/lscpu.txt" 2>&1
"$APP/.venv/bin/python" -c "import numpy as np; np.show_runtime(); print(np.__version__)" > "$P/numpy_runtime.txt" 2>&1
"$APP/.venv/bin/python" -c "import sys; sys.path.insert(0, '$MODEL_EXT'); import slope_mult as sm; print(sm.blas_info())" > "$P/blas_info.txt" 2>&1
CMD=("$APP/.venv/bin/python" -W ignore "$MODEL_EXT/nonneg_start.py" --run "$R/base" --d2 "$D2/fit_restart.npz"
     --ft-templates-from "$D2/fit_manifest.json" --out "$OUT" --threads "$THREADS")
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; printf '%q ' "${CMD[@]}"; echo; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname) ($(grep -m1 'Model name' "$P/lscpu.txt" | sed 's/  */ /g')): THREADS=$THREADS -> $OUT"
"${CMD[@]}" > "$P/nonneg_start.log" 2>&1; rc=$?
grep -E '^nonneg_start: (quality|b0)|Traceback|Error' "$P/nonneg_start.log" | cut -c1-300
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
