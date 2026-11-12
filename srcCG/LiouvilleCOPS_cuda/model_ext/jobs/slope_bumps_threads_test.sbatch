#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_bumps_threads_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_bumps_threads_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0a9ad16dbc51aa70de5f1a67bfd0ac25003b45958c09db378e514104b7c04b79
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test). Timing / bitwise test of further L_THREADS values (default 64) for the bump + slope-multiplier short
#   fit, compared with an existing short fit (default the 32-thread run of 123984); written 2026-11-10 because the
#   64-thread run of 123984 was stopped by the 1-h limit.
# RUN: sbatch --export=ALL,MODEL_EXT=<snapshot>,TAG=<tag> slope_bumps_threads_test.sbatch. NEVER SUBMITTED (64 threads
#   not timed; L_THREADS=32 kept).
# ==== end of header ====
# Timing / bitwise test of further L_THREADS values for the bump + slope-multiplier fit (option 3a), one node: the short
# fit of jobs/slope_bumps_test.sbatch STEP=fits (START=nnls from D2's grain, NFEV 5, CHECK 1) with L_THREADS = each of
# THREADS (default "64"), then test_slope_mult.py --compare-runs REF <the new runs> (REF = an existing short fit's
# segment, default the 32-thread run of job 123984). Written 2026-11-10 because the 64-thread run of job 123984 was
# stopped by gpu_fit's wall-time guard (four 14-min fits in the 1-h test limit) before its first evaluation.
# No floor-0 check (job 123984 passed it with the same snapshot).
# Settings: MODEL_EXT (required, frozen snapshot), TAG (required), THREADS, REF.
#SBATCH --job-name=cops12m_bumps_threads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=00:40:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_bumps_threads-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
MODEL_EXT=$(readlink -f "${MODEL_EXT:?set MODEL_EXT}"); TAG=${TAG:?set TAG}
THREADS=${THREADS:-64}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs; T_=$R/slope_mult_tests
REF=$(readlink -f "${REF:-$T_/g_bumps_short_l32_20261110T201212Z/seg01}")
[ -f "$REF/evaluations.txt" ] || { echo "no reference run $REF" >&2; exit 2; }
BUMPS=$R/slope_bumps/bundle/bumps32.npz
PY="$APP/.venv/bin/python -W ignore"
D=$T_/bumps_threads_test_$TAG; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$D"; P=$D/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/slope_bumps_threads_test.sbatch"; cp -p "$J/slope_mult_fit.sbatch" "$P/"
cp -p "$MODEL_EXT"/*.py "$MODEL_EXT/SHA256SUMS" "$P/" 2>/dev/null
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ)"; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): THREADS=$THREADS REF=$REF MODEL_EXT=$MODEL_EXT TAG=$TAG"
rc=0
SEGS=()
for N in $THREADS; do
    label=short_l$N
    echo "$(date -u +%FT%TZ) $label"
    echo "env TESTNAME=g_bumps_short_l${N}_$TAG BUMPS=$BUMPS START=nnls CHECK=1 NFEV=5 L_THREADS=$N MODEL_EXT=$MODEL_EXT bash $J/slope_mult_fit.sbatch # $label" >> "$P/commands.sh"
    env TESTNAME=g_bumps_short_l${N}_$TAG BUMPS=$BUMPS START=nnls CHECK=1 NFEV=5 L_THREADS=$N MODEL_EXT=$MODEL_EXT \
        bash "$J/slope_mult_fit.sbatch" > "$P/$label.log" 2>&1 || rc=1
    echo "  exit $?"
    grep -E "slope_mult: (start|L = Y a)|start check \(start|^lambda|gpu_fit: [0-9]|no time|complete" "$P/$label.log" | cut -c1-420
    SEGS+=("$T_/g_bumps_short_l${N}_$TAG/seg01")
done
CCMD=($PY "$MODEL_EXT/test_slope_mult.py" --compare-runs "$REF" "${SEGS[@]}" --json-out "$D/compare.json")
printf '%q ' "${CCMD[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
"${CCMD[@]}" > "$P/compare.log" 2>&1 || rc=1
grep -E "PASS|FAIL|SUMMARY|Error" "$P/compare.log" | cut -c1-900
echo "$(date -u +%FT%TZ) exit $rc"
exit $rc
