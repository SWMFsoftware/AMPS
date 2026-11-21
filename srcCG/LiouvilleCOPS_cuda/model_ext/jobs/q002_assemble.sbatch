#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_base_20261118/jobs/q002_assemble.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_base_20261118/jobs/q002_assemble.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 99a251cc4fd85c1f994d30acbfc6c66e336d8d71b2f29ff931a05f90424d49eb
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18, Phase 0). assemble base_q002 (v3 tool; requires the full C++ run complete and the v3 report
#   passed).
# 2026-11-18: comment for v4; resubmitted as job 136533 (136522 failed before writing anything)
# ==== end of header ====
#SBATCH --job-name=q002_assemble
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/q002_base_20261118/logs/%x-%j.out
# Phase 0 step 2 (2026-11-18): assemble runs/shap8_012m_lategrid_dz1.5/base_q002 from the full new-law C++ run of job
# 136506 (cpp_full_q002, WriteKernels off) and hard links of the existing base, with nightfix_base_q002.py v4 assemble (v4: a WriteKernels = off run has no elements.bin; no run_manifest.json in the base)
# (identical to v1's). The tool refuses unless sample_report.json passed, every geometry / data cell of the full table
# equals the base's and elements.bin is byte-identical; this script also requires the C++ run to have completed.
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; W=$COPS_OUTPUT_ROOT/tests/q002_base_20261118
CODE=${CODE:?CODE required}; PY=$APP/.venv/bin/python
P=$W/provenance/job-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cp -p "$0" "$P/q002_assemble.sbatch.used"; cp -p "$CODE/nightfix_base_q002.py" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
grep -q '^Done: ' "$W/cpp_full_q002/build.log" || { echo "the full C++ run has not completed" >&2; exit 2; }
"$PY" -c "import json,sys; m=json.load(open(sys.argv[1])); sys.exit(0 if m['states']['process']=='complete' else 1)" \
    "$W/cpp_full_q002/run_manifest.json" || { echo "the full run's manifest is not complete" >&2; exit 2; }
cd "$APP" || exit 1
CMD=("$PY" -W ignore "$CODE/nightfix_base_q002.py" assemble --old-base "$R/base" --cpp-run "$W/cpp_full_q002" \
     --cpp-input "$W/inputs/full_q002.in" --sample-report "$W/sample_report_v3.json" --out "$R/base_q002")
printf '%q ' "${CMD[@]}" > "$P/commands.sh"; echo >> "$P/commands.sh"
"${CMD[@]}"; rc=$?
echo "$(date -u +%FT%TZ) assemble rc=$rc"
exit $rc
