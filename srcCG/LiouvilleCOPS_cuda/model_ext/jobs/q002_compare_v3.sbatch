#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_base_20261118/jobs/q002_compare_v3.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_base_20261118/jobs/q002_compare_v3.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6814c605d8c4789959180e99dffdef75bc85b4439dee9fd83d2aefa37b2f5e57
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18, Phase 0). compare-sample v3 -> tests/q002_base_20261118/sample_report_v3.json (job 136521).
# ==== end of header ====
#SBATCH --job-name=q002_compare
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/q002_base_20261118/logs/%x-%j.out
# Phase 0 (2026-11-18): compare-sample of nightfix_base_q002.py (v2 fixed the reload; v3: NGBelowDeckFraction reported, not gated) on the sample re-runs A-D of job 136506,
# after v1 crashed there (IndexError: it reloaded only the 21 used rows). Report: tests/q002_base_20261118/sample_report.json
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; W=$COPS_OUTPUT_ROOT/tests/q002_base_20261118
CODE=${CODE:?CODE required}; PY=$APP/.venv/bin/python
P=$W/provenance/job-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cp -p "$0" "$P/q002_compare.sbatch.used"; cp -p "$CODE/nightfix_base_q002.py" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
"$PY" -m pip list > "$P/pip_list.txt" 2>&1; lscpu | grep 'Model name' > "$P/cpu.txt"
cd "$APP" || exit 1
CMD=("$PY" -W ignore "$CODE/nightfix_base_q002.py" compare-sample --old-base "$R/base" --new "$W/cpp_sample_A" \
     --new "$W/cpp_sample_B" --new "$W/cpp_sample_C" --new "$W/cpp_sample_D" --out "$W/sample_report_v3.json")
printf '%q ' "${CMD[@]}" > "$P/commands.sh"; echo >> "$P/commands.sh"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "${CMD[@]}"; rc=$?
echo "$(date -u +%FT%TZ) compare-sample rc=$rc"
exit $rc
