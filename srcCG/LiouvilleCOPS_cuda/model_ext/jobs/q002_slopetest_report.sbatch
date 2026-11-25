#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_slopetest_report.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_slopetest_report.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: c9b56c6ad559f84b9c66319e32b3c20656bd405a0b3e39d4d0752ea98b24f709
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job of figures/q002_slopetest_report.py
# RUN: sbatch --export=ALL,RUNROOT=<run root> jobs/q002_slopetest_report.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_slope3c_report
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_slope3c_report-%j.out
# Table of the finished slope-effect control fits (decision 3c): sbatch --export=ALL,RUNROOT=<run root>[,MULTISTART=<refits>] <this file>
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
A=/data001/highresolu_vt/AMPS
PY=$A/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_slopetest_report.py
RUNROOT=${RUNROOT:?RUNROOT is required}
O=$RUNROOT/report_$SLURM_JOB_ID
P=$RUNROOT/provenance-report-$SLURM_JOB_ID
mkdir "$P"
cp -p "$0" "$P/q002_slopetest_report.sbatch.used"
cp -p "$S" "$P/"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
MS=()
[ -n "${MULTISTART:-}" ] && MS=(--multistart "$MULTISTART")
exec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy $PY "$S" --runroot "$RUNROOT" --out "$O" "${MS[@]}"
