#!/bin/bash
# ==== version-controlled copy: header added 2026-11-23, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/plausibility_q002_v2.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/plausibility_q002_v2.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: ec2255b6bf4a321b0109a427678214ac21eb5864dc70cb78828f4eda3e7b1126
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job of figures/plausibility_q002_v2.py
# RUN: sbatch jobs/plausibility_q002_v2.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_plaus_v2
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_plausibility/logs/%x_%j.log
# Physical plausibility of every new-law fit with one method (figures_cluster/scripts/plausibility_q002_v2.py; user
# 2026-11-23). Evaluations only.
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
D=$R/figures_cluster/12M/q002_plausibility
S=$R/figures_cluster/scripts
PY=$A/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
P=$D/provenance_v2_$SLURM_JOB_ID
mkdir "$P"
cp -p "$0" "$P/plausibility_q002_v2.sbatch.used"
cp -p $S/plausibility_q002_v2.py $S/q002_slope_activity.py "$P/"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
cd $S
exec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy $PY plausibility_q002_v2.py --out "$D/v2_$SLURM_JOB_ID"
