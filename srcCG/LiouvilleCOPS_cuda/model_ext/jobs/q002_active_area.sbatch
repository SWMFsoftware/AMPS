#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_active_area.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_active_area.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0c2951917860a1027decc28102069ea54bed82adb91f4ef6aff12cf3bf4f377a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job of figures/q002_active_area.py
# RUN: sbatch jobs/q002_active_area.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_active
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121/logs/%x_%j.log
# Switched-off share of the seen surface, user decision 1b of 2026-11-21 (figures_cluster/scripts/q002_active_area.py).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
D=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts
A=/data001/highresolu_vt/AMPS
PY=$A/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
P=$D/provenance_$SLURM_JOB_ID
mkdir "$P"
cp -p "$0" "$P/q002_active_area.sbatch.used"
cp -p $S/q002_active_area.py $S/q002_slope_activity.py "$P/"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; $PY -c "import numpy; print('numpy', numpy.__version__)"; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
cd $S
exec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy $PY q002_active_area.py --out "$D"
