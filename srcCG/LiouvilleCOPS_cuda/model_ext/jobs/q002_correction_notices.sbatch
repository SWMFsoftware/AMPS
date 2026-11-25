#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_correction_notices.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_correction_notices.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 21efe504c803986db17c119e2e3522deb4534a1b0d7d5ba869f4ef84ec4d3ad9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job of figures/q002_correction_notices.py + schema check
# RUN: sbatch jobs/q002_correction_notices.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_notices
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121/logs/%x_%j.log
# Correction notices on the older documents and figure sets (decision 2b, 2026-11-21; q002_correction_notices.py),
# then the OOXML schema check of the three documents.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
R=/data001/highresolu_vt/cops_output
F=$R/figures_cluster/12M
A=/data001/highresolu_vt/AMPS
P=$R/tools/docx_venv/bin/python
D=/home/vtenishev/.claude/skills/synced/7064e95a-1953-4569-8bdf-a09fe2695c3b_de6ecb44-8ed6-4b14-b43c-146225d8337e/docx/scripts/office
V=$F/q002_decisions_20261121/provenance_$SLURM_JOB_ID
mkdir "$V"
cp -p "$0" "$V/q002_correction_notices.sbatch.used"
cp -p $R/figures_cluster/scripts/q002_correction_notices.py "$V/"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$V/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; $P --version; } > "$V/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$V/module_list.txt" 2>&1 || echo "no module command" >> "$V/module_list.txt"
$P -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$V/pip_list.txt"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy \
  $P -I $R/figures_cluster/scripts/q002_correction_notices.py --root "$F" --stamp "$(date -u +%Y%m%dT%H%M%SZ)"
cd $D
for f in stage2_q002_comparison/12M_stage2_q002_comparison.docx full_sets_stage2_comparison/12M_stage2_full_sets_comparison.docx cliff_effect_stage2/12M_stage2_cliff_effect.docx; do
  echo "schema: $f"
  $P validate.py "$F/$f" | tail -1
done
