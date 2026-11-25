#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_all_fits_summary_docx.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_all_fits_summary_docx.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2992f9a9d466715c4791c04f5d3095f6e9ca416eb759941948a160ee295cbe58
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: build and checks of the all-fits summary document
# RUN: move the output folder aside, then sbatch
# ==== end of header ====
#SBATCH --job-name=cops_allfitsdoc
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_allfitsdoc-%j.out
# Summary document of all fits of the new-law round (figures_cluster/scripts/q002_all_fits_summary_docx.py): build,
# validate, OOXML schema check, plain-text dump. The output folder must be new: move an existing one aside first.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
O=$R/figures_cluster/12M/q002_all_fits_summary
S=$R/figures_cluster/scripts
OBS=$S/data/q002_all_fits_observations_20261121.json
P=$R/tools/docx_venv/bin/python
D=/home/vtenishev/.claude/skills/synced/7064e95a-1953-4569-8bdf-a09fe2695c3b_de6ecb44-8ed6-4b14-b43c-146225d8337e/docx/scripts/office
[ ! -e "$O" ] || { echo "$O exists: move it aside first" >&2; exit 2; }
mkdir "$O"
V=$O/provenance-$SLURM_JOB_ID
mkdir "$V"
cp -p "$0" "$V/q002_all_fits_summary_docx.sbatch.used"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$V/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $P --version; } > "$V/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$V/module_list.txt" 2>&1 || echo "no module command" >> "$V/module_list.txt"
$P -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$V/pip_list.txt"
cd "$O"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy \
  $P -I $S/q002_all_fits_summary_docx.py build --out "$O" --observations "$OBS" --jobs 8 2>&1 | tee "$O/build.log"
$P -I $S/q002_all_fits_summary_docx.py validate --out "$O" 2>&1 | tee -a "$O/build.log"
cd $D
$P validate.py "$O/12M_q002_all_fits_summary.docx" | tee "$O/validation/xsd_validate.txt"
$P -I -c "
import docx,sys
d=docx.Document(sys.argv[1]); out=open(sys.argv[2],'w')
for p in d.paragraphs: out.write(p.text+'\n')
for t in d.tables:
    out.write('--- table\n')
    for r in t.rows: out.write(' | '.join(c.text for c in r.cells)+'\n')
" "$O/12M_q002_all_fits_summary.docx" "$O/validation/text_dump.txt"
echo "end=$(date -Is)"
