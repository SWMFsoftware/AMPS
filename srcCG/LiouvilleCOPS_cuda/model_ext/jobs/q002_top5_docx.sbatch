#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/q002_top5_docx.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/q002_top5_docx.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: b05577982733514caa48b35e53b01a37772900fd6af583b65b382c1d4c9ca732
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm CPU job: build, schema check and text dump of 12M_q002_top5_models.docx (2026-11-21). Job 137815
#   waited on QOSMaxCpuPerUserLimit and was cancelled; the same commands ran on the login node (logs/doc_build_*)
# RUN: sbatch jobs/q002_top5_docx.sbatch (move the document, images/ and validation/ aside first)
# ==== end of header ====
#SBATCH --job-name=cops_top5doc
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/logs/%x_%j.log
# Word document on the five most plausible new-law fits (figures_cluster/scripts/q002_top5_docx.py), 2026-11-21:
# build, reopen check, OOXML schema check (docx skill validate.py), plain-text dump for proofreading.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4
O=/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models
P=/data001/highresolu_vt/cops_output/tools/docx_venv/bin/python
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts
D=/home/vtenishev/.claude/skills/synced/7064e95a-1953-4569-8bdf-a09fe2695c3b_de6ecb44-8ed6-4b14-b43c-146225d8337e/docx/scripts/office
echo "host=$(hostname) start=$(date -Is) job=$SLURM_JOB_ID"
{ echo "job $SLURM_JOB_ID $(date -Is) $(hostname)"; $P --version; $P -c "import docx, PIL; print('python-docx', docx.__version__, 'Pillow', PIL.__version__)";
  sha256sum $S/q002_top5_docx.py $S/q002_top5_content.py; } > $O/provenance_doc_$SLURM_JOB_ID.txt
mkdir -p $O/validation
$P -I $S/q002_top5_docx.py --out $O
cd $D && $P validate.py $O/12M_q002_top5_models.docx | tee $O/validation/xsd_validate.txt
$P -I -c "
import docx,sys
d=docx.Document(sys.argv[1]); out=open(sys.argv[2],'w')
for p in d.paragraphs: out.write(p.text+'\n')
for t in d.tables:
    out.write('--- table\n')
    for r in t.rows: out.write(' | '.join(c.text for c in r.cells)+'\n')
" $O/12M_q002_top5_models.docx $O/validation/text_dump.txt
echo "end=$(date -Is)"
