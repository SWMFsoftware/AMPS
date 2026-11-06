#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_44m_c30.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_44m_c30.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: e4a2c2e641a241a65306f0250c5e185d39bf26955895c083bb2ee52512d84c25
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 44M with the best 12M fit's conditions.
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/44M/model_vs_observations_grain_best12m_conditions/
# ==== end of header ====
# 44M rerun with exactly the conditions of the best 12M grain fit (C30_n2000_grain_cross10; job 122078, user
# 2026-11-05: "rerun 44M exactly for the same conditions as those of the best fit of 12M"): the 44M fit and its fine
# tuning against the 12M fit and its fine tuning, and the COPS observations; writes HOW-TO-RERUN.md.
# Login node, no GPU:   figures_44m_c30.sh [<job id>]
set -u
JOB=${1:-122078}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R=$COPS_OUTPUT_ROOT/runs/shap8_044m_flat_earlier
O=$COPS_OUTPUT_ROOT/figures_cluster/44M/model_vs_observations_grain_best12m_conditions; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
FIT=C30_n2000_grain_cross10
OF=$R12/grain/fits/$FIT; OFT=$R12/grain/finetune/${FIT}_m5_i5_cpu
F=$R/grain/fits/$FIT;    FT=$R/grain/finetune/${FIT}_m5_i5_cpu
.venv/bin/python $S/finetune_table.py "$OFT" "$OF" "$TB/12M_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$FT" "$F" "$TB/44M_finetuned.dat" || exit 1
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/44M_vs_12M_cliff60_free_T" --ref "$F" \
    --title 'SHAP8, cliff $\geq$60$^\circ$ template, free temperatures, dust grain in NG + crossing constraint: 44M vs 12M, same fit settings' \
    --fit "12M fit=$OF:dashed" --fit "12M fit, fine-tuned=$TB/12M_finetuned.dat:dotted" \
    --fit "44M fit=$F:solid" --fit "44M fit, fine-tuned=$TB/44M_finetuned.dat:dashdot" \
    && echo "wrote $O/44M_vs_12M_cliff60_free_T_*"
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 44M fit with the conditions of the best 12M grain fit (C30_n2000_grain_cross10)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_44m_c30.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

The 12M command recorded in \`$OF/fit_manifest.json\` (default start, free temperatures, --tbounds 30 400,
cliff_ge_60deg, grain + crossing weight 10, mixed, max_nfev 2000) with the 44M base, record store and template
bundle; then grain_finetune.py (5 min, 5 intervals, cpu). Job $JOB, provenance \`$R/grain/provenance-$JOB/\`.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_44m_c30.sh $JOB
\`\`\`

## A1. The 44M fit and fine tuning (one 8-GPU job, about 2 h 25 min)

Move \`$F\` and \`$FT\` aside first (as \`*.incomplete_<UTC>\`), then

\`\`\`bash
sbatch $J/grain_44m_c30.sbatch
\`\`\`

## A2. The 44M base, store, gravity and bundle from scratch (genus 6 gravity, user 2026-11-04)

\`sbatch $J/prepare_44m.sbatch\` (8 h limit; 6 GPUs for the gravity step, see 44M STATUS.md), then A1.

## B. The 12M fit itself

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_repro--2369abeacf15807199ce76130dcea79c4bfcf68d/HOW-TO-RERUN.md\`.
MD
echo "wrote $O/HOW-TO-RERUN.md"
