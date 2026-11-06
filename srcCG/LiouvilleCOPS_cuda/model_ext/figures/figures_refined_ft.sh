#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_refined_ft.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_refined_ft.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 63a0ec3130e752240777688d9f8960e6c474424dd0663a940154e70c29d04a39
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M refinement test.
# RUN: on the login node, no GPU: <this> [<job id>]
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_refined_late_mesh/
# ==== end of header ====
# Figures of the 12M refinement test (jobs/refine_12m.sbatch; user request 2026-11-04: the 12M elements seen by NG and
# RG after the gap refined by a factor of 2 (edges halved), the two best 12M grain fits carried to the refined mesh,
# their fine tuning rerun there), against the COPS observations (model_obs_figures.py), and HOW-TO-RERUN.md.
# Login node, no GPU:
#     figures_refined_ft.sh [<job id of the refine_12m.sbatch run that made the results>]
set -u
JOB=${1:-122009}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; X=$COPS_OUTPUT_ROOT/model_ext; J=$COPS_OUTPUT_ROOT/jobs
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R=$COPS_OUTPUT_ROOT/runs/shap8_012m_refined_late
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_refined_late_mesh; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R12/bundle/restart_templates.npz --overwrite"
GR="dust grain in NG + crossing constraint"
declare -A D=([C30]='cliff $\geq$60$^\circ$ template, free temperatures'
              [NS31]='late NG-only + sin(slope) templates, monotonic penalty + ordered temperatures')
declare -A N=([C30]=cliff60_free_T [NS31]=late_NG_sin_slope_monotonic_ordered_T)
BEST=()
for FIT in C30_n2000_grain_cross10 NS31_n2000_mono10_ordered_grain_cross10; do
    M=${FIT%%_*}
    OD=$R12/grain/fits/$FIT; OFT=$R12/grain/finetune/${FIT}_m5_i5_cpu
    RFT=$R/grain/translated/finetune/${FIT}_refined_m5_i5_cpu
    T12=$TB/12M_${N[$M]}_finetuned.dat; TR=$TB/refined_${N[$M]}_finetuned.dat
    .venv/bin/python $S/finetune_table.py "$OFT" "$OD" "$T12" || exit 1
    .venv/bin/python $S/finetune_table.py "$RFT" "$R/base" "$TR" || exit 1
    $P --out "$O/refined_vs_12M_${N[$M]}" --ref "$OD" \
        --title "SHAP8 12M, ${D[$M]}, $GR: fine tuning on 12M and on the refined mesh" \
        --fit "12M fit=$OD:dashed" --fit "12M fit, fine-tuned on 12M=$T12:dotted" \
        --fit "12M fit, fine-tuned on the refined mesh=$TR:solid" \
        && echo "wrote refined_vs_12M_${N[$M]}_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
done
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 12M refinement test (late-seen elements refined by 2, fine tuning on the refined mesh)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_refined_ft.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

What was done (job $JOB; provenance \`$R/provenance/refine-$JOB/\`; the first attempt, job 122002, is in
\`$R/*.incomplete_20261105T002215Z\` and \`$R/grain/translated.incomplete_20261105T002215Z\`):
1. \`model_ext/refine_mesh.py\`: the 265 12M elements first seen in the 5-min fine-tuning window plus all within 10 m
   split 1 -> 4 (edge midpoints, geometry unchanged), conforming closure (1 -> 4 or 1 -> 2):
   12,364,290 -> 12,468,264 elements (\`$R/mesh/refine.json\`).
2. Runner prepare (base + store) on the refined mesh, 12M settings (\`$R/pipeline.json\`).
3. \`model_ext/refined_bundle.py\`: templates (late NG/RG-only, cliff, sin(slope)) carried from the parent elements.
4. Per fit: \`model_ext/translate_grain.py\` (the 12M fit's restart -> the refined base), then
   \`model_ext/grain_finetune.py\` (5 min, 5 intervals, cpu, mixed, max_nfev 2000).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_refined_ft.sh $JOB
\`\`\`

## A1. Translation and fine tuning from the 12M fits' restart files (one 8-GPU job, about 25 min)

The job skips the steps whose products exist (mesh, base + store, bundle, translated restart). To redo the
translations and fine tunings, move \`$R/grain/translated\` aside first; to redo the bundle, also
\`$R/bundle_parent\`.

\`\`\`bash
sbatch $J/refine_12m.sbatch
\`\`\`

## A2. From scratch on the refined mesh (mesh, base, store: about 75 min more)

Move \`$R\` aside (as \`$R.incomplete_<UTC>\`), then \`sbatch $J/refine_12m.sbatch\`.

## B. The 12M fits themselves, from their restart files or from scratch

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_repro--2369abeacf15807199ce76130dcea79c4bfcf68d/HOW-TO-RERUN.md\`
(C30_n2000_grain_cross10, NS31_n2000_mono10_ordered_grain_cross10), then A1, A0.
MD
echo "wrote $O/HOW-TO-RERUN.md"
