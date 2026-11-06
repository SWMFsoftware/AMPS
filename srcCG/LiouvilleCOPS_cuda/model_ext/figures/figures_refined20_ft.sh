#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_refined20_ft.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_refined20_ft.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 060b29b0d62eb97439253f402dadd84fcd27419def730fc68db115f6a75561a5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 20M refinement test.
# RUN: on the login node, no GPU: <this> [<job id>]
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/20M/model_vs_observations_grain_refined_seen_after_gap/
# ==== end of header ====
# Figures of the 20M local refinement test (jobs/refine_20m.sbatch; user request 2026-11-04: the 20M surface seen by NG
# or RG after the gap refined one level (edges halved), the 12M fit with T75 = Tmin = 133 K carried to it and fine-tuned
# there), against the COPS observations (model_obs_figures.py), with the unrefined 20M and the 12M fine tunings, and
# HOW-TO-RERUN.md. Login node, no GPU:
#     figures_refined20_ft.sh [<job id of the refine_20m.sbatch run>]
set -u
JOB=${1:-122011}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R20=$COPS_OUTPUT_ROOT/runs/shap8_020m_flat_earlier
R=$COPS_OUTPUT_ROOT/runs/shap8_020m_refined_late
SRC=C30_n2000_grain_cross10_T75_133K_Tmin_133K; N=${SRC}_from_12M
SF=$R12/grain/fixed_T/$SRC; SFT=$R12/grain/fixed_T/finetune/${SRC}_m5_i5_cpu
UFT=$R20/grain/translated/finetune/${N}_m5_i5_cpu; RFT=$R/grain/translated/finetune/${N}_m5_i5_cpu
O=$COPS_OUTPUT_ROOT/figures_cluster/20M/model_vs_observations_grain_refined_seen_after_gap; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
.venv/bin/python $S/finetune_table.py "$SFT" "$SF" "$TB/12M_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$UFT" "$R20/base" "$TB/20M_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$RFT" "$R/base" "$TB/20M_refined_finetuned.dat" || exit 1
L='cliff $\geq$60$^\circ$ template, T75 and Tmin fixed at 133 K, dust grain in NG + crossing constraint'
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R20/bundle/restart_templates.npz --overwrite \
    --out "$O/refined_20M_vs_20M_vs_12M" --ref "$SF" --title "$L: 12M fit fine-tuned on 12M, 20M and locally refined 20M" \
    --fit "12M fit=$SF:dashed" --fit "fine-tuned on 12M=$TB/12M_finetuned.dat:dotted" \
    --fit "fine-tuned on 20M=$TB/20M_finetuned.dat:dashdot" \
    --fit "fine-tuned on 20M refined after the gap=$TB/20M_refined_finetuned.dat:solid" \
    && echo "wrote refined_20M_vs_20M_vs_12M_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 20M refined over the surface seen after the gap, the 12M fixed-T fit fine-tuned on it

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_refined20_ft.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

What was done (job $JOB; provenance \`$R/provenance/refine-$JOB/\`):
0. \`model_ext/seen_after_gap.py\` (CPU, run with srun before the job): the 20M facets with a record at any of the 95
   points from 10:33:10 UTC in the 20M record store -> \`$R/seen_after_gap.npz\` (1,466,454 facets).
1. \`model_ext/refine_mesh.py --seed-mask ... --margin 0\`: those facets split 1 -> 4 (edge midpoints, geometry
   unchanged), conforming closure: 19,999,504 -> 24,530,882 elements (\`$R/mesh/refine.json\`).
2. Runner prepare (base + store) on the refined mesh, 20M settings (\`$R/pipeline.json\`).
3. \`model_ext/refined_bundle.py\`: the 20M templates carried from the parent elements (gravity provenance kept).
4. \`model_ext/translate_grain.py\` (the 12M fixed-T fit -> the refined base), then \`model_ext/grain_finetune.py\`
   (5 min, 5 intervals, cpu, mixed, max_nfev 2000).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_refined20_ft.sh $JOB
\`\`\`

## A1. Translation and fine tuning from the 12M fit's restart file (one 8-GPU job, about 30 min)

Move \`$R/grain/translated\` aside first (the job skips steps whose products exist), then:

\`\`\`bash
sbatch $J/refine_20m.sbatch
\`\`\`

## A2. From scratch (mesh, base, store: about 2.5 h more)

Move \`$R/mesh\`, \`$R/base\`, \`$R/store\`, \`$R/stamps\`, \`$R/pipeline.json\`, \`$R/bundle_parent\`, \`$R/grain\` aside
(keep \`$R/seen_after_gap.npz\`, or rebuild it:
\`srun -p gpu -A general -N1 -n1 -c 8 --mem=64G --time=00:30:00 $COPS_OUTPUT_ROOT/../AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -W ignore $COPS_OUTPUT_ROOT/model_ext/seen_after_gap.py --run $R20 --out $R/seen_after_gap.npz\`),
then \`sbatch $J/refine_20m.sbatch\`.

## B. The 12M fit and the unrefined 20M fine tuning

\`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_fixed_T/HOW-TO-RERUN.md\` and
\`$COPS_OUTPUT_ROOT/figures_cluster/20M/model_vs_observations_grain_translated_finetune/HOW-TO-RERUN.md\`.
MD
echo "wrote $O/HOW-TO-RERUN.md"
