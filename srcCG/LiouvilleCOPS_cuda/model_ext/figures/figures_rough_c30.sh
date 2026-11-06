#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_rough_c30.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_rough_c30.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8e06dabe19b40428e4e5d405a7b15415ace9b9a8df8ce45f6164c244eb4267f7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M rough-surface C30 fit.
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_rough_surface/
# ==== end of header ====
# Figures of the rough-surface-BC test (user 2026-11-05: "prepare the rough base for 12M on gpu002", "run the best
# grain fit on the rough base"): C30_n2000_grain_cross10 rerun from scratch on the 12M base built with
# --surface-bc rough (job 122033 base, job 122034 fit + fine tuning), against the flat-BC C30 and the COPS
# observations (model_obs_figures.py); writes HOW-TO-RERUN.md. Login node, no GPU:
#     figures_rough_c30.sh [<job id of grain_rough.sbatch>]
set -u
JOB=${1:-122034}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R=$COPS_OUTPUT_ROOT/runs/shap8_012m_rough_earlier
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_rough_surface; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
FIT=C30_n2000_grain_cross10
OF=$R12/grain/fits/$FIT; OFT=$R12/grain/finetune/${FIT}_m5_i5_cpu
RF=$R/grain/fits/$FIT;   RFT=$R/grain/finetune/${FIT}_m5_i5_cpu
.venv/bin/python $S/finetune_table.py "$OFT" "$OF" "$TB/flat_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$RFT" "$RF" "$TB/rough_finetuned.dat" || exit 1
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R12/bundle/restart_templates.npz --overwrite \
    --out "$O/rough_vs_flat_cliff60_free_T" --ref "$OF" \
    --title 'SHAP8 12M, cliff $\geq$60$^\circ$ template, free temperatures, dust grain in NG + crossing constraint: flat vs rough surface emission' \
    --fit "flat surface=$OF:dashed" --fit "flat surface, fine-tuned=$TB/flat_finetuned.dat:dotted" \
    --fit "rough surface=$RF:solid" --fit "rough surface, fine-tuned=$TB/rough_finetuned.dat:dashdot" \
    && echo "wrote $O/rough_vs_flat_cliff60_free_T_*"
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 12M rough-surface emission test (C30 grain fit + fine tuning on a rough-BC base)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_rough_c30.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

Flat surface (all earlier fits): emission cos(theta)/pi (Lambertian). Rough surface (\`SurfaceBC = rough\`,
../LiouvilleCOPS/src/Liouville.h): F/(2 pi), uniform over the hemisphere. Same 12M mesh, same fit command.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_rough_c30.sh $JOB
\`\`\`

## A1. The fit and fine tuning on the existing rough base (one 8-GPU job, about 50 min)

Move \`$R/grain\` aside first (as \`$R/grain.incomplete_<UTC>\`), then

\`\`\`bash
sbatch $J/grain_rough.sbatch
\`\`\`

The fit command is in \`$RF/fit_manifest.json\`; provenance in \`$R/grain/provenance-$JOB/\`.

## A2. From scratch (rough base + record store: about 55 min more on 8 GPUs)

Move \`$R\` aside (as \`$R.incomplete_<UTC>\`), then \`sbatch $J/prepare_12m_rough.sbatch\` (job 122033), then A1.

## B. The flat-surface C30 fit itself

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_repro--2369abeacf15807199ce76130dcea79c4bfcf68d/HOW-TO-RERUN.md\`.
MD
echo "wrote $O/HOW-TO-RERUN.md"
