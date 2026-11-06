#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_cycles.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_cycles.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 18d729b4456a53abde556b0abd507c731c6f9b0dde3cd3163af4e276ae1f4a5a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M fit/fine-tuning cycles 1-3.
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles/
# ==== end of header ====
# 12M fit / fine-tuning cycles from the best 12M grain fit (jobs 122080-122082, jobs/grain_ft_cycle.sbatch; user
# 2026-11-05: "restart after the fine tuning, so 4000 iteration fitting, then fine tuning ... 3 times"; corrections
# as templates held at 1, accumulated; grain start fine-tuned), against the COPS observations; writes HOW-TO-RERUN.md.
# Login node, no GPU:   figures_12m_cycles.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; CY=$R/grain/cycles
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
F0=$R/grain/fits/C30_n2000_grain_cross10; FT0=$R/grain/finetune/C30_n2000_grain_cross10_m5_i5_cpu
.venv/bin/python $S/finetune_table.py "$FT0" "$F0" "$TB/c0_finetuned.dat" || exit 1
for k in 1 2 3; do
    N=C30_n4000_grain_cross10_ftc$k
    .venv/bin/python $S/finetune_table.py "$CY/c$k/finetune/${N}_m5_i5_cpu" "$CY/c$k/$N" "$TB/c${k}_finetuned.dat" || exit 1
done
F3=$CY/c3/C30_n4000_grain_cross10_ftc3
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_ft_cycles_finetuned" --ref "$F3" \
    --title 'SHAP8 12M, cliff $\geq$60$^\circ$ template, free temperatures, dust grain in NG + crossing constraint: fit / fine-tuning cycles' \
    --fit "best 12M fit, fine-tuned=$TB/c0_finetuned.dat:dotted" \
    --fit "cycle 1, fine-tuned=$TB/c1_finetuned.dat:dashdot" \
    --fit "cycle 2, fine-tuned=$TB/c2_finetuned.dat:dashed" \
    --fit "cycle 3, fine-tuned=$TB/c3_finetuned.dat:solid" \
    && echo "wrote $O/12M_ft_cycles_finetuned_*"
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_ft_cycles_fits" --ref "$F3" \
    --title 'SHAP8 12M, cliff $\geq$60$^\circ$ template, free temperatures, dust grain in NG + crossing constraint: fits of the cycles' \
    --fit "best 12M fit=$F0:dotted" --fit "cycle 1 fit=$CY/c1/C30_n4000_grain_cross10_ftc1:dashdot" \
    --fit "cycle 2 fit=$CY/c2/C30_n4000_grain_cross10_ftc2:dashed" --fit "cycle 3 fit=$F3:solid" \
    && echo "wrote $O/12M_ft_cycles_fits_*"
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 12M fit / fine-tuning cycles from the best 12M grain fit

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_cycles.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

Cycle k (\`jobs/grain_ft_cycle.sbatch\`, CYCLE=k; jobs 122080, 122081, 122082; test 122079 in
\`$R/grain/cycles_test-122079\`): the previous fine tuning's activity correction becomes template ft_c<k>
(model_ext/ft_template.py, same mesh), held at amplitude 1 with ft_c1..ft_c<k-1>; the previous fit's restart with
the fine-tuned grain start (model_ext/restart_grain_start.py); the C30 fit command with these templates,
--restart-from, --template-start 0, max_nfev 4000 (model_ext/grain_fit_fixed.py); then grain_finetune.py (5 min,
5 intervals, cpu). Cycle 1 starts from \`$F0\` and \`$FT0\`. Provenance: \`$CY/c<k>/provenance-<job>/\`.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_cycles.sh
\`\`\`

## A1. The cycles (one 8-GPU job each, 25-80 min)

Move \`$CY\` aside first (as \`$CY.incomplete_<UTC>\`), then

\`\`\`bash
c1=\$(sbatch --parsable --export=ALL,CYCLE=1 $J/grain_ft_cycle.sbatch)
c2=\$(sbatch --parsable --dependency=afterok:\$c1 --export=ALL,CYCLE=2 $J/grain_ft_cycle.sbatch)
sbatch --dependency=afterok:\$c2 --export=ALL,CYCLE=3 $J/grain_ft_cycle.sbatch
\`\`\`

## B. The best 12M fit and its fine tuning (the start of cycle 1)

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_repro--2369abeacf15807199ce76130dcea79c4bfcf68d/HOW-TO-RERUN.md\`.
MD
echo "wrote $O/HOW-TO-RERUN.md"
