#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_cycles_Tmin130.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_cycles_Tmin130.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: f73ef9e884befc717f6821985cd4e2e2fa4f630ffba2a4a8367f8af464a85ce4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M cycles 4-6 with Tmin = 130 K, Tmax held.
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles_Tmin130/
# ==== end of header ====
# 12M fit / fine-tuning cycles 4-6 with Tmin held at 130 K and Tmax held at the best fit's 231.982 K (jobs 122117,
# 122140, 122141, jobs/grain_ft_cycle_Tmin130.sbatch; user 2026-11-05: "130 should be the minimum temperature, not a
# temperature everywhere; set the max temperature from the max temperature of the original model"), against the COPS
# observations, with cycle 3 (30 K floor, the start) and the cancelled all-at-130 K cycle 4 (job 122096); writes
# HOW-TO-RERUN.md.   Login node, no GPU:   figures_12m_cycles_Tmin130.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; CY=$R/grain/cycles; CT=$R/grain/cycles_T130; CM=$R/grain/cycles_Tmin130
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles_Tmin130; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
F3=$CY/c3/C30_n4000_grain_cross10_ftc3; FT4=$CT/c4/C30_n4000_T130_grain_cross10_ftc4
.venv/bin/python $S/finetune_table.py "$CY/c3/finetune/C30_n4000_grain_cross10_ftc3_m5_i5_cpu" "$F3" "$TB/c3_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$CT/c4/finetune/C30_n4000_T130_grain_cross10_ftc4_m5_i5_cpu" "$FT4" "$TB/T130_c4_finetuned.dat" || exit 1
for k in 4 5 6; do
    N=C30_n4000_Tmin130_grain_cross10_ftc$k
    .venv/bin/python $S/finetune_table.py "$CM/c$k/finetune/${N}_m5_i5_cpu" "$CM/c$k/$N" "$TB/Tmin130_c${k}_finetuned.dat" || exit 1
done
F6=$CM/c6/C30_n4000_Tmin130_grain_cross10_ftc6
T='SHAP8 12M, cliff $\geq$60$^\circ$ template, dust grain in NG + crossing constraint'
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_ft_cycles_Tmin130_finetuned" --ref "$F6" \
    --title "$T: cycles with Tmin 130 K, Tmax = 232 K held (fine-tuned)" \
    --fit "cycle 3, Tmin $\geq$ 30 K, fine-tuned=$TB/c3_finetuned.dat:dotted" \
    --fit "cycle 4, all T at 130 K, fine-tuned=$TB/T130_c4_finetuned.dat:dashdot" \
    --fit "cycle 4, Tmin 130 K, fine-tuned=$TB/Tmin130_c4_finetuned.dat:dashed" \
    --fit "cycle 5, Tmin 130 K, fine-tuned=$TB/Tmin130_c5_finetuned.dat:dashed" \
    --fit "cycle 6, Tmin 130 K, fine-tuned=$TB/Tmin130_c6_finetuned.dat:solid" \
    && echo "wrote $O/12M_ft_cycles_Tmin130_finetuned_*"
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_ft_cycles_Tmin130_fits" --ref "$F6" \
    --title "$T: fits of the cycles with Tmin 130 K, Tmax = 232 K held" \
    --fit "cycle 3 fit, Tmin $\geq$ 30 K=$F3:dotted" --fit "cycle 4 fit, all T at 130 K=$FT4:dashdot" \
    --fit "cycle 4 fit, Tmin 130 K=$CM/c4/C30_n4000_Tmin130_grain_cross10_ftc4:dashed" \
    --fit "cycle 5 fit, Tmin 130 K=$CM/c5/C30_n4000_Tmin130_grain_cross10_ftc5:dashed" \
    --fit "cycle 6 fit, Tmin 130 K=$F6:solid" \
    && echo "wrote $O/12M_ft_cycles_Tmin130_fits_*"
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 12M fit / fine-tuning cycles 4-6 with Tmin = 130 K and Tmax = 231.982 K held

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_cycles_Tmin130.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

Cycle k = 4, 5, 6 (\`jobs/grain_ft_cycle_Tmin130.sbatch\`, CYCLE=k; jobs 122117, 122140, 122141; test 122116 in
\`$R/grain/cycles_Tmin130_test-122116\`): the scheme of \`jobs/grain_ft_cycle.sbatch\` (previous fine tuning's
correction -> template ft_c<k>, ft_c1..ft_c<k> held at 1; grain start from the fine tuning; max_nfev 4000; fine
tuning 5 min, 5 intervals, cpu) with \`--tbounds 30 400 --fix Tmin=130 --fix Tmax=231.9820451565314\` (Tmax of
\`$R/grain/fits/C30_n2000_grain_cross10\`); T75 free. Cycle 4 starts from \`$F3\` (cycle 3, 30 K floor) and its
fine tuning. Provenance: \`$CM/c<k>/provenance-<job>/\`.
Cycle 5 needs \`model_ext/fixed_parameters.py\` as fixed on 2026-11-06 (the first attempt, job 122118, failed in
reg_jacobian; in \`$CM/c5.incomplete_20261106T003556Z\`).
Shown for comparison: cycle 3 (\`figures_12m_cycles.sh\`) and the all-temperatures-130 K cycle 4
(\`jobs/grain_ft_cycle_T130.sbatch\`, CYCLE=4, job 122096; its cycles 5-6 were cancelled).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_cycles_Tmin130.sh
\`\`\`

## A1. The cycles (one 8-GPU job each, about 80-90 min)

Move \`$CM\` aside first (as \`$CM.incomplete_<UTC>\`), then

\`\`\`bash
c4=\$(sbatch --parsable --export=ALL,CYCLE=4 $J/grain_ft_cycle_Tmin130.sbatch)
c5=\$(sbatch --parsable --dependency=afterok:\$c4 --export=ALL,CYCLE=5 $J/grain_ft_cycle_Tmin130.sbatch)
sbatch --dependency=afterok:\$c5 --export=ALL,CYCLE=6 $J/grain_ft_cycle_Tmin130.sbatch
\`\`\`

## B. Cycles 1-3 and the best 12M fit (the start of cycle 4)

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles/HOW-TO-RERUN.md\`.
The all-130 K cycle 4: \`sbatch --export=ALL,CYCLE=4 $J/grain_ft_cycle_T130.sbatch\` (after moving \`$CT/c4\` aside).
MD
echo "wrote $O/HOW-TO-RERUN.md"
