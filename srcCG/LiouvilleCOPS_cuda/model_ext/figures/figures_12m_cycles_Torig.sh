#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_cycles_Torig.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_cycles_Torig.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 9d014cd894c8902da4cdd3e63f1024d18c11fb09aec8da971a6cbff08bad2980
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M cycles 4-6 with the original temperature law held (also writes code/).
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles_Torig/
# ==== end of header ====
# 12M fit / fine-tuning cycles 4-6 with the original model's temperature law held (Tmax 182.09922, T75 136.11149,
# Tmin 133 K; jobs 122219-122221, jobs/grain_ft_cycle_Torig.sbatch; user 2026-11-06), against the COPS observations,
# with cycle 3 (30 K floor, free T: the start) and the Tmin = 130 K cycle 6; writes HOW-TO-RERUN.md (with the AMPS
# commit of each fit) and code/ (a copy of the code that ran; user 2026-11-06).   Login node, no GPU:
#     figures_12m_cycles_Torig.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; CY=$R/grain/cycles; CO=$R/grain/cycles_Torig; CM=$R/grain/cycles_Tmin130
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles_Torig; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
F3=$CY/c3/C30_n4000_grain_cross10_ftc3; FM6=$CM/c6/C30_n4000_Tmin130_grain_cross10_ftc6
.venv/bin/python $S/finetune_table.py "$CY/c3/finetune/C30_n4000_grain_cross10_ftc3_m5_i5_cpu" "$F3" "$TB/c3_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$CM/c6/finetune/C30_n4000_Tmin130_grain_cross10_ftc6_m5_i5_cpu" "$FM6" "$TB/Tmin130_c6_finetuned.dat" || exit 1
for k in 4 5 6; do
    N=C30_n4000_Torig_grain_cross10_ftc$k
    .venv/bin/python $S/finetune_table.py "$CO/c$k/finetune/${N}_m5_i5_cpu" "$CO/c$k/$N" "$TB/Torig_c${k}_finetuned.dat" || exit 1
done
F6=$CO/c6/C30_n4000_Torig_grain_cross10_ftc6
T='SHAP8 12M, cliff $\geq$60$^\circ$ template, dust grain in NG + crossing constraint'
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_ft_cycles_Torig_finetuned" --ref "$F6" \
    --title "$T: cycles with the original temperature law held (fine-tuned)" \
    --fit "cycle 3, free T (Tmin $\geq$ 30 K), fine-tuned=$TB/c3_finetuned.dat:dotted" \
    --fit "cycle 6, Tmin 130 K + Tmax 232 K held, fine-tuned=$TB/Tmin130_c6_finetuned.dat:dashdot" \
    --fit "cycle 4, original T, fine-tuned=$TB/Torig_c4_finetuned.dat:dashed" \
    --fit "cycle 5, original T, fine-tuned=$TB/Torig_c5_finetuned.dat:dashed" \
    --fit "cycle 6, original T, fine-tuned=$TB/Torig_c6_finetuned.dat:solid" \
    && echo "wrote $O/12M_ft_cycles_Torig_finetuned_*"
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_ft_cycles_Torig_fits" --ref "$F6" \
    --title "$T: fits of the cycles with the original temperature law held" \
    --fit "cycle 3 fit, free T (Tmin $\geq$ 30 K)=$F3:dotted" --fit "cycle 6 fit, Tmin 130 K + Tmax 232 K held=$FM6:dashdot" \
    --fit "cycle 4 fit, original T=$CO/c4/C30_n4000_Torig_grain_cross10_ftc4:dashed" \
    --fit "cycle 5 fit, original T=$CO/c5/C30_n4000_Torig_grain_cross10_ftc5:dashed" \
    --fit "cycle 6 fit, original T=$F6:solid" \
    && echo "wrote $O/12M_ft_cycles_Torig_fits_*"
# code copy (provenance copies of what ran) and figure scripts
rows=""
for d in $CO/c4/provenance-122219 $CO/c5/provenance-122220 $CO/c6/provenance-122221 $CY/c3/provenance-122082 $CM/c6/provenance-122141; do
    job=${d##*-}; lab=$(basename $(dirname $d)); case $d in *cycles_Torig*) set=Torig;; *cycles_Tmin130*) set=Tmin130;; *) set=cycles;; esac
    mkdir -p "$C/${set}_${lab}-$job"; cp -p "$d"/*.sbatch "$d"/*.py "$d"/model_ext.sha256 "$d"/git.txt "$C/${set}_${lab}-$job/"
    rows+="| ${set} ${lab} | $job | \`$(head -1 $d/git.txt)\` | $(tail -n +2 $d/git.txt | grep -v '^??' | tr '\n' ' ') |"$'\n'
done
mkdir -p "$C/figures"; cp -p $S/figures_12m_cycles_Torig.sh $S/model_obs_figures.py $S/finetune_table.py "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 12M fit / fine-tuning cycles 4-6 with the original temperature law held

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_cycles_Torig.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

## Code version

| Plotted case | Job | AMPS commit | Tracked local changes |
|---|---|---|---|
$rows
The tracked LiouvilleCOPS/LiouvilleCOPS_cuda code is the same in 2369abe and dc5cbce (dc5cbce, committed by the user
2026-11-06 21:56 UTC during cycle 6, adds only the model_ext copies and srcMars changes: \`git diff --stat 2369abe dc5cbce\`).

Untracked code that ran (\`$COPS_OUTPUT_ROOT/model_ext\`, \`jobs\`) is copied in \`code/\` (from each job's provenance, with
\`model_ext.sha256\` and \`git.txt\`; the figure scripts; \`code/SHA256SUMS\`). Version-controlled copies (with header
comments): \`srcCG/LiouvilleCOPS_cuda/model_ext/\` (README.md, RUN_COPIES.sha256, check_copies.sh).

## What was done

Temperature law of the original model (\`srcCG/LiouvilleCOPS/src/Source.cpp:159-161\` defaults, TemperatureModel =
analytical): lit (SZA <= 75 deg) T = Tmax + KT (Tmax - T75)(1/cos SZA - 1), KT = 1/(1 - 1/cos 75 deg) = -0.349,
floored at Tmin; SZA > 75 deg, night, shadow: Tmin; Tmax 182.09922, T75 136.11149, Tmin 133 K, all three held
(\`--fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0\`, --tbounds 30 400).
Cycle k = 4, 5, 6 (\`jobs/grain_ft_cycle_Torig.sbatch\`, CYCLE=k; test 122218 in \`$R/grain/cycles_Torig_test-122218\`):
the scheme of \`jobs/grain_ft_cycle.sbatch\` (previous fine tuning's correction -> template ft_c<k>, ft_c1..ft_c<k>
held at 1; grain start from the fine tuning; max_nfev 4000; fine tuning 5 min, 5 intervals, cpu). Cycle 4 starts from
\`$F3\` (cycle 3) and its fine tuning. Provenance: \`$CO/c<k>/provenance-<job>/\`.
Shown for comparison: cycle 3 and the Tmin = 130 K + Tmax = 231.98 K cycle 6.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_cycles_Torig.sh
\`\`\`

## A1. The cycles (one 8-GPU job each, 72-90 min)

Move \`$CO\` aside first (as \`$CO.incomplete_<UTC>\`), then

\`\`\`bash
c4=\$(sbatch --parsable --export=ALL,CYCLE=4 $J/grain_ft_cycle_Torig.sbatch)
c5=\$(sbatch --parsable --dependency=afterok:\$c4 --export=ALL,CYCLE=5 $J/grain_ft_cycle_Torig.sbatch)
sbatch --dependency=afterok:\$c5 --export=ALL,CYCLE=6 $J/grain_ft_cycle_Torig.sbatch
\`\`\`

## B. Cycles 1-3 and the best 12M fit (the start of cycle 4)

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles/HOW-TO-RERUN.md\`, then A1, A0.
The Tmin = 130 K cycles: \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles_Tmin130/HOW-TO-RERUN.md\`.
MD
echo "wrote $O/HOW-TO-RERUN.md"
