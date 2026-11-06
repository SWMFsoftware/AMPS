#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_44m_resolutions.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_44m_resolutions.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: c04923a60323df0972e2d59f69daff3f78aad066f8b4d9073d688b1a592da31e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: fixed-T fits fine-tuned on 12M, 20M, 20M refined, 44M.
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/44M/model_vs_observations_grain_resolutions/
# ==== end of header ====
# Cross-resolution comparison of the fixed-T (T75 = Tmin = 133 K) grain fits after fine tuning, against the COPS
# observations (user 2026-11-05: "prepare comparison plots for 44M case"): fine-tuned on 12M, on 20M, on 20M refined
# after the gap, and the two 44M template cases. Reads the tables written by figures_refined20_ft.sh and
# figures_44m_cases.sh (run those first). Login node, no GPU:   figures_44m_resolutions.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
F=$COPS_OUTPUT_ROOT/figures_cluster; S=$F/scripts
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
T20=$F/20M/model_vs_observations_grain_refined_seen_after_gap/tables; T44=$F/44M/model_vs_observations_grain_cases/tables
O=$F/44M/model_vs_observations_grain_resolutions; mkdir -p "$O"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
for f in $T20/12M_finetuned.dat $T20/20M_finetuned.dat $T20/20M_refined_finetuned.dat \
         $T44/44M_cliff60_sin_slope_finetuned.dat $T44/44M_cliff40_sin_slope_finetuned.dat; do
    [ -s "$f" ] || { echo "missing $f (run figures_refined20_ft.sh and figures_44m_cases.sh first)"; exit 1; }
done
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R12/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_20M_44M_finetuned" --ref "$R12/grain/fixed_T/C30_n2000_grain_cross10_T75_133K_Tmin_133K" \
    --title 'SHAP8, dust grain in NG + crossing constraint, T75 = Tmin = 133 K: fine-tuned fits on 12M, 20M and 44M' \
    --fit "12M, cliff \$\\geq\$60\$^\\circ\$=$T20/12M_finetuned.dat:dotted" \
    --fit "20M, cliff \$\\geq\$60\$^\\circ\$=$T20/20M_finetuned.dat:dashdot" \
    --fit "20M refined after the gap, cliff \$\\geq\$60\$^\\circ\$=$T20/20M_refined_finetuned.dat:dashed" \
    --fit "44M, cliff \$\\geq\$60\$^\\circ\$ + sin(slope)=$T44/44M_cliff60_sin_slope_finetuned.dat:solid" \
    --fit "44M, cliff \$\\geq\$40\$^\\circ\$ + sin(slope)=$T44/44M_cliff40_sin_slope_finetuned.dat:solid" \
    && echo "wrote $O/12M_20M_44M_finetuned_*"
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: fine-tuned fits on 12M, 20M, 20M refined and 44M (T75 = Tmin = 133 K)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_44m_resolutions.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_refined20_ft.sh      # writes the 12M, 20M and 20M-refined tables
$S/figures_44m_cases.sh 121999 122000   # writes the 44M tables
$S/figures_44m_resolutions.sh
\`\`\`

## A1. The fits and fine tunings themselves

* 12M fixed-T fit and fine tuning, 20M and 20M refined (job 122011):
  \`$F/20M/model_vs_observations_grain_refined_seen_after_gap/HOW-TO-RERUN.md\`
* 44M cases (preparation 122020, translation 121998, cases 121999 and 122000):
  \`$F/44M/model_vs_observations_grain_cases/HOW-TO-RERUN.md\`
MD
echo "wrote $O/HOW-TO-RERUN.md"
