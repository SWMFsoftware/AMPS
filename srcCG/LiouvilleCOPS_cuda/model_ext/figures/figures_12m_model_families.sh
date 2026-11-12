#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_model_families.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_model_families.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 95ad9d6ce329fdb049a2565003180a893259436c90e3947b63ea2a7ab732cbc4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The four 12M activity-model families of 2026-11-08 ... 11-10 (clipped refit D2; positive ratio
#   1e-5 lambda 1e-4 after 3 cycles; slope multipliers x harmonics; bumps x slope multipliers after 3 cycles) against
#   the COPS observations (user, 2026-11-11 04:05 UTC: "make comparison plot for those modes"); summary.txt,
#   HOW-TO-RERUN.md (pointers to each family's own), code/.
# RUN: figures_12m_model_families.sh   (login node, no GPU) -> figures_cluster/12M/model_vs_observations_model_families/
# ==== end of header ====
# The four 12M activity-model families of 2026-11-08 ... 11-10 against the COPS observations in one plot set (user
# 2026-11-10: "make comparison plot for those modes"), lategrid base (gauge +1.5 m, W10x5), best result of each family:
#   clipped refit D2 (cliff + 6 templates, A = max(0, ...))                     runs/.../clip_refit/D2/final
#   positive activity, ratio 1e-5 + CO2 floor, lambda 1e-4, after 3 cycles      runs/.../positive_r1e-5/lam1e-4/c3/ft
#   slope multipliers x 25 harmonics, floor 1e-5 A_orig (srcCG start)          runs/.../slope_mult/fit/final
#   32 bumps x slope multipliers, floor 1e-5 A_orig, after 3 cycles (3a)        runs/.../slope_bumps/chain_.../c3/ft/final
# Writes figures_cluster/12M/model_vs_observations_model_families/: 12M_model_families_*, summary.txt (cost, rms,
# termination, kappa, grain of each), HOW-TO-RERUN.md (pointers to each model's own HOW-TO-RERUN), code/.
#     Login node, no GPU:   figures_12m_model_families.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
F12=$COPS_OUTPUT_ROOT/figures_cluster/12M
D2=$RUN/clip_refit/D2/final
PB=$RUN/positive_r1e-5/lam1e-4/c3/ft;            PBF=$RUN/positive_r1e-5/lam1e-4/c3/fit/final
SM=$RUN/slope_mult/fit/final
CH=$RUN/slope_bumps/chain_slopemul_20261110T201212Z; BF=$CH/c3/ft/final; BFF=$CH/c3/fit/final
O=$F12/model_vs_observations_model_families; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C/figures"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
$PY $S/positive_table.py "$D2" "$TB/clip_D2.dat" || exit 1
$PY $S/positive_table.py "$PB" "$TB/positive_r1e-5_lam1e-4_c3.dat" || exit 1
$PY $S/positive_table.py "$SM" "$TB/slope_mult_harmonics.dat" || exit 1
$PY $S/positive_table.py "$BF" "$TB/slope_bumps_c3.dat" || exit 1
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_model_families" \
    --ref "$TB/clip_D2.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): best fit of each activity-model family" \
    --fit "clipped refit D2, data cost 4.89e-3=$TB/clip_D2.dat:solid" \
    --fit "positive, ratio 1e-5, lambda 1e-4, 3 cycles, 8.54e-3=$TB/positive_r1e-5_lam1e-4_c3.dat:dashed" \
    --fit "slope multipliers x harmonics, 8.72e-3=$TB/slope_mult_harmonics.dat:dotted" \
    --fit "bumps x slope multipliers, 3 cycles, 1.41e-2=$TB/slope_bumps_c3.dat:dashdot" || exit 3
fitrow() {  # name, fit dir (fit_result.txt)
    local f=$2/fit_result.txt
    printf "%-34s %-13s %-8s %-8s %-9s %-13s %-10s\n" "$1" $(awk '/^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) \
        "$(grep -oE 'ftol|maximum number' $f | head -1 | sed 's/maximum number/max_nfev/')/$(grep -oE 'nfev [0-9]+' $f | head -1 | cut -d' ' -f2)" \
        "$(awk '/^kappa /{print $3; exit}' $f)" "$(awk '/^grain_A_uPa /{print $3; exit}' $f)"
}
ftrow() {   # name, fine-tuning dir (finetune_result.txt), its fit dir (kappa)
    local f=$2/finetune_result.txt
    printf "%-34s %-13s %-8.4f %-8.4f %-9s %-13s %-10s\n" "$1" $(awk '/^cost_data/{print $5}' $f) $(awk '/^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $f) \
        "$(awk '/^stop /{print $2}' $f)/$(awk '/^nfev_total/{print $2}' $f)" "$(awk '/^kappa /{print $3; exit}' $3/fit_result.txt)" \
        "$(awk '/^grain_A_uPa/{print $5; exit}' $f)"
}
{
echo "# best 12M fit of each activity-model family, lategrid base (gauge +1.5 m, W10x5); compare by cost, termination, active bounds"
echo "# kappa = CO2 scale of the last fit of the family (lower bound 0); grain A in uPa (fine tunings: refitted value)"
printf "%-34s %-13s %-8s %-8s %-9s %-13s %-10s\n" "# model" data_cost rms_NG rms_RG stop/nfev kappa grain_A
fitrow "clipped_refit_D2" "$D2"
ftrow  "positive_r1e-5_lam1e-4_c3ft" "$PB" "$PBF"
fitrow "slope_mult_x_harmonics" "$SM"
ftrow  "bumps_x_slope_mult_c3ft" "$BF" "$BFF"
} > "$O/summary.txt"
cp -p $S/figures_12m_model_families.sh $S/positive_table.py $S/model_obs_figures.py "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 12M model-family comparison (figures only)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_model_families.sh\` (login node, no GPU):

\`\`\`bash
$S/figures_12m_model_families.sh
\`\`\`

The four fits and the commands that produced them (each with its own HOW-TO-RERUN.md, code/ and AMPS commit):

| model | products | how to rerun |
|---|---|---|
| clipped refit D2 | \`$D2\` | \`$F12/model_vs_observations_clip_vs_positive/HOW-TO-RERUN.md\` |
| positive, ratio 1e-5, lambda 1e-4, 3 cycles | \`$PB\` | \`$F12/model_vs_observations_positive_cycles_r1e-5/HOW-TO-RERUN.md\` |
| slope multipliers x harmonics | \`$SM\` | \`$F12/model_vs_observations_slope_mult/HOW-TO-RERUN.md\` |
| bumps x slope multipliers, 3 cycles | \`$BF\` | \`$F12/model_vs_observations_slope_bumps/HOW-TO-RERUN.md\` |

Figure code that ran: \`code/figures/\` (SHA256SUMS).
MD
echo "wrote $O"; cat "$O/summary.txt"
