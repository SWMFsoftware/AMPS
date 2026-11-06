#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_44m_cases.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_44m_cases.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5a423e8a90369f8d8985d381d67ae094b4bedb4c0c174e81cf5dc39abe4bbf6e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 44M grain cases cliff >=40/60 + sin(slope).
# RUN: on the login node, no GPU: <this> [<job id cliff40> <job id cliff60>]
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/44M/model_vs_observations_grain_cases/
# ==== end of header ====
# Figures of the 44M grain cases (model_ext/grain_44m_case.sbatch; user request 2026-11-04: "separate cases for >=40,
# >=60 plus sin; add"; the 12M fit with T75 = Tmin = 133 K translated to 44M, refitted with the cliff and sin(slope)
# templates, then fine-tuned), against the COPS observations (model_obs_figures.py), with the 12M fine tuning of the
# translated fit as reference, and HOW-TO-RERUN.md. Login node, no GPU:
#     figures_44m_cases.sh [<job id cliff40> <job id cliff60>]
set -u
J40=${1:-121999}; J60=${2:-122000}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; X=$COPS_OUTPUT_ROOT/model_ext
R12=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; R=$COPS_OUTPUT_ROOT/runs/shap8_044m_flat_earlier
SRC=C30_n2000_grain_cross10_T75_133K_Tmin_133K; SF=$R12/grain/fixed_T/$SRC; SFT=$R12/grain/fixed_T/finetune/${SRC}_m5_i5_cpu
XL=$R/grain/translated/${SRC}_from_12M; C=$R/grain/cases
O=$COPS_OUTPUT_ROOT/figures_cluster/44M/model_vs_observations_grain_cases; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
.venv/bin/python $S/finetune_table.py "$SFT" "$SF" "$TB/12M_finetuned.dat" || exit 1
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite"
GR='T75 and Tmin fixed at 133 K, dust grain in NG + crossing constraint'
declare -A D=([C40S]='cliff $\geq$40$^\circ$ + sin(slope) templates' [C60S]='cliff $\geq$60$^\circ$ + sin(slope) templates')
declare -A N=([C40S]=cliff40_sin_slope [C60S]=cliff60_sin_slope)
BEST=()
for M in C40S C60S; do
    F=$C/${M}_n2000_grain_cross10_T75_133K_Tmin_133K_from12M; FT=$C/finetune/$(basename "$F")_m5_i5_cpu
    T=$TB/44M_${N[$M]}_finetuned.dat
    .venv/bin/python $S/finetune_table.py "$FT" "$F" "$T" || exit 1
    $P --out "$O/44M_${N[$M]}" --ref "$F" --title "SHAP8 44M, ${D[$M]}, $GR: fitted and fine-tuned" \
        --fit "44M fit=$F:dashed" --fit "44M fit, fine-tuned=$T:solid" \
        --fit "12M fit (cliff \$\\geq\$60\$^\\circ\$), fine-tuned on 12M=$TB/12M_finetuned.dat:dotted" \
        && echo "wrote 44M_${N[$M]}_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
    BEST+=("${D[$M]}, fine-tuned=$T:$( [ ${#BEST[@]} -eq 0 ] && echo solid || echo dashed )")
done
$P --out "$O/44M_cases_finetuned" --ref "$C/C40S_n2000_grain_cross10_T75_133K_Tmin_133K_from12M" \
    --title "SHAP8 44M, $GR: the two template cases, fine-tuned" "${BEST[@]/#/--fit=}" \
    --fit "12M fit (cliff \$\\geq\$60\$^\\circ\$), fine-tuned on 12M=$TB/12M_finetuned.dat:dotted" \
    && echo "wrote 44M_cases_finetuned_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 44M grain cases (cliff >=40 / >=60 deg + sin(slope), T75 = Tmin = 133 K, from the 12M fit)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_44m_cases.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

Jobs: $J40 (cliff >=40 + sin(slope)), $J60 (cliff >=60 + sin(slope)); provenance \`$C/provenance-<job>/\`.
Each: \`model_ext/grain_fit_fixed.py --fix T75=133 --fix Tmin=133\` (crossing weight 10, mixed, 2000 nfev,
--restart-from the translated 44M state \`$XL/fit_restart.npz\`, template start 1), then \`model_ext/grain_finetune.py\`
(5 min, 5 intervals, cpu). The exact fit commands: \`fit_manifest.json\` ["model_extension"]["command"] of each fit dir.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_44m_cases.sh $J40 $J60
\`\`\`

## A1. A case from the translated restart (one 8-GPU job, about 3 h)

Move the case's fit and fine-tuning directories (\`$C/<case>\`, \`$C/finetune/<case>_m5_i5_cpu\`) aside first.

\`\`\`bash
sbatch --export=ALL,COPS_RUN=$R,XL=$XL,CASE=cliff40 $X/grain_44m_case.sbatch
sbatch --export=ALL,COPS_RUN=$R,XL=$XL,CASE=cliff60 $X/grain_44m_case.sbatch
\`\`\`

## A2. The translation (12M restart -> 44M, about 20 min), then A1

Move \`$XL\` aside, then
\`sbatch --export=ALL,SRC_FIT=$SF,COPS_RUN=$R,NAME=${SRC}_from_12M,NO_FT=1 $X/grain_translate_finetune.sbatch\`.

## B. From scratch

The 44M base, store, visibility, gravity (genus 6, 6 GPUs), cache and bundles:
\`sbatch $COPS_OUTPUT_ROOT/jobs/prepare_44m.sbatch\` (move \`$R\` aside first; about 7 h). The 12M fit:
\`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_fixed_T/HOW-TO-RERUN.md\`. Then A2, A1, A0.
MD
echo "wrote $O/HOW-TO-RERUN.md"
