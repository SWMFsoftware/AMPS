#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_late_weights_ramp.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_late_weights_ramp.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: ac2c5313d2e2424f1e27c19a03cb3db7fc177abcb244a057d8e6c982556e99ef
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M fits with the last 6 points weighted 1 -> 20, from cycle 3 and from the original-T cycle 6 (HOW-TO-RERUN.md with full commands via rerun_commands.py; code/).
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_late_weights_ramp/
# ==== end of header ====
# 12M fits with the last 6 points weighted linearly from 1 (6th from the end) to 20 (last), NG and RG (user
# 2026-11-06): from cycle 3 (free T; job 122248) and from the original-temperature cycle 6 (Tmax/T75/Tmin held; job
# 122250), with their starting fits, against the COPS observations. Writes HOW-TO-RERUN.md with every command and all
# arguments (rerun_commands.py), the AMPS commit and the commits of the code, and code/ (the code that ran).
# Login node, no GPU:   figures_12m_late_weights_ramp.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; CY=$R/grain/cycles; CO=$R/grain/cycles_Torig; LW=$R/grain/late_weights
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_late_weights_ramp; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
F3=$CY/c3/C30_n4000_grain_cross10_ftc3; FT3=$CY/c3/finetune/C30_n4000_grain_cross10_ftc3_m5_i5_cpu
F6=$CO/c6/C30_n4000_Torig_grain_cross10_ftc6; FT6=$CO/c6/finetune/C30_n4000_Torig_grain_cross10_ftc6_m5_i5_cpu
NA=C30_n4000_grain_cross10_c3_last6_w1to20_linear_ngrg; NB=C30_n4000_grain_cross10_torigc6_last6_w1to20_linear_ngrg
.venv/bin/python $S/finetune_table.py "$FT3" "$F3" "$TB/c3_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$FT6" "$F6" "$TB/torig_c6_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$LW/$NA/finetune/${NA}_m5_i5_cpu" "$LW/$NA/fit" "$TB/c3_ramp_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$LW/$NB/finetune/${NB}_m5_i5_cpu" "$LW/$NB/fit" "$TB/torig_c6_ramp_finetuned.dat" || exit 1
T='SHAP8 12M, cliff $\geq$60$^\circ$, grain + crossing: last 6 points weighted 1$\rightarrow$20 (NG, RG)'
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_last6_ramp_finetuned" --ref "$LW/$NA/fit" --title "$T, fine-tuned" \
    --fit "cycle 3 (free T), fine-tuned=$TB/c3_finetuned.dat:dotted" \
    --fit "cycle 3 + ramp weights, fine-tuned=$TB/c3_ramp_finetuned.dat:solid" \
    --fit "original T cycle 6, fine-tuned=$TB/torig_c6_finetuned.dat:dashdot" \
    --fit "original T cycle 6 + ramp weights, fine-tuned=$TB/torig_c6_ramp_finetuned.dat:dashed" \
    && echo "wrote $O/12M_last6_ramp_finetuned_*"
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_last6_ramp_fits" --ref "$LW/$NA/fit" --title "$T, fits" \
    --fit "cycle 3 fit (free T)=$F3:dotted" --fit "cycle 3 + ramp weights, fit=$LW/$NA/fit:solid" \
    --fit "original T cycle 6 fit=$F6:dashdot" --fit "original T cycle 6 + ramp weights, fit=$LW/$NB/fit:dashed" \
    && echo "wrote $O/12M_last6_ramp_fits_*"
# code copy: the job scripts and modules as they ran (provenance copies), the figure scripts
for d in $LW/$NA/provenance-122248 $LW/$NB/provenance-122250 $CY/c3/provenance-122082 $CO/c6/provenance-122221; do
    k=$(basename $(dirname $d))-${d##*-}; mkdir -p "$C/$k"
    cp -p "$d"/*.sbatch "$d"/*.py "$d"/model_ext.sha256 "$d"/git.txt "$C/$k/"
done
mkdir -p "$C/figures"; cp -p $S/figures_12m_late_weights_ramp.sh $S/model_obs_figures.py $S/finetune_table.py $S/rerun_commands.py "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
RC=".venv/bin/python $S/rerun_commands.py"
{
cat <<MD
# How to rerun: 12M fits with the last 6 points weighted linearly 1 -> 20 (NG and RG)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_late_weights_ramp.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

Weights: the cost terms of the last 6 used NG and RG points (10:39:06-10:39:26 UTC) multiplied by W = 1, 4.8, 8.6,
12.4, 16.2, 20 (residuals by sqrt W); \`--weight-last 6 --weight-profile linear --weight-start 1 --weight-factor 20\`.
Reported costs, rms and restart files are unweighted. Each fit starts from its start fit's restart with the start
fine tuning's grain (restart_grain_start.py); the correction templates of the start fit are held at 1.

The commands below are complete (every argument as run), so they do not depend on the current job scripts. Code:
the "commit" column names the git commit (in \`srcCG/LiouvilleCOPS_cuda/model_ext\`) whose version-controlled copy
is identical to the code that ran ("not in git": only in \`code/\` here); \`code/\` has the files as they ran, with
\`SHA256SUMS\`. Job scripts for the same runs: \`$J/grain_late_weights.sbatch\` (START=c3 or torig_c6, WLAST=6,
WFACTOR=20, WPROFILE=linear, WSTART=1), not needed for the commands below.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_late_weights_ramp.sh
\`\`\`

## A1. The two weighted fits and their fine tunings (each on one GPU node; inside a Slurm job)

MD
$RC --fit "$LW/$NA/fit" --finetune "$LW/$NA/finetune/${NA}_m5_i5_cpu" --provenance "$LW/$NA/provenance-122248" \
    --restart-start "$LW/$NA/restart_start.npz" --label "From cycle 3, free T (job 122248, gpu002, 3 GPUs)"
$RC --fit "$LW/$NB/fit" --finetune "$LW/$NB/finetune/${NB}_m5_i5_cpu" --provenance "$LW/$NB/provenance-122250" \
    --restart-start "$LW/$NB/restart_start.npz" --label "From the original-temperature cycle 6, Tmax/T75/Tmin held (job 122250, gpu001, 8 GPUs)"
cat <<MD
## A2. The start fits (cycle 3, original-temperature cycle 6)

MD
$RC --fit "$F3" --finetune "$FT3" --provenance "$CY/c3/provenance-122082" --restart-start "$CY/c3/restart_start.npz" \
    --template "ft_c3=$CY/c2/finetune/C30_n4000_grain_cross10_ftc2_m5_i5_cpu,$CY/c2/C30_n4000_grain_cross10_ftc2/fit_restart.npz" \
    --label "Cycle 3 (job 122082)"
$RC --fit "$F6" --finetune "$FT6" --provenance "$CO/c6/provenance-122221" --restart-start "$CO/c6/restart_start.npz" \
    --template "ft_c6=$CO/c5/finetune/C30_n4000_Torig_grain_cross10_ftc5_m5_i5_cpu,$CO/c5/C30_n4000_Torig_grain_cross10_ftc5/fit_restart.npz" \
    --label "Original-temperature cycle 6 (job 122221)"
cat <<MD
Earlier cycles: \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles/HOW-TO-RERUN.md\` (1-3),
\`.../model_vs_observations_grain_ft_cycles_Torig/HOW-TO-RERUN.md\` (original temperature, 4-6).
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O/HOW-TO-RERUN.md"
