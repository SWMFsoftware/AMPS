#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_grain_fixT.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_grain_fixT.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: ac8652ad1ea509801ce8fc9b1cb4aba80e853778c01eb3460ceb512ffde6f7e5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M T75 = Tmin = 133 K fit.
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_fixed_T/
# ==== end of header ====
# Figures of the grain fit with T75 = Tmin = 133 K fixed (model_ext/grain_fixT.sbatch, user request 2026-11-04)
# against the original cliff >=60 deg free-temperature grain fit; the same model_obs_figures.py as figures_grain.sh
# (observations only, NG / RG / NG+RG with crossings, time and altitude, .eps + .png 300 dpi, late window shaded).
# Run on the login node (no GPU):   figures_grain_fixT.sh <resolution label, e.g. 12M> <runner restart dir>
# Output: $COPS_OUTPUT_ROOT/figures_cluster/<label>/model_vs_observations_grain_fixed_T/ (+ tables/).
set -u
RES=${1:?resolution label}; R=${2:?runner restart dir}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
O=$COPS_OUTPUT_ROOT/figures_cluster/$RES/model_vs_observations_grain_fixed_T; TB=$O/tables
mkdir -p "$O" "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite"
OF=$R/grain/fits/C30_n2000_grain_cross10;               OFT=$R/grain/finetune/C30_n2000_grain_cross10_m5_i5_cpu
XF=$R/grain/fixed_T/C30_n2000_grain_cross10_T75_133K_Tmin_133K; XFT=$R/grain/fixed_T/finetune/$(basename $XF)_m5_i5_cpu
for p in "$OFT:$OF" "$XFT:$XF"; do
    .venv/bin/python $S/finetune_table.py "${p%%:*}" "${p#*:}" "$TB/$(basename ${p%%:*}).dat" || exit 1
done
TO=$TB/$(basename $OFT).dat; TX=$TB/$(basename $XFT).dat
C='cliff $\geq$60$^\circ$ template'; GR="dust grain in NG + crossing constraint"
FREE='T75, Tmin free (fit: 30 K)'; FIX='T75 and Tmin fixed at 133 K'
plot() { $P --out "$O/$1" --title "$2" --ref "$3" "${@:4}" && echo "wrote ${1}_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"; }
plot grain_fixed_T_vs_free_T_cliff60_fitted_and_finetuned "SHAP8 $RES, $C, $GR: T75 = Tmin = 133 K fixed vs free temperatures" "$XF" \
    --fit "$FREE=$OF:dashed" --fit "$FREE, fine-tuned=$TO:solid" --fit "$FIX=$XF:dotted" --fit "$FIX, fine-tuned=$TX:dashdot"
plot grain_fixed_T_cliff60_fitted_and_finetuned "SHAP8 $RES, $C, $GR, T75 = Tmin = 133 K fixed" "$XF" \
    --fit "fitted=$XF:dashed" --fit "fine-tuned=$TX:solid"
# exact commands to rerun every plotted case from the restart files and from scratch (user request 2026-11-03)
.venv/bin/python -W ignore $S/how_to_rerun.py --res "$RES" --run "$R" --figdir "$O" \
    --figure-command "$S/$(basename "$0") $RES $R" "$OF" "$OFT" "$XF" "$XFT"
