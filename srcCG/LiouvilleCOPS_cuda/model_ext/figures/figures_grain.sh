#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_grain.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_grain.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 73edc4cef6431b6ffd6a556ebbf9d4a91b84e565c5e5fa4166f48f446439ef53
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: all grain fits of a resolution.
# RUN: on the login node, no GPU: <this> <RES> <run>
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/<RES>/model_vs_observations_grain/
# ==== end of header ====
# Fits with the dust grain in NG and the NG = RG crossing constraint (model_ext, user request 2026-10-27) against the
# COPS observations; the same model_obs_figures.py as figures_model_vs_obs.sh (observations only, NG / RG / NG+RG with
# crossings, time and altitude, .eps + .png 300 dpi, late window shaded). Run on the login node (no GPU):
#     figures_grain.sh <resolution label, e.g. 12M> <runner restart dir>
# Output: $COPS_OUTPUT_ROOT/figures_cluster/<label>/model_vs_observations_grain/ (+ tables/ of the fine-tuned models).
# The NG of a grain fit includes the grain (surface + grain); fine-tuned = after grain_finetune.py.
set -u
RES=${1:?resolution label}; R=${2:?runner restart dir}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
O=$COPS_OUTPUT_ROOT/figures_cluster/$RES/model_vs_observations_grain; TB=$O/tables
mkdir -p "$O" "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite"
G=$R/grain/fits; F=$R/grain/finetune
gf() { echo "$G/$1_n2000_grain_cross10"; }                    # free temperatures + grain
gm() { echo "$G/$1_n2000_mono10_ordered_grain_cross10"; }     # monotonic + ordered + grain
nf() { echo "$R/fits/$1_n2000/seg01"; }                       # free temperatures, no grain
nm() { echo "$R/c16/fits/$1_n2000_mono10_ordered"; }          # monotonic + ordered, no grain
tf() { local d=$(gf $1); .venv/bin/python $S/finetune_table.py "$F/$(basename $d)_m5_i5_cpu" "$d" "$TB/$(basename $d)_m5_i5_cpu.dat" && echo "$TB/$(basename $d)_m5_i5_cpu.dat"; }
tm() { local d=$(gm $1); .venv/bin/python $S/finetune_table.py "$F/$(basename $d)_m5_i5_cpu" "$d" "$TB/$(basename $d)_m5_i5_cpu.dat" && echo "$TB/$(basename $d)_m5_i5_cpu.dat"; }
declare -A D=([B29]="no template" [N30]="late NG-only template" [C30]='cliff $\geq$60$^\circ$ template'
              [S30]='sin(slope) template' [NC31]='late NG-only + cliff $\geq$60$^\circ$ templates'
              [NS31]='late NG-only + sin(slope) templates')
declare -A N=([B29]=no_template [N30]=late_NG_only [C30]=cliff60 [S30]=sin_slope [NC31]=late_NG_cliff60
              [NS31]=late_NG_sin_slope)   # file names: templates in words, no model codes (user 2026-10-27)
RERUN=()   # every plotted fit / fine tuning, for HOW-TO-RERUN.md (user 2026-11-03)
SF="free temperatures"; SM="monotonic penalty + ordered temperatures"; GR="dust grain in NG + crossing constraint"
plot() {   # plot <name> <title> <ref dir> <label=path:style>...
    local name=$1 title=$2 ref=$3; shift 3
    local args=() d
    for s in "$@"; do
        args+=(--fit "$s"); d=${s#*=}; d=${d%:*}
        case $d in "$TB"/*.dat) d=$F/$(basename "$d" .dat);; esac   # a fine-tuned table: its fine-tuning directory
        RERUN+=("$d")
    done
    $P --out "$O/$name" --title "$title" --ref "$ref" "${args[@]}" && echo "wrote ${name}_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
}
# 1. all templates, each fit setting, with the grain
plot grain_free_T_all_templates "SHAP8 $RES, $SF, $GR: fitted models vs COPS observations" "$(gf B29)" \
    "${D[B29]}=$(gf B29):dotted" "${D[N30]}=$(gf N30):dashdot" "${D[C30]}=$(gf C30):solid" "${D[NC31]}=$(gf NC31):dashed"
plot grain_monotonic_ordered_T_all_templates "SHAP8 $RES, $SM, $GR: fitted models vs COPS observations" "$(gm B29)" \
    "${D[B29]}=$(gm B29):dotted" "${D[N30]}=$(gm N30):dashdot" "${D[S30]}=$(gm S30):solid" "${D[NS31]}=$(gm NS31):dashed"
# 2. with vs without the grain, every fit
for M in B29 N30 C30 NC31; do
    plot "with_vs_without_grain_free_T_${N[$M]}" "SHAP8 $RES, $SF, ${D[$M]}: with and without the dust grain" "$(gf $M)" \
        "without grain=$(nf $M):dashed" "with grain + crossing=$(gf $M):solid"
done
for M in B29 N30 S30 NS31; do
    plot "with_vs_without_grain_monotonic_ordered_T_${N[$M]}" "SHAP8 $RES, $SM, ${D[$M]}: with and without the dust grain" "$(gm $M)" \
        "without grain=$(nm $M):dashed" "with grain + crossing=$(gm $M):solid"
done
# 3. fine-tuned (grain refitted with the late groups)
plot grain_finetuned_free_T_all_templates "SHAP8 $RES, $SF, $GR, fine-tuned: models vs COPS observations" "$(gf B29)" \
    "${D[B29]}=$(tf B29):dotted" "${D[N30]}=$(tf N30):dashdot" "${D[C30]}=$(tf C30):solid" "${D[NC31]}=$(tf NC31):dashed"
plot grain_finetuned_monotonic_ordered_T_all_templates "SHAP8 $RES, $SM, $GR, fine-tuned: models vs COPS observations" "$(gm B29)" \
    "${D[B29]}=$(tm B29):dotted" "${D[N30]}=$(tm N30):dashdot" "${D[S30]}=$(tm S30):solid" "${D[NS31]}=$(tm NS31):dashed"
# 4. the lowest-cost grain fit of each setting (12M: free T cliff >=60, monotonic NS31), fitted and fine-tuned
plot grain_best_fitted_and_finetuned "SHAP8 $RES, $GR: lowest-cost fit of each setting, fitted and fine-tuned" "$(gf C30)" \
    "${D[C30]}, $SF=$(gf C30):dashed" "${D[C30]}, $SF, fine-tuned=$(tf C30):solid" \
    "${D[NS31]}, $SM=$(gm NS31):dotted" "${D[NS31]}, $SM, fine-tuned=$(tm NS31):dashdot"
# exact commands to rerun every plotted case from the restart files and from scratch (user request 2026-11-03)
[ ${#RERUN[@]} -gt 0 ] && .venv/bin/python -W ignore $COPS_OUTPUT_ROOT/figures_cluster/scripts/how_to_rerun.py --res "$RES" --run "$R" \
    --figdir "$O" --figure-command "$COPS_OUTPUT_ROOT/figures_cluster/scripts/$(basename "$0") $RES $R" "${RERUN[@]}"
