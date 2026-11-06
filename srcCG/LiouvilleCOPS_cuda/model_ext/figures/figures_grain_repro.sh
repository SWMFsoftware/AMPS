#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_grain_repro.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_grain_repro.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: ca3b8b3889dc522fcd825f44d5fc38b7783b306e7b29d3cda5b87b19d37660c0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: reproductions of the best grain fits.
# RUN: on the login node, no GPU: <this> <RES> <run> <repro job ids>
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_repro/
# ==== end of header ====
# Figures from the reproductions of the best grain fits (model_ext/grain_repro.sbatch, user request 2026-11-03:
# "start a new fit that reproduces that fit for further minor tweaks"; figures requested the same day). The same
# model_obs_figures.py as figures_grain.sh (observations only, NG / RG / NG+RG with crossings, time and altitude,
# .eps + .png 300 dpi, late window shaded). Run on the login node (no GPU):
#     figures_grain_repro.sh <resolution label, e.g. 12M> <runner restart dir> <repro job id>...
# Output: $COPS_OUTPUT_ROOT/figures_cluster/<label>/model_vs_observations_grain_repro/ (+ tables/).
# The NG of a grain fit includes the grain (surface + grain); fine-tuned = after grain_finetune.py.
set -u
RES=${1:?resolution label}; R=${2:?runner restart dir}; shift 2
[ $# -gt 0 ] || { echo "give the reproduction job ids" >&2; exit 2; }
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
O=$COPS_OUTPUT_ROOT/figures_cluster/$RES/model_vs_observations_grain_repro; TB=$O/tables
mkdir -p "$O" "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite"
declare -A D=([B29]="no template" [N30]="late NG-only template" [C30]='cliff $\geq$60$^\circ$ template'
              [S30]='sin(slope) template' [NC31]='late NG-only + cliff $\geq$60$^\circ$ templates'
              [NS31]='late NG-only + sin(slope) templates')
declare -A N=([B29]=no_template [N30]=late_NG_only [C30]=cliff60 [S30]=sin_slope [NC31]=late_NG_cliff60
              [NS31]=late_NG_sin_slope)   # file names: templates in words, no model codes (user 2026-10-27)
declare -A ST=([grain_cross10]="free temperatures" [mono10_ordered_grain_cross10]="monotonic penalty + ordered temperatures")
GR="dust grain in NG + crossing constraint"
declare -A FT=()   # table -> fine-tuning directory
RERUN=()           # every plotted fit / fine tuning, for HOW-TO-RERUN.md (user 2026-11-03)
table() {   # table <fine-tuning dir> <fit dir>: echo the table path
    local t=$TB/$(basename "$1").dat
    .venv/bin/python $S/finetune_table.py "$1" "$2" "$t" && FT[$t]=$1 && echo "$t"
}
plot() {   # plot <name> <title> <ref dir> <label=path:style>...
    local name=$1 title=$2 ref=$3; shift 3
    local args=() d
    for s in "$@"; do
        args+=(--fit "$s"); d=${s#*=}; d=${d%:*}
        [ -n "${FT[$d]:-}" ] && d=${FT[$d]}
        RERUN+=("$d")
    done
    $P --out "$O/$name" --title "$title" --ref "$ref" "${args[@]}" && echo "wrote ${name}_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
}
BEST=()
for J in "$@"; do
    RD=$(find "$R/grain/repro" -maxdepth 1 -name "*_repro_$J" -type d | head -1)
    [ -n "$RD" ] || { echo "no reproduction of job $J in $R/grain/repro" >&2; exit 2; }
    SRC=$(basename "$RD"); SRC=${SRC%_repro_$J}           # e.g. C30_n2000_grain_cross10
    M=${SRC%%_*}; SET=${SRC#*_n2000_}                      # model code, fit setting
    OD=$R/grain/fits/$SRC
    # tables are built here in the main shell (FT must survive), so no $(...) around table()
    table "$R/grain/repro/finetune/$(basename "$RD")_m5_i5_cpu" "$RD" > /dev/null || exit 1
    TR=$TB/$(basename "$RD")_m5_i5_cpu.dat
    table "$R/grain/finetune/${SRC}_m5_i5_cpu" "$OD" > /dev/null || exit 1
    TO=$TB/${SRC}_m5_i5_cpu.dat
    # 1. reproduction vs the original fit, fitted and fine-tuned (identical curves are expected to coincide)
    plot "grain_repro_vs_original_${N[$M]}_$( [ "$SET" = grain_cross10 ] && echo free_T || echo monotonic_ordered_T )" \
        "SHAP8 $RES, ${D[$M]}, ${ST[$SET]}, $GR: reproduction (job $J) vs original" "$OD" \
        "original fit=$OD:dashed" "reproduction fit=$RD:dotted" \
        "original, fine-tuned=$TO:solid" "reproduction, fine-tuned=$TR:dashdot"
    BEST+=("${D[$M]}, ${ST[$SET]} (reproduction)=$RD:$( [ ${#BEST[@]} -eq 0 ] && echo dashed || echo dotted )"
           "${D[$M]}, ${ST[$SET]} (reproduction), fine-tuned=$TR:$( [ ${#BEST[@]} -eq 0 ] && echo solid || echo dashdot )")
    [ ${#BEST[@]} -eq 2 ] && REF=$RD
done
# 2. the reproduced fits together, fitted and fine-tuned (as grain_best_fitted_and_finetuned of figures_grain.sh)
plot grain_repro_best_fitted_and_finetuned "SHAP8 $RES, $GR: reproduced best fits, fitted and fine-tuned" "$REF" "${BEST[@]}"
# exact commands to rerun every plotted case from the restart files and from scratch (user request 2026-11-03)
[ ${#RERUN[@]} -gt 0 ] && .venv/bin/python -W ignore $S/how_to_rerun.py --res "$RES" --run "$R" \
    --figdir "$O" --figure-command "$S/$(basename "$0") $RES $R $*" "${RERUN[@]}"
