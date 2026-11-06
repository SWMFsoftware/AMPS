#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_model_vs_obs.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_model_vs_obs.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 15c796e560ed1304e59ab6579dd4191d1fd758e6a586d1c631c9af71c4fe7159
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: the runner and C16 fits (no grain) of a resolution.
# RUN: on the login node, no GPU: <this> <RES> <run>
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/<RES>/model_vs_observations/
# ==== end of header ====
# Fitted model pressures against the COPS observations (NG, RG), for one resolution (user request 2026-10-27: plots
# compare the model with the observations; no passion comparison). Publication format .eps + .png, 300 dpi, late
# window shaded and labelled (scripts/model_obs_figures.py; observations only, no ratio panel). Run on the login node (no GPU):
#     figures_model_vs_obs.sh <resolution label, e.g. 12M> <runner restart dir>
# Output: $COPS_OUTPUT_ROOT/figures_cluster/<label>/model_vs_observations/ and NOMENCLATURE.txt there.
# Fit sets: C14 = <run>/fits/<M>_n2000/seg01 (runner, the "earlier" set), C16 = <run>/c16/fits/<M>_n2000_mono10_ordered.
# A plot whose fits are missing (run still going) is skipped with a note.
set -u
RES=${1:?resolution label}; R=${2:?runner restart dir}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
O=$COPS_OUTPUT_ROOT/figures_cluster/$RES/model_vs_observations
mkdir -p "$O"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1

P=".venv/bin/python -W ignore $COPS_OUTPUT_ROOT/figures_cluster/scripts/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite"
C14=$R/fits; C16=$R/c16/fits
c14() { echo "$C14/$1_n2000/seg01"; }
c16() { echo "$C16/$1_n2000_mono10_ordered"; }
# labels must not contain "=" (pub_figures splits LABEL=DIR there). The figures name the templates and the fit
# settings in words (user request 2026-10-27); the model codes appear only in NOMENCLATURE.txt and the fit directories.
declare -A D=([B29]="no template" [N30]="late NG-only template" [C30]='cliff $\geq$60$^\circ$ template'
              [S30]='sin(slope) template' [NC31]='late NG-only + cliff $\geq$60$^\circ$ templates'
              [NS31]='late NG-only + sin(slope) templates')
RERUN=()   # every plotted fit, for HOW-TO-RERUN.md (user 2026-11-03)
S14="free temperatures"
S16="monotonic penalty + ordered temperatures"
ok() { for d in "$@"; do [ -f "$d/fit_result.txt" ] || return 1; done; }
plot() {   # plot <name> <title> <ref dir> <label=dir:style>...
    local name=$1 title=$2 ref=$3; shift 3
    local args=() dirs=()
    for s in "$@"; do args+=(--fit "$s"); d=${s#*=}; dirs+=("${d%:*}"); done
    if ! ok "${dirs[@]}"; then echo "skip $name (fit not finished)"; return; fi
    RERUN+=("${dirs[@]}")
    # panels NG, RG, NG and RG together (crossings marked), vs time and altitude; separate altitude figures
    $P --out "$O/$name" --title "$title" --ref "$ref" "${args[@]}" && echo "wrote ${name}_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
}

# one plot per fit setting, all templates
plot free_T_all_templates "SHAP8 $RES, $S14: fitted models vs COPS observations" "$(c14 B29)" \
    "${D[B29]}=$(c14 B29):dotted" "${D[N30]}=$(c14 N30):dashdot" \
    "${D[C30]}=$(c14 C30):solid" "${D[NC31]}=$(c14 NC31):dashed"
plot monotonic_ordered_T_all_templates "SHAP8 $RES, $S16: fitted models vs COPS observations" "$(c16 B29)" \
    "${D[B29]}=$(c16 B29):dotted" "${D[N30]}=$(c16 N30):dashdot" \
    "${D[S30]}=$(c16 S30):solid" "${D[NS31]}=$(c16 NS31):dashed"
# each template choice with free temperatures vs its counterpart with the penalty and ordered temperatures
pair() {   # pair <file name> <C14 model> <C16 model>
    local n=$1 a=$2 b=$3 t
    if [ "$a" = "$b" ]; then t="${D[$a]}"; else t="templates and fit settings compared"; fi   # the legend names templates and settings
    plot "$n" "SHAP8 $RES, $t: fitted models vs COPS observations" "$(c16 $b)" \
        "${D[$a]}, $S14=$(c14 $a):dashed" "${D[$b]}, $S16=$(c16 $b):solid"
}
pair no_template_free_vs_monotonic_ordered_T B29 B29
pair late_NG_only_free_vs_monotonic_ordered_T N30 N30
pair cliff60_free_T_vs_sin_slope_monotonic_ordered_T C30 S30
pair late_NG_cliff60_free_T_vs_late_NG_sin_slope_monotonic_ordered_T NC31 NS31
# the lowest-cost fit of each setting (12M: C14 NC31, C16 S30)
plot best_late_NG_cliff60_free_T_vs_sin_slope_monotonic_ordered_T \
    "SHAP8 $RES vs COPS observations: lowest-cost fit of each setting" "$(c16 S30)" \
    "${D[NC31]}, $S14=$(c14 NC31):dashed" "${D[S30]}, $S16=$(c16 S30):solid"

cat > "$O/NOMENCLATURE.txt" <<EOF
SHAP8 $RES fits on the cluster (run $R), fitted model pressures against the COPS observations.

Model code (fit directory names only; the figures use the template names) = activity templates + number of fitted parameters
(25 H2O spherical-harmonic coefficients a00..a24, 3 gas temperatures, CO2 scale kappa, plus one amplitude per template):
  B29   no template (baseline)
  N30   + late_ng_not_rg   facets seen by the nude gauge but not the ram gauge in the late window
  C30   + cliff_ge_60deg   binary mask, effective slope >= 60 deg
  S30   + slope_sin        continuous, sin(min(effective slope, 90 deg)); replaces the cliff mask in C16
  NC31  N + C              NS31  N + S

Fit set = the stage in STATUS.md that defined the settings:
  C14   the "earlier" set (runner fits B29, N30, C30, NC31): free Tmax, T75, Tmin, no penalty
  C16   B29, N30, S30, NS31 with --monotonic-weight 10 (penalty on modelled pressures rising along the descent)
        and --ordered-temperature (Tmin <= T75 <= Tmax)
  both: mixed precision, --tbounds 30 400, max_nfev 2000 (none converged), default start
There is no "N31": the 31-parameter fits are NC31 (C14) and NS31 (C16).
EOF
echo "wrote NOMENCLATURE.txt"
# exact commands to rerun every plotted case from the restart files and from scratch (user request 2026-11-03)
[ ${#RERUN[@]} -gt 0 ] && .venv/bin/python -W ignore $COPS_OUTPUT_ROOT/figures_cluster/scripts/how_to_rerun.py --res "$RES" --run "$R" \
    --figdir "$O" --figure-command "$COPS_OUTPUT_ROOT/figures_cluster/scripts/$(basename "$0") $RES $R" "${RERUN[@]}"
