#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_area_scaling.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_area_scaling.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: cf410a5dc5e832e1b762308ba3111efc57d10e6aa0b91f765525d7dfde27a7ec
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: the area-scaling test.
# RUN: on the login node, no GPU: <this> <RES> <run> <area_scaling output dir> <figure folder name>
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_fixed_T_area_scaling/
# ==== end of header ====
# Figures of an area-scaling test (model_ext/area_scaling.py; user request 2026-11-04): the fine-tuned model and the
# same model with the activity of the target areas scaled, against the COPS observations (model_obs_figures.py,
# as figures_grain.sh), plus the region map, and HOW-TO-RERUN.md. Login node, no GPU:
#     figures_area_scaling.sh <resolution label> <runner restart dir> <area_scaling output dir> <figure folder name>
set -u
RES=${1:?}; R=${2:?}; A=${3:?}; NAME=${4:?}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; X=$COPS_OUTPUT_ROOT/model_ext
O=$COPS_OUTPUT_ROOT/figures_cluster/$RES/$NAME
mkdir -p "$O"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
eval "$(.venv/bin/python - "$A/area_scaling.json" <<'PY'
import json, shlex, sys
s = json.load(open(sys.argv[1]))
t = s["targets"]
rg = [x for x in t if x["gauge"] == "RG"]; ng = [x for x in t if x["gauge"] == "NG"]
hm = lambda x: x["UTC"].split("T")[1]
print("FIT=%s FT=%s" % (shlex.quote(s["fit"]), shlex.quote(s["finetune"])))
lab = ("RG area (%d%% of the H2O pressure at %s) x%g, NG area (%s-%s) x%g, both x%g" % (
    round(100 * s["fraction"]), hm(rg[0]), s["rg_factor"], hm(ng[0]), hm(ng[-1]), s["ng_factor"], s["overlap_factor"]))
print("LAB=%s" % shlex.quote(lab))
PY
)"
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite"
$P --out "$O/finetuned_vs_area_scaled" --ref "$FIT" \
    --title "SHAP8 $RES, cliff \$\\geq\$60\$^\\circ\$ template, T75 and Tmin fixed at 133 K, dust grain: fine-tuned and area-scaled" \
    --fit "fine-tuned=$A/tables/finetuned.dat:dashed" --fit "fine-tuned, $LAB=$A/tables/scaled.dat:solid" \
    && echo "wrote finetuned_vs_area_scaled_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
cp -p "$A/region_map.png" "$A/region_map.eps" "$A/area_scaling.json" "$O/"
J=$(ls -d "$A".provenance-* | tail -1)
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: SHAP8 $RES area-scaling figures in this folder

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_area_scaling.sh\`. No restart file was
written or changed; the scaled model is a recomputation of the fine-tuned state (\`model_ext/area_scaling.py\`).
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_area_scaling.sh $RES $R $A $NAME
\`\`\`

## A1. The area scaling from the fit's restart file and its fine tuning (GPU job, about 5 min)

Exact command of this result: \`$J/command.txt\`. Resubmit with a new output directory:

\`\`\`bash
sbatch --export=ALL,COPS_RUN=$R,FIT=$FIT,FT=$FT,OUT=${A}_rerun $X/area_scaling.sbatch
\`\`\`

Other settings: \`ARGS="--fraction F --rg-last K --rg-factor V --ng-from-end K0 K1 --ng-factor V"\` in the
\`--export\` list (defaults 0.9, 1, 1.2, 4 6, 1.5).

## B. The fit and its fine tuning, from the restart files or from scratch

See \`$COPS_OUTPUT_ROOT/figures_cluster/$RES/model_vs_observations_grain_fixed_T/HOW-TO-RERUN.md\` (the fit
\`$(basename $FIT)\` and its fine tuning \`$(basename $FT)\`), then A1 and A0.
MD
echo "wrote $O/HOW-TO-RERUN.md"
