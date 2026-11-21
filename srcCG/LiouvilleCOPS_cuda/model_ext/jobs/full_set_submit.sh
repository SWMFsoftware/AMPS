#!/bin/bash
# ==== version-controlled copy: header added 2026-11-20, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/full_set_submit.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/full_set_submit.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2f70a23b9487f5e8be8e4de0da2c7205e8b788e3181728c7a3e76d13728dae13
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: submits the figure jobs of one full set (the c9 set's job list, parameterized; 2026-11-17); 2026-11-20 (full
#   sets of the new-law fits): MEM_CAP, SKIP; BASE / D2REF / MCODE / MODEL_TEXT passed through (defaults unchanged,
#   tested with a stub sbatch)
# RUN: see model_ext/README.md (full sets of the stage-2 states, 2026-11-17; new-law states 2026-11-20 via
#   figures/full_set_q002.py)
# ==== end of header ====
# Submit the figure jobs of ONE full figure set (the jobs of the c9 set figures_cluster/12M/full_set_softplus_smm1_9cycles,
# unchanged; 2026-11-17, lead: full sets of the six cases of 12M_stage2_reg_best). Prints the job ids and logs every
# submission in model_ext/logs/submissions.txt. The maps / contributions of the state must exist or be submitted with
# MAPS_JOB / CONTRIB_JOB (dependencies). Run from anywhere (login node).
# Environment (all required unless a default is given):
#   KEY        short name (folder full_set_<KEY>, file tags)          CHAIN   root of the state (FT = CHAIN/c$NCYC/ft/final)
#   NCYC       last cycle (default 5)        FT / FIT  final fine tuning / fit (default CHAIN/c$NCYC/{ft,fit}/final)
#   MAPS       maps folder of the state      GEOM      geometry of its maps job       CONTRIB  contributions folder
#   XC         total-rate cross-check folder (default CHAIN/v2_total_crosscheck_c$NCYC)
#   PHASE      contrib_phase_data output (default CHAIN/contrib_phase_c$NCYC)
#   MX         frozen nonneg_maps dir of the maps (default model_ext_frozen/nonneg_maps_20261117T224340Z)
#   CX         frozen contrib code of the state (default model_ext_frozen/nonneg_contrib_20261117T224340Z)
#   ADAPTER    contrib_phase_data --state-adapter (option C / D2 exports; default none); SBINS slope_bins.npz of the
#              contributions (default CONTRIB/slope_bins.npz)
#   PUB, PUBSTEM, PUBTABLE  publication NG/RG figure folder, stem, table (made beforehand with pub_ng_rg_altitude.py)
#   MODEL_NAME, MODEL_LABEL, MODEL_SHORT, TITLE_PRESSURES  caption / legend texts
#   MAPS_JOB, CONTRIB_JOB  job ids to wait for (optional)        NODE  GPU node (default gpu002)
#   2026-11-20 (full sets of the new-law fits; defaults unchanged): MEM_CAP  --mem of every GPU job (default: the values
#   below); SKIP  comma list of figure jobs not submitted (activity, h2o: products copied from an identical earlier run);
#   BASE, D2REF, MCODE, MODEL_TEXT  passed through to the jobs (export them; see the job scripts)
set -euo pipefail
J=/data001/highresolu_vt/cops_output/jobs; L=/data001/highresolu_vt/cops_output/model_ext/logs/submissions.txt
F=/data001/highresolu_vt/cops_output/figures_cluster/12M
: "${KEY:?}" "${CHAIN:?}" "${MAPS:?}" "${GEOM:?}" "${CONTRIB:?}" "${PUB:?}" "${PUBSTEM:?}" "${PUBTABLE:?}" "${MODEL_NAME:?}"
: "${MODEL_LABEL:?}" "${MODEL_SHORT:?}" "${TITLE_PRESSURES:?}"
export NCYC=${NCYC:-5}
export FT=${FT:-$CHAIN/c$NCYC/ft/final} FIT=${FIT:-$CHAIN/c$NCYC/fit/final}
export XC=${XC:-$CHAIN/v2_total_crosscheck_c$NCYC} PHASE=${PHASE:-$CHAIN/contrib_phase_c$NCYC}
export MX=${MX:-/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_maps_20261117T224340Z}
export CX=${CX:-/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261117T224340Z}
export ADAPTER=${ADAPTER:-} SBINS=${SBINS:-$CONTRIB/slope_bins.npz}
export RUN=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5 OUT=$F/full_set_$KEY
export CHAIN MAPS GEOM CONTRIB PUB PUBSTEM PUBTABLE MODEL_NAME MODEL_LABEL MODEL_SHORT TITLE_PRESSURES
NODE=${NODE:-gpu002}
D="--nodelist=$NODE --kill-on-invalid-dep=yes"
DM=${MAPS_JOB:+--dependency=afterok:$MAPS_JOB}; DC=${CONTRIB_JOB:+--dependency=afterok:$CONTRIB_JOB}
cd "$J"
log() { echo "$(date -u +%FT%TZ) $1 FULL SET $KEY: $2 (CHAIN=$CHAIN NCYC=$NCYC MAPS=$MAPS GEOM=$GEOM CONTRIB=$CONTRIB MX=$MX CX=$CX OUT=$OUT)" >> "$L"; }
m() { echo "--mem=${MEM_CAP:-$1}"; }
a=skipped; b=skipped; c=skipped
if [[ ",${SKIP:-}," != *,activity,* ]]; then a=$(QUANTITY=activity sbatch --parsable $D $DM $(m 24G) --export=ALL pub_maps_v2.sbatch); log $a "pub_maps_v2 activity"; fi
if [[ ",${SKIP:-}," != *,h2o,* ]]; then
b=$(QUANTITY=h2o TIME=10:33:10 sbatch --parsable $D $DM $(m 24G) --export=ALL pub_maps_v2.sbatch); log $b "pub_maps_v2 h2o 10:33:10"
c=$(QUANTITY=h2o TIME=10:39:26 sbatch --parsable $D $DM $(m 24G) --export=ALL pub_maps_v2.sbatch); log $c "pub_maps_v2 h2o 10:39:26"
fi
d=$(QUANTITY=total TIME=10:33:10 OUT=$XC sbatch --parsable $D $DM $(m 24G) --export=ALL pub_maps_v2.sbatch); log $d "pub_maps_v2 total 10:33:10 OUT=$XC"
e=$(QUANTITY=total TIME=10:39:26 OUT=$XC sbatch --parsable $D $DM $(m 24G) --export=ALL pub_maps_v2.sbatch); log $e "pub_maps_v2 total 10:39:26 OUT=$XC"
DCM=${CONTRIB_JOB:+--dependency=afterok:$CONTRIB_JOB}      # the contributions job itself waits for the maps job
f=$(STEPS=sources sbatch --parsable $D $DCM $(m 32G) --export=ALL pub_source_areas.sbatch); log $f "pub_source_areas sources"
g=$(STEPS=cliffs V2=$XC CLIFF_GROUPS=all sbatch --parsable $D $(m 16G) --dependency=afterok:$d:$e${CONTRIB_JOB:+:$CONTRIB_JOB} --export=ALL pub_source_areas.sbatch); log $g "pub_source_areas cliffs V2=$XC CLIFF_GROUPS=all"
h=$(TAG=$KEY STEPS=set1,set4,pressures sbatch --parsable $D $DM $(m 40G) --export=ALL pub_maps_figures.sbatch); log $h "pub_maps_figures set1,set4,pressures PUB=$PUB"
i=$(OUTD=$PHASE sbatch --parsable --kill-on-invalid-dep=yes $DCM --export=ALL contrib_phase_data.sbatch); log $i "contrib_phase_data OUTD=$PHASE ADAPTER=$ADAPTER"
k=$(DATA=$PHASE OUTF=$OUT sbatch --parsable $D ${MEM_CAP:+--mem=$MEM_CAP} --dependency=afterok:$i --export=ALL pub_contrib_density.sbatch); log $k "pub_contrib_density DATA=$PHASE"
echo "$KEY: activity $a h2o $b $c total $d $e sources $f cliffs $g set1/set4/pressures $h contrib_phase $i contrib_density $k"
