#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_cliff_submit.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_cliff_submit.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: a64c587247887308edc8d72bf336161b4f8c69cd6a74d239a5f090de585b7484
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: plan / submission of the cliff test
# RUN: localnn_cliff_submit.sh [all|cliffs] [--submit]
# ==== end of header ====
# The cliff test (user 2026-11-22: "prepare the cliff test"; plan figures_cluster/12M/q002_decisions_20261121/PLAN_cliff.md).
#   localnn_cliff_submit.sh [all|cliffs]             show what would be submitted; submits NOTHING
#   localnn_cliff_submit.sh [all|cliffs] --submit    submit (only after the user's OK of the plan)
#   all = cliffs + the moderate-slope control (proposed); cliffs = without the control (no verdict on specificity)
set -euo pipefail
R=/data001/highresolu_vt/cops_output
TOOLS=$R/model_ext_frozen/cliff_20261122T200453Z
L=$R/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn
CONTROLS_C=$L/slopetest_c_q002_f10_pm0.01_tploff_20261122T195607Z
SET=all SUBMIT=0
for x in "$@"; do
  case $x in all|cliffs) SET=$x;; --submit) SUBMIT=1;; *) echo "unknown argument $x" >&2; exit 2;; esac
done
(cd "$TOOLS" && sha256sum -c --quiet SHA256SUMS)
(cd "$TOOLS/code" && sha256sum -c --quiet SHA256SUMS)
BASE_FIT=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["fit"])' "$CONTROLS_C/chain_pick.json")
FITS="cliffs_off_from_held cliffs_off_from_d2"
[ $SET = all ] && FITS="$FITS moderate_off_from_held moderate_off_from_d2"
RUNROOT=$L/clifftest_q002_f10_tploff_$(date -u +%Y%m%dT%H%M%SZ)
echo "runroot  $RUNROOT"
echo "tools    $TOOLS"
echo "baseline $BASE_FIT"
echo "every job: normal128cores, default account, 16 CPUs, 32 GB, 12:00:00 (the end check: 1:00:00)"
n=0
for f in $FITS; do
  case $f in cliffs_*) d="60-180 deg held at 0.01, 40-60 at 1";; *) d="40-60 deg held at 0.01, 60-180 at 1 (control)";; esac
  printf '  %-24s %-46s 3 jobs (afterany)\n' "$f" "$d"
  n=$((n + 1))
done
echo "  check                    evaluations only, the plan's reading rules             1 job"
echo "at most $((n * 16)) CPUs at once"
if [ $SUBMIT -ne 1 ]; then echo "plan only: nothing submitted"; exit 0; fi

mkdir "$RUNROOT"
cp -p "$0" "$R/jobs/localnn_cliff_fit.sbatch" "$R/jobs/localnn_cliff_check.sbatch" "$RUNROOT/"
E=ALL,RUNROOT=$RUNROOT,TOOLS=$TOOLS,CONTROLS_C=$CONTROLS_C
last=
for f in $FITS; do
  ids= dep=
  for k in 1 2 3; do
    j=$(sbatch --parsable --job-name="lnn_cliff_${f}_$k" ${dep:+--dependency=afterany:$dep} --export="$E,FIT=$f" "$R/jobs/localnn_cliff_fit.sbatch")
    ids="$ids $j" dep=$j
  done
  echo "$f$ids" | tee -a "$RUNROOT/submission.txt"
  last=$last:$dep
done
k=$(sbatch --parsable --job-name=lnn_cliff_check --dependency=afterany$last --export="$E" "$R/jobs/localnn_cliff_check.sbatch")
echo "check $k" | tee -a "$RUNROOT/submission.txt"
