#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_slope_test_submit.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_slope_test_submit.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 07d1fefc0f072a528b52fb7f49ce82a0aaa4ae383c1af097934943522547a465
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: plan / submission of the eight control fits (decision 3c)
# RUN: localnn_slope_test_submit.sh --min L [--submit]
# ==== end of header ====
# The slope-effect control set (user 2026-11-21, decision 3c): 8 fits with jobs/localnn_slope_test.sbatch, each as two
# jobs (the second continues after a wall-time stop, or only confirms a finished fit).
#   localnn_slope_test_submit.sh --min L            show the plan; submits NOTHING
#   localnn_slope_test_submit.sh --min L --submit   submit (only after the user's OK of the plan and of L)
#   L = the lower bound of the six multipliers in the free_min fit, 0 < L < 1. No default: the user's decision.
set -euo pipefail
R=/data001/highresolu_vt/cops_output
JOB=$R/jobs/localnn_slope_test.sbatch
CODE=$R/model_ext_frozen/localnn_slopetest_20261122T033334Z
L= SUBMIT=0
while [ $# -gt 0 ]; do
  case $1 in
    --min) L=${2:?--min needs a value}; shift 2;;
    --submit) SUBMIT=1; shift;;
    *) echo "unknown option $1" >&2; exit 2;;
  esac
done
[ -n "$L" ] || { echo "--min L is required (no default: the user's decision)" >&2; exit 2; }
python3 -I -c "import sys; v=float(sys.argv[1]); sys.exit(0 if 0<v<1 else 1)" "$L" || { echo "--min must be between 0 and 1" >&2; exit 2; }
(cd "$CODE" && sha256sum -c --quiet SHA256SUMS)
RUNROOT=$R/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/slopetest_q002_f10_pm0.01_tploff_$(date -u +%Y%m%dT%H%M%SZ)
V="none only_40_50 only_50_60 only_60_70 only_70_80 only_80_90 only_90_180 free_min"
echo "code    $CODE"
echo "runroot $RUNROOT"
echo "8 fits x 2 jobs (second: afterany); each job normal128cores, 16 CPUs, 32 GB, 12:00:00; 128 CPUs at once"
for v in $V; do
  case $v in none) d="all six held at 1";; free_min) d="all six free, >= $L";; *) d="only ${v#only_} free (>= 1), five held at 1";; esac
  printf '  %-12s %s\n' "$v" "$d"
done
if [ $SUBMIT -ne 1 ]; then echo "plan only: nothing submitted"; exit 0; fi
mkdir "$RUNROOT"
cp -p "$0" "$JOB" "$RUNROOT/"
echo "lower bound of free_min: $L" > "$RUNROOT/settings.txt"
for v in $V; do
  E=ALL,RUNROOT=$RUNROOT,VARIANT=$v,MODEL_EXT=$CODE
  [ $v = free_min ] && E=$E,SLOPE_MIN=$L
  a=$(sbatch --parsable --job-name=lnn_3c_$v --export=$E "$JOB")
  b=$(sbatch --parsable --job-name=lnn_3c_${v}_2 --dependency=afterany:$a --export=$E "$JOB")
  echo "$v $a $b" | tee -a "$RUNROOT/submission.txt"
done
