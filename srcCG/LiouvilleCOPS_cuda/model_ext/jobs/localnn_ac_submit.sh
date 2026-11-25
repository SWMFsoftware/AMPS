#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_ac_submit.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_ac_submit.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: acc849032d8dda60bba14524f3e6c1d1aae1fa1af22ced6f904decd9b58d216d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: plan / submission of the follow-up a + c
# RUN: localnn_ac_submit.sh [a|c|ac] [--submit]
# ==== end of header ====
# Follow-up 'a' and 'c' of the slope test (user 2026-11-22: "as recommended: a, c"; plan figures_cluster/12M/
# q002_decisions_20261121/PLAN_ac.md).
#   localnn_ac_submit.sh [a|c|ac]            show what would be submitted; submits NOTHING
#   localnn_ac_submit.sh [a|c|ac] --submit   submit (only after the user's OK of the plan)
set -euo pipefail
R=/data001/highresolu_vt/cops_output
TOOLS=$R/model_ext_frozen/ac_20261122T183034Z
CODE=$R/model_ext_frozen/localnn_slopetest_20261122T033334Z
PART=ac SUBMIT=0
for x in "$@"; do
  case $x in a|c|ac) PART=$x;; --submit) SUBMIT=1;; *) echo "unknown argument $x" >&2; exit 2;; esac
done
(cd "$TOOLS" && sha256sum -c --quiet SHA256SUMS)
(cd "$CODE" && sha256sum -c --quiet SHA256SUMS)
A_FITS="a_from_40_50 a_from_50_60 a_from_60_70 a_from_70_80 a_from_80_90"
C_FITS="c_from_40_50 c_from_50_60 c_from_60_70 c_from_70_80 c_from_80_90 c_from_uniform"
case $PART in a) FITS=$A_FITS;; c) FITS=$C_FITS;; ac) FITS="$A_FITS $C_FITS";; esac
RUNROOT=$R/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/slopetest_${PART}_q002_f10_pm0.01_tploff_$(date -u +%Y%m%dT%H%M%SZ)
echo "runroot $RUNROOT"
echo "code    $CODE"
echo "tools   $TOOLS"
echo "every job: normal128cores, default account, 16 CPUs, 32 GB, 12:00:00 (the end check: 1:00:00)"
n=0
for f in $FITS; do
  case $f in a_*) d="six free, >= 0.01";; *) d="six held at 1";; esac
  printf '  %-16s %-18s 3 jobs (afterany)\n' "$f" "$d"
  n=$((n + 1))
done
if [ "$PART" != a ]; then printf '  %-16s %-18s 2 jobs, after the six c fits\n' chain "5 cycles, held at 1"; fi
echo "  check            evaluations only    1 job, after everything"
echo "at most $((n * 16)) CPUs at once"
if [ $SUBMIT -ne 1 ]; then echo "plan only: nothing submitted"; exit 0; fi

mkdir "$RUNROOT"
cp -p "$0" "$R/jobs/localnn_ac_fit.sbatch" "$R/jobs/localnn_ac_chain.sbatch" "$R/jobs/localnn_ac_check.sbatch" "$RUNROOT/"
E=ALL,RUNROOT=$RUNROOT,TOOLS=$TOOLS
last= clast=
for f in $FITS; do
  ids= dep=
  for k in 1 2 3; do
    j=$(sbatch --parsable --job-name="lnn_ac_${f}_$k" ${dep:+--dependency=afterany:$dep} --export="$E,FIT=$f" "$R/jobs/localnn_ac_fit.sbatch")
    ids="$ids $j" dep=$j
  done
  echo "$f$ids" | tee -a "$RUNROOT/submission.txt"
  last=$last:$dep
  case $f in c_*) clast=$clast:$dep;; esac
done
if [ -n "$clast" ]; then
  a=$(sbatch --parsable --job-name=lnn_ac_chain_1 --dependency=afterany$clast --export="$E" "$R/jobs/localnn_ac_chain.sbatch")
  b=$(sbatch --parsable --job-name=lnn_ac_chain_2 --dependency=afterany:$a --export="$E" "$R/jobs/localnn_ac_chain.sbatch")
  echo "chain $a $b" | tee -a "$RUNROOT/submission.txt"
  last=$last:$b
fi
k=$(sbatch --parsable --job-name=lnn_ac_check --dependency=afterany$last --export="$E" "$R/jobs/localnn_ac_check.sbatch")
echo "check $k" | tee -a "$RUNROOT/submission.txt"
