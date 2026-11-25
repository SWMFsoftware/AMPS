#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_multistart_submit.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_multistart_submit.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2387c89a2c7f7bf87d507c9f732ffc431c9919b93485c03342a9f21bee85e6fa
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: plan / submission of the multi-start refits (up to five)
# RUN: localnn_multistart_submit.sh [--submit]
# ==== end of header ====
# Multi-start refits of the option C reference problem (user 2026-11-22: "multi-start when ready"; up to five fits).
#   localnn_multistart_submit.sh            show what would be submitted; submits NOTHING
#   localnn_multistart_submit.sh --submit   submit
# One pair of jobs per possible parent. A pair waits for its parent (afterany), then the gate in
# jobs/localnn_multistart.sbatch runs or skips it. only_40_50 waits for all parents: it runs only if fewer than five
# of the others qualify.
set -euo pipefail
R=/data001/highresolu_vt/cops_output
L=$R/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn
CONTROLS=$L/slopetest_q002_f10_pm0.01_tploff_20261122T044019Z
GATE=$R/model_ext_frozen/multistart_20261122T064943Z/localnn_multistart.py
JOB=$R/jobs/localnn_multistart.sbatch
SUBMIT=0
case "${1:-}" in "") ;; --submit) SUBMIT=1;; *) echo "unknown option $1" >&2; exit 2;; esac
(cd "$(dirname "$GATE")" && sha256sum -c --quiet SHA256SUMS)
queued=$(squeue -u "$USER" -h -o %i)
last() { awk -v v="$1" '$1 == v {print $3}' "$CONTROLS/submission.txt"; }     # the parent's second job
wait_for() {
  local dep=
  for p in "$@"; do
    local j; j=$(last "$p")
    [ -n "$j" ] || { echo "no job of $p in submission.txt" >&2; exit 2; }
    if grep -qx "$j" <<< "$queued"; then dep=$dep:$j; fi
  done
  if [ -n "$dep" ]; then echo "--dependency=afterany$dep"; fi
}
RUNROOT=$L/multistart_q002_f10_pm0.01_tploff_$(date -u +%Y%m%dT%H%M%SZ)
FIVE="only_60_70 only_70_80 only_50_60 only_80_90 free_min"
echo "runroot $RUNROOT"
echo "each pair: 2 x (normal128cores, 16 CPUs, 32 GB, 12:00:00); at most five fits"
if [ $SUBMIT -eq 1 ]; then
  mkdir "$RUNROOT"
  cp -p "$0" "$JOB" "$RUNROOT/"
fi
for v in $FIVE only_40_50; do
  if [ "$v" = only_40_50 ]; then D=$(wait_for $FIVE only_40_50); else D=$(wait_for "$v"); fi
  echo "from_$v ${D:-(parent finished: no dependency)}"
  [ $SUBMIT -eq 1 ] || continue
  E=ALL,RUNROOT=$RUNROOT,CONTROLS=$CONTROLS,PARENT=$v,GATE=$GATE
  a=$(sbatch --parsable --job-name=lnn_ms_$v ${D:+$D} --export=$E "$JOB")
  b=$(sbatch --parsable --job-name=lnn_ms_${v}_2 --dependency=afterany:$a --export=$E "$JOB")
  echo "from_$v $a $b" >> "$RUNROOT/submission.txt"
  echo "  jobs $a $b"
done
[ $SUBMIT -eq 1 ] || echo "plan only: nothing submitted"
