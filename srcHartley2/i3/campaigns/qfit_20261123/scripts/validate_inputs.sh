#!/bin/bash
# Check that an executable accepts the campaign deck for all six cases with both OH-loss variants:
# AMPS --validate-runtime-input (2 MPI ranks; it resolves the case, checks the complete comet, chemistry,
# collision, cooling, brightness and radial input chain, and exits before mesh or particle allocation).
#
#   validate_inputs.sh --amps PATH --work DIR
#
# For each case it requires "PASS comet runtime input validation" and the selected variant
# (Chemistry.Variant.OH_LOSS=SA1982_VH from the deck, STATE_RESOLVED from --chemistry-variant).
# MPI launcher: $MPIRUN (default "mpiexec"). Logs: DIR/<chemistry>/caseN/validate.log.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CAMPAIGN=$(dirname "$HERE")
REPO=$(cd "$CAMPAIGN/../../../.." && pwd)
AMPS="" WORK=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --amps) AMPS=$2; shift 2 ;;
    --work) WORK=$2; shift 2 ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option $1 (see --help)"; exit 2 ;;
  esac
done
[[ -x $AMPS && -n $WORK ]] || { echo "--amps PATH (executable) and --work DIR are required"; exit 2; }
[[ ! -e $WORK ]] || { echo "refusing to use existing $WORK"; exit 2; }
AMPS=$(cd "$(dirname "$AMPS")" && pwd)/$(basename "$AMPS")
mkdir -p "$WORK"
cp -p "$CAMPAIGN/i3.in" "$REPO/srcHartley2/i3/reference/oh_lifetime_sa1982.csv" "$WORK/"
fail=0
for CHEM in sa1982_vh state_resolved; do
  want=${CHEM^^}
  extra=()
  [[ $CHEM == state_resolved ]] && extra=(--chemistry-variant OH_LOSS=STATE_RESOLVED)
  for C in 1 2 3 4 5 6; do
    d=$WORK/$CHEM/case$C
    mkdir -p "$d"
    if ( cd "$d" && ${MPIRUN:-mpiexec} -n 2 "$AMPS" --validate-runtime-input -i ../../i3.in -case $C \
           ${extra[@]+"${extra[@]}"} > validate.log 2>&1 ) \
       && grep -aq 'PASS comet runtime input validation' "$d/validate.log" \
       && grep -aq "Chemistry.Variant.OH_LOSS=$want" "$d/validate.log"; then
      echo "PASS $CHEM case $C"
    else
      echo "FAIL $CHEM case $C (log: $d/validate.log)"; fail=1
    fi
  done
done
exit $fail
