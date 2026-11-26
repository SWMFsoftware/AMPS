#!/bin/bash
# Set up and run one case of the I3 campaign qfit_20261123 with an executable built by compile_amps.sh.
#
#   run_case.sh --amps PATH --case N --out DIR [--np NP] [--chemistry sa1982_vh|state_resolved]
#               [--short ITERATIONS] [--setup-only]
#
#   --amps PATH       the executable (compile_amps.sh: <dest>/amps)
#   --case N          1-6
#   --out DIR         campaign output directory; the run directory is DIR/<chemistry>/caseN_np<NP>/production
#                     (the layout that srcHartley2/i3/tools/make_qfit_figures.py --campaign DIR reads)
#   --np NP           MPI ranks (default 128, as in the campaign)
#   --chemistry       sa1982_vh (default: the deck's OH loss, the campaign) or state_resolved (the earlier
#                     default OH chemistry, selected on the command line: --chemistry-variant
#                     OH_LOSS=STATE_RESOLVED; same deck and Q(H2O), which were fitted for SA1982_VH)
#   --short K         test run: max_iterations, sampling_length and weight_adaptation_interval set to K in
#                     the run's copy of the deck (the run directory name then ends in _short<K>)
#   --setup-only      prepare the run directory and print the command, do not run
#
# MPI launcher: $MPIRUN (default "mpiexec"), called as  $MPIRUN -n NP ./amps -i i3.in -case N [...].
# The run directory receives amps, i3.in, i3.conf, oh_lifetime_sa1982.csv (the deck names it relative to
# itself), SHA256SUMS and provenance.txt; it must not exist yet. The production run converges after
# 70000-100000 iterations (1.7-3.1 days with 128 ranks in the campaign; about 1.8 GB per rank).
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CAMPAIGN=$(dirname "$HERE")
REPO=$(cd "$CAMPAIGN/../../../.." && pwd)                # AMPS root of this checkout
TABLE=$REPO/srcHartley2/i3/reference/oh_lifetime_sa1982.csv
TABLE_SHA=d3eb13f8522c44e2c284435dccae232d5c604c05e2b038adc2c49f2f455b57d0
AMPS="" CASE="" OUT="" NP=128 CHEM=sa1982_vh SHORT="" SETUP_ONLY=0
while [[ $# -gt 0 ]]; do
  case $1 in
    --amps) AMPS=$2; shift 2 ;;
    --case) CASE=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --np) NP=$2; shift 2 ;;
    --chemistry) CHEM=$2; shift 2 ;;
    --short) SHORT=$2; shift 2 ;;
    --setup-only) SETUP_ONLY=1; shift ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option $1 (see --help)"; exit 2 ;;
  esac
done
[[ -x $AMPS && -n $OUT ]] || { echo "--amps PATH (executable) and --out DIR are required"; exit 2; }
[[ $CASE =~ ^[1-6]$ ]] || { echo "--case must be 1-6"; exit 2; }
case $CHEM in
  sa1982_vh) VARIANT=() ;;
  state_resolved) VARIANT=(--chemistry-variant OH_LOSS=STATE_RESOLVED) ;;
  *) echo "--chemistry must be sa1982_vh or state_resolved"; exit 2 ;;
esac
( cd "$CAMPAIGN" && sha256sum -c --quiet SHA256SUMS ) || { echo "campaign files differ from SHA256SUMS"; exit 3; }
[[ $(sha256sum "$TABLE" | cut -d' ' -f1) == "$TABLE_SHA" ]] || { echo "$TABLE is not the campaign's table"; exit 3; }
AMPS=$(cd "$(dirname "$AMPS")" && pwd)/$(basename "$AMPS")
BUILD=$(dirname "$AMPS")

RUN=$OUT/$CHEM/case${CASE}_np${NP}${SHORT:+_short$SHORT}/production
[[ ! -e $RUN ]] || { echo "refusing to reuse existing run directory $RUN"; exit 2; }
mkdir -p "$RUN"
RUN=$(cd "$RUN" && pwd)
cp -p "$AMPS" "$RUN/amps"
cp -p "$CAMPAIGN/i3.conf" "$RUN/i3.conf"
cp -p "$TABLE" "$RUN/oh_lifetime_sa1982.csv"
if [[ -n $SHORT ]]; then
  sed -e "s/^\( *simulation.max_iterations *=\).*/\1 $SHORT/" \
      -e "s/^\( *simulation.sampling_length *=\).*/\1 $SHORT/" \
      -e "s/^\( *simulation.weight_adaptation_interval *=\).*/\1 $SHORT/" "$CAMPAIGN/i3.in" > "$RUN/i3.in"
else
  cp -p "$CAMPAIGN/i3.in" "$RUN/i3.in"
fi
CMD="${MPIRUN:-mpiexec} -n $NP ./amps -i i3.in -case $CASE ${VARIANT[*]:-}"
{
  echo "case=$CASE chemistry=$CHEM np=$NP short=${SHORT:-no} host=$(hostname) date=$(date -Is)"
  echo "campaign=$CAMPAIGN checkout_commit=$(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "amps=$AMPS"
  [[ -f $BUILD/build_provenance.txt ]] && sed 's/^/build: /' "$BUILD/build_provenance.txt"
  echo "command: $CMD"
} > "$RUN/provenance.txt"
( cd "$RUN" && sha256sum amps i3.in i3.conf oh_lifetime_sa1982.csv > SHA256SUMS )
echo "run directory: $RUN"
echo "command: $CMD"
[[ $SETUP_ONLY == 1 ]] && exit 0
cd "$RUN"
set +e
${MPIRUN:-mpiexec} -n "$NP" ./amps -i i3.in -case "$CASE" ${VARIANT[@]+"${VARIANT[@]}"} > "case${CASE}.run.log" 2>&1
rc=$?
set -e
echo "end=$(date -Is) rc=$rc"
tail -3 "case${CASE}.run.log"
exit $rc
