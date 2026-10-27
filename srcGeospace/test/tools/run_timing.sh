#!/bin/bash
# run_timing.sh -- timing comparison of several amps executables on the same tests.
#
# Usage:
#   srcGeospace/test/tools/run_timing.sh <results-dir> <label>=<amps> [<label>=<amps>...]
#
#   e.g. run_timing.sh $GEOSPACE_RUNS/results/timing_2026-10-30 \
#            cpu=$GEOSPACE_RUNS/builds/cpu_double/amps gpu=$GEOSPACE_RUNS/builds/gpu_double/amps
#
# Runs list_timing (12 of srcEarth's tests: gridless cutoff and density, directional
# transmission, the C8/C12/C13/C17 campaigns, Mode3D MESH and a gridded C6 deck) for
# every executable with srcEarth's runner, in two rounds of opposite order (first..last,
# then last..first), so that slow changes of the machine state do not favour one
# executable.  Each run goes into its own run root <results-dir>/timing_<label>_<round>/
# (run_srcearth_tests.sh).  A srcEarth executable can be timed the same way.
# The machine load before each run is printed; for meaningful numbers the machine
# should otherwise be idle.  Summarize with
#   summarize.py --root <label>1=<results-dir>/timing_<label>_1 ...
set -euo pipefail
[ $# -ge 2 ] || { sed -n 2,17p "$0"; exit 2; }
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT="$1"; shift
mkdir -p "$OUT"
LABELS=(); EXES=()
for spec in "$@"; do
  LABELS+=("${spec%%=*}"); EXES+=("${spec#*=}")
done
run() {   # run <index> <round>
  local label=${LABELS[$1]} exe=${EXES[$1]}
  echo "=== $(date '+%F %T') round $2: $label (load: $(cut -d' ' -f1-3 /proc/loadavg))"
  "$HERE/run_srcearth_tests.sh" --run-root "$OUT/timing_${label}_$2" --amps "$exe" \
      --list "$HERE/list_timing" -- -j 1 --no-start-progress \
      > "$OUT/timing_${label}_$2.log" 2>&1 || true
  grep -E '^\[(OK|MISMATCH)' "$OUT/timing_${label}_$2.log" || true
}
n=${#LABELS[@]}
for ((i=0; i<n; i++)); do run $i 1; done
for ((i=n-1; i>=0; i--)); do run $i 2; done
echo "=== $(date '+%F %T') done"
