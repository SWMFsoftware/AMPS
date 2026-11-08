#!/usr/bin/env bash
# Exercise the production AMPS boundary-flux sampler for the three directly
# injected I3 parent species.  This runner does not configure or build AMPS:
# EXE must already be the DIM1 executable compiled from input/i3.conf.  The
# physical model still comes from the runtime file and explicit case selector,
# so the test exercises the same two-layer contract as a production launch.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
executable=${1:-"$source_root/amps"}
runtime_input=${I3_RUNTIME_INPUT:-"$source_root/input/i3.in"}
case_index=${I3_CASE:-1}
ranks=${I3_INJECTION_MPI_RANKS:-1}
particles=${I3_INJECTION_PARTICLES:-200000}
log=${I3_INJECTION_LOG:-"$source_root/i3.surface-injection.validation.log"}
launcher=${MPIEXEC:-mpiexec}

if [[ ! -x "$executable" ]]; then
  printf 'I3 surface-injection executable is not executable: %s\n' "$executable" >&2
  exit 2
fi
if [[ ! -f "$runtime_input" ]]; then
  printf 'I3 runtime input does not exist: %s\n' "$runtime_input" >&2
  exit 2
fi
if [[ ! "$case_index" =~ ^[1-9][0-9]*$ ]]; then
  printf 'I3_CASE must be a positive integer: %s\n' "$case_index" >&2
  exit 2
fi
if [[ ! "$ranks" =~ ^[1-9][0-9]*$ ]]; then
  printf 'I3_INJECTION_MPI_RANKS must be a positive integer: %s\n' "$ranks" >&2
  exit 2
fi
if [[ ! "$particles" =~ ^[1-9][0-9]*$ ]] || (( particles < 1000 )); then
  printf 'I3_INJECTION_PARTICLES must be an integer of at least 1000: %s\n' "$particles" >&2
  exit 2
fi

"$launcher" -n "$ranks" "$executable" \
  -i "$runtime_input" -case "$case_index" \
  --validate-surface-injection \
  --validation-particles "$particles" 2>&1 | tee "$log"

if ! grep -Fq 'PASS I3 surface-injection validation' "$log"; then
  printf 'I3 surface-injection PASS marker is missing; log: %s\n' "$log" >&2
  exit 1
fi

printf 'PASS I3 production surface-injection sampler case=%s input=%s; log: %s\n' \
  "$case_index" "$runtime_input" "$log"
