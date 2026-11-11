#!/usr/bin/env bash
# Parse archived Stage-3 logs without silently launching MPI jobs.  Generate
# all three logs with one frozen executable and the same input/seed, then pass
# their absolute paths through these environment variables.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
required=(HARTLEY2_STAGE3_LOG_1 HARTLEY2_STAGE3_LOG_4 HARTLEY2_STAGE3_LOG_12)
missing=()
for name in "${required[@]}"; do
  if [[ -z "${!name:-}" ]]; then
    missing+=("$name")
  fi
done
if (( ${#missing[@]} != 0 )); then
  printf 'SKIP Stage-3 integrated logs: set'
  printf ' %s' "${missing[@]}"
  printf '\n'
  exit 0
fi

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage3_integrated_logs.py" \
  --log 1 "$HARTLEY2_STAGE3_LOG_1" \
  --log 4 "$HARTLEY2_STAGE3_LOG_4" \
  --log 12 "$HARTLEY2_STAGE3_LOG_12"
