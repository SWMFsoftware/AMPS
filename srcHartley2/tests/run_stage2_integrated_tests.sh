#!/usr/bin/env bash
# Validate the Stage-2 runtime transaction at the required MPI decompositions.
# The runner consumes saved logs so a source checkout can rerun the strict
# parser without launching an expensive executable.  Producing those logs is
# intentionally separate and must use one clean, identical executable/input.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

required=(
  HARTLEY2_STAGE2_LOG_1
  HARTLEY2_STAGE2_LOG_4
  HARTLEY2_STAGE2_LOG_12
)
missing=()
for name in "${required[@]}"; do
  if [[ -z "${!name:-}" ]]; then
    missing+=("$name")
  fi
done

if (( ${#missing[@]} != 0 )); then
  printf 'SKIP Stage-2 integrated logs: set'
  printf ' %s' "${missing[@]}"
  printf '\n'
  exit 0
fi

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage2_integrated_logs.py" \
  --log 1 "$HARTLEY2_STAGE2_LOG_1" \
  --log 4 "$HARTLEY2_STAGE2_LOG_4" \
  --log 12 "$HARTLEY2_STAGE2_LOG_12"
