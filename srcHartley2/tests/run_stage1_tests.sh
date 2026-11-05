#!/usr/bin/env bash
# Focused Stage-1 parser tests always run. Runtime policy/cadence and legacy
# output checks are enabled by explicit archive/log paths so the source suite
# remains usable in a checkout without large validation artifacts.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)

perl "$source_root/input/test_comet_population_policy.pl"
perl "$source_root/input/test_comet_reaction_network.pl"
perl -c "$source_root/input/hartley2.pl"

if [[ -n "${HARTLEY2_STAGE1_POLICY_LOG:-}" ]]; then
  runtime_arguments=("$HARTLEY2_STAGE1_POLICY_LOG")
  if [[ -n "${HARTLEY2_STAGE1_ALTERNATE_LOG:-}" ]]; then
    runtime_arguments+=(--alternate-log "$HARTLEY2_STAGE1_ALTERNATE_LOG")
  fi
  if [[ -n "${HARTLEY2_STAGE1_CADENCE_LOG:-}" ]]; then
    runtime_arguments+=(--cadence-log "$HARTLEY2_STAGE1_CADENCE_LOG")
  fi
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage1_population_policy_runtime.py" \
    "${runtime_arguments[@]}"
else
  echo "SKIP Stage-1 runtime table: set HARTLEY2_STAGE1_POLICY_LOG"
fi

if [[ -n "${HARTLEY2_STAGE1_LEGACY_BASELINE:-}" &&
      -n "${HARTLEY2_STAGE1_LEGACY_CANDIDATE:-}" ]]; then
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage1_legacy_equivalence.py" \
    "$HARTLEY2_STAGE1_LEGACY_BASELINE" \
    "$HARTLEY2_STAGE1_LEGACY_CANDIDATE"
else
  echo "SKIP Stage-1 legacy output: set both legacy plot directories"
fi
