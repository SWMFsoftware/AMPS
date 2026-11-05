#!/usr/bin/env bash
# Run the Stage-7 input policy (when present) and pure channel-importance
# planner without configuring or modifying the AMPS build tree.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
compiler=${CXX:-c++}
temporary=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage7.XXXXXX")
trap 'rm -rf -- "$temporary"' EXIT

policy_test="$source_root/input/test_comet_channel_importance_policy.pl"
if [[ -f "$policy_test" ]]; then
  perl "$policy_test"
else
  printf 'SKIP Stage-7 input policy test (not present): %s\n' "$policy_test"
fi

"$compiler" -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I"$source_root" \
  "$script_dir/test_stage7_channel_importance_planner.cpp" \
  -o "$temporary/test_stage7_channel_importance_planner"
"$temporary/test_stage7_channel_importance_planner"
