#!/usr/bin/env bash
# Run the Stage-8 input policy and pure deterministic forcing planner without
# configuring AMPS or touching the shared build tree.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
compiler=${CXX:-c++}
temporary=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage8.XXXXXX")
trap 'rm -rf -- "$temporary"' EXIT

perl "$source_root/input/test_comet_total_reaction_forcing_policy.pl"
perl -c "$source_root/input/hartley2.pl"

"$compiler" -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I"$source_root" \
  "$script_dir/test_stage8_total_reaction_forcing_planner.cpp" \
  -o "$temporary/test_stage8_total_reaction_forcing_planner"
"$temporary/test_stage8_total_reaction_forcing_planner"
