#!/usr/bin/env bash
# Run Stage-6 policy and pure packet-transaction tests without configuring or
# modifying the AMPS build tree.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
compiler=${CXX:-c++}
temporary=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage6.XXXXXX")
trap 'rm -rf -- "$temporary"' EXIT

perl "$source_root/input/test_comet_daughter_packet_policy.pl"

"$compiler" -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I"$source_root" \
  "$script_dir/test_stage6_daughter_packet_planner.cpp" \
  -o "$temporary/test_stage6_daughter_packet_planner"
"$temporary/test_stage6_daughter_packet_planner"
