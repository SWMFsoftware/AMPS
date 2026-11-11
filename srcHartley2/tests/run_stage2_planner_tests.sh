#!/usr/bin/env bash
# Compile and run the geometry-independent part of the Stage-2 component
# carrier planner.  No AMPS build tree is needed, and the temporary executable
# is removed automatically so the source tree never retains test binaries.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
stage2_tmp=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage2.XXXXXX")
trap 'rm -rf -- "$stage2_tmp"' EXIT

${CXX:-c++} -std=c++11 -O2 -Wall -Wextra -Werror \
  -I"$source_root" \
  "$script_dir/test_stage2_component_carrier_planner.cpp" \
  -o "$stage2_tmp/test_stage2_component_carrier_planner"

"$stage2_tmp/test_stage2_component_carrier_planner"
