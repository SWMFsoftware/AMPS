#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${TMPDIR:-/tmp}/sep_in_geospace_trajectory_core_test"
mkdir -p "$build_dir"

"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic \
  -I"$script_dir" -I"$script_dir/../.." \
  "$script_dir/test_trajectory_core.cpp" \
  "$script_dir/../../gridless/GridlessParticleMovers.cpp" \
  -o "$build_dir/test_trajectory_core"
"$build_dir/test_trajectory_core"
