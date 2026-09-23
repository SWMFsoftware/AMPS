#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${TMPDIR:-/tmp}/sep_in_geospace_trajectory_core_test"
mkdir -p "$build_dir"

# Production Earth sources use C++17.  This syntax-only gate includes both public
# backend headers together and catches contract alias/overload drift without needing
# to link MPI, Geopack, SPICE, or AMPS.
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -pedantic \
  -I"$script_dir" -I"$script_dir/../.." \
  -fsyntax-only "$script_dir/test_public_headers.cpp"

# The common contract and mover reference test intentionally retain the stricter
# dependency-free C++11 compilation path.
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic \
  -I"$script_dir" -I"$script_dir/../.." \
  "$script_dir/test_trajectory_core.cpp" \
  "$script_dir/../../gridless/GridlessParticleMovers.cpp" \
  -o "$build_dir/test_trajectory_core"
"$build_dir/test_trajectory_core"
