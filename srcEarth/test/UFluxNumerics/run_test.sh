#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${TMPDIR:-/tmp}/sep_in_geospace_flux_numerics_test"
mkdir -p "$build_dir"

"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic \
  "$script_dir/test_flux_numerics.cpp" -o "$build_dir/test_flux_numerics"
"$build_dir/test_flux_numerics"
