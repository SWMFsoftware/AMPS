#!/bin/sh
set -eu
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-background-derivatives.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Wpedantic -Werror -O1 -g \
  -I"$src_root/util" -I"$src_root/.." \
  "$src_root/util/sep_background_derivative_validation.cpp" \
  "$src_root/util/sep_parker_core.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_background_derivatives.cpp" \
  -o "$build_dir/test_background_derivatives"
"$build_dir/test_background_derivatives"
