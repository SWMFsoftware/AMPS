#!/bin/sh
set -eu
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-flux-state.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
${CXX:-c++} -std=c++17 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_flux_tube_geometry_core.cpp" \
  "$src_dir/util/sep_flux_tube_state.cpp" "$src_dir/util/sep_flux_tube_state_validation.cpp" \
  "$src_dir/test/test_flux_tube_state.cpp" -o "$build_dir/test_flux_tube_state"
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_flux_tube_state)
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/flux-state-results.json" "$build_dir/flux-state-results.xml" "$SRCSEP_REPORT_DIR/"
fi
