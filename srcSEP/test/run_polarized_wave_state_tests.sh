#!/bin/sh
set -eu
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-polarized-wave-state.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Wpedantic -Werror -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/.." -I"$src_root/util" \
  "$src_root/util/sep_polarized_wave_state.cpp" \
  "$src_root/util/sep_polarized_wave_state_validation.cpp" \
  "$src_root/util/sep_resonant_collision.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_polarized_wave_state.cpp" \
  -o "$build_dir/test_polarized_wave_state"
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_polarized_wave_state"
