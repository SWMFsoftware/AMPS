#!/bin/sh
set -eu
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-conservative-resampling.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
${CXX:-c++} -std=c++11 -O1 -g -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/../src/pic" \
  -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_population_control.cpp" \
  "$src_dir/util/sep_conservative_resampling.cpp" \
  "$src_dir/util/sep_conservation_resampling_validation.cpp" \
  "$src_dir/test/test_conservation_resampling.cpp" \
  -o "$build_dir/test_conservation_resampling"
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_conservation_resampling)
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/conservative-resampling.json" "$build_dir/conservative-resampling.xml" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir"/seeded-campaign-*.json "$build_dir"/seeded-campaign-*.xml "$SRCSEP_REPORT_DIR/"
fi
