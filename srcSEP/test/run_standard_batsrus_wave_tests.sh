#!/bin/sh
set -eu
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-batsrus-waves.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
# This contract suite has no PIC/imported-array dependency. The native factory
# separately verifies actual allocated vertex storage and read-only ownership.
${CXX:-c++} -std=c++17 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_standard_batsrus_waves.cpp" "$src_dir/util/sep_standard_batsrus_wave_validation.cpp" \
  "$src_dir/test/test_standard_batsrus_waves.cpp" -o "$build_dir/test_standard_batsrus_waves"
test_status=0
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_standard_batsrus_waves) || test_status=$?
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/standard-batsrus-results.json" "$build_dir/standard-batsrus-results.xml" "$SRCSEP_REPORT_DIR/"
fi
exit "$test_status"
