#!/bin/sh
set -eu
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-birth-window.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
${CXX:-c++} -std=c++17 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_birth_window.cpp" \
  "$src_dir/util/sep_birth_window_validation.cpp" \
  "$src_dir/test/test_birth_window.cpp" -o "$build_dir/test_birth_window"
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_birth_window)
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/birth-window-results.json" "$build_dir/birth-window-results.xml" "$SRCSEP_REPORT_DIR/"
fi
