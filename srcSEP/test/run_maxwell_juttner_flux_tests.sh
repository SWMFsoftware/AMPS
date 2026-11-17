#!/bin/sh
set -eu
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-maxwell-juttner.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
# Identical callbacks are linked by the native registry. This bounded build
# proves the standalone mathematical kernel under strict warnings and memory/
# undefined-behavior sanitizers, not the as-yet gated native source bridge.
${CXX:-c++} -std=c++17 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/test" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_maxwell_juttner_flux.cpp" "$src_dir/util/sep_maxwell_juttner_flux_validation.cpp" \
  "$src_dir/test/test_maxwell_juttner_flux.cpp" -o "$build_dir/test_maxwell_juttner_flux"
test_status=0
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_maxwell_juttner_flux) || test_status=$?
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/maxwell-juttner-results.json" "$build_dir/maxwell-juttner-results.xml" \
    "$build_dir/seeded-campaign-A.json" "$build_dir/seeded-campaign-reordered.json" \
    "$build_dir/seeded-campaign-B.json" "$build_dir/seeded-campaign-A.xml" "$SRCSEP_REPORT_DIR/"
fi
exit "$test_status"
