#!/bin/sh
set -eu
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-juttner-volume.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
# Compile the exact production sampler and native-registry callbacks. These
# independent distribution/measure tests do not qualify SWMF thermodynamics
# or the linked particle source; those remain separate native gates.
${CXX:-c++} -std=c++17 -O1 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir/util" -I"$src_dir/test" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_maxwell_juttner_volume.cpp" \
  "$src_dir/util/sep_maxwell_juttner_volume_validation.cpp" \
  "$src_dir/test/test_maxwell_juttner_volume.cpp" -o "$build_dir/test_volume"
status=0
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_volume) || status=$?
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/maxwell-juttner-volume-results.json" "$build_dir/maxwell-juttner-volume-results.xml" "$SRCSEP_REPORT_DIR/"
fi
exit "$status"
