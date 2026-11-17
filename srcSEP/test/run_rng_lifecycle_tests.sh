#!/bin/sh
set -eu

src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-rng-lifecycle.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_parker_core.cpp" \
  "$src_dir/util/sep_focused_transport_core.cpp" \
  "$src_dir/util/sep_focused_transport_mfp_core.cpp" \
  "$src_dir/util/sep_source_lifecycle.cpp" \
  "$src_dir/util/sep_rng_validation.cpp" \
  "$src_dir/test/rng/test_rng_lifecycle.cpp" \
  -o "$build_dir/test_rng_lifecycle"

(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_rng_lifecycle)

for mover in parker_mover.cpp focused_transport_dmumu.cpp focused_transport_mfp.cpp; do
  grep -q 'ResumeRandomStream' "$src_dir/$mover"
  grep -q 'StageRandomStream' "$src_dir/$mover"
done
echo 'PASS RNG-SOURCE: every production mover resumes and stages complete persistent state'

if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/rng-lifecycle-results.json" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir/rng-lifecycle-results.xml" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir"/seeded-campaign-*.json "$build_dir"/seeded-campaign-*.xml "$SRCSEP_REPORT_DIR/"
fi
