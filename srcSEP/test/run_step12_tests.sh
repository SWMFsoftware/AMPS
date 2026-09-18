#!/bin/sh
set -eu

# Step 12 runs the same physical contributions through different synthetic
# worker/rank layouts.  Equality of the canonical evidence hash is the declared
# bitwise reproducibility criterion for reduced particle-to-wave accumulators.
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step12.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir/util" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/util/sep_reproducible_reduction.cpp" \
  "$src_dir/test/step12/test_reproducibility.cpp" \
  -o "$build_dir/test_reproducibility"

ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_reproducibility"

grep -q 'Rank' "$src_dir/util/sep_reproducible_reduction.h"
grep -q 'thread number are deliberately absent' \
  "$src_dir/util/sep_reproducible_reduction.h"
grep -q 'CanonicalPartitionReduction' \
  "$src_dir/util/sep_reproducible_reduction.h"
echo "PASS PAR-SOURCE: physical keys and gather-then-canonical-reduce policy present"
