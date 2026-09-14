#!/bin/sh
set -eu

# Build the exact dependency-free shared transport kernels.  A temporary build
# directory keeps sanitizer and object output out of the source handoff.
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step6.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" \
  "$src_dir/util/sep_transport_common.cpp" \
  "$src_dir/test/step6_8/test_transport_numerics.cpp" \
  -o "$build_dir/test_transport_numerics"

# LeakSanitizer cannot enumerate /proc tasks in the managed execution profile;
# AddressSanitizer and UndefinedBehaviorSanitizer remain active for the run.
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_transport_numerics" step6

# The source contract complements the dependency-free numerical checks: all
# three production shells must enter through the common PIC adapter, and the
# common implementation must own the only segment-list insertion sequence.
for mover in parker_mover.cpp focused_transport_dmumu.cpp focused_transport_mfp.cpp; do
  grep -q 'transport_common.h' "$src_dir/$mover" || {
    echo "FAIL CORE-SOURCE: $mover bypasses the common transport boundary" >&2
    exit 1
  }
done
grep -q 'CommitAndAttach' "$src_dir/parker_mover.cpp"
grep -q 'CommitAndAttach' "$src_dir/focused_transport_dmumu.cpp"
grep -q 'CommitAndAttach' "$src_dir/focused_transport_mfp.cpp"
grep -q 'AccountAdiabaticCoolingFlag' "$src_dir/parker_mover.cpp"
grep -q 'AccountAdiabaticCoolingFlag' "$src_dir/focused_transport_dmumu.cpp"
grep -q 'AccountAdiabaticCoolingFlag' "$src_dir/focused_transport_mfp.cpp"
echo "PASS CORE-SOURCE: all production shells use common validation/attachment"
