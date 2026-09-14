#!/bin/sh
set -eu

# Compile the exact dependency-light event-driven production core. Keeping the
# build in a disposable directory prevents sanitizer products from entering the
# source handoff and makes the test usable without AMPS, PIC, or MPI headers.
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step9.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" \
  "$src_dir/QLT.cpp" \
  "$src_dir/util/sep_transport_common.cpp" \
  "$src_dir/util/sep_parker_core.cpp" \
  "$src_dir/util/sep_focused_transport_core.cpp" \
  "$src_dir/util/sep_focused_transport_mfp_core.cpp" \
  "$src_dir/util/sep_test_registry.cpp" \
  "$src_dir/util/sep_mover_validation.cpp" \
  "$src_dir/test/controlled/test_registered_movers.cpp" \
  -o "$build_dir/test_registered_movers"

(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 \
  ./test_registered_movers fte-mfp)
for id in FTEM01 FTEM02 FTEM03 FTEM04 FTEM05 FTEM06 FTEM07 FTEM08; do
  grep -q "$id" "$build_dir/controlled-fte-mfp-results.json"
  grep -q "$id" "$build_dir/controlled-fte-mfp-results.xml"
done

# Retain structured evidence only when the orchestration layer requests it.
# This preserves the source-only target's clean-tree behavior and ensures that
# neither a sanitizer executable nor an intermediate object enters a handoff.
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/controlled-fte-mfp-results.json" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir/controlled-fte-mfp-results.xml" "$SRCSEP_REPORT_DIR/"
fi

# The focused executable and production component-test CLI deliberately share
# this factory/object pair. These checks guard against future registry drift.
grep -q 'ControlledMoverDescriptors()' "$src_dir/component_tests.cpp"
grep -q 'util/sep_mover_validation.o' "$src_dir/makefile"

# The public registry must dispatch the new thin adapter, while the old large
# shell must be absent rather than lingering as unlinked migration material.
grep -q 'return SEP::ParticleMover_FocusedTransport_EventDriven;' \
  "$src_dir/production_mover_runtime.cpp"
grep -q 'focused_transport_mfp.o' "$src_dir/makefile"
if grep -q 'fte_mover.o' "$src_dir/makefile"; then
  echo "FAIL FTEM-SOURCE: legacy fte_mover.cpp is still a production object" >&2
  exit 1
fi
if test -e "$src_dir/fte_mover.cpp"; then
  echo "FAIL FTEM-SOURCE: retired legacy source still exists" >&2
  exit 1
fi
echo "PASS FTEM-SOURCE: the registry uses the canonical event-driven MFP adapter"
