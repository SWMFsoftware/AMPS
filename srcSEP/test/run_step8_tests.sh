#!/bin/sh
set -eu

src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step8.XXXXXX")
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
  ./test_registered_movers fte-dmumu)
for id in FTED01 FTED02 FTED03 FTED04 FTED05 FTED06 FTED07 FTED08; do
  grep -q "$id" "$build_dir/controlled-fte-dmumu-results.json"
  grep -q "$id" "$build_dir/controlled-fte-dmumu-results.xml"
done

# Evidence retention is opt-in so direct Make invocations leave no generated
# files in the source tree.  The destination is selected by the Python runner,
# while this script copies only reports that passed the ID/schema smoke checks
# above; transient sanitizer binaries remain inside build_dir.
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/controlled-fte-dmumu-results.json" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir/controlled-fte-dmumu-results.xml" "$SRCSEP_REPORT_DIR/"
fi

# Prove that this exact descriptor catalog is reachable from the linked CLI,
# not only from the focused sanitizer executable built above.
grep -q 'ControlledMoverDescriptors()' "$src_dir/component_tests.cpp"
grep -q 'util/sep_mover_validation.o' "$src_dir/makefile"

grep -q 'return SEP::ParticleMover_FocusedTransport_Dmumu;' \
  "$src_dir/production_mover_runtime.cpp" || {
  echo "FAIL FTED-SOURCE: registry does not dispatch the coefficient-driven mover" >&2
  exit 1
}
if grep -En 'ParkerSpiral|QLT::calculateDmuMu|rnd\(|CellIntegratedWaveEnergy' \
    "$src_dir/focused_transport_dmumu.cpp" >/dev/null; then
  echo "FAIL FTED-SOURCE: canonical Dmumu mover embeds a background, RNG, or wave-state update" >&2
  exit 1
fi
grep -q 'QueueWaveContribution' "$src_dir/focused_transport_dmumu.cpp"
grep -q 'turbulenceStateIdentity' "$src_dir/focused_transport_dmumu.cpp"
grep -q 'CouplingEventType::BoundaryExit' "$src_dir/focused_transport_dmumu.cpp"
echo "PASS FTED-SOURCE: Dmumu mover is provider-driven with escape-safe coupling"
