#!/bin/sh
set -eu

src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step7.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/QLT.cpp" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/util/sep_parker_core.cpp" \
  "$src_dir/util/sep_focused_transport_core.cpp" \
  "$src_dir/util/sep_focused_transport_mfp_core.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_mover_validation.cpp" \
  "$src_dir/test/controlled/test_registered_movers.cpp" \
  -o "$build_dir/test_registered_movers"

(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 \
  ./test_registered_movers parker)
for id in PARK01 PARK02 PARK03 PARK04 PARK05 PARK06 PARK07; do
  grep -q "$id" "$build_dir/controlled-parker-results.json"
  grep -q "$id" "$build_dir/controlled-parker-results.xml"
done

# The ordinary focused target remains artifact-free because build_dir is
# disposable.  The Python campaign runner opts into evidence retention with an
# explicit directory; copy only the already-validated machine-readable files
# and never expose sanitizer executables or temporary compiler products.
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/controlled-parker-results.json" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir/controlled-parker-results.xml" "$SRCSEP_REPORT_DIR/"
fi

# Keep the focused executable and the linked production CLI on the same
# descriptor factory.  These source checks complement the runtime checks above:
# a later edit cannot silently copy the cases back into a private test-only
# implementation or omit the factory object from the production archive.
grep -q 'ControlledMoverDescriptors()' "$src_dir/component_tests.cpp"
grep -q 'util/sep_mover_validation.o' "$src_dir/makefile"

grep -q 'return SEP::ParticleMover_Parker;' \
  "$src_dir/production_mover_runtime.cpp" || {
  echo "FAIL PARK-SOURCE: registry does not dispatch the canonical Parker mover" >&2
  exit 1
}
if grep -En 'rnd\(|Vector3D::Distribution::Normal|ParticleMover_ParkerEquation' \
    "$src_dir/parker_mover.cpp" >/dev/null; then
  echo "FAIL PARK-SOURCE: canonical Parker mover contains hidden RNG or legacy dispatch" >&2
  exit 1
fi
echo "PASS PARK-SOURCE: canonical Parker dispatch uses the explicit-stream core"
