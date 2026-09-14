#!/bin/sh
set -eu

# Build the exact production validation descriptors in a disposable directory.
# ASan/UBSan cover both the stochastic mover campaigns and the independent
# finite-volume reference implementation without requiring AMPS/PIC headers.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step15-numerical.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -O1 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" \
  "$src_root/util/sep_transport_common.cpp" \
  "$src_root/util/sep_parker_core.cpp" \
  "$src_root/util/sep_focused_transport_core.cpp" \
  "$src_root/util/sep_focused_transport_mfp_core.cpp" \
  "$src_root/util/sep_coefficient_physics.cpp" \
  "$src_root/util/sep_coefficient_registry.cpp" \
  "$src_root/util/sep_test_registry.cpp" \
  "$src_root/util/sep_scientific_validation.cpp" \
  "$src_root/test/step15/test_numerical_validation.cpp" \
  -o "$build_dir/test_numerical_validation"

(cd "$build_dir" && ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} \
  ./test_numerical_validation)
for report in "$build_dir/step15-numerical-results.json" \
              "$build_dir/step15-numerical-results.xml"; do
  test -s "$report" || { echo "FAIL VAL-REPORT: missing $report" >&2; exit 1; }
  for id in VAL01 VAL02 VAL03; do
    grep -q "$id" "$report" || {
      echo "FAIL VAL-REPORT: $report omits $id" >&2
      exit 1
    }
  done
done

# The default focused run remains artifact-free.  The campaign orchestrator may
# opt in to copies in its own disposable evidence directory, keeping report
# ownership explicit and avoiding generated files in a source distribution.
if test -n "${STEP15_REPORT_DIR:-}"; then
  mkdir -p "$STEP15_REPORT_DIR"
  cp "$build_dir/step15-numerical-results.json" \
     "$build_dir/step15-numerical-results.xml" "$STEP15_REPORT_DIR/"
fi
echo "PASS VAL-REPORT: numerical, cross-mover, and cross-model evidence retained"
