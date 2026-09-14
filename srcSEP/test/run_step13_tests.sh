#!/bin/sh
set -eu

# Step 13 must remain runnable from a source archive that does not contain the
# enclosing AMPS/PIC checkout.  Compile the exact production descriptors and
# transport/background kernels with both ASan and UBSan, then verify that both
# machine-readable report formats contain every stable acceptance ID.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step13.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" \
  "$src_root/util/sep_transport_common.cpp" \
  "$src_root/util/sep_focused_transport_core.cpp" \
  "$src_root/util/sep_focused_transport_mfp_core.cpp" \
  "$src_root/util/sep_coefficient_physics.cpp" \
  "$src_root/util/sep_coefficient_registry.cpp" \
  "$src_root/util/sep_background_snapshot.cpp" \
  "$src_root/util/sep_test_registry.cpp" \
  "$src_root/util/sep_acceptance_cases.cpp" \
  "$src_root/test/step13/test_acceptance_cases.cpp" \
  -pthread -o "$build_dir/test_acceptance_cases"

(cd "$build_dir" && ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} \
  ./test_acceptance_cases)
for report in "$build_dir/acceptance-results.json" \
              "$build_dir/acceptance-results.xml"; do
  test -s "$report" || { echo "FAIL REPORT: missing $report" >&2; exit 1; }
  for id in BG01 BG02 CROSS01 CROSS02; do
    grep -q "$id" "$report" || {
      echo "FAIL REPORT: $report omits $id" >&2
      exit 1
    }
  done
done
echo "PASS STEP13-REPORTS: JSON and JUnit retain all acceptance evidence"
