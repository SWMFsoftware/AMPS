#!/bin/sh
set -eu

# Compile the actual SWCME 1-D SEP adapter and srcSEP transport/background
# components together.  This target is deliberately independent of AMPS so a
# source archive can verify the coupling boundary under ASan/UBSan.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
amps_root=$(CDPATH= cd -- "$src_root/.." && pwd)

# Match the production makefile and srcSEP3D runner: the normal installed path
# is AMPS/src/models/swcme, while SWCME_DIR supports an intentionally detached
# source checkout. Never fall back to srcSEP/swcme, because that would hide an
# incomplete relocation behind a private provider copy.
swcme_dir=${SWCME_DIR:-"$amps_root/src/models/swcme"}
test -d "$swcme_dir" || {
  echo "ERROR: canonical SWCME directory is absent: $swcme_dir" >&2
  echo "       Install it at AMPS/src/models/swcme or set SWCME_DIR." >&2
  exit 2
}
swcme_dir=$(CDPATH= cd -- "$swcme_dir" && pwd)

build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step15-swcme.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++17 -O1 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" -I"$swcme_dir" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/util/sep_parker_core.cpp" \
  "$src_root/../src/models/sep_common/sep_background_snapshot.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/step15/test_swcme_srcsep_integration.cpp" \
  -pthread -o "$build_dir/test_swcme_srcsep_integration"

(cd "$build_dir" && ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} \
  ./test_swcme_srcsep_integration)
for report in "$build_dir/step15-swcme-results.json" \
              "$build_dir/step15-swcme-results.xml"; do
  test -s "$report" || { echo "FAIL VAL04-REPORT: missing $report" >&2; exit 1; }
  grep -q "VAL04" "$report" || {
    echo "FAIL VAL04-REPORT: $report omits VAL04" >&2
    exit 1
  }
done

# An enclosing campaign can request durable copies without making ordinary
# focused tests leave generated evidence inside the source tree.
if test -n "${STEP15_REPORT_DIR:-}"; then
  mkdir -p "$STEP15_REPORT_DIR"
  cp "$build_dir/step15-swcme-results.json" \
     "$build_dir/step15-swcme-results.xml" "$STEP15_REPORT_DIR/"
fi
echo "PASS VAL04-REPORT: real SWCME states consumed by srcSEP transport"
