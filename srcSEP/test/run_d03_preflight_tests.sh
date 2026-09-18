#!/bin/sh
set -eu

# D03PRE is the dependency-light half of D03: it verifies the exact native
# registry descriptor, linked three-mover catalog, canonical fingerprint, and
# monotonically refreshed SWCME state identity.  The Python D03 runner executes
# afterward and owns the external build/MPI/restart campaign contract.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
amps_root=$(CDPATH= cd -- "$src_root/.." && pwd)
build_dir=${TMPDIR:-/tmp}/srcsep-d03-preflight-tests
mkdir -p "$build_dir"

# The focused executable below proves descriptor behavior.  This source audit
# separately proves that the production ComponentTestRegistry consumes the
# descriptor factory; without it D03PRE could pass in isolation yet remain
# invisible to `amps --list-tests` and Python `--all`.
grep -q 'SwcmeImprovementDescriptors()' "$src_root/component_tests.cpp" || {
  echo "ERROR: ComponentTestRegistry does not include D01/D02/D03PRE" >&2
  exit 1
}
grep -q 'util/sep_swcme_validation.o' "$src_root/makefile" || {
  echo "ERROR: production archive omits sep_swcme_validation.o" >&2
  exit 1
}

${CXX:-c++} -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror \
  -I"$src_root/adapters" -I"$src_root/util" \
  -I"$amps_root/src/models/sep_common" \
  -I"$amps_root/src/models/swcme" \
  "$src_root/adapters/swcme1d_adapter.cpp" \
  "$src_root/util/sep_production_mover.cpp" \
  "$src_root/util/sep_swcme_validation.cpp" \
  "$amps_root/src/models/sep_common/sep_background_snapshot.cpp" \
  "$amps_root/src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_d03_native_preflight.cpp" \
  -o "$build_dir/test_d03_native_preflight"

"$build_dir/test_d03_native_preflight"
