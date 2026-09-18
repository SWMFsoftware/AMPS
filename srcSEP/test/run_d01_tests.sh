#!/bin/sh
set -eu

# D01 is deliberately compiled outside AMPS. The test links the real adapter,
# canonical SWCME header implementation, canonical background snapshot store,
# and the same descriptor factory inserted into ComponentTestRegistry().  It
# therefore exercises the production status/policy path *and* native registry
# callback without requiring a generated PIC tree or an MPI launcher.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
amps_root=$(CDPATH= cd -- "$src_root/.." && pwd)
build_dir=${TMPDIR:-/tmp}/srcsep-d01-tests
mkdir -p "$build_dir"

${CXX:-c++} -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror \
  -I"$src_root/adapters" -I"$src_root/util" \
  -I"$amps_root/src/models/sep_common" \
  -I"$amps_root/src/models/swcme" \
  "$src_root/adapters/swcme1d_adapter.cpp" \
  "$src_root/util/sep_production_mover.cpp" \
  "$src_root/util/sep_swcme_validation.cpp" \
  "$amps_root/src/models/sep_common/sep_background_snapshot.cpp" \
  "$amps_root/src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_d01_swcme_fail_closed.cpp" \
  -o "$build_dir/test_d01_swcme_fail_closed"

"$build_dir/test_d01_swcme_fail_closed"
