#!/bin/sh
set -eu

# D02 compiles the production adapter, canonical header-only resolver, and the
# native descriptor factory also used by `amps --test D02`. It needs neither
# generated PIC headers nor MPI; native application closure remains the D03
# release gate.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
amps_root=$(CDPATH= cd -- "$src_root/.." && pwd)
build_dir=${TMPDIR:-/tmp}/srcsep-d02-tests
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
  "$src_root/test/test_d02_swcme_configuration.cpp" \
  -o "$build_dir/test_d02_swcme_configuration"

"$build_dir/test_d02_swcme_configuration"
