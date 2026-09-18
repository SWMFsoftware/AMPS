#!/bin/sh
set -eu

# Compile the exact registry and production CLI parser without PIC/MPI.  This
# keeps canonical selection, capabilities, and rejection policy testable from a
# source archive while the runtime adapter is covered by the native AMPS gate.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step4.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -DSEP_CLI_PARSE_ONLY \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  "$src_root/test/step4/test_production_mover_cli.cpp" \
  "$src_root/util/sep_cli.cpp" \
  "$src_root/util/sep_configuration_matrix.cpp" \
  "$src_root/util/sep_production_mover.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/util/sep_turbulence_core.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_physics.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_registry.cpp" \
  -o "$build_dir/test_production_mover_cli"

"$build_dir/test_production_mover_cli"

# This source assertion protects the main-loop dispatch regression that the
# dependency-light binary cannot instantiate without PIC: policy must query
# MoverCapabilities rather than compare raw implementation addresses.
if grep -En 'ParticleMoverPtr[[:space:]]*[!=]=' "$src_root/main.cpp" >/dev/null; then
  echo "FAIL MOVCLI03: main loop compares ParticleMoverPtr addresses" >&2
  exit 1
fi
echo "PASS MOVCLI03: main-loop policy uses mover capabilities"
