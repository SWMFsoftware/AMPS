#!/bin/sh
set -eu

# Compile into a disposable directory so documentation generation never leaves
# a binary in the source archive.  The emitted Markdown is deterministic and
# comes directly from the same classifier used by CLI/run preflight.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-config-matrix.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
${CXX:-c++} -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_background_snapshot.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_physics.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_registry.cpp" \
  "$src_root/util/sep_production_mover.cpp" \
  "$src_root/util/sep_turbulence_core.cpp" \
  "$src_root/util/sep_configuration_matrix.cpp" \
  "$src_root/test/print_configuration_matrix.cpp" \
  -o "$build_dir/print_configuration_matrix"
"$build_dir/print_configuration_matrix"
