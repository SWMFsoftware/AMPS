#!/bin/sh
set -eu

# PDB01-PDB05: compile the production srcSEP schema-6 input path: the
# parallel-diffusion binding and the application.sep.transport key
# (physics: parallel_diffusion block records, adapter, canonical registry)
# together with the shared library objects, without AMPS or MPI.  The library,
# adapter and the srcSEP input reader (generic AMPS runtime input,
# src/INPUT_FILE.md; header-only reader in src/general/input) are C++17, as in
# the production AMPS build.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
common="$source_root/../src/models/sep_common"
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-parallel-diffusion.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cxx=${CXX:-c++}
"$cxx" -std=c++17 -Wall -Wextra -Werror -pedantic \
  -I"$common" -I"$source_root/util" \
  -I"$source_root/../src/general/input" \
  -I"$source_root/../src/models/swcme" \
  "$source_root/adapters/parallel_diffusion_adapter.cpp" \
  "$source_root/util/sep_initialization.cpp" \
  "$source_root/util/sep_production_mover.cpp" \
  "$source_root/util/sep_turbulence_core.cpp" \
  "$common/sep_transport_common.cpp" \
  "$common/sep_coefficient_physics.cpp" \
  "$common/sep_coefficient_registry.cpp" \
  "$common/sep_background_snapshot.cpp" \
  "$common/sep_test_registry.cpp" \
  "$common/parallel_diffusion/parallel_diffusion.cpp" \
  "$common/parallel_diffusion/parallel_diffusion_advanced.cpp" \
  "$script_dir/parallel_diffusion_binding/test_parallel_diffusion_binding.cpp" \
  -o "$build_dir/test_parallel_diffusion_binding"

(cd "$build_dir" && ./test_parallel_diffusion_binding \
  "$source_root/examples/sep_parker_mesh_parallel_diffusion.in")
