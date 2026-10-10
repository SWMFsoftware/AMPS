#!/bin/sh
set -eu

# PDB01-PDB05: compile the production srcSEP schema-4 input path: the
# parallel-diffusion binding and the run.particle_mover key
# (schema-4 section location, adapter, canonical registry) together with the
# shared library objects, without AMPS or MPI.  The library and adapter are
# C++17; the srcSEP INI reader is compiled with the same C++17 flag here, as in
# the production AMPS build.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
common="$source_root/../src/models/sep_common"
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-parallel-diffusion.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cxx=${CXX:-c++}
"$cxx" -std=c++17 -Wall -Wextra -Werror -pedantic \
  -I"$common" -I"$source_root/util" \
  "$source_root/adapters/parallel_diffusion_adapter.cpp" \
  "$source_root/util/sep_initialization.cpp" \
  "$source_root/util/sep_production_mover.cpp" \
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
