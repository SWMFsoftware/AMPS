#!/bin/sh
set -eu

# Compile the exact production parser/geometry and native callbacks without
# AMPS or MPI.  The linked executable registers the same INIT01-INIT03 records,
# so this fast gate and `amps --all-tests` cannot diverge semantically.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-initialization.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cxx=${CXX:-c++}
"$cxx" -std=c++11 -Wall -Wextra -Werror \
  -I"$source_root/util" -I"$source_root/../src/models/sep_common" \
  "$source_root/util/sep_initialization.cpp" \
  "$source_root/util/sep_initialization_validation.cpp" \
  "$source_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$source_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$script_dir/test_initialization.cpp" \
  -o "$build_dir/test_initialization"

(cd "$build_dir" && ./test_initialization \
  "$source_root/examples/sep_parker_mesh.in")
