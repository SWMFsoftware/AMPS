#!/bin/sh
set -eu

# Build in a disposable directory so a focused parser/registry check never
# leaves objects or executables in the source tree.  SEP_CLI_PARSE_ONLY excludes
# only the function that writes production globals; parsing and help text are
# compiled from the exact production CLI source.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step1.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cxx=${CXX:-c++}
"$cxx" -std=c++11 -Wall -Wextra -Werror -DSEP_CLI_PARSE_ONLY \
  -I"$source_root/util" -I"$source_root/../src/models/sep_common" \
  "$source_root/util/sep_cli.cpp" \
  "$source_root/util/sep_configuration_matrix.cpp" \
  "$source_root/util/sep_production_mover.cpp" \
  "$source_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$source_root/util/sep_turbulence_core.cpp" \
  "$source_root/../src/models/sep_common/sep_coefficient_physics.cpp" \
  "$source_root/../src/models/sep_common/sep_coefficient_registry.cpp" \
  "$source_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$script_dir/step1/test_registry_cli.cpp" \
  -o "$build_dir/test_registry_cli"

(cd "$build_dir" && ./test_registry_cli)
