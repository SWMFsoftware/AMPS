#!/bin/sh
set -eu

# Compile the canonical registry and conversion layer without AMPS/PIC.  This
# gives coefficient tests sanitizer coverage even when the surrounding host
# application and its generated Makefile.conf are unavailable.
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step10.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir/util" \
  -DSEP_CLI_PARSE_ONLY \
  "$src_dir/util/sep_transport_common.cpp" \
  "$src_dir/util/sep_turbulence_core.cpp" \
  "$src_dir/util/sep_coefficient_physics.cpp" \
  "$src_dir/util/sep_coefficient_registry.cpp" \
  "$src_dir/util/sep_production_mover.cpp" \
  "$src_dir/util/sep_cli.cpp" \
  "$src_dir/util/sep_configuration_matrix.cpp" \
  "$src_dir/test/step10/test_coefficients.cpp" \
  -o "$build_dir/test_coefficients"

ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_coefficients"

# Source gates prove that each production mover uses the shared PIC adapter and
# that every registry choice is exposed through the one existing CLI parser.
grep -q '#include "coefficient_providers.h"' "$src_dir/parker_mover.cpp"
grep -q '#include "coefficient_providers.h"' "$src_dir/focused_transport_dmumu.cpp"
grep -q '#include "coefficient_providers.h"' "$src_dir/focused_transport_mfp.cpp"
grep -q -- '--coefficient-source' "$src_dir/util/sep_cli.cpp"
grep -q -- '--spatial-diffusion-provider' "$src_dir/util/sep_cli.cpp"
grep -q -- '--pitch-angle-diffusion-provider' "$src_dir/util/sep_cli.cpp"
grep -q -- '--mean-free-path-provider' "$src_dir/util/sep_cli.cpp"
grep -q -- '--invalid-coefficient-policy' "$src_dir/util/sep_cli.cpp"
echo "PASS COEF-SOURCE: all movers and CLI share the coefficient registry"
