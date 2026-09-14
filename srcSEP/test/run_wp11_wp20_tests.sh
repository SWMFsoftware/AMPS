#!/bin/sh
set -eu

# Build the exact dependency-light transport and coefficient kernels under the
# project's strict C++11 warning policy and ASan/UBSan.  The disposable build
# directory ensures this gate never leaves binaries or evidence in the source
# archive.  Native PIC-adapter execution remains a separate AMPS gate.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-wp11-wp20.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -O1 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" \
  "$src_root/util/sep_transport_common.cpp" \
  "$src_root/util/sep_coefficient_physics.cpp" \
  "$src_root/util/sep_species_source.cpp" \
  "$src_root/util/sep_coefficient_registry.cpp" \
  "$src_root/util/sep_focused_transport_core.cpp" \
  "$src_root/test/test_wp11_wp20.cpp" \
  -o "$build_dir/test_wp11_wp20"

ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} "$build_dir/test_wp11_wp20"

# Source-level gates complement the dependency-light executable by proving the
# repaired kernels are wired into production adapters and legacy compatibility
# entry points.  These are integration claims, not substitutes for a linked
# native AMPS run.
grep -q 'EvaluateAdaptiveKappa' "$src_root/coefficient_providers.cpp"
grep -q 'BuildSourceBoundInput' "$src_root/coefficient_providers.cpp"
grep -q 'EvaluateConstantDmumu' "$src_root/diffusion.cpp"
grep -q 'EvaluateJokipiiSlab' "$src_root/diffusion.cpp"
grep -q 'D=D_mu_mu' "$src_root/diffusion.cpp"
grep -q 'P_s_minus(-Omega/t1' "$src_root/diffusion.cpp"
grep -q 'ActiveNumericalTolerances' "$src_root/parker_mover.cpp"
grep -q 'ActiveNumericalTolerances' "$src_root/focused_transport_dmumu.cpp"
grep -q 'ActiveNumericalTolerances' "$src_root/focused_transport_mfp.cpp"
if grep -Eq 'QLT1::B[[:space:]]*=' "$src_root/sampling.cpp"; then
  echo "FAIL WP19-SOURCE: sampling still uses a hard-coded QLT1 field" >&2
  exit 1
fi
echo "PASS WP11-WP20-SOURCE: production adapters use repaired source-bound contracts"
