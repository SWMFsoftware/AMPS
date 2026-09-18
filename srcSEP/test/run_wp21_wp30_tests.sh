#!/bin/sh
set -eu

# WP21--WP30 are deliberately exercised through dependency-free SI kernels.
# ASan/UBSan and the project's strict warning policy turn undefined behavior,
# implicit conversions, and unused compatibility remnants into gate failures.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-wp21-wp30.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -O1 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_physics.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_registry.cpp" \
  "$src_root/util/sep_production_mover.cpp" \
  "$src_root/util/sep_turbulence_core.cpp" \
  "$src_root/../src/models/sep_common/sep_injection_spectrum.cpp" \
  "$src_root/util/sep_shock_source_core.cpp" \
  "$src_root/util/sep_flux_tube_geometry_core.cpp" \
  "$src_root/util/sep_sampling_products.cpp" \
  "$src_root/util/sep_transactional_output.cpp" \
  "$src_root/util/sep_run_configuration.cpp" \
  "$src_root/util/sep_population_control.cpp" \
  "$src_root/util/sep_configuration_matrix.cpp" \
  "$src_root/test/test_wp21_wp30.cpp" \
  -o "$build_dir/test_wp21_wp30"

ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} "$build_dir/test_wp21_wp30"

# These source assertions verify that the pure contracts are reached by the
# production adapters.  They complement rather than impersonate a linked AMPS
# integration run, which remains a separately reported native gate.
grep -q 'SEP::Shock::IntersectSphere' "$src_root/shock_analytical_model.cpp"
grep -q 'QueueShockContribution' "$src_root/shock_analytical_model.cpp"
grep -q 'pendingShockPlusJ' "$src_root/turbulence_production_adapter.cpp"
grep -q 'SEP::Injection::InverseCdf' "$src_root/field_line.cpp"
grep -q 'GetMagneticFluxWb' "$src_root/flux_tube_geometry.cpp"
grep -q 'luminal/superluminal particle excluded' "$src_root/sampling.cpp"
if grep -q 'maxSamplingSpeed' "$src_root/sampling.cpp" "$src_root/sep.h"; then
  echo "FAIL WP27-SOURCE: a sampling path still caps invalid particle speed" >&2
  exit 1
fi
grep -q 'WriteTransactional' "$src_root/sampling_output.cpp"
grep -q 'RunConfiguration fingerprint=' "$src_root/main.cpp"
if grep -Eq '\bsystem[[:space:]]*\(' "$src_root/sampling_output.cpp" "$src_root/sep.h"; then
  echo "FAIL WP29-SOURCE: sampling output still invokes a shell" >&2
  exit 1
fi
if grep -Eq '\bsprintf[[:space:]]*\(' "$src_root/sampling_output.cpp"; then
  echo "FAIL WP29-SOURCE: sampling output still uses unchecked sprintf" >&2
  exit 1
fi
echo "PASS WP21-WP30-SOURCE: production adapters reach the new contracts"
