#!/bin/sh
set -eu

# WP31--WP41 share a strict dependency-light gate for deterministic numerical
# contracts.  External native/SWMF/observational/performance campaigns remain
# separate fail-closed targets because this runner cannot manufacture their
# executables, inputs, decompositions, or hardware identity.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-wp31-wp41.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -O1 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_background_snapshot.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_physics.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_registry.cpp" \
  "$src_root/util/sep_production_mover.cpp" \
  "$src_root/util/sep_turbulence_core.cpp" \
  "$src_root/util/sep_population_control.cpp" \
  "$src_root/util/sep_configuration_matrix.cpp" \
  "$src_root/util/sep_evidence.cpp" \
  "$src_root/util/sep_system_ledger.cpp" \
  "$src_root/util/sep_validation_tools.cpp" \
  "$src_root/util/sep_observation_forward_model.cpp" \
  "$src_root/test/test_wp31_wp41.cpp" \
  -o "$build_dir/test_wp31_wp41"

ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} "$build_dir/test_wp31_wp41"

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
"$build_dir/print_configuration_matrix" > "$build_dir/configuration-matrix.md"
[ "$(grep -c '^| ' "$build_dir/configuration-matrix.md")" -eq 92 ]
echo "PASS WP35-DOC: generated matrix contains header, separator, and 90 registry rows"

# Source checks prove the dependency-light contracts reach the production
# configuration and driver seams.  They are SourceIntegration evidence only;
# the native harness target below still requires a linked AMPS executable.
grep -q 'PlanAdvance' "$src_root/turbulence_production_adapter.cpp" || \
  grep -q 'Turbulence::Advance' "$src_root/turbulence_production_adapter.cpp"
grep -q 'populationControl.spatialBins' "$src_root/main_lib.cpp"
grep -q 'ConfigurationMatrix::Preflight' "$src_root/util/sep_cli.cpp"
grep -q 'PICAdapter::Advance' "$src_root/main.cpp"
grep -q 'PICAdapter::Advance' "$src_root/main_lib.cpp"
flush_owners=$(grep -R 'FlushWaveContributions[[:space:]]*(' \
  "$src_root"/*.cpp "$src_root"/util/*.cpp | \
  grep -v 'void[[:space:]]\+FlushWaveContributions' | wc -l | tr -d ' ')
if [ "$flush_owners" -ne 1 ]; then
  echo "FAIL WP41-SOURCE: expected one queue-flush owner, found $flush_owners" >&2
  exit 1
fi
echo "PASS WP31-WP41-SOURCE: production configuration and single-owner seams verified"
