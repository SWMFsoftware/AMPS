#!/bin/sh
set -eu

# WP42--WP64 use one strict dependency-light executable for the mathematical,
# ownership, transaction, restart, and evidence contracts.  Native AMPS/SWMF,
# MPI, observational, and multi-node performance evidence remains in explicit
# fail-closed runners below; this source suite cannot promote itself to those
# higher evidence levels.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-wp42-wp64.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -O1 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" \
  "$src_root/util/sep_transport_common.cpp" \
  "$src_root/util/sep_turbulence_core.cpp" \
  "$src_root/util/sep_configuration_matrix.cpp" \
  "$src_root/util/sep_coefficient_registry.cpp" \
  "$src_root/util/sep_coefficient_physics.cpp" \
  "$src_root/util/sep_production_mover.cpp" \
  "$src_root/util/sep_focused_transport_core.cpp" \
  "$src_root/util/sep_focused_transport_mfp_core.cpp" \
  "$src_root/util/sep_evidence.cpp" \
  "$src_root/util/sep_runtime_contracts.cpp" \
  "$src_root/util/sep_physics_extensions.cpp" \
  "$src_root/util/sep_numerical_extensions.cpp" \
  "$src_root/util/sep_validation_extensions.cpp" \
  "$src_root/test/test_wp42_wp64.cpp" \
  -o "$build_dir/test_wp42_wp64"

ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} "$build_dir/test_wp42_wp64"

# Production seams are checked separately from pure numerical tests.  These
# assertions prevent the two P0 defects identified by the review from being
# reintroduced: per-step State construction and direct legacy wave mutation.
grep -q 'TurbulenceRuntimeStore' "$src_root/turbulence_production_adapter.cpp"
if grep -q 'State state;' "$src_root/turbulence_production_adapter.cpp"; then
  echo "FAIL WP42-SOURCE: production adapter still constructs per-step state" >&2
  exit 1
fi
if grep -q 'WaveParticleCouplingManager' "$src_root/turbulence_production_adapter.cpp"; then
  echo "FAIL WP43-SOURCE: production adapter still calls legacy coupling manager" >&2
  exit 1
fi
grep -q 'DrainWaveContributions' "$src_root/turbulence_production_adapter.cpp"
grep -q 'CampaignRandomSeed' "$src_root/field_line.cpp"
grep -q 'StableSourceIdentity' "$src_root/field_line.cpp"
grep -q 'relative_normal_speed_m_s' "$src_root/field_line.cpp"
if grep -q 'n_sw_end=AMPS2SWMF::ShockData.*DownStreamDensity' \
    "$src_root/field_line.cpp"; then
  echo "FAIL WP45-SOURCE: downstream density still normalizes upstream flux" >&2
  exit 1
fi
grep -q 'emitted.preWaveMomentumKgMPerS' \
  "$src_root/focused_transport_mfp.cpp"
grep -q 'VelocityGradientInput' "$src_root/transport_common.cpp"
if grep -Eq '^[[:space:]]*goto[[:space:]]+end' "$src_root/sampling.cpp"; then
  echo "FAIL SAMPLING-SOURCE: root-only output still jumps across C++ locals" >&2
  exit 1
fi
grep -q 'PIC::ThisThread==0' "$src_root/sampling.cpp"
echo "PASS WP42-WP46-SOURCE: persistent owner, common coupling, derivatives, and source identity seams"

# Higher evidence gates are executable contracts, not silent skips.  The
# optional native runner returns BLOCKED when no linked AMPS command is supplied.
"$src_root/test/run_wp59_wp64_native_gates.sh"
