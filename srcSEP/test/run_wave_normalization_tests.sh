#!/bin/sh
set -eu

src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-wave-normalization.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir" -I"$src_dir/util" -I"$src_dir/../src/models/sep_common" \
  "$src_dir/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_dir/../src/models/sep_common/sep_coefficient_physics.cpp" \
  "$src_dir/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_dir/util/sep_turbulence_core.cpp" \
  "$src_dir/util/sep_focused_transport_core.cpp" \
  "$src_dir/util/sep_focused_transport_mfp_core.cpp" \
  "$src_dir/util/sep_wave_normalization_validation.cpp" \
  "$src_dir/test/rng/test_wave_normalization.cpp" \
  -o "$build_dir/test_wave_normalization"

(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_wave_normalization)

# Routing signatures supplement the numerical callbacks; the SBWN native
# accessor fixtures are the separate runtime evidence. A count of two direct
# conversions in coefficient_providers became invalid when raw SWMF totals
# acquired their explicit read-only half-energy adapter.
# The lookups use POSIX grep with -F (the symbol is plain text, not a regular
# expression), not ripgrep: rg is not installed on every test machine, and a
# missing tool used to stop this script with exit 127 before WN-SOURCE was
# reported (every nightly from 2026-11-17 to 2026-11-27). A symbol that is
# absent prints which routing check failed and returns 1.
require_symbol() {
  if ! grep -qF -- "$1" "$2"; then
    echo "FAIL WN-SOURCE: '$1' not found in $2" >&2
    exit 1
  fi
}
require_symbol 'CP::MagneticVarianceFromEnergyDensity' "$src_dir/coefficient_providers.cpp"
require_symbol 'ReadStandardBatsrusWaveView' "$src_dir/coefficient_providers.cpp"
require_symbol 'InterpolateMagneticView' "$src_dir/standard_batsrus_wave_adapter.cpp"
require_symbol 'CP::MagneticVarianceFromEnergyDensity' "$src_dir/util/sep_standard_batsrus_waves.cpp"
require_symbol 'CoefficientPhysics::MagneticVarianceFromEnergyDensity' "$src_dir/util/sep_turbulence_core.cpp"
echo 'PASS WN-SOURCE: local magnetic conversion and typed SWMF-total half-energy routing use the common magnetic formula'

if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/wave-normalization-results.json" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir/wave-normalization-results.xml" "$SRCSEP_REPORT_DIR/"
fi
