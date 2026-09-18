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
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/util/sep_turbulence_core.cpp" \
  "$src_root/util/sep_configuration_matrix.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_registry.cpp" \
  "$src_root/../src/models/sep_common/sep_coefficient_physics.cpp" \
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

# B05 does not promote prototype sources merely to preserve an obsolete
# overlay.  The manifest checker requires every WP42--WP64 item to be either an
# experimental component contract or an external gate, keeps the four
# prototype objects out of the production archive, verifies unified-runner
# discovery, and rejects historical documents that claim production status.
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$src_root/test/check_wp42_wp64_disposition.py" --root "$src_root" --self-test

# Two supporting data contracts are accepted into current code: the canonical
# Parker measure vocabulary and branch-local wave-frame event metadata.  The
# latter must be copied from each emitted event, never from enclosing shell
# endpoints that include unrelated focusing/cooling work.
grep -q 'enum class ParkerMeasure' \
  "$src_root/../src/models/sep_common/sep_transport_common.h"
grep -q 'emitted.preWaveMomentumKgMPerS' "$src_root/focused_transport_mfp.cpp"
grep -q 'emitted.postWaveMomentumKgMPerS' "$src_root/focused_transport_mfp.cpp"
grep -q 'emitted.resonantBranch' "$src_root/focused_transport_mfp.cpp"
if grep -Eq '^[[:space:]]*goto[[:space:]]+end' "$src_root/sampling.cpp"; then
  echo "FAIL SAMPLING-SOURCE: root-only output still jumps across C++ locals" >&2
  exit 1
fi
grep -q 'PIC::ThisThread==0' "$src_root/sampling.cpp"
echo "PASS B05-SOURCE: accepted data seams and experimental ownership boundary"
