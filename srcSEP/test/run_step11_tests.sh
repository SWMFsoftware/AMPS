#!/bin/sh
set -eu

# Step 11's physics/ownership tests are dependency-light by design.  Running
# them under both address and undefined-behavior sanitizers catches storage and
# restart errors without requiring the enclosing AMPS application checkout.
src_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step11.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_dir/util" \
  "$src_dir/util/sep_transport_common.cpp" \
  "$src_dir/util/sep_focused_transport_core.cpp" \
  "$src_dir/util/sep_focused_transport_mfp_core.cpp" \
  "$src_dir/util/sep_turbulence_core.cpp" \
  "$src_dir/util/sep_test_registry.cpp" \
  "$src_dir/util/sep_turbulence_validation.cpp" \
  "$src_dir/test/step11/test_turbulence.cpp" \
  -o "$build_dir/test_turbulence"

(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_turbulence)
for id in TURBOWN01 TURB02 TURB03 TURB04 TURB05 TURB06 TURB07 TURB08 \
          TURB09 TURB10 TURB11 TURB12 TURB13 TURB14 TURB15 TURB16 \
          TURB17 TURB18 TURB19 TURB20 TURB21 TURB22 TURB23; do
  grep -q "$id" "$build_dir/controlled-turbulence-results.json"
  grep -q "$id" "$build_dir/controlled-turbulence-results.xml"
done

# The analytical plotting runner needs the registry metrics after this
# disposable directory is removed.  Retention is deliberately opt-in and
# limited to validated JSON/JUnit evidence; build and sanitizer products are
# never copied to the requested report directory.
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/controlled-turbulence-results.json" "$SRCSEP_REPORT_DIR/"
  cp "$build_dir/controlled-turbulence-results.xml" "$SRCSEP_REPORT_DIR/"
fi

# The dependency-light runner and native --test CLI must consume one catalog,
# rather than two implementations that happen to reuse the same case IDs.
grep -q 'ControlledTurbulenceDescriptors()' "$src_dir/component_tests.cpp"
grep -q 'util/sep_turbulence_validation.o' "$src_dir/makefile"

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -DSEP_CLI_PARSE_ONLY -I"$src_dir/util" \
  "$src_dir/util/sep_transport_common.cpp" \
  "$src_dir/util/sep_turbulence_core.cpp" \
  "$src_dir/util/sep_coefficient_physics.cpp" \
  "$src_dir/util/sep_coefficient_registry.cpp" \
  "$src_dir/util/sep_production_mover.cpp" \
  "$src_dir/util/sep_cli.cpp" \
  "$src_dir/util/sep_configuration_matrix.cpp" \
  "$src_dir/test/step11/test_turbulence_cli.cpp" \
  -o "$build_dir/test_turbulence_cli"
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_turbulence_cli"
echo "PASS TURB-CLI: source, representation, boundary, spectral, and tolerance controls"

# Source gates prevent later refactors from silently reintroducing a fourth
# mover or a second mutable wave-energy authority.
grep -q 'enum class Source' "$src_dir/util/sep_turbulence_core.h"
grep -q 'SwmfReadOnly' "$src_dir/util/sep_turbulence_core.h"
grep -q 'EnergyLedger' "$src_dir/util/sep_turbulence_core.h"
grep -q 'RemapConservatively' "$src_dir/util/sep_turbulence_core.h"
grep -q 'SerializeCheckpoint' "$src_dir/util/sep_turbulence_core.h"
grep -q 'EvolvesLocally' "$src_dir/main.cpp"
grep -q 'ActiveConfiguration().cascadeCoefficient' "$src_dir/main.cpp"
echo "PASS TURB-SOURCE: shared registry, ownership, ledger, remap, and restart contracts present"
