#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/sep_in_geospace_directional_access.XXXXXX")"
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

# Compile the exact production headers with strict warnings.  The suite is deliberately
# dependency-free: failures identify the Step-5 numerics/data contract rather than an
# MPI, AMPS, SPICE, or empirical-field installation problem.
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic \
  -I"$script_dir/../.." \
  "$script_dir/test_directional_access.cpp" \
  -o "$build_dir/test_directional_access"

"$build_dir/test_directional_access"

# Verify that the observation-facing consumer enforces the same complete saved-product
# schema before testing optional command-line overrides.
python3 "$script_dir/test_c19_step5_reader.py"

# Compile the real command-line parser with only the AMPS fatal-exit hook replaced by
# an exception.  This exercises option spelling, value validation, defaults, and help
# text without duplicating parser logic in the test.
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic \
  -I"$script_dir/stubs" -I"$script_dir/../.." \
  "$script_dir/test_step5_cli.cpp" \
  "$script_dir/../../util/cutoff_cli.cpp" \
  -o "$build_dir/test_step5_cli"
"$build_dir/test_step5_cli"

# Run the production C8 parser/auditor regression as U-F18.  This exercises both the
# corrected append-only 45-column product and already-generated files from the brief
# historical Step-5 column order, then proves that incomplete cubes are stopped before
# any blank-fraction reduction or comparison is attempted.
python3 "$script_dir/../C8/run_C8.py" --self-test
echo "UDirectionalAccess C8 schema/auditor: PASS (U-F18)"
