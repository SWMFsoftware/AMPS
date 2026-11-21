#!/usr/bin/env bash
# Focused tests of the chemistry variants (srcHartley2/CHEMISTRY_VARIANTS.md):
# the rate-table module against an independent Python evaluation (test 1),
# and the build-input preprocessor checks.  This runner does not configure
# AMPS, mutate build/, or launch a simulation.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
test_root=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-chemistry-variants.XXXXXX")
trap 'rm -rf -- "$test_root"' EXIT

${CXX:-g++} -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I"$source_root/srcHartley2" \
  "$script_dir/test_chemistry_variants.cpp" \
  "$source_root/srcHartley2/CometChemistryVariants.cpp" \
  -o "$test_root/test_chemistry_variants"

PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_chemistry_variants.py" \
  "$test_root/test_chemistry_variants"
perl -I"$source_root/input" "$source_root/input/test_comet_chemistry_variants.pl"

printf '%s\n' 'PASS chemistry-variants focused test group'
