#!/usr/bin/env bash
# Dependency-free focused tests for the DIM1 brightness geometry, units, and
# strict input policy.  This runner does not configure AMPS, mutate build/, or
# launch a science simulation.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
test_root=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-brightness.XXXXXX")
trap 'rm -rf -- "$test_root"' EXIT

${CXX:-g++} -std=c++11 -O2 -Wall -Wextra -pedantic \
  -I"$source_root" \
  "$script_dir/test_comet_brightness_math.cpp" \
  -o "$test_root/test_comet_brightness_math"
"$test_root/test_comet_brightness_math"

make -C "$source_root/src/models/resonance_fluorescence" test

perl -I"$source_root/input" -c "$source_root/input/hartley2.pl"
perl "$source_root/input/test_comet_brightness_config.pl"
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_comet_brightness_parser.py"

printf '%s\n' 'PASS comet-brightness focused test group tests=7'
