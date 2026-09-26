#!/usr/bin/env bash
# Roadmap Step 10 dependency-free numerical and production-wiring tests.
#
# This runner intentionally compiles only the shared contract.  It therefore catches
# cadence, naming, Shue-boundary, and status defects even on a machine without SWMF,
# MPI, SPICE, or a configured AMPS executable.  The linked 1x1/2x8/8x16 campaign is a
# separate acceptance action documented in README.md; it is not simulated here.

set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/amps_step10.XXXXXXXX")
cleanup() {
  # build_dir is created by mktemp immediately above and is never supplied by a user.
  rm -rf -- "$build_dir"
}
trap cleanup EXIT HUP INT TERM

cxx=${CXX:-c++}
"$cxx" -std=c++17 -O2 -Wall -Wextra -pedantic \
  "$script_dir/test_swmf_coupled_access.cpp" \
  -o "$build_dir/test_swmf_coupled_access"

"$build_dir/test_swmf_coupled_access"
PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_compare_cutoff_access.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_legacy_cutoff_consumers.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_step10_source_contract.py"

echo "RESULT: PASS"
