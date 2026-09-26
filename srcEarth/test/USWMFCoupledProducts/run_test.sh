#!/usr/bin/env bash
# Roadmap Step 11 dependency-free numerical, manifest, replay, and wiring tests.
#
# This fast suite deliberately does not pretend to execute a linked SWMF component.
# The required live/replay, MPI/OpenMP/scheduler, F6/F7/F13/F17, and O3 campaigns are
# specified in README.md.  Here we test the exact production contracts plus analytic
# product references on any host with a C++17 compiler and Python 3.

set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/amps_step11.XXXXXXXX")
cleanup() {
  # build_dir is a private path returned by mktemp immediately above.
  rm -rf -- "$build_dir"
}
trap cleanup EXIT HUP INT TERM

cxx=${CXX:-c++}
"$cxx" -std=c++17 -O2 -Wall -Wextra -Werror -pedantic \
  -I"$script_dir/../UBoundaryProducts" \
  "$script_dir/test_swmf_coupled_products.cpp" \
  -o "$build_dir/test_swmf_coupled_products"

"$build_dir/test_swmf_coupled_products"
PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_compare_flux_products.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_step11_source_contract.py"

echo "RESULT: PASS"
