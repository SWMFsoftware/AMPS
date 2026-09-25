#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build_dir="${TMPDIR:-/tmp}/earth_ustandaloneproducts_$$"
mkdir -p "$build_dir"
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

# Warnings are errors.  This target compiles the exact production header used by
# main.cpp and Mode3D.cpp and deliberately has no MPI/PIC/SPICE/Geopack dependency.
${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  "$here/test_standalone_products.cpp" \
  -o "$build_dir/test_standalone_products"

"$build_dir/test_standalone_products"
python3 "$here/test_step7_source_contract.py"

