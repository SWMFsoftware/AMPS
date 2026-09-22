#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build_dir="${TMPDIR:-/tmp}/earth_ufieldprovider_$$"
mkdir -p "$build_dir"
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

# Warnings are errors.  The test links the production dipole implementation and uses
# the production FieldProvider header; it has no MPI/PIC/SPICE/Geopack dependency.
${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic -pthread \
  -I"$here/../.." \
  "$here/test_field_provider.cpp" \
  "$here/../../gridless/DipoleInterface.cpp" \
  -o "$build_dir/test_field_provider"

"$build_dir/test_field_provider"
