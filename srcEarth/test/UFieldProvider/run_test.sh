#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build_dir="${TMPDIR:-/tmp}/earth_ufieldprovider_$$"
mkdir -p "$build_dir"
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic -pthread \
  -I"$here/../.." \
  "$here/test_field_provider.cpp" \
  "$here/../../gridless/DipoleInterface.cpp" \
  -o "$build_dir/test_field_provider"

"$build_dir/test_field_provider"

# Compile the production Boris mover and trap/boundary headers against the same small
# generated-constant shim used by UTrajectoryCore.  This second executable protects
# the exact 5-location x 32-energy F4 setup without requiring a configured AMPS tree.
${CXX:-c++} -std=c++11 -O2 -Wall -Wextra -Werror -pedantic -pthread \
  -I"$here/../UTrajectoryCore" \
  -I"$here/../.." \
  "$here/test_f4_dipole_probe.cpp" \
  "$here/../../gridless/GridlessParticleMovers.cpp" \
  "$here/../../gridless/DipoleInterface.cpp" \
  -o "$build_dir/test_f4_dipole_probe"

"$build_dir/test_f4_dipole_probe"
