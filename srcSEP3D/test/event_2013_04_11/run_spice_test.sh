#!/bin/sh
# SPICE ephemeris check of the 2013-04-11 3-D test (E3D07,
# test_spice_ephemeris.cpp): the SPICE body-observer path reproduces the
# checked-in Earth table.  Needs a CSPICE toolkit and the generic kernels:
#
#   run_spice_test.sh CSPICE_DIR LEAPSECONDS_KERNEL PLANETARY_SPK
#   e.g. run_spice_test.sh /home/vtenishe/SPICE/cspice \
#          /home/vtenishe/SPICE/Kernels/NAIF/naif0012.tls \
#          /home/vtenishe/SPICE/Kernels/NAIF/de430.bsp
#
#   CSPICE_DIR          a CSPICE toolkit (include/SpiceUsr.h, lib/cspice.a)
#   LEAPSECONDS_KERNEL  an LSK, for the UTC -> ephemeris time conversion
#   PLANETARY_SPK       an SPK with the Earth and Sun states
#
# The test prints one PASS/FAIL line with the largest position and velocity
# differences and exits non-zero on failure (README.md section 5).  CXX
# selects the compiler (default c++).  Nothing outside a temporary build
# directory is written.
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
sep3d=$(CDPATH= cd -- "$here/../.." && pwd)
cspice=$1; lsk=$2; spk=$3
# Shared model sources: the observer table checksum uses the SHA-256 of
# src/models/sep_coronal_cme.
models="$sep3d/../src/models"
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep3d-spice.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
# A meta-kernel that loads the two kernels given on the command line.
cat > "$build_dir/event.tm" <<EOF
\begindata
KERNELS_TO_LOAD = ( '$lsk', '$spk' )
\begintext
EOF
# Build the body-observer sources with the SPICE path compiled in
# (-DSEP3D_WITH_SPICE, as srcSEP3D/makefile does for an AMPS configuration
# that names a SPICE toolkit) and link them with the toolkit's cspice.a.
${CXX:-c++} -std=c++17 -O1 -Wall -Wextra -DSEP3D_WITH_SPICE \
  -I"$sep3d" -I"$sep3d/core" -I"$sep3d/runtime" -I"$sep3d/mesh" \
  -I"$sep3d/output" -I"$cspice/include" -I"$models/sep_coronal_cme/include" \
  "$sep3d/output/observer_trajectory.cpp" "$sep3d/output/spice_ephemeris.cpp" \
  "$models/sep_coronal_cme/src/sha256.cpp" "$here/test_spice_ephemeris.cpp" \
  "$cspice/lib/cspice.a" -lm -o "$build_dir/test_spice_ephemeris"
# Compare the SPICE Earth state with the checked-in table, using the
# ecliptic -> model rotation recorded in deck3d_summary.json.
"$build_dir/test_spice_ephemeris" "$build_dir/event.tm" "$here/ephemeris/earth.csv" \
  "$here/deck3d_summary.json"
