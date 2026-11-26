#!/bin/sh
# Dependency-light unit tests of the srcSEP3D additions for the 2013-04-11
# 3-D real-event test (test_event3d_units.cpp, E3D01-E3D06): computation
# frame, focused-mover inertial terms, ephemeris interpolation, field-line
# tracing and the observer flux tube, the absorbing shock, and the
# 1-D-equivalent tube-line source.  The canonical sources are compiled
# directly with sanitizers; no AMPS/PIC build is needed.  A failing check
# exits non-zero.
#
# Usage: sh run_unit_tests.sh   (CXX selects the compiler, default c++)
# Each check prints PASS or FAIL (README.md section 5).  The binary is built
# in a temporary directory that is removed on exit; nothing else is written.
set -eu

# Source roots: srcSEP3D itself, the shared sep_common kernels, and the
# coronal-CME model (local-compression DSA law, focused escape law, SHA-256).
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
sep3d=$(CDPATH= cd -- "$here/../.." && pwd)
models="$sep3d/../src/models"
common="$models/sep_common"
cme="$models/sep_coronal_cme"
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep3d-event-2013.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

# One translation unit list: the srcSEP3D sources under test, the shared
# sources they call, and the test driver.  Warnings are errors.
# spice_ephemeris.cpp is compiled without SEP3D_WITH_SPICE, i.e. as the
# typed refusal of a nospice build; E3D07 covers the SPICE path.
${CXX:-c++} -std=c++17 -O1 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$sep3d" -I"$sep3d/core" -I"$sep3d/runtime" -I"$sep3d/mesh" \
  -I"$sep3d/transport" -I"$sep3d/adapters" -I"$sep3d/output" \
  -I"$sep3d/background" -I"$common" -I"$cme/include" -I"$cme/src" \
  -I"$models/swcme" -I"$models/sep_corona_swcme/include" \
  "$sep3d/transport/focused_transport.cpp" \
  "$sep3d/transport/keyed_random.cpp" \
  "$sep3d/transport/perpendicular_transport.cpp" \
  "$sep3d/mesh/observer_flux_tube.cpp" \
  "$sep3d/output/observer_trajectory.cpp" \
  "$sep3d/output/spice_ephemeris.cpp" \
  "$sep3d/adapters/shock_return.cpp" \
  "$sep3d/adapters/tube_line_source.cpp" \
  "$common/sep_transport_common.cpp" \
  "$common/sep_injection_spectrum.cpp" \
  "$cme/src/particle_source.cpp" \
  "$cme/src/sha256.cpp" \
  "$here/test_event3d_units.cpp" \
  -o "$build_dir/test_event3d_units"
# Run under the sanitizers with LeakSanitizer off (detect_leaks=0).
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_event3d_units"
