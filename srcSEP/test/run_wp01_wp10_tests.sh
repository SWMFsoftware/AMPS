#!/bin/sh
set -eu

# This bounded runner validates the dependency-free numerical contracts.  The
# production adapter itself requires a configured AMPS/MPI build and is checked
# structurally below; a native executable remains a separate validation gate.
cxx=${CXX:-c++}
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

"$cxx" -std=c++11 -Wall -Wextra -Werror -pedantic \
  -I. -Iutil \
  test/test_wp01_wp10.cpp \
  util/sep_transport_common.cpp \
  util/sep_background_snapshot.cpp \
  util/sep_focused_transport_core.cpp \
  util/sep_focused_transport_mfp_core.cpp \
  util/sep_reproducible_reduction.cpp \
  -o "$build_dir/wp01_wp10"
"$build_dir/wp01_wp10"

grep -q 'Turbulence::PICAdapter::Advance' main.cpp
grep -q 'Turbulence::PICAdapter::Advance' main_lib.cpp
grep -q 'if (false && SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag)' main.cpp
grep -q 'QueueWaveContribution' parker_mover.cpp
grep -q 'QueueWaveContribution' focused_transport_dmumu.cpp
grep -q 'QueueWaveContribution' focused_transport_mfp.cpp
! grep -q 'particlePointer' transport_common.cpp
! grep -q 'densityEvolutionIntervalS' transport_common.cpp
printf '%s\n' 'PASS WP01/WP03 production entry and escape-safe source contracts'
