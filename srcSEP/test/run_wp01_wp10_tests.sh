#!/bin/sh
set -eu

# This bounded runner validates the dependency-free numerical contracts.  The
# production adapter itself requires a configured AMPS/MPI build and is checked
# structurally below; a native executable remains a separate validation gate.
cxx=${CXX:-c++}
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

"$cxx" -std=c++11 -Wall -Wextra -Werror -pedantic \
  -I. -Iutil -I../src/models/sep_common \
  test/test_wp01_wp10.cpp \
  ../src/models/sep_common/sep_transport_common.cpp \
  ../src/models/sep_common/sep_background_snapshot.cpp \
  util/sep_focused_transport_core.cpp \
  util/sep_focused_transport_mfp_core.cpp \
  util/sep_reproducible_reduction.cpp \
  -o "$build_dir/wp01_wp10"
"$build_dir/wp01_wp10"

# Both standalone and coupled callers reach the same fixed physical-step
# owner in main_lib.cpp. Requiring a second driver call would reintroduce a
# double wave update. Reuse the stricter common audit: one executable library
# call inside amps_time_step(), none in main.cpp, plus duplicate/missing-owner
# negative controls. The numerical groups above and remaining source-record
# checks below are still required and propagate failures unchanged.
python3 -B test/check_stage1_stage2_contracts.py
grep -q 'if (false && SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag)' main.cpp
grep -q 'QueueWaveContribution' parker_mover.cpp
grep -q 'QueueWaveContribution' focused_transport_dmumu.cpp
grep -q 'QueueWaveContribution' focused_transport_mfp.cpp
! grep -q 'particlePointer' transport_common.cpp
! grep -q 'densityEvolutionIntervalS' transport_common.cpp
printf '%s\n' 'PASS WP01/WP03 production entry and escape-safe source contracts'
