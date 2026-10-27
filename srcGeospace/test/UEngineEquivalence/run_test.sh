#!/usr/bin/env bash
# UEngineEquivalence -- bitwise equivalence of batched-engine components with the
# legacy (per-trajectory) implementations they replace:
#
#   * gpu/EngineTrap.h (fixed-capacity arena) vs util/TrajectoryTrapDetector.h
#     (std::vector state), 440 synthetic bounce/drift configurations;
#   * gpu/EngineMovers.h vs gridless/GridlessParticleMovers.cpp (StepParticleChecked)
#     for BORIS, HC4, RK2/4/6, GC2/4/6 and HYBRID in a tilted dipole.
#
# Builds into $TMPDIR (or /tmp); nothing is written into the source tree.  Compiled
# with -O3 without -march (as AMPS) and -ffp-contract=off, i.e. the same contraction
# rule as the AMPS -std=c++17 build.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src_dir="$(cd "$script_dir/../.." && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/geospace_engine_equivalence.XXXXXX")"
trap 'rm -rf "$build_dir"' EXIT
CXX="${CXX:-g++}"
FLAGS="-std=c++17 -O3 -ffp-contract=off -Wall"

echo "== trap detector (gpu/EngineTrap.h vs util/TrajectoryTrapDetector.h)"
$CXX $FLAGS -I"$src_dir" "$script_dir/test_trap_equivalence.cpp" -o "$build_dir/trap"
"$build_dir/trap"

echo "== movers (gpu/EngineMovers.h vs gridless/GridlessParticleMovers.cpp)"
$CXX $FLAGS -I"$src_dir" -I"$src_dir/gridless" -I"$src_dir/../src/general" \
    "$script_dir/test_mover_equivalence.cpp" "$src_dir/gridless/GridlessParticleMovers.cpp" \
    -o "$build_dir/movers"
"$build_dir/movers"
echo "UEngineEquivalence: ALL TESTS PASSED"
