#!/bin/sh
set -eu

# Build outside the source tree so the focused Step 2 gate never leaves object,
# coverage, or executable artifacts in the deliverable.  This unit target uses
# the exact production snapshot implementation and exercises its concurrent
# read/publication contract without requiring the enclosing AMPS/MPI checkout.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step2.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cxx=${CXX:-c++}
"$cxx" -std=c++11 -Wall -Wextra -Werror -pthread \
  -I"$source_root/../src/models/sep_common" \
  "$source_root/../src/models/sep_common/sep_background_snapshot.cpp" \
  "$script_dir/step2/test_background_snapshot.cpp" \
  -o "$build_dir/test_background_snapshot"

"$build_dir/test_background_snapshot"

# The second half of the focused gate verifies that no independent elapsed-time
# variable or direct access to PIC's mutable TimeCounter has returned to the
# production files touched by Step 2.  PIC::SimulationTime::Get() is permitted
# only in the single runtime adapter.
if grep -nE 'double[[:space:]]+SimulationTime|t_since_launch' \
    "$source_root/main.cpp"; then
  echo "FAIL: a second standalone simulation clock remains in main.cpp" >&2
  exit 1
fi

if grep -n 'PIC::SimulationTime::TimeCounter' \
    "$source_root/main.cpp" "$source_root/main_lib.cpp" \
    "$source_root/sep.cpp" "$source_root/shock_analytical_model.cpp"; then
  echo "FAIL: production code bypasses the authoritative clock adapter" >&2
  exit 1
fi

clock_reads=$(grep -R 'PIC::SimulationTime::Get()' \
    "$source_root/main.cpp" "$source_root/main_lib.cpp" \
    "$source_root/sep.cpp" "$source_root/shock_analytical_model.cpp" \
    "$source_root/util/sep_background_runtime.cpp" | wc -l | tr -d ' ')
if [ "$clock_reads" -ne 1 ]; then
  echo "FAIL: expected exactly one PIC clock read adapter, found $clock_reads" >&2
  exit 1
fi

echo "Step 2 single-clock source contract: PASS"
