#!/bin/bash
#======================================================================================
# make_run_root.sh -- create a run root (outside srcGeospace) in which an individual
# test runner can be started by hand against a chosen amps executable.
#
# srcGeospace's end-to-end tests are srcEarth's tests.  Their runners expect to be
# started from a directory that contains the application trees and the executable as
# ./amps: they call ./amps, refer to srcEarth/test/... (srcGeospace/test/... for
# srcGeospace's own unit test) by relative path, and write test_output/ into the current
# directory.  A run root reproduces that layout with symbolic links:
#
#   <root>/srcEarth    -> the repository srcEarth   (tests and runners, read only)
#   <root>/srcGeospace -> this srcGeospace tree     (its own unit test, read only)
#   <root>/amps        -> <executable under test>
#
# test/tools/run_srcearth_tests.sh creates the same layout for a whole list and starts
# srcEarth's test_runner.py there; this script is for running one runner directly:
#
#   srcGeospace/scripts/make_run_root.sh $GEOSPACE_RUNS/results/c1_gpu $GEOSPACE_RUNS/builds/gpu_double/amps
#   cd $GEOSPACE_RUNS/results/c1_gpu
#   TMPDIR=$PWD/tmp PYTHONDONTWRITEBYTECODE=1 python3 srcEarth/test/C1/run_C1.py --mode gridless -np 4 -nt 16
#
# Usage: make_run_root.sh <root> <amps-executable>
#======================================================================================
set -euo pipefail
[ $# -eq 2 ] || { echo "usage: $0 <root> <amps-executable>" >&2; exit 2; }
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
REPO=$(cd "$SRC_DIR/.." && pwd)
ROOT="$1"; EXE=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
[ -x "$EXE" ] || { echo "not executable: $EXE" >&2; exit 2; }
mkdir -p "$ROOT"
ROOT=$(cd "$ROOT" && pwd)
case "$ROOT/" in
  "$REPO/"*) echo "refusing to create a run root inside the repository $REPO" >&2; exit 2;;
esac
for f in srcEarth srcGeospace amps; do
  if [ -e "$ROOT/$f" ] && [ ! -L "$ROOT/$f" ]; then
    echo "refusing to replace $ROOT/$f (not a symbolic link)" >&2; exit 2
  fi
done
ln -sfn "$REPO/srcEarth" "$ROOT/srcEarth"
ln -sfn "$SRC_DIR" "$ROOT/srcGeospace"
ln -sfn "$EXE" "$ROOT/amps"
mkdir -p "$ROOT/tmp"
echo "run root $ROOT: amps -> $EXE"
