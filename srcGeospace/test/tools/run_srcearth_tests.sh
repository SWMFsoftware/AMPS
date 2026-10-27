#!/bin/bash
#======================================================================================
# run_srcearth_tests.sh -- run a test list with srcEarth's test runner against a chosen
# amps executable (normally a srcGeospace build), in a run root outside the source tree.
#
# srcGeospace has no copy of srcEarth's tests: its end-to-end tests are srcEarth's
# (srcEarth/test/list and the per-test runners next to it), run with srcEarth's runner
# srcEarth/test/test_runner.py.  The runners call ./amps, refer to srcEarth/test/... (and
# srcGeospace/test/... for srcGeospace's own tests) by relative path and write
# test_output/, logs and reports into the current directory.  This script creates that
# layout as a run root with symbolic links and starts the runner there:
#
#   <run-root>/srcEarth    -> <repository>/srcEarth      (tests and runner, read only)
#   <run-root>/srcGeospace -> <repository>/srcGeospace   (its own tests, read only)
#   <run-root>/amps        -> <executable under test>
#   <run-root>/tmp/                                      (TMPDIR of every test)
#
# Usage:
#   srcGeospace/test/tools/run_srcearth_tests.sh --run-root DIR --amps EXE
#       [--list LIST] [--only REGEX]... [--exclude REGEX]...
#       [--np N] [--nt M] [--np-nt-exclude REGEX] [-- RUNNER-OPTIONS...]
#
#   --list LIST      test list (default: srcEarth/test/list; srcGeospace/test/list holds
#                    srcGeospace's own unit test)
#   --only/--exclude select commands (test/tools/sublist.py)
#   --np/--nt        replace the values of existing -np/-nt tokens (sublist.py; C13 and
#                    C18 keep theirs unless --np-nt-exclude says otherwise), e.g.
#                    --np 1 --nt 32 for one GPU
#   RUNNER-OPTIONS   passed to srcEarth/test/test_runner.py, e.g. -j 6 --no-start-progress
#
# With --only/--exclude/--np/--nt the selected list is written to <run-root>/list.
# The runner always gets --commit-id (the run root is not a git work tree) and
# --last-pass-results-file inside the run root; --update-last-pass is refused (it would
# rewrite srcEarth/test/list).  At the end the runner prints how to apply the saved
# results with --commit-last-pass: do not do that for srcEarth's list from a srcGeospace
# run (the last-pass lines of srcEarth/test/list record srcEarth's own passes).
#
# Examples (from the repository root):
#   R=$GEOSPACE_RUNS/results/<run>
#   srcGeospace/test/tools/run_srcearth_tests.sh --run-root $R/cpu_double \
#       --amps $GEOSPACE_RUNS/builds/cpu_double/amps -- -j 6
#   GEOSPACE_AMPS=$GEOSPACE_RUNS/builds/gpu_double/amps GEOSPACE_GPUS=0 \
#   srcGeospace/test/tools/run_srcearth_tests.sh --run-root $R/gpu1_double \
#       --amps srcGeospace/scripts/amps_gpu_select.sh --exclude 'test/SA0[123]/' \
#       --np 1 --nt 32 -- -j 1
#   srcGeospace/test/tools/run_srcearth_tests.sh --run-root $R/unit \
#       --amps $GEOSPACE_RUNS/builds/cpu_double/amps --list srcGeospace/test/list
#======================================================================================
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GEO=$(cd "$HERE/../.." && pwd)
REPO=$(cd "$GEO/.." && pwd)

ROOT=""; EXE=""; LIST="$REPO/srcEarth/test/list"
SEL=(); RUNNER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --run-root) ROOT="$2"; shift 2;;
    --amps) EXE="$2"; shift 2;;
    --list) LIST="$2"; shift 2;;
    --only|--exclude|--np|--nt|--np-nt-exclude) SEL+=("$1" "$2"); shift 2;;
    --) shift; RUNNER=("$@"); break;;
    -h|--help) sed -n 2,46p "$0"; exit 0;;
    *) echo "unknown argument: $1 (runner options go after --)" >&2; exit 2;;
  esac
done
[ -n "$ROOT" ] && [ -n "$EXE" ] || { sed -n 2,46p "$0"; exit 2; }
for a in ${RUNNER[@]+"${RUNNER[@]}"}; do
  case "$a" in
    --update-last-pass|--commit-last-pass)
      echo "refusing $a: it would rewrite srcEarth/test/list" >&2; exit 2;;
  esac
done

EXE=$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")
[ -x "$EXE" ] || { echo "not executable: $EXE" >&2; exit 2; }
LIST=$(cd "$(dirname "$LIST")" && pwd)/$(basename "$LIST")
[ -f "$LIST" ] || { echo "no test list: $LIST" >&2; exit 2; }

mkdir -p "$ROOT"
ROOT=$(cd "$ROOT" && pwd)
case "$ROOT/" in
  "$REPO/"*) echo "refusing to use $ROOT as a run root: it is inside the repository" >&2; exit 2;;
esac
for f in srcEarth srcGeospace amps; do
  if [ -e "$ROOT/$f" ] && [ ! -L "$ROOT/$f" ]; then
    echo "refusing to replace $ROOT/$f (not a symbolic link)" >&2; exit 2
  fi
done
ln -sfn "$REPO/srcEarth" "$ROOT/srcEarth"
ln -sfn "$GEO" "$ROOT/srcGeospace"
ln -sfn "$EXE" "$ROOT/amps"
mkdir -p "$ROOT/tmp"

if [ ${#SEL[@]} -gt 0 ]; then
  python3 -B "$HERE/sublist.py" "$LIST" "${SEL[@]}" -o "$ROOT/list"
  LIST="$ROOT/list"
fi

COMMIT=$(git -C "$REPO" rev-parse HEAD)
echo "[run_srcearth_tests] run root $ROOT: amps -> $EXE; list $LIST; commit $COMMIT"
cd "$ROOT"
export TMPDIR="$ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
exec python3 -B srcEarth/test/test_runner.py "$LIST" --commit-id "$COMMIT" \
    --last-pass-results-file "$ROOT/.last_pass_results.json" ${RUNNER[@]+"${RUNNER[@]}"}
