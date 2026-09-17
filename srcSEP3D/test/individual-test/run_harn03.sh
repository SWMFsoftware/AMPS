#!/bin/bash
# ==============================================================================
# test/individual-test/run_harn03.sh
#
# Verifies that the Stage-1 runner exits ZERO when the only selected test
# returns Skip.  (Skip is not a failure — exit code must be 0.)
#
# Run from the srcSEP3D/ directory:
#   bash test/individual-test/run_harn03.sh
#
# Exit codes:
#   0  the runner correctly exited 0 for a Skip-only run
#   1  the runner exited non-zero (wrong behaviour)
#   2  the test binary could not be found
# ==============================================================================
# NOTE: do NOT use 'set -e' here.
# The binary itself exits 0 for a Skip result, but if there is any startup
# problem the binary might exit non-zero, and set -e would abort this script
# before RC=$? is captured.  We capture and test the exit code explicitly.

RUNNER="./test/stage1"

if [ ! -x "$RUNNER" ]; then
  echo "ERROR: $RUNNER not found or not executable." >&2
  echo "       Build it first: python3 test/run_tests.py --rebuild" >&2
  exit 2
fi

# Run the binary with --test HARN_SKIP_BEACON only.
# HARN03 always returns Skip, so the runner should exit 0.
# Redirect stdout and stderr to /dev/null; we only care about the exit code.
"$RUNNER" --test HARN_SKIP_BEACON >/dev/null 2>&1
RC=$?

if [ "$RC" -eq 0 ]; then
  echo "PASS: runner exited 0 when HARN_SKIP_BEACON (Skip) was the only test."
  exit 0
else
  echo "FAIL: runner exited $RC instead of 0 when HARN_SKIP_BEACON (Skip) was the only test." >&2
  echo "      This means Skip is being treated as a failure, which is wrong." >&2
  echo "      Run  ./test/stage1 --test HARN_SKIP_BEACON  directly to see the output." >&2
  exit 1
fi
