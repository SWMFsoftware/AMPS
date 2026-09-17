#!/bin/bash
# ==============================================================================
# test/individual-test/run_harn02.sh
#
# Verifies that the Stage-1 runner exits NON-ZERO when a test returns Fail.
# HARN02 is a deliberately failing test; running it alone must produce exit 1.
#
# Run from the srcSEP3D/ directory:
#   bash test/individual-test/run_harn02.sh
#
# Exit codes:
#   0  the runner correctly exited 1 for a Fail-only run
#   1  the runner exited unexpectedly (wrong exit code)
#   2  the test binary could not be found
# ==============================================================================
# NOTE: do NOT use 'set -e' here.
# We deliberately run a test that exits non-zero.  set -e would abort the
# script before RC=$? is captured.

RUNNER="./test/stage1"

if [ ! -x "$RUNNER" ]; then
  echo "ERROR: $RUNNER not found or not executable." >&2
  echo "       Build it first: python3 test/run_tests.py --rebuild" >&2
  exit 2
fi

# Run the binary with --test HARN_FAIL_BEACON only.
# HARN02 always returns Fail, so the runner should exit 1.
# Redirect stdout and stderr to /dev/null; we only care about the exit code.
"$RUNNER" --test HARN_FAIL_BEACON >/dev/null 2>&1
RC=$?

if [ "$RC" -eq 1 ]; then
  echo "PASS: runner exited 1 (non-zero) when HARN_FAIL_BEACON (Fail) was the only test."
  exit 0
else
  echo "FAIL: runner exited $RC instead of 1 when HARN_FAIL_BEACON (Fail) was the only test." >&2
  echo "      This means Fail is not being reported as a non-zero exit code." >&2
  echo "      Run  ./test/stage1 --test HARN_FAIL_BEACON  directly to see the output." >&2
  exit 1
fi
