#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

# These tests use temporary, hash-pinned evidence bundles.  Bytecode is disabled so
# the validation gate never modifies a shared source checkout.
if PYTHONDONTWRITEBYTECODE=1 python3 "$here/test_release_validation.py" && \
   PYTHONDONTWRITEBYTECODE=1 python3 "$here/test_step12_source_contract.py"; then
  echo "RESULT: PASS"
else
  status=$?
  echo "RESULT: FAIL" >&2
  exit "$status"
fi
