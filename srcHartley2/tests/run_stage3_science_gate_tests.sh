#!/usr/bin/env bash
# Run the permanent Stage-3 science-gate tests.  The suite creates complete
# synthetic eight-run archives under a temporary directory; it neither starts
# AMPS/MPI nor reads or modifies a production validation archive.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage3_science_gate.py"
