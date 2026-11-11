#!/usr/bin/env bash
# The comparison tests use an authenticated synthetic seed matrix and a tiny
# synthetic paper table. They validate tooling only and never constitute AMPS
# or Rubin science evidence.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/stage9_comparison_selftest.py"
printf '%s\n' 'PASS Stage-9 ensemble comparison tooling'
