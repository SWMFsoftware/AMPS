#!/usr/bin/env bash
# Exercise only the dependency-free Stage-9 archive and summary contracts.
# No AMPS executable or MPI launcher is invoked by this test.

set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/stage9_tooling_selftest.py"
printf '%s\n' 'PASS Stage-9 production archive tooling'
