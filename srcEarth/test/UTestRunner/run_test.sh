#!/usr/bin/env bash
# Run the dependency-free tests for test_runner.py scheduling metadata and pool
# exclusivity.  The script deliberately does not require AMPS, MPI, or compiler
# toolchains, so the main validation list can execute it as an independent gate.
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 "$script_dir/test_exclusive_runner.py"
python3 "$script_dir/test_deferred_last_pass.py"
