#!/usr/bin/env bash
# Run the dependency-free tests for test_runner.py scheduling metadata and pool
# exclusivity.  The script deliberately does not require AMPS, MPI, or compiler
# toolchains, so the main validation list can execute it as an independent gate.
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 "$script_dir/test_exclusive_runner.py"
python3 "$script_dir/test_interrupted_runner.py"
python3 "$script_dir/test_deferred_last_pass.py"
python3 "$script_dir/earth_test_runtime_test.py"
python3 "$script_dir/test_failure_reporting.py"
# Storm-extremum timing criteria of C9 (minimum cutoff latitude) and C10
# (maximum boundary expansion), both built on
# earth_test_runtime.evaluate_extremum_time().  Offline fixtures only: no AMPS,
# MPI, aacgmv2, or matplotlib.
python3 "$script_dir/../C9/test_minimum_time.py"
python3 "$script_dir/../C10/tests/test_expansion_time.py"
