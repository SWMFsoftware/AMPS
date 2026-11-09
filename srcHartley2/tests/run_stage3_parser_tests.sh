#!/usr/bin/env bash
# Exercise the Stage-3 random-seed and residual-sampling input records in
# source-only temporary trees. No generated build/, assembled input, AMPS
# executable, or machine-local configuration is needed.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_stage3_parser.py"
