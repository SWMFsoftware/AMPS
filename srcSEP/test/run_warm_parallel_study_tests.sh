#!/bin/sh
set -eu
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PYTHONDONTWRITEBYTECODE=1 python3 "$src_root/test/test_warm_parallel_maxwellian.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$src_root/test/test_warm_study_protocol.py"
