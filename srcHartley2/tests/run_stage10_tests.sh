#!/usr/bin/env bash
# Run the dependency-free Stage-10 release-evidence mutation suites.
#
# These tests exercise JSON receipt serialization and the fail-closed,
# read-only merge-readiness collector.  They use temporary files and temporary
# Git repositories only.  A PASS is tooling evidence; it is neither Stage-9
# comet science acceptance nor permission to modify master.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)

for required_tool in python3 git sha256sum; do
  if ! command -v -- "$required_tool" >/dev/null 2>&1; then
    printf 'required Stage-10 test tool is not on PATH: %s\n' \
      "$required_tool" >&2
    exit 2
  fi
done

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage10_all_stage_result.py"

"$source_root/srcHartley2/validation/run_stage10_merge_readiness_tests.sh"

printf '%s\n' \
  'PASS Stage-10 receipt and merge-readiness tooling (not science or merge authorization)'
