#!/usr/bin/env bash
# Run the deterministic Stage-10 release-boundary mutation suite.  The Python
# test creates and destroys only its own temporary repository in a unique
# stage10-readiness-selftest-* child of the AMPS root; it never changes
# the caller's branch, refs, index, worktree, or validation evidence.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/stage10_merge_readiness_selftest.py"

printf '%s\n' 'PASS Stage-10 merge-readiness tooling'
