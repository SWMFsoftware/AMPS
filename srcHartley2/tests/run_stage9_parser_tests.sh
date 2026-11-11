#!/usr/bin/env bash
# Run the Stage-9 diagnostics input-policy tests against the tracked parser.
# The Python fixture constructs source-only temporary working trees; it never
# compiles AMPS, launches MPI, edits build/, or invokes the generated root
# hartley2.pl copy.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)

usage() {
  cat <<'EOF'
usage: run_stage9_parser_tests.sh

Tests Stage9ScienceDiagnostics default/OFF/ON parsing, invalid values, and
every prerequisite of the ON mode. No arguments are accepted and no AMPS
compile or MPI run is performed.
EOF
}

if (( $# != 0 )); then
  usage >&2
  exit 2
fi

perl -c "$source_root/input/hartley2.pl"
PYTHONDONTWRITEBYTECODE=1 python3 \
  "$script_dir/test_stage9_science_diagnostics_parser.py"

printf '%s\n' 'PASS Stage-9 science-diagnostics parser policy'
