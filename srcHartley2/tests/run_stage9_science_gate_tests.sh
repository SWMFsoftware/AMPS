#!/usr/bin/env bash
# Run only the dependency-free Stage-9 parser/statistics mutation suite.  This
# is a tooling test, not a substitute for the eight transported MPI runs.  A
# real archive is evaluated only when --manifest is supplied explicitly.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
usage: run_stage9_science_gate_tests.sh [--manifest RUN_MANIFEST.json] [--output RESULT.json]

Runs deterministic synthetic mutation tests for the Stage-9 independent-seed
science gate. If --manifest is supplied, also evaluates that checksum-
protected ANALOG/FORCED archive. --output is valid only with --manifest.

A synthetic PASS validates the gate implementation only. It is never evidence
that a comet simulation is scientifically accepted.
EOF
}

manifest=
output=
while (( $# > 0 )); do
  case "$1" in
    --manifest)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      manifest=$2
      shift 2
      ;;
    --output)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      output=$2
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      printf 'unknown Stage-9 test option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ -z "$output" || -n "$manifest" ]] || {
  printf '%s\n' '--output requires --manifest' >&2
  exit 2
}

PYTHONDONTWRITEBYTECODE=1 python3 "$script_dir/test_stage9_science_gate_selftest.py"

if [[ -n "$manifest" ]]; then
  command=(python3 "$script_dir/test_stage9_science_gate.py" "$manifest")
  if [[ -n "$output" ]]; then
    command+=(--output "$output")
  fi
  PYTHONDONTWRITEBYTECODE=1 "${command[@]}"
fi

printf 'PASS Stage-9 science-gate tooling'
if [[ -n "$manifest" ]]; then
  printf ' manifest=%s' "$manifest"
fi
printf '\n'
