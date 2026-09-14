#!/bin/sh
set -eu

# Re-run every dependency-light stochastic transport/reduction campaign with
# its recorded keyed seeds.  The numerical evidence must be byte-identical,
# while wall-clock duration is deliberately normalized because scheduler load
# is not part of the physics.  The anchored expression changes only the
# registry result header; configuration entries such as duration_s=1 remain in
# the evidence and must still compare exactly.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-repeat.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

for step in 7 8 9 12; do
  "$src_root/test/run_step${step}_tests.sh" |
    sed -E '/^\[[A-Z0-9]+\] (PASS|FAIL|SKIP|ERROR) duration_s=/ s/duration_s=[0-9.]+/duration_s=<runtime>/' \
    >"$build_dir/step${step}.first"
  "$src_root/test/run_step${step}_tests.sh" |
    sed -E '/^\[[A-Z0-9]+\] (PASS|FAIL|SKIP|ERROR) duration_s=/ s/duration_s=[0-9.]+/duration_s=<runtime>/' \
    >"$build_dir/step${step}.second"
  if ! cmp -s "$build_dir/step${step}.first" "$build_dir/step${step}.second"; then
    echo "FAIL REPEAT${step}: fixed-seed evidence changed between runs" >&2
    diff -u "$build_dir/step${step}.first" \
            "$build_dir/step${step}.second" >&2 || true
    exit 1
  fi
  echo "PASS REPEAT${step}: fixed-seed evidence is byte-identical"
done
