#!/bin/sh
set -eu
# Bounded reference/protocol checks only. The broad study is opt-in and is
# never started by this runner. Concurrent source edits discard the verdict.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-finite-amplitude-tests.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
hash_inputs() {
  python3 -B - "$src_root" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
relative=["validation/physics/finite_amplitude_response.py","validation/physics/run_finite_amplitude_study.py","validation/physics/plot_finite_amplitude_study.py",
 "validation/physics/controlled_finite_amplitude_pilot.json","validation/physics/controlled_ep_maxwellian_baseline.json",
 "validation/physics/warm_parallel_maxwellian.py","test/test_finite_amplitude_response.py",
 "test/run_finite_amplitude_study_tests.sh"]
print(json.dumps({"schema":"srcsep-finite-amplitude-test-inputs-v1","inputs":[
 {"path":str(root/path),"sha256":hashlib.sha256((root/path).read_bytes()).hexdigest()} for path in relative]},sort_keys=True))
PY
}
hash_inputs > "$build_dir/before.json"
cat "$build_dir/before.json"
set +e
python3 -B "$src_root/test/test_finite_amplitude_response.py" -v
run_status=$?
set -e
hash_inputs > "$build_dir/after.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after.json"; then
  echo 'ERROR finite-amplitude source changed during execution; all verdicts are discarded' >&2
  exit 2
fi
exit "$run_status"
