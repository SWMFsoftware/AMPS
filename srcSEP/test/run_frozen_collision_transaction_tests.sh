#!/bin/sh
set -eu
# Compile only neutral value components into an owned temporary tree. Native
# AMPS artifacts and the user's configured build are never cleaned or changed.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-frozen-collision-transaction.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

# Record every source-owned translation unit and transitive local header used
# by this deliberately small neutral build. Start/end equality is a validity
# gate, not just metadata: a live peer edit during compilation discards this
# run before execution, and an execution-time edit discards its verdict too.
# Keep this explicit list synchronized with new local includes below.
hash_inputs() {
  python3 -B - "$src_root" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
relative=["util/sep_frozen_collision_transaction.h","util/sep_frozen_collision_transaction.cpp",
 "util/sep_frozen_collision_transaction_validation.h","util/sep_frozen_collision_transaction_validation.cpp",
 "util/sep_polarized_wave_state.h","util/sep_polarized_wave_state.cpp",
 "util/sep_resonant_collision.h","util/sep_resonant_collision.cpp","util/sep_common_header_path.h",
 "../src/models/sep_common/sep_transport_common.h","../src/models/sep_common/sep_transport_common.cpp",
 "../src/models/sep_common/sep_test_registry.h","../src/models/sep_common/sep_test_registry.cpp",
 "test/test_frozen_collision_transaction.cpp","test/run_frozen_collision_transaction_tests.sh"]
records=[{"path":str((root/path).resolve()),"sha256":hashlib.sha256((root/path).read_bytes()).hexdigest()}
         for path in relative]
print(json.dumps({"schema":"srcsep-fct-source-inputs-v1","inputs":records},sort_keys=True))
PY
}
hash_inputs > "$build_dir/before.json"
cat "$build_dir/before.json"
set +e
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Wpedantic -Werror -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/.." -I"$src_root/util" \
  "$src_root/util/sep_frozen_collision_transaction.cpp" \
  "$src_root/util/sep_frozen_collision_transaction_validation.cpp" \
  "$src_root/util/sep_polarized_wave_state.cpp" \
  "$src_root/util/sep_resonant_collision.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_frozen_collision_transaction.cpp" \
  -o "$build_dir/test_frozen_collision_transaction"
compile_status=$?
set -e
hash_inputs > "$build_dir/after-compile.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-compile.json"; then
  echo 'ERROR FCT source changed during compilation; this run is discarded' >&2
  exit 2
fi
if [ "$compile_status" -ne 0 ]; then exit "$compile_status"; fi
# LeakSanitizer cannot enumerate /proc tasks in the managed test environment;
# memory access and undefined-behavior instrumentation remain active.
set +e
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_frozen_collision_transaction"
run_status=$?
set -e
hash_inputs > "$build_dir/after-run.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-run.json"; then
  echo 'ERROR FCT source changed during execution; all verdicts from this run are discarded' >&2
  exit 2
fi
exit "$run_status"
