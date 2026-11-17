#!/bin/sh
set -eu
# Compile in an owned temporary directory. A source edit during compilation or
# execution invalidates the verdict, so a parallel development session cannot
# accidentally report a mixed source revision as verified.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-material-tube.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
hash_inputs() {
  python3 -B - "$src_root" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
relative=["util/sep_material_tube_kinematics.h","util/sep_material_tube_kinematics.cpp",
 "util/sep_material_tube_kinematics_validation.h","util/sep_material_tube_kinematics_validation.cpp",
 "util/sep_common_header_path.h","../src/models/sep_common/sep_transport_common.h",
 "../src/models/sep_common/sep_transport_common.cpp","../src/models/sep_common/sep_test_registry.h",
 "../src/models/sep_common/sep_test_registry.cpp","test/test_material_tube_kinematics.cpp",
 "test/run_material_tube_kinematics_tests.sh"]
records=[{"path":str((root/path).resolve()),"sha256":hashlib.sha256((root/path).read_bytes()).hexdigest()}
         for path in relative]
print(json.dumps({"schema":"srcsep-material-tube-inputs-v1","inputs":records},sort_keys=True))
PY
}
hash_inputs > "$build_dir/before.json"
cat "$build_dir/before.json"
set +e
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Wpedantic -Werror -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/.." -I"$src_root/util" \
  "$src_root/util/sep_material_tube_kinematics.cpp" \
  "$src_root/util/sep_material_tube_kinematics_validation.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_material_tube_kinematics.cpp" -o "$build_dir/test_material_tube"
compile_status=$?
set -e
hash_inputs > "$build_dir/after-compile.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-compile.json"; then
  echo 'ERROR material source changed during compilation; this run is discarded' >&2
  exit 2
fi
if [ "$compile_status" -ne 0 ]; then exit "$compile_status"; fi
# Leak detection is disabled. Evidence claims ASan memory-access and UBSan
# instrumentation only; no unobserved environment explanation is inferred.
set +e
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_material_tube"
run_status=$?
set -e
hash_inputs > "$build_dir/after-run.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-run.json"; then
  echo 'ERROR material source changed during execution; all verdicts are discarded' >&2
  exit 2
fi
exit "$run_status"
