#!/bin/sh
set -eu
# Build in an owned temporary directory; source/native AMPS artifacts are never
# cleaned. Pin the complete source-owned closure across compile and execution:
# a concurrent edit invalidates a verdict instead of reporting a mixed revision.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-polarized-wkb.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
hash_inputs() {
  python3 -B - "$src_root" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
relative=["util/sep_polarized_wkb_transport.h","util/sep_polarized_wkb_transport.cpp",
 "util/sep_exact_nonnegative_sum.h","util/sep_exact_nonnegative_sum.cpp",
 "util/sep_polarized_wkb_validation.h","util/sep_polarized_wkb_validation.cpp",
 "util/sep_polarized_wave_state.h","util/sep_polarized_wave_state.cpp",
 "util/sep_resonant_collision.h","util/sep_resonant_collision.cpp","util/sep_common_header_path.h",
 "../src/models/sep_common/sep_transport_common.h","../src/models/sep_common/sep_transport_common.cpp",
 "../src/models/sep_common/sep_test_registry.h","../src/models/sep_common/sep_test_registry.cpp",
 "test/test_polarized_wkb.cpp","test/run_polarized_wkb_tests.sh"]
records=[{"path":str((root/path).resolve()),"sha256":hashlib.sha256((root/path).read_bytes()).hexdigest()}
         for path in relative]
print(json.dumps({"schema":"srcsep-wkb-source-inputs-v1","inputs":records},sort_keys=True))
PY
}
hash_inputs > "$build_dir/before.json"
cat "$build_dir/before.json"
set +e
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Wpedantic -Werror -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/.." -I"$src_root/util" \
  "$src_root/util/sep_polarized_wkb_transport.cpp" \
  "$src_root/util/sep_exact_nonnegative_sum.cpp" \
  "$src_root/util/sep_polarized_wkb_validation.cpp" \
  "$src_root/util/sep_polarized_wave_state.cpp" \
  "$src_root/util/sep_resonant_collision.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_polarized_wkb.cpp" -o "$build_dir/test_polarized_wkb"
compile_status=$?
set -e
hash_inputs > "$build_dir/after-compile.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-compile.json"; then
  echo 'ERROR WKB source changed during compilation; this run is discarded' >&2
  exit 2
fi
if [ "$compile_status" -ne 0 ]; then exit "$compile_status"; fi
# Leak detection is disabled for this runner. The claimed sanitizer scope is
# AddressSanitizer memory-access checks and UndefinedBehaviorSanitizer only.
set +e
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_polarized_wkb"
run_status=$?
set -e
hash_inputs > "$build_dir/after-run.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-run.json"; then
  echo 'ERROR WKB source changed during execution; all verdicts are discarded' >&2
  exit 2
fi
exit "$run_status"
