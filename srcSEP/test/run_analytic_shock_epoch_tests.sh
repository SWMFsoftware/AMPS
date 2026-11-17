#!/bin/sh
set -eu
# Build the exact trajectory/view callbacks in an owned temporary directory.
# This source gate checks explicit-time mathematics; actual host-clock/source
# ordering must pass the separate fresh-process HOST campaign.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-analytic-shock-epoch.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
hash_inputs() {
  # Refuse mixed-revision evidence even if another agent edits a dependency
  # while this bounded source gate compiles or runs. The native host campaign
  # separately records its staged executable and cannot inherit this result.
  python3 -B - "$src_root" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
relative=["util/sep_analytic_shock_epoch.h","util/sep_analytic_shock_epoch.cpp",
 "util/sep_analytic_shock_epoch_validation.h","util/sep_analytic_shock_epoch_validation.cpp",
 "util/sep_shock_source_core.h","util/sep_shock_source_core.cpp",
 "util/sep_common_header_path.h","../src/models/sep_common/sep_transport_common.h",
 "../src/models/sep_common/sep_transport_common.cpp","../src/models/sep_common/sep_test_registry.h",
 "../src/models/sep_common/sep_test_registry.cpp","test/test_analytic_shock_epoch.cpp",
 "test/run_analytic_shock_epoch_tests.sh"]
records=[{"path":str((root/path).resolve()),"sha256":hashlib.sha256((root/path).read_bytes()).hexdigest()} for path in relative]
print(json.dumps({"schema":"srcsep-analytic-shock-epoch-source-inputs-v1","inputs":records},sort_keys=True))
PY
}
hash_inputs > "$build_dir/before.json"
cat "$build_dir/before.json"
set +e
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Wpedantic -Werror -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  "$src_root/util/sep_analytic_shock_epoch.cpp" "$src_root/util/sep_analytic_shock_epoch_validation.cpp" \
  "$src_root/util/sep_shock_source_core.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/test/test_analytic_shock_epoch.cpp" -o "$build_dir/test_analytic_shock_epoch"
compile_status=$?
set -e
hash_inputs > "$build_dir/after-compile.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-compile.json"; then
  echo 'ERROR analytic shock epoch source changed during compilation; evidence discarded' >&2
  exit 2
fi
if [ "$compile_status" -ne 0 ]; then exit "$compile_status"; fi
set +e
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_analytic_shock_epoch)
run_status=$?
set -e
hash_inputs > "$build_dir/after-run.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-run.json"; then
  echo 'ERROR analytic shock epoch source changed during execution; all verdicts discarded' >&2
  exit 2
fi
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/analytic-shock-epoch-results.json" "$build_dir/analytic-shock-epoch-results.xml" \
     "$build_dir/before.json" "$SRCSEP_REPORT_DIR/"
fi
exit "$run_status"
