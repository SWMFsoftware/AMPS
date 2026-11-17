#!/bin/sh
set -eu
# Compile the exact explicit radial kernel and public geometry builder. Native
# field-line primitives/history are a separate ARBN gate. Build/run only in an
# owned directory and discard mixed-revision evidence after any source edit.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-analytic-radial.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
hash_inputs() {
  python3 -B - "$src_root" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root=Path(sys.argv[1])
relative=["util/sep_analytic_radial_background.h","util/sep_analytic_radial_background.cpp",
 "util/sep_analytic_radial_validation.h","util/sep_analytic_radial_validation.cpp",
 "util/sep_background_derivatives.h","util/sep_initialization.h","util/sep_initialization.cpp",
 "util/sep_production_mover.h","util/sep_production_mover.cpp","util/sep_turbulence_core.h","util/sep_turbulence_core.cpp",
 "util/sep_common_header_path.h","../src/models/sep_common/sep_transport_common.h",
 "../src/models/sep_common/sep_transport_common.cpp","../src/models/sep_common/sep_test_registry.h",
 "../src/models/sep_common/sep_test_registry.cpp","../src/models/sep_common/sep_input_layout.h",
 "../src/models/sep_common/sep_background_snapshot.h","../src/models/sep_common/sep_background_snapshot.cpp",
 "../src/general/input/amps_input.h","../src/models/swcme/swcme1d_input.hpp",
 "test/test_analytic_radial_background.cpp","test/run_analytic_radial_background_tests.sh"]
records=[{"path":str((root/path).resolve()),"sha256":hashlib.sha256((root/path).read_bytes()).hexdigest()} for path in relative]
print(json.dumps({"schema":"srcsep-analytic-radial-source-inputs-v1","inputs":records},sort_keys=True))
PY
}
hash_inputs > "$build_dir/before.json"
cat "$build_dir/before.json"
set +e
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Wpedantic -Werror -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  -I"$src_root/../src/general/input" -I"$src_root/../src/models/swcme" \
  "$src_root/util/sep_analytic_radial_background.cpp" "$src_root/util/sep_analytic_radial_validation.cpp" \
  "$src_root/util/sep_initialization.cpp" "$src_root/util/sep_production_mover.cpp" \
  "$src_root/util/sep_turbulence_core.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/../src/models/sep_common/sep_test_registry.cpp" \
  "$src_root/../src/models/sep_common/sep_background_snapshot.cpp" \
  "$src_root/test/test_analytic_radial_background.cpp" -o "$build_dir/test_analytic_radial"
compile_status=$?
set -e
hash_inputs > "$build_dir/after-compile.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-compile.json"; then
  echo 'ERROR analytic radial source changed during compilation; evidence discarded' >&2
  exit 2
fi
if [ "$compile_status" -ne 0 ]; then exit "$compile_status"; fi
set +e
(cd "$build_dir" && ASAN_OPTIONS=detect_leaks=0 ./test_analytic_radial)
run_status=$?
set -e
hash_inputs > "$build_dir/after-run.json"
if ! cmp -s "$build_dir/before.json" "$build_dir/after-run.json"; then
  echo 'ERROR analytic radial source changed during execution; all verdicts discarded' >&2
  exit 2
fi
if [ -n "${SRCSEP_REPORT_DIR:-}" ]; then
  mkdir -p "$SRCSEP_REPORT_DIR"
  cp "$build_dir/analytic-radial-results.json" "$build_dir/analytic-radial-results.xml" \
     "$build_dir/before.json" "$SRCSEP_REPORT_DIR/"
fi
exit "$run_status"
