#!/bin/sh
set -eu

# This target is intentionally independent of MPI/PIC.  It validates the exact
# SI geometry/source core used by production while remaining runnable from the
# source archive before an enclosing AMPS application has been configured.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-step3.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$src_root/util" -I"$src_root/../src/models/sep_common" \
  "$src_root/test/step3/test_flux_tube_geometry.cpp" \
  "$src_root/../src/models/sep_common/sep_transport_common.cpp" \
  "$src_root/util/sep_flux_tube_geometry_core.cpp" \
  -o "$build_dir/test_flux_tube_geometry"

# LeakSanitizer cannot inspect /proc while this browser sandbox traces child
# processes.  AddressSanitizer and UndefinedBehaviorSanitizer remain active;
# native CI may override ASAN_OPTIONS to enable leak detection as well.
ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} "$build_dir/test_flux_tube_geometry"
