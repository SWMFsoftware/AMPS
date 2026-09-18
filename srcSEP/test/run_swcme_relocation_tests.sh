#!/bin/sh
set -eu

# SWCME-R1 relocation and copied-build path gate.
#
# This test mirrors srcSEP3D's BLDL3D05 contract. It verifies both filesystem
# layouts used by AMPS:
#
#   AMPS/srcSEP/makefile   (source and dependency-light tests)
#   AMPS/build/main/makefile (production application copy)
#
# Both must resolve the same AMPS root, Makefile.conf, and canonical SWCME
# directory without consulting the process working directory. It also rejects
# any application-local SWCME source/header tree, which would create a second
# apparent owner and make builds depend on which -I option happened to win.

src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
amps_root=$(CDPATH= cd -- "$src_root/.." && pwd)
swcme_dir=${SWCME_DIR:-"$amps_root/src/models/swcme"}

test -d "$swcme_dir" || {
  echo "FAIL SWCME-R1: canonical directory is absent: $swcme_dir" >&2
  echo "Set SWCME_DIR only when testing a detached srcSEP source tree." >&2
  exit 2
}
swcme_dir=$(CDPATH= cd -- "$swcme_dir" && pwd)

for header in swcme1d.hpp swcme_sep_interface.hpp; do
  test -f "$swcme_dir/$header" || {
    echo "FAIL SWCME-R1: missing canonical public header: $header" >&2
    exit 1
  }
done

for obsolete in "$src_root/swcme" "$src_root/swcme----moved-out" \
                "$src_root/shared_model_compat" \
                "$src_root/swcme1d.hpp"; do
  test ! -e "$obsolete" || {
    echo "FAIL SWCME-R1: application-local SWCME tree remains: $obsolete" >&2
    exit 1
  }
done

test ! -e "$src_root/demo1d.cpp" || {
  echo "FAIL SWCME-R1: obsolete application-local SWCME demo remains" >&2
  exit 1
}

# Production and focused-validation consumers use public header names. A path
# containing swcme/ would again bind them to a directory beneath srcSEP.
if grep -R -n -E '#include[[:space:]]*[<"][^>"]*swcme/' \
    "$src_root" --include='*.h' --include='*.hpp' --include='*.cpp'; then
  echo "FAIL SWCME-R1: application-local SWCME include path detected" >&2
  exit 1
fi

# Mirror srcSEP3D's dependency boundary: the umbrella header is provider-free,
# and only the private adapter implementation sees the canonical SWCME header.
if grep -n -E 'swcme1d::|#include[[:space:]]*[<"]swcme' \
    "$src_root/sep.h"; then
  echo "FAIL SWCME-R1: sep.h leaks the SWCME implementation dependency" >&2
  exit 1
fi
grep -q '#include "swcme1d.hpp"' \
  "$src_root/adapters/swcme1d_adapter.cpp"
if grep -n -E 'swcme1d::|#include[[:space:]]*[<"]swcme' \
    "$src_root/adapters/swcme1d_adapter.h"; then
  echo "FAIL SWCME-R1: adapter API leaks a canonical SWCME type/header" >&2
  exit 1
fi
grep -q '#include "swcme_sep_interface.hpp"' \
  "$src_root/test/step15/test_swcme_srcsep_integration.cpp"
for token in SWCME_ARCHIVE SWCME_OBJECTS SWCME_ADAPTER_CXXFLAGS \
             'adapters/%.o' audit-production-swcme; do
  grep -q "$token" "$src_root/makefile" || {
    echo "FAIL SWCME-R1: makefile omits canonical contract token $token" >&2
    exit 1
  }
done

fixture=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-swcme-layout.XXXXXX")
trap 'rm -rf "$fixture"' EXIT HUP INT TERM
mkdir -p "$fixture/srcSEP" "$fixture/build/main" \
         "$fixture/src/models/swcme"
: > "$fixture/Makefile.conf"
cp "$src_root/makefile" "$fixture/srcSEP/makefile"
cp "$src_root/makefile" "$fixture/build/main/makefile"

expected_root=$(CDPATH= cd -- "$fixture" && pwd)
expected_swcme="$expected_root/src/models/swcme"
for makefile in "$fixture/srcSEP/makefile" \
                "$fixture/build/main/makefile"; do
  # The outer SWCME_DIR may intentionally point at a detached canonical tree
  # for the header compile below. Do not allow that override to mask the path
  # discovery being tested inside the synthetic installed layout.
  output=$(cd / && env -u AMPS_ROOT -u AMPS_CONFIG -u SWCME_DIR \
    -u MAKEFLAGS -u MAKEOVERRIDES -u MFLAGS -u MAKELEVEL \
    make --no-print-directory -f "$makefile" print-layout-paths)
  printf '%s\n' "$output" | grep -Fx "AMPS_ROOT=$expected_root" >/dev/null || {
    echo "FAIL SWCME-R1: incorrect AMPS_ROOT from $makefile" >&2
    printf '%s\n' "$output" >&2
    exit 1
  }
  printf '%s\n' "$output" | grep -Fx \
    "AMPS_CONFIG=$expected_root/Makefile.conf" >/dev/null || {
    echo "FAIL SWCME-R1: incorrect AMPS_CONFIG from $makefile" >&2
    printf '%s\n' "$output" >&2
    exit 1
  }
  printf '%s\n' "$output" | grep -Fx "SWCME_DIR=$expected_swcme" >/dev/null || {
    echo "FAIL SWCME-R1: incorrect SWCME_DIR from $makefile" >&2
    printf '%s\n' "$output" >&2
    exit 1
  }
  printf '%s\n' "$output" | grep -Fx \
    "MODEL_INCLUDE_FLAGS=-I$expected_swcme" >/dev/null || {
    echo "FAIL SWCME-R1: canonical include flag is absent from $makefile" >&2
    printf '%s\n' "$output" >&2
    exit 1
  }
done

# Finally compile a tiny consumer through only the canonical -I root. This
# catches a forwarding header or implicit current-directory dependency that a
# source-text path check cannot detect.
${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$swcme_dir" "$src_root/test/step15/test_swcme_header_path.cpp" \
  -o "$fixture/test_swcme_header_path"
"$fixture/test_swcme_header_path"

# The public adapter header must compile with no SWCME include path at all.
printf '%s\n' \
  '#include "swcme1d_adapter.h"' \
  'int main() {' \
  '  (void)SEP::SW1DAdapter::Scenario::Fast;' \
  '  return 0;' \
  '}' > "$fixture/test_adapter_boundary.cpp"
${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$src_root/adapters" "$fixture/test_adapter_boundary.cpp" \
  -o "$fixture/test_adapter_boundary"
"$fixture/test_adapter_boundary"

# Conversely, the one provider-facing implementation compiles only when the
# canonical model root is supplied, just like srcSEP3D's adapters/%.o rule.
${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$src_root/adapters" -I"$src_root/util" -I"$swcme_dir" \
  -c "$src_root/adapters/swcme1d_adapter.cpp" \
  -o "$fixture/swcme1d_adapter.o"

echo "SWCME-R1 PASS: srcSEP3D-style private adapter uses canonical src/models/swcme"
