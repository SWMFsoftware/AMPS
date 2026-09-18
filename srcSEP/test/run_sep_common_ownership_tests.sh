#!/bin/sh
set -eu

# B02 canonical-kernel ownership and relocation gate.
#
# The test intentionally checks source, archive, and installed-layout
# boundaries.  A successful compile alone is insufficient: an application-
# local copy can compile cleanly while still changing physics with include
# order or introducing a second strong definition at the final AMPS link.

src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
amps_root=$(CDPATH= cd -- "$src_root/.." && pwd)
sep_common_dir=${SEP_COMMON_DIR:-"$amps_root/src/models/sep_common"}
swcme_dir=${SWCME_DIR:-"$amps_root/src/models/swcme"}

test -d "$sep_common_dir" || {
  echo "FAIL B02: canonical sep_common directory is absent: $sep_common_dir" >&2
  exit 2
}
sep_common_dir=$(CDPATH= cd -- "$sep_common_dir" && pwd)

# These names are the exact archive contract. coefficient_registry was born in
# the shared directory; the other six formerly had private srcSEP/util copies.
kernels="sep_transport_common sep_coefficient_physics sep_coefficient_registry sep_background_snapshot sep_test_registry sep_injection_spectrum sep_species_source"
for kernel in $kernels; do
  for suffix in h cpp; do
    test -f "$sep_common_dir/$kernel.$suffix" || {
      echo "FAIL B02: missing canonical source $kernel.$suffix" >&2
      exit 1
    }
    test ! -e "$src_root/util/$kernel.$suffix" || {
      echo "FAIL B02: private shared-kernel copy remains: util/$kernel.$suffix" >&2
      exit 1
    }
  done
done

# Direct path-qualified includes remain forbidden.  Application translation
# units use -I$(SEP_COMMON_DIR); AMPS-wide public headers use the single
# sep_common_header_path.h resolver tested below.  Keeping the fallback in one
# named resolver prevents individual headers from inventing relative paths or
# reaching a retired util copy.
if grep -R -n -E '#include[[:space:]]*[<"][^>"]+/sep_(transport_common|coefficient_physics|coefficient_registry|background_snapshot|test_registry|injection_spectrum|species_source)\.h[>"]' \
    "$src_root" --include='*.h' --include='*.hpp' --include='*.cpp' \
    >/dev/null; then
  echo "FAIL B02: a shared-kernel include still names an application-relative path" >&2
  exit 1
fi

# Production .cpp files are compiled after srcSEP has been copied to
# build/main.  They therefore have the same AMPS-root-only constraint as
# public headers included by sibling libraries.  Focused tests are excluded
# here deliberately: their short includes verify the documented detached
# component contract and their commands explicitly provide SEP_COMMON_DIR.
if grep -R -n -E '#include[[:space:]]*[<"]sep_(transport_common|coefficient_physics|coefficient_registry|background_snapshot|test_registry|injection_spectrum|species_source)\.h[>"]' \
    "$src_root" --include='*.cpp' --exclude-dir=test >/dev/null; then
  echo "FAIL B02: a production source bypasses sep_common_header_path.h" >&2
  exit 1
fi
grep -q '#include "util/sep_common_header_path.h"' "$src_root/sep.h" || {
  echo "FAIL B02: public sep.h does not load the canonical-header resolver" >&2
  exit 1
}
grep -q '#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)' \
  "$src_root/sep.h" || {
  echo "FAIL B02: public sep.h bypasses the canonical transport-header resolver" >&2
  exit 1
}
if grep -R -n 'src/models/sep_common/' "$src_root" \
    --include='*.h' --include='*.hpp' \
    --exclude='sep_common_header_path.h' >/dev/null; then
  echo "FAIL B02: canonical root fallback escaped sep_common_header_path.h" >&2
  exit 1
fi

for token in SEP_COMMON_DIR SEP_COMMON_ARCHIVE SEP_COMMON_OBJECTS \
             audit-production-shared-models; do
  grep -q "$token" "$src_root/makefile" || {
    echo "FAIL B02: srcSEP makefile omits $token" >&2
    exit 1
  }
done

# Build the sole archive owner, then require exact membership and unique
# strong symbols. Recreating the archive in its makefile prevents stale ar
# members from surviving after a source-list change.
make -C "$sep_common_dir" --no-print-directory verify
archive="$sep_common_dir/sep_common.a"
members=$(ar t "$archive" | sort)
expected=$(for kernel in $kernels; do printf '%s.o\n' "$kernel"; done | sort)
test "$members" = "$expected" || {
  echo "FAIL B02: sep_common.a has unexpected members" >&2
  printf 'archive:\n%s\nexpected:\n%s\n' "$members" "$expected" >&2
  exit 1
}
duplicates=$(nm -g --defined-only "$archive" | \
  awk '$2 ~ /^[TDBRS]$/ { count[$3]++ } END { for (s in count) if (count[s] > 1) print s }')
test -z "$duplicates" || {
  echo "FAIL B02: duplicate strong definitions in sep_common.a" >&2
  printf '%s\n' "$duplicates" >&2
  exit 1
}

fixture=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-common-layout.XXXXXX")
trap 'rm -rf "$fixture"' EXIT HUP INT TERM
mkdir -p "$fixture/srcSEP" "$fixture/build/main" \
         "$fixture/src/models/sep_common" "$fixture/src/models/swcme"
: > "$fixture/Makefile.conf"
cp "$src_root/makefile" "$fixture/srcSEP/makefile"
cp "$src_root/makefile" "$fixture/build/main/makefile"

expected_root=$(CDPATH= cd -- "$fixture" && pwd)
expected_common="$expected_root/src/models/sep_common"
expected_swcme="$expected_root/src/models/swcme"
for makefile in "$fixture/srcSEP/makefile" "$fixture/build/main/makefile"; do
  output=$(cd / && env -u AMPS_ROOT -u AMPS_CONFIG -u SEP_COMMON_DIR \
    -u SWCME_DIR -u MAKEFLAGS -u MAKEOVERRIDES -u MFLAGS -u MAKELEVEL \
    make --no-print-directory -f "$makefile" print-layout-paths)
  printf '%s\n' "$output" | grep -Fx "AMPS_ROOT=$expected_root" >/dev/null
  printf '%s\n' "$output" | grep -Fx \
    "AMPS_CONFIG=$expected_root/Makefile.conf" >/dev/null
  printf '%s\n' "$output" | grep -Fx \
    "SEP_COMMON_DIR=$expected_common" >/dev/null
  printf '%s\n' "$output" | grep -Fx \
    "MODEL_INCLUDE_FLAGS=-I$expected_common -I$expected_swcme" >/dev/null
done

# Link a real symbol through the archive.  This catches an archive that merely
# has correctly named members but cannot satisfy a canonical public consumer.
cat > "$fixture/consumer.cpp" <<'EOF'
#include "sep_transport_common.h"

int main() {
  const SEP::Transport::NumericalTolerances tolerances;
  return SEP::Transport::ValidateNumericalTolerances(tolerances).ok() ? 0 : 1;
}
EOF
${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$sep_common_dir" "$fixture/consumer.cpp" "$archive" \
  -o "$fixture/consumer"
"$fixture/consumer"

# Reproduce the failing AMPS sibling-submake topology.  meshAMR, PIC, and other
# libraries see -I<AMPS_ROOT> while processing build/main/sep.h, but variables
# appended by build/main/makefile cannot propagate back to those sibling
# submakes.  This consumer deliberately omits -I$(SEP_COMMON_DIR); successful
# compilation proves the public header resolves the same canonical source via
# the AMPS-root include already present in every production command.
cat > "$fixture/amps_root_consumer.cpp" <<'EOF'
#include "srcSEP/util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_species_source.h)
#include SRCSEP_SEP_COMMON_HEADER(sep_injection_spectrum.h)
#include "srcSEP/util/sep_transport_coefficients.h"

int main() {
  SEP::Transport::SpatialDiffusionSample sample;
  return sample.status.ok() ? 0 : 0;
}
EOF
${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$amps_root" -fsyntax-only "$fixture/amps_root_consumer.cpp"

# A detached component checkout may not have an AMPS-root include at all. Its
# documented SEP_COMMON_DIR must therefore remain a complete secondary path.
cat > "$fixture/detached_consumer.cpp" <<'EOF'
#include "sep_transport_coefficients.h"

int main() {
  SEP::Transport::MeanFreePathSample sample;
  return sample.status.ok() ? 0 : 0;
}
EOF
${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$src_root/util" -I"$sep_common_dir" -fsyntax-only \
  "$fixture/detached_consumer.cpp"

# The same resolver must work after AMPS copies the application to build/main.
# Only headers needed by this focused probe are staged; the canonical model is
# intentionally not copied, linked, or forwarded into the application tree.
mkdir -p "$fixture/build/main/util"
cp "$src_root/util/sep_common_header_path.h" \
   "$src_root/util/sep_transport_coefficients.h" \
   "$fixture/build/main/util/"
cat > "$fixture/copied_header_consumer.cpp" <<'EOF'
#include "sep_transport_coefficients.h"

int main() {
  SEP::Transport::PitchAngleDiffusionSample sample;
  return sample.status.ok() ? 0 : 0;
}
EOF
${CXX:-c++} -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$fixture/build/main/util" -I"$amps_root" -fsyntax-only \
  "$fixture/copied_header_consumer.cpp"

echo "B02 PASS: canonical archive and AMPS-root public-header/production-source fallback are relocation-safe"
