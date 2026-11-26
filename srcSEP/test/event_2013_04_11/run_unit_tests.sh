#!/bin/sh
# Dependency-light unit tests for the production changes made for the
# 2013-04-11 real-event nightly test (see README.md in this directory):
#   SRC01-SRC03  log-uniform importance sampling of the DSA shock source
#   COEF07-COEF10 radial-rigidity power-law lambda and isotropic closure
#                (run through the existing Step 10 coefficient suite)
# Both compile the canonical sources directly with sanitizers; no AMPS/PIC
# build is needed.  A failing check exits non-zero.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
src_dir=$(CDPATH= cd -- "$here/../.." && pwd)
common="$src_dir/../src/models/sep_common"
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-event-2013.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$common" \
  "$common/sep_transport_common.cpp" \
  "$common/sep_injection_spectrum.cpp" \
  "$here/test_dsa_source_sampling.cpp" \
  -o "$build_dir/test_dsa_source_sampling"
ASAN_OPTIONS=detect_leaks=0 "$build_dir/test_dsa_source_sampling"

sh "$src_dir/test/run_step10_tests.sh"
