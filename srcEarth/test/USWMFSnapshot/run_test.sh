#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build_dir="${TMPDIR:-/tmp}/earth_uswmfsnapshot_$$"
mkdir -p "$build_dir"
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

# This suite compiles the exact dependency-free production contract.  Warnings are
# errors so format/hash/status changes cannot hide behind compiler diagnostics.
${CXX:-c++} -std=c++11 -Wall -Wextra -Werror -pedantic \
  "$here/test_swmf_snapshot.cpp" \
  -o "$build_dir/test_swmf_snapshot"

"$build_dir/test_swmf_snapshot" \
  "$build_dir/roundtrip-a.csv" "$build_dir/roundtrip-b.csv"
PYTHONDONTWRITEBYTECODE=1 python3 "$here/test_step9_source_contract.py"
