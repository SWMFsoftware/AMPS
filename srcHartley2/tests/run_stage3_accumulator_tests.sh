#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
stage3_tmp=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage3.XXXXXX")
trap 'rm -rf -- "$stage3_tmp"' EXIT

${CXX:-g++} -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I"$source_root" \
  "$script_dir/test_stage3_residual_accumulator.cpp" \
  -o "$stage3_tmp/test_stage3_residual_accumulator"

"$stage3_tmp/test_stage3_residual_accumulator"
