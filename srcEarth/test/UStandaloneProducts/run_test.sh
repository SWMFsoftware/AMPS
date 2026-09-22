#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/amps-standalone-products.XXXXXX")"
trap 'rm -rf "$build_dir"' EXIT

"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic \
  "$script_dir/test_standalone_products.cpp" -o "$build_dir/test_standalone_products"
"$build_dir/test_standalone_products"
