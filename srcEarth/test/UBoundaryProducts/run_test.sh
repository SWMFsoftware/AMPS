#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../.." && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/amps-boundary-products.XXXXXX")"
trap 'rm -rf "$build_dir"' EXIT

"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic \
  -I"$script_dir" \
  "$script_dir/test_boundary_products.cpp" -o "$build_dir/test_boundary_products"
"$build_dir/test_boundary_products"

# Compile and execute the production anisotropy adapter as a separate focused
# target.  This guards declaration/order errors in AnisotropicSpectrum.cpp that
# header-only BoundaryProducts tests cannot detect.
"${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -pedantic \
  -I"$source_root" \
  "$script_dir/test_anisotropic_spectrum.cpp" \
  "$source_root/gridless/AnisotropicSpectrum.cpp" \
  -o "$build_dir/test_anisotropic_spectrum"
"$build_dir/test_anisotropic_spectrum"
