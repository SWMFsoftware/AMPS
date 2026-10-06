#!/bin/sh
set -eu

# This is deliberately a provider/adapter test, not native AMPS qualification.
# It compiles the production srcSEP boundary against the maintained shared
# archives in a disposable directory and therefore detects API drift without
# creating another application binary in the source tree.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
src_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
repository_root=$(CDPATH= cd -- "$src_root/.." && pwd)
model_root="$repository_root/src/models"
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-reduced-shock.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

make -C "$model_root/sep_coronal_cme" -j16 lib
make -C "$model_root/sep_corona_swcme" -j16 lib

cxx=${CXX:-c++}
"$cxx" -std=c++17 -Wall -Wextra -Wpedantic -Werror \
  -I"$src_root/adapters" \
  -I"$model_root/sep_corona_swcme/include" \
  -I"$model_root/sep_corona_swcme/shock_front" \
  -I"$model_root/sep_coronal_cme/include" \
  -I"$model_root/sep_coronal_cme/src" \
  -I"$model_root/sep_common" \
  "$src_root/adapters/reduced_shock_background_adapter.cpp" \
  "$script_dir/test_reduced_shock_background_adapter.cpp" \
  -Wl,--start-group \
  "$model_root/sep_corona_swcme/build/libsep_corona_swcme.a" \
  "$model_root/sep_coronal_cme/build/libsep_coronal_cme.a" \
  -Wl,--end-group \
  -o "$build_dir/test_reduced_shock_background_adapter"

"$build_dir/test_reduced_shock_background_adapter" \
  "$src_root/examples/reduced-shock/positive_1au.event"
