#!/usr/bin/env bash
# Compile and execute the AMPS/MPI-independent Stage-5 correction-window
# planner test.  The temporary executable is deleted on every exit, so this
# runner cannot leave a binary or generated source in the repository.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)
stage5_tmp=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-stage5-window.XXXXXX")
trap 'rm -rf -- "$stage5_tmp"' EXIT

# Validate the generic input selector and its generated-policy contract before
# compiling the numerical planner.  Keeping both halves in one Stage-5 runner
# prevents a parser-valid window from drifting away from the C++ invariants.
perl "$source_root/input/test_comet_correction_window_policy.pl"

${CXX:-c++} -std=c++11 -O2 -Wall -Wextra -Werror -pedantic \
  -I"$source_root" \
  "$script_dir/test_stage5_correction_window_planner.cpp" \
  -o "$stage5_tmp/test_stage5_correction_window_planner"

"$stage5_tmp/test_stage5_correction_window_planner"
