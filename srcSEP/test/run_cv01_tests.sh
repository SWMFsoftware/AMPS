#!/bin/sh
set -eu

# This target always performs a source-level compile/contract check. When a
# linked srcSEP/AMPS executable is available it additionally runs the complete
# application-level CV01 campaign. A missing binary is reported as SKIP rather
# than substituting the former standalone model and mislabelling that result as
# validation of the linked application.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
output_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-cv01.XXXXXX")
trap 'rm -rf "$output_dir"' EXIT HUP INT TERM

# Deliberately omit a srcSEP/util include flag.  Add only the repository root,
# which is the canonical production include root for
# <src/models/sep_common/...>.  This makes the detached compile hermetic: it
# neither relies on the caller's current directory nor on stale headers in an
# AMPS build tree, while still exercising the same public path used by sibling
# AMPS libraries.
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic -O1 -g \
  -I"$src_root/.." \
  -c "$src_root/validation/cases/CV01/cv01_model.cpp" \
  -o "$output_dir/cv01_model.o"
test -s "$output_dir/cv01_model.o"

amps_executable=${SEP_EXECUTABLE:-"$src_root/../amps"}
if [ ! -x "$amps_executable" ]; then
  echo "PASS CV01-UNIT: linked model source compiles with strict warnings"
  echo "SKIP CV01-LINKED: set SEP_EXECUTABLE to the linked srcSEP/AMPS application"
  exit 0
fi

PYTHONDONTWRITEBYTECODE=1 python3 "$src_root/test/run_tests.py" \
  --amps "$amps_executable" --validation-case CV01 \
  --output-dir "$output_dir" --formats png,eps

report="$output_dir/srcsep-tests.json"
test -s "$report"
test -s "$output_dir/srcsep-tests.xml"
grep -q '"id": "CV01"' "$report"
grep -q '"status": "PASS"' "$report"
grep -q 'negative_control_error_over_line_length' "$report"
grep -q 'boundary_time_error_over_final_time' "$report"

for artifact in \
  "$output_dir/CV01/CV01_solution.csv" \
  "$output_dir/CV01/CV01_packet_moments.csv" \
  "$output_dir/CV01/CV01_final_profiles.csv" \
  "$output_dir/CV01/CV01_four_panel.png" \
  "$output_dir/CV01/CV01_four_panel.eps" \
  "$output_dir/plots/CV01_comparison.png" \
  "$output_dir/plots/CV01_comparison.eps" \
  "$output_dir/CV01/provenance.json" \
  "$output_dir/validation-run-manifest.json"; do
  test -s "$artifact" || {
    echo "FAIL CV01-ARTIFACT: missing or empty $artifact" >&2
    exit 1
  }
done

# These source assertions prevent a future refactor from replacing the linked
# model with the independent reference or from silently enabling scattering.
grep -q 'AdvanceFocusedTransportDmumu' \
  "$src_root/validation/cases/CV01/cv01_model.cpp"
grep -q 'dMuMuPerS = 0.0' \
  "$src_root/validation/cases/CV01/cv01_model.cpp"
grep -q 's(t) = s0' \
  "$src_root/validation/cases/CV01/reference_solution.py"
grep -q -- '--test CV01' "$output_dir/CV01/CV01_run.log"
echo "PASS CV01-LINKED: AMPS model/reference, boundaries, metrics, and figures"
