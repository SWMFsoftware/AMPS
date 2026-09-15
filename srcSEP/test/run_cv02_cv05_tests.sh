#!/bin/sh
set -eu

# Always verify the nested validation source against the production header
# layout and strict warnings.  Full scientific evidence is generated only by a
# supplied linked srcSEP/AMPS application; no standalone harness is promoted to
# application validation when that prerequisite is unavailable.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
output_dir=$(mktemp -d "${TMPDIR:-/tmp}/srcsep-cv02-cv05.XXXXXX")
trap 'rm -rf "$output_dir"' EXIT HUP INT TERM

"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic -O1 -g \
  -c "$src_root/validation/cases/controlled_transport_models.cpp" \
  -o "$output_dir/controlled_transport_models.o"
test -s "$output_dir/controlled_transport_models.o"

# Syntax-check every orchestration/reference module without leaving pycache in
# the source tree.  Registry loading below additionally validates all paths and
# versioned JSON descriptors.
PYTHONPYCACHEPREFIX="$output_dir/pycache" python3 -m py_compile \
  "$src_root/validation/cases/linked_case_common.py" \
  "$src_root/validation/cases/controlled_case_runner.py" \
  "$src_root/validation/cases/CV02/case.py" \
  "$src_root/validation/cases/CV02/reference_solution.py" \
  "$src_root/validation/cases/CV03/case.py" \
  "$src_root/validation/cases/CV03/reference_solution.py" \
  "$src_root/validation/cases/CV04/case.py" \
  "$src_root/validation/cases/CV04/reference_solution.py" \
  "$src_root/validation/cases/CV05/case.py" \
  "$src_root/validation/cases/CV05/reference_solution.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$src_root/validation/run_case.py" --list \
  | grep -q 'CV05'

amps_executable=${SEP_EXECUTABLE:-"$src_root/../amps"}
if [ ! -x "$amps_executable" ]; then
  echo "PASS CV02-CV05-UNIT: native source, Python modules, and registry compile"
  echo "SKIP CV02-CV05-LINKED: set SEP_EXECUTABLE to the linked srcSEP/AMPS application"
  exit 0
fi

PYTHONDONTWRITEBYTECODE=1 python3 "$src_root/test/run_tests.py" \
  --amps "$amps_executable" \
  --validation-case CV02 --validation-case CV03 \
  --validation-case CV04 --validation-case CV05 \
  --output-dir "$output_dir" --formats png,eps

for case_id in CV02 CV03 CV04 CV05; do
  grep -q "\"id\": \"$case_id\"" "$output_dir/srcsep-tests.json"
  for artifact in \
    "$output_dir/$case_id/${case_id}_solution.csv" \
    "$output_dir/$case_id/${case_id}_reference.csv" \
    "$output_dir/$case_id/${case_id}_comparison.png" \
    "$output_dir/$case_id/${case_id}_comparison.eps" \
    "$output_dir/$case_id/native/${case_id}_model.csv" \
    "$output_dir/$case_id/provenance.json" \
    "$output_dir/plots/${case_id}_comparison.png" \
    "$output_dir/plots/${case_id}_comparison.eps"; do
    test -s "$artifact" || {
      echo "FAIL $case_id-ARTIFACT: missing or empty $artifact" >&2
      exit 1
    }
  done
done
grep -q '"failed": 0' "$output_dir/srcsep-tests.json"
grep -q '"errors": 0' "$output_dir/srcsep-tests.json"
echo "PASS CV02-CV05-LINKED: application models, references, metrics, and figures"
