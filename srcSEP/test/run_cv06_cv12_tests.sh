#!/bin/sh
# Source/registry gate for CV06-CV12, plus optional linked-application science.
# A missing application is an explicit SKIP: this script never substitutes the
# standalone model harness for the executable users intend to validate.
set -eu
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
output_dir=${SEP_TEST_OUTPUT_DIR:-/tmp/srcsep-cv06-cv12-test}
mkdir -p "$output_dir"

g++ -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -c "$src_root/validation/cases/advanced_validation_models.cpp" \
  -o "$output_dir/advanced_validation_models.o"
# py_compile intentionally writes bytecode even when PYTHONDONTWRITEBYTECODE is
# set. Redirect that cache into the external test-output directory so running
# this source gate cannot make a clean source/archive fail DOC03 afterward.
PYTHONPYCACHEPREFIX="$output_dir/pycache" python3 -m py_compile \
  "$src_root/validation/cases/advanced_case_runner.py" \
  "$src_root"/validation/cases/CV0[6-9]/case.py \
  "$src_root"/validation/cases/CV1[0-2]/case.py \
  "$src_root"/validation/cases/CV0[6-9]/reference_solution.py \
  "$src_root"/validation/cases/CV1[0-2]/reference_solution.py
python3 -m json.tool "$src_root/validation/case_registry.json" >/dev/null
python3 "$src_root/validation/run_case.py" --list | grep -q CV12
grep -q 'RunAdvancedValidationModel' "$src_root/component_tests.cpp"
grep -q '"CV12", "Controlled particle-wave' "$src_root/component_tests.cpp"
grep -q 'validation/cases/advanced_validation_models.o' "$src_root/makefile"
echo "PASS CV06-CV12-UNIT: native source, Python modules, and registry compile"

if [ -z "${SEP_EXECUTABLE:-}" ]; then
  echo "SKIP CV06-CV12-LINKED: set SEP_EXECUTABLE to the linked srcSEP/AMPS application"
  exit 0
fi

python3 "$src_root/test/run_tests.py" --amps "$SEP_EXECUTABLE" \
  --validation-case CV06 --validation-case CV07 --validation-case CV08 \
  --validation-case CV09 --validation-case CV10 --validation-case CV11 \
  --validation-case CV12 --output-dir "$output_dir/linked"
for case_id in CV06 CV07 CV08 CV09 CV10 CV11 CV12; do
  test -s "$output_dir/linked/$case_id/${case_id}_solution.csv"
  test -s "$output_dir/linked/$case_id/${case_id}_comparison.png"
  test -s "$output_dir/linked/$case_id/${case_id}_comparison.eps"
done
echo "PASS CV06-CV12-LINKED: linked models, references, metrics, and figures"
