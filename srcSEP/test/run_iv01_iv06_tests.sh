#!/bin/sh
# Strict source gate and optional linked-application execution for IV01-IV06.
# The source-only phase checks compilation/orchestration but reports no physics
# PASS. Scientific acceptance requires SEP_EXECUTABLE and native registry rows.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
out=${SEP_TEST_OUTPUT_DIR:-/tmp/srcsep-iv01-iv06-test}
mkdir -p "$out"
# Keep interpreter bytecode beside the disposable test evidence.  This makes
# the gate safe to run immediately before packaging and prevents validation
# imports from contaminating the source tree with __pycache__ directories.
export PYTHONPYCACHEPREFIX="$out/pycache"
# The integrated native source reaches canonical shared headers through
# srcSEP's public resolver.  Supplying the repository root explicitly makes
# that path independent of both the process CWD and generated AMPS includes.
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -I"$root/.." \
  -c "$root/validation/cases/integrated_validation_models.cpp" \
  -o "$out/integrated_validation_models.o"
python3 -m py_compile \
  "$root/validation/cases/integrated_case_runner.py" \
  "$root"/validation/cases/IV0[1-6]/case.py \
  "$root"/validation/cases/IV0[1-6]/reference_solution.py
python3 -m json.tool "$root/validation/case_registry.json" >/dev/null
python3 "$root/validation/run_case.py" --list | grep -q IV06
grep -q 'RunIntegratedValidationModel' "$root/component_tests.cpp"
grep -q 'integrated_validation_models.o' "$root/makefile"
echo "PASS IV01-IV06-UNIT: integrated native source, Python, and registry compile"
if [ -z "${SEP_EXECUTABLE:-}" ]; then
  echo "SKIP IV01-IV06-LINKED: set SEP_EXECUTABLE to linked srcSEP/AMPS"
  exit 0
fi
python3 "$root/test/run_tests.py" --amps "$SEP_EXECUTABLE" \
  --validation-case IV01 --validation-case IV02 --validation-case IV03 \
  --validation-case IV04 --validation-case IV05 --validation-case IV06 \
  --output-dir "$out/linked"
for id in IV01 IV02 IV03 IV04 IV05 IV06; do
  test -s "$out/linked/$id/${id}_solution.csv"
  test -s "$out/linked/$id/${id}_comparison.png"
  test -s "$out/linked/$id/${id}_comparison.eps"
done
echo "PASS IV01-IV06-LINKED: native models, references, metrics, and figures"
