#!/bin/sh
# Dependency-light integrity gate plus optional linked OV01-OV05 campaign.
# This script never promotes the standalone harness to observational evidence:
# the harness catches callback/schema regressions, while only SEP_EXECUTABLE
# can create a scientific PASS/FAIL record for the linked application.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
out=${SEP_TEST_OUTPUT_DIR:-/tmp/srcsep-ov01-ov05-test}
mkdir -p "$out"
export PYTHONPYCACHEPREFIX="$out/pycache"

g++ -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -DSRCSEP_CROSS_MODEL_STANDALONE_TEST_HARNESS \
  "$root/validation/cases/cross_model_validation_models.cpp" \
  "$root/util/sep_focused_transport_core.cpp" \
  "$root/util/sep_parker_core.cpp" \
  "$root/util/sep_transport_common.cpp" \
  -o "$out/ov_native_harness"

# Exercise the generic observational callback with the minimum accepted
# ensemble.  The resulting 24 hourly bins at two energies plus one header
# verify the output schema without pretending to validate the event.
"$out/ov_native_harness" OV02 "$out/OV02_model.csv" \
  --observer PSP --output-mode time-profile --source-radius-mode inner-boundary \
  --energies-mev 2.2,12.3 \
  --source-history-csv "$root/validation/cases/OV02/input/coronal_shock_source.csv" \
  --injection-radius-solar-radii 2.5 --observer-radius-au 0.33 \
  --solar-wind-speed-m-per-s 337000 --shock-speed-m-per-s 337000 \
  --solar-rotation-rate-rad-per-s 2.8653290846e-6 \
  --launch-offset-s 0 --connection-delay-s 0 \
  --mfp-normalization-au 0.0465 --mfp-radial-exponent 1 \
  --mfp-rigidity-exponent 0.3333333333333333 \
  --injection-momentum-index 6.16 --particles-per-energy 100 \
  --time-step-s 120 --duration-s 86400 --output-cadence-s 3600 \
  --campaign-seed 40202
test "$(wc -l < "$out/OV02_model.csv")" -eq 49
grep -q '^PSP,PSP_2.2MeV,' "$out/OV02_model.csv"
grep -q '^PSP,PSP_12.3MeV,' "$out/OV02_model.csv"

python3 -m py_compile "$root/validation/cases/observational_case_runner.py" \
  "$root"/validation/cases/OV0[1-5]/case.py
python3 -m json.tool "$root/validation/case_registry.json" >/dev/null
for case_id in OV01 OV02 OV03 OV04 OV05; do
  python3 -m json.tool "$root/validation/cases/$case_id/input.json" >/dev/null
  python3 -m json.tool \
    "$root/validation/cases/$case_id/publication_input.json" >/dev/null
  grep -q '"publication_input_only"[[:space:]]*:[[:space:]]*true' \
    "$root/validation/cases/$case_id/input.json"
done
python3 "$root/validation/run_case.py" --list | grep -q 'OV05 | observational-validation'
echo "PASS OV01-OV05-UNIT: native callback, fixed inputs, references, and registry"

if [ -z "${SEP_EXECUTABLE:-}" ]; then
  echo "SKIP OV01-OV05-LINKED: set SEP_EXECUTABLE to linked srcSEP/AMPS"
  exit 0
fi
python3 "$root/test/run_tests.py" --amps "$SEP_EXECUTABLE" \
  --validation-case OV01 --validation-case OV02 --validation-case OV03 \
  --validation-case OV04 --validation-case OV05 --output-dir "$out/linked"
for case_id in OV01 OV02 OV03 OV04 OV05; do
  test -s "$out/linked/$case_id/provenance.json"
done
test -s "$out/linked/OV01/OV01_earth_observation_comparison.png"
test -s "$out/linked/OV05/OV05_observational_comparison.eps"
echo "PASS OV01-OV05-LINKED-RUNNER: linked observational campaign completed"
