#!/bin/sh
# Dependency-light source gate plus optional linked XM01-XM03 execution.
# A source compile verifies registration/build integrity. Scientific PASS is
# possible only through SEP_EXECUTABLE. XM02 and XM03 both generate their model
# results inside that application; neither accepts an external result table.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
out=${SEP_TEST_OUTPUT_DIR:-/tmp/srcsep-xm01-xm03-test}
mkdir -p "$out"
export PYTHONPYCACHEPREFIX="$out/pycache"

g++ -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -c "$root/validation/cases/cross_model_validation_models.cpp" \
  -o "$out/cross_model_validation_models.o"

# Build the dependency-light standalone form of the same native callback and
# execute XM02 with the registered parameters.  This catches runtime/schema
# regressions even on a workstation where the complete AMPS link is absent;
# scientific evidence still requires SEP_EXECUTABLE in the block below.
g++ -std=c++11 -Wall -Wextra -Wpedantic -Werror \
  -DSRCSEP_CROSS_MODEL_STANDALONE_TEST_HARNESS \
  "$root/validation/cases/cross_model_validation_models.cpp" \
  "$root/util/sep_focused_transport_core.cpp" \
  "$root/util/sep_parker_core.cpp" \
  "$root/util/sep_transport_common.cpp" \
  -o "$out/xm_native_harness"
"$out/xm_native_harness" XM02 "$out/XM02_native_model.csv" \
  --injection-radius-solar-radii 2.5 --observer-radius-au 1.0 \
  --particle-energy-mev 10.1 --mfp-0-au 0.05 --mfp-1-au 0.3 \
  --mfp-2-au 1.0 --plasma-advection-m-per-s 0.0 \
  --dlnb-ds-per-m 0.0 --source-rise-time-s 1800 \
  --source-decay-time-s 10800 --particles 6000 --time-step-s 120 \
  --duration-s 158400 --output-cadence-s 7200 --campaign-seed 30202
test "$(wc -l < "$out/XM02_native_model.csv")" -eq 67
grep -q '^1,mfp_0.05au_integral_gt10mev,' "$out/XM02_native_model.csv"
grep -q '^43,mfp_1.0au_integral_gt10mev,' "$out/XM02_native_model.csv"

# Use the minimum accepted ensemble for a fast XM03 callback/schema smoke test.
# Scientific acceptance is evaluated only by the linked run below with the
# registered 2500 particles per energy bin.
"$out/xm_native_harness" XM03 "$out/XM03_native_model.csv" \
  --source-history-csv "$root/validation/cases/XM03/input/earth_shock_thermal_source.csv" \
  --injection-radius-solar-radii 2.5 --observer-radius-au 1 \
  --solar-wind-speed-m-per-s 363000 --shock-speed-m-per-s 675000 \
  --solar-rotation-rate-rad-per-s 2.8653290846e-6 \
  --launch-offset-s 5040 --connection-delay-s 900 \
  --mfp-normalization-au 0.3 --mfp-radial-exponent 1 \
  --mfp-rigidity-exponent 0.3333333333333333 \
  --injection-min-energy-mev 0.01 --injection-max-energy-mev 200 \
  --injection-momentum-index 5 --injection-flux-factor 1.2 \
  --energy-bins 16 --particles-per-energy 100 --time-step-s 120 \
  --duration-s 158400 --snapshot-window-s 7200 --campaign-seed 30303
test "$(wc -l < "$out/XM03_native_model.csv")" -eq 49
grep -q '^4,' "$out/XM03_native_model.csv"
grep -q '^36,' "$out/XM03_native_model.csv"
python3 -m py_compile \
  "$root/validation/cases/cross_model_case_runner.py" \
  "$root/validation/cases/XM01/reference_solution.py" \
  "$root"/validation/cases/XM0[1-3]/case.py
python3 -m json.tool "$root/validation/case_registry.json" >/dev/null

# The literature reconstructions are scientific inputs to the comparison
# workflow even though they are intentionally not executable SWMF decks.  Parse
# both here so malformed hand edits fail before a costly linked application run.
python3 -m json.tool "$root/validation/cases/XM02/publication_input.json" >/dev/null
python3 -m json.tool "$root/validation/cases/XM03/publication_input.json" >/dev/null
python3 -m json.tool "$root/validation/cases/XM03/reference/provenance.json" >/dev/null
grep -q '"reproduction_status": "partial"' \
  "$root/validation/cases/XM02/publication_input.json"
grep -q '"reproduction_status": "partial"' \
  "$root/validation/cases/XM03/publication_input.json"
# The reviewed vector extraction is immutable test evidence.  The line counts
# include one header and catch truncated or accidentally raster-derived tables.
test "$(wc -l < "$root/validation/cases/XM03/reference/liu_figure12_earth_observations.csv")" -eq 81
test "$(wc -l < "$root/validation/cases/XM03/input/earth_shock_thermal_source.csv")" -eq 98
grep -q '"reference_kind": "digitized-earth-observations"' \
  "$root/validation/cases/XM03/input.json"
if grep -q 'model_source_csv' "$root/validation/cases/XM03/input.json"; then
  echo "ERROR XM03 must not accept an external model result" >&2
  exit 1
fi

# Exercise the publication renderer without requiring the complete AMPS link.
# The smoke-model CSV was produced by the native XM03 callback above, while the
# observation CSV is the immutable reviewed Figure-12 extraction.  Each time
# must produce its own PNG/EPS pair; a combined three-panel figure no longer
# satisfies the XM03 publication-output contract.
python3 - "$root" "$out" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve() / "XM03_plot_contract"
output.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(root))
sys.path.insert(0, str(root / "validation/cases"))
from validation.cases import cross_model_case_runner as runner

case = json.loads((root / "validation/cases/XM03/input.json").read_text(
    encoding="utf-8"))
model = runner.read_csv(Path(sys.argv[2]) / "XM03_native_model.csv")
reference = runner.read_csv(
    root / "validation/cases/XM03/reference/liu_figure12_earth_observations.csv")
_, scale, _ = runner._xm03_score(case, model, reference)
paths = runner._xm03_plot(
    case, output, model, reference, scale, case["plot"]["formats"])
if len(paths) != 6:
    raise SystemExit(f"XM03 renderer returned {len(paths)} artifacts, expected 6")
PY
for time in 04h 12h 36h; do
  test -s "$out/XM03_plot_contract/XM03_earth_observation_comparison_${time}.png"
  test -s "$out/XM03_plot_contract/XM03_earth_observation_comparison_${time}.eps"
done
python3 "$root/validation/run_case.py" --list | grep -q 'XM03 | cross-model'
grep -q 'RunCrossModelValidationModel' "$root/component_tests.cpp"
grep -q 'cross_model_validation_models.o' "$root/makefile"

# Run the independent solver in isolation to catch stability/schema failures;
# this is a software PASS only and is not labelled cross-model validation.
python3 "$root/validation/cases/XM01/reference_solution.py" \
  --input "$root/validation/cases/XM01/input.json" \
  --output "$out/XM01_reference.csv"
test -s "$out/XM01_reference.csv"
echo "PASS XM01-XM03-UNIT: native source, references, Python, and registry compile"

if [ -z "${SEP_EXECUTABLE:-}" ]; then
  echo "SKIP XM01-XM03-LINKED: set SEP_EXECUTABLE to linked srcSEP/AMPS"
  exit 0
fi
python3 "$root/test/run_tests.py" --amps "$SEP_EXECUTABLE" \
  --validation-case XM01 --validation-case XM02 --validation-case XM03 \
  --output-dir "$out/linked"
test -s "$out/linked/XM01/XM01_comparison.png"
test -s "$out/linked/XM01/XM01_comparison.eps"
test -s "$out/linked/XM02/XM02_comparison.png"
for time in 04h 12h 36h; do
  test -s "$out/linked/XM03/XM03_earth_observation_comparison_${time}.png"
  test -s "$out/linked/XM03/XM03_earth_observation_comparison_${time}.eps"
done
# The runner must preserve the exact literature-input interpretation next to
# each reference plot.  These artifacts make a reference-only SKIP auditable.
test -s "$out/linked/XM02/XM02_publication_input.json"
test -s "$out/linked/XM03/XM03_publication_input.json"
echo "PASS XM01-XM03-LINKED-RUNNER: all linked cross-model/observation cases completed"
