# Registered validation-case layout

The end-to-end validation portfolio uses one registry and one directory
contract. `validation/case_registry.json` is the authoritative catalog;
`validation/run_case.py` resolves selections and aggregates standard result
records; `test/run_tests.py` is the user-facing orchestrator that also renders
saved numerical/reference comparisons.

Each new test from the validation plan should use this layout:

```text
validation/cases/<ID>/
  README.md              physics, assumptions, commands, and failure diagnosis
  input.json             reviewed production-style configuration
  case.py                model/reference/scoring adapter with run_case()
  <model adapter>        calls a production kernel or linked executable
  <reference solver>     mathematically independent expected solution
  reference/*.csv        immutable observation/digitized publication points
  reference/provenance.json publication, figure, extraction, hash, uncertainty
  publication_input.json reported inputs, reduced assumptions, missing inputs
```

The case descriptor must provide a stable ID, group, runtime class, entrypoint,
default input, and description. A case entrypoint receives `source_root`,
`input_path`, `output_dir`, `executable`, and `timeout`; it returns the existing
`srcsep-component-tests-v1` result fields. `executable` is the resolved linked
srcSEP/AMPS application, never a case-specific substitute. The entrypoint must
write comparisons from saved model/reference files, record effective physics
flags and checksums, and return `PASS`, `FAIL`, `SKIP`, or `ERROR` without
terminating the parent runner.

This structure intentionally mirrors production campaign setup: a versioned
input is resolved once; the selected physics executable consumes it; raw model
output is immutable; a separate reference/observation stage produces matched
quantities; scoring and plots read saved artifacts; and a provenance manifest
links inputs, code, commands, and results. Later analytical, cross-model, and
observational cases may use different executables or data acquisition, but they
must preserve these lifecycle and evidence contracts.

OV01-OV05 use the same lifecycle with the `observational-validation` group. EV01-EV02 add the `campaign-evidence` group and use real NASA CCMC/GOES observations with a frozen train/validation/holdout split.
OV01/OV02 are release-gating comparisons. OV03-OV05 are diagnostic-only because
their compound or wide-longitude structure is outside a one-field-line model;
diagnostic metrics are retained with `gating=false`, never discarded or
misrepresented as acceptance. Evidence coverage and linked execution remain
hard gates in every case.

To add a case:

1. Copy the directory structure, choose the next stable plan ID, and add one
   descriptor to `case_registry.json`.
2. Express units and every enabled/disabled operator in `input.json`; never
   inherit an undocumented production default.
3. Register the numerical stage in the linked application's native test
   registry and invoke that exact executable through `--test <ID>`.
4. Implement the reference independently and add a negative control capable of
   detecting the principal sign, factor, or normalization failure.
5. Emit metrics with frozen comparisons/tolerances, CSV comparison series,
   PNG/EPS figures, a run log, resolved input, and provenance.
6. Add a bounded `make test-<id>-unit` gate and document which SWMF or
   observational evidence remains outside the linked-application result.

The production CLI carries case-specific paths with `--test-input` and
`--test-output-dir`. These options must be supplied together and only with one
explicit `--test` selector, preventing a reviewed input from leaking into a
group or `--all-tests` run. A generated narrow native protocol is allowed, but
the reviewed public input and resolved snapshot remain portable, unit-bearing
JSON.

List or execute registered cases with:

```sh
python3 validation/run_case.py --list
python3 test/run_tests.py --amps /path/to/amps --validation-case CV01 \
  --output-dir test_output/CV01
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case CV02 --validation-case CV03 \
  --validation-case CV04 --validation-case CV05 \
  --output-dir test_output/CV02-CV05
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case CV06 --validation-case CV07 \
  --validation-case CV08 --validation-case CV09 \
  --validation-case CV10 --validation-case CV11 \
  --validation-case CV12 --output-dir test_output/CV06-CV12
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case IV01 --validation-case IV02 \
  --validation-case IV03 --validation-case IV04 \
  --validation-case IV05 --validation-case IV06 \
  --output-dir test_output/IV01-IV06
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case XM01 --validation-case XM02 --validation-case XM03 \
  --output-dir test_output/XM01-XM03
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case EV01 --validation-case EV02 --output-dir test_output/EV01-EV02
python3 test/run_tests.py --amps /path/to/amps --validation-all \
  --output-dir test_output/validation-all
```

CV02-CV12 share `linked_case_common.py` for strict native process/report
verification. CV02-CV05 use `controlled_case_runner.py`; CV06-CV12 use
`advanced_case_runner.py`; IV01-IV06 use `integrated_case_runner.py`; XM01-XM03
use `cross_model_case_runner.py`. Their
expected physics remains in separate case-local `reference_solution.py`
programs. The shared C++ adapter calls the production Parker and focused-
transport/turbulence cores and is compiled into the requested AMPS application;
it is not a substitute executable. CV09 explicitly remains a controlled
shock-cycle test rather than a full resolved heliospheric shock campaign.

XM02/XM03 each use one registry-owned paper reconstruction and reject input
overrides. XM02 runs a controlled production-core first-passage ensemble for
the three published MFPs. XM03 runs production Parker transport on a reduced
Earth-connected Parker spiral, uses the Earth shock thermal-energy trace from
Liu et al. Figure 12(d) as source timing, and scores Earth observations from
Figure 12(a–c). Both obtain every model row from the selected linked executable
and require no external model CSV. See each XM README for the precise physics
scope, assumptions, columns, provenance, and provisional gates.
