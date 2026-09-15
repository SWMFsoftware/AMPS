# EV01-EV02 change manifest

Implemented campaign-level evidence cases EV01 and EV02 on top of the supplied
OV01-OV05 codebase.

## New code and data

- `validation/cases/campaign_evidence_runner.py`
- `validation/cases/EV01/{case.py,input.json,README.md}`
- `validation/cases/EV02/{case.py,input.json,README.md}`
- `validation/cases/EV01/reference/ccmc_2021_goes_observations.csv`
- `validation/cases/EV02/reference/ccmc_2021_goes_observations.csv`
- reference provenance JSON for both cases
- `test/run_ev01_ev02_tests.sh`
- `EV01_EV02_VALIDATION_IMPLEMENTATION.md`

## Integrated files

- `validation/case_registry.json`
- `validation/run_case.py`
- `validation/cases/cross_model_validation_models.cpp`
- `component_tests.cpp`
- `test/run_tests.py`
- `makefile`
- relevant validation/test READMEs

## Evidence policy

The expected/reference quantities are direct NASA CCMC/GOES observed event
values.  No synthetic, theoretical, or assumed expected solution is generated.
Model-side source histories are explicitly identified as inputs rather than
reference data.

## Verification performed

- `./test/run_ev01_ev02_tests.sh` — PASS
- `python3 test/test_python_test_runner.py` — PASS (16 tests)
- strict C++11 compile of the native linked-validation callback with
  `-Wall -Wextra -Wpedantic -Werror` — PASS
- `./test/run_ov01_ov05_tests.sh` — unit checks PASS; linked execution skipped
  because no linked srcSEP/AMPS executable was included in the supplied archive.

A full EV01/EV02 end-to-end scientific run still requires the linked production
srcSEP/AMPS executable selected with `--amps`.

## Fix 2 — CCMC 2012-07-12 equal peak/end source timestamp

- Fixed `campaign_evidence_runner._write_source()` for the observed CCMC
  2012-07-12 record, whose flare peak and flare end are both 16:49 UTC.
  The previous triangular serialization emitted duplicate elapsed times and was
  rejected by the native source reader before particle transport began.
- Equal observed peak/end times are now encoded as a right-hand zero-rate sample
  at a portable IEEE-754 next-representable-value helper. This satisfies the native strict-monotonic
  source contract without inventing a finite physical decay timescale.
- Added validation of `onset <= peak <= end`, campaign-window coverage, detailed
  comments, and a regression check that generates all nine CCMC event sources
  and verifies strict monotonicity.
- Updated EV01/EV02 READMEs to document this observational edge case.


## Fix 3 — Python < 3.9 compatibility

- Removed the dependency on `math.nextafter`, which is unavailable in Python
  versions earlier than 3.9 and caused EV01/EV02 to stop during source-file
  generation on some NASA/HPC systems using `/usr/bin/python3`.
- Added `_next_float_toward_positive_infinity()`, a dependency-free IEEE-754
  binary64 helper implemented with Python's standard `struct` module. It has
  exactly the required next-representable-timestamp semantics and does not
  introduce an arbitrary physical source duration.
- Extended the EV unit gate to reject accidental reintroduction of
  `math.nextafter` and to verify the compatibility helper directly.
