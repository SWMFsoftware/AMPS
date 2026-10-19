# Stage A6 validation: scope, tests, and acceptance thresholds

## Purpose and limits

Stage A6 verifies implementation consistency before controlled flight-data
experiments. It tests the shared selection logic, Python/C++ harmonic
convention, coefficient exchange format, fitter derivatives and recovery,
regularization limits, the C++ forward application, and Python/C++ pressure
reproduction.

The synthetic tests use the production pressure and Jacobian algebra with a
small deterministic transport operator. They are same-model recovery tests.
They do **not** independently validate coma physics, COPS calibration, the
SHAP5 mesh, the real descent geometry, or the uniqueness of a source inferred
from flight data. Passing A6 is therefore required before Stage A8, but it is
not evidence for a localized physical source.

## One-command gate

From `srcCG/LiouvilleCOPS`:

```bash
# Fast tests plus the C++ analytic/transport self-test. Solver integrations
# that need the configured external SHAP5/SPICE/observation data are skipped.
python3 tools/run_a6_validation.py --threads 2

# Complete gate, including isolated temporary-directory solver runs.
python3 tools/run_a6_validation.py --integration --threads 2 \
  --report a6_validation/validation_20261018.json
```

Relative report paths are placed under `<AMPS>/cops_output`; an existing
report is refused. The JSON record contains the exact commands, executable and
A6-test hashes, thread environment, return codes, elapsed times, test counts,
output hashes/tails, and a statement of what the gate does not claim. The
integration tests write only to temporary directories and remove their
products when finished.

The driver fixes `OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1`. Its
`--threads` value controls OpenMP for the C++ solver and self-test. This avoids
multiplying BLAS and OpenMP workers during validation.

## Requirement-to-test map

| A6 requirement | Automated evidence |
|---|---|
| 1. Zero-enhancement regression | `ProductIntegrationTests.test_zero_enhancement_forward_regression` runs the real SHAP5 solver with the complete all-zero degree-5 file. Numeric model rows, `elements.bin`, and every kernel must be byte-identical to the disabled baseline. Headers are intentionally different because provenance records the enabled component. |
| 2. Python/C++ basis equality through degree 8 | `tests/test_a2_activity_basis.py` compares every mode at poles, axes, near-poles, and generic coordinates. Maximum allowed absolute difference is `2e-14`. |
| 3. Ordering and file round trip | `CoefficientRoundTripTests` writes all 56 degree-5--8 coefficients with the production Python writer, verifies exact `(l,m)` ordering and binary64 values, and evaluates the file with the strict C++ reader. |
| 4. Selection counts and reasons | `tests/test_a1_selection_profiles.py` checks the three audited counts, closed endpoints, saturation, run/profile exclusions, missing data, inside-shape points, reason codes, and fatal count drift. |
| 5. Finite-difference Jacobians | A4 covers smooth and clipping-boundary derivatives. A6 checks all 24 degree-5--6 columns at a nonzero multi-mode state using central differences. |
| 6. One injected degree-5 mode | The overdetermined synthetic geometry must have rank 11/11, recover the complete coefficient vector, and reproduce the exact injected pressures. |
| 7. Several modes with controlled noise | Five injected degree-5--6 modes use independent deterministic multiplicative NG/RG noise with log standard deviation `2e-4`; recovery is judged in facet field space. |
| 8. Signal primarily on unseen facets | A north-pole-only geometry injects `(l,m)=(5,+1)`, which is exactly zero on every contributing facet. The test requires rank deficiency and non-recovery while demonstrating a nonzero field on unseen facets. |
| 9. Regularization limits | `lambda=0` recovers a supported exact signal; `lambda=1e10` must shrink its coefficient norm by at least `1e4`. The integration manifest must explicitly warn that `lambda=0` is unregularized. |
| 10. C++ verification of fitted pressures | The opt-in degree-5 Python fit is replayed by the real C++ solver and registered by `validate_cops_run.py`; measured NG and RG maximum relative differences must each be at most `1e-6`. |
| Legacy workflows remain green | The complete normal test discovery and `liouville_cops selftest.in` are part of the driver. With `--integration`, the existing collection, plotting, boundary, overwrite-protection, and small fit workflows also run. |

## Predeclared numerical thresholds

The following thresholds are test contracts, not values selected after seeing
flight-data results:

| Quantity | Required threshold |
|---|---|
| Python/C++ harmonic values | maximum absolute difference `<= 2e-14` |
| Smooth nonzero synthetic Jacobian | `rtol=4e-6`, `atol=3e-9` columnwise |
| Exact single-mode coefficients | `rtol=2e-8`, `atol=2e-10` |
| Exact single-mode weighted residual norm | `< 1e-10` |
| Noisy multi-mode facet field | relative L2 error `< 0.08`; correlation `> 0.997` |
| Unsupported north-cap signal | observed-facet field `< 1e-14`; Jacobian rank `< parameter count`; fitted norm `< 1e-12` |
| Strong-regularization limit | `||b(lambda=1e10)|| < 1e-4 ||b(lambda=0)||` |
| Python/C++ fitted pressures | maximum relative NG and RG difference `<= 1e-6` |

The noisy recovery threshold is intentionally applied to the activity field,
not coefficient-by-coefficient error: correlated harmonic coefficients can
represent nearly the same field over the observed geometry. The unseen test
is the converse safeguard—it must expose, rather than hide, an unsupported
coefficient direction.

## Measured result for the current implementation

On 2026-10-18, the complete `--integration --threads 2` gate passed 66/66
Python tests with no skips, followed by `SelfTest: ALL PASSED (0 failures)`.
The deterministic A6 synthetic suite measured:

- Python-writer/C++-reader residual difference: `1.943e-16` for 56 modes.
- Exact degree-5 recovery: rank `11/11`, maximum coefficient error
  `7.013e-10`, weighted residual norm `1.801e-11`.
- Nonzero degree-5--6 Jacobian: 24 columns, worst relative difference
  `1.087e-7`.
- Strong- to zero-lambda coefficient-norm ratio: `3.256e-12`.
- Noisy five-mode recovery: facet-field relative error `2.833e-2` and
  correlation `0.999598733`.
- Unsupported north-cap experiment: rank `1/11`, observed-facet maximum
  `0`, unseen-facet RMS `5.561e-2`, fitted coefficient norm `0`.
- Fitted Python/C++ pressure replay: maximum relative difference
  `5.360e-8` for NG and `5.368e-8` for RG.

These values are deterministic for the checked-in random seeds and current
implementation on the tested platform; optimizer/library versions can change
last digits while remaining within the predeclared thresholds. The
authoritative pass/fail evidence for another machine is the report produced
there by `run_a6_validation.py`.

## Interpretation of failures

- Basis, round-trip, zero-enhancement, or C++ pressure-replay failures are
  implementation blockers. Do not proceed to flight-data fitting.
- Exact or noisy recovery failure with a full-rank Jacobian is a fitter or
  derivative blocker.
- Rank deficiency in the controlled north-cap case is expected and required.
  Recovery of its injected non-zonal coefficient would indicate leakage or a
  faulty test geometry.
- Rank deficiency in a real flight-data fit is not automatically a software
  defect. It is an identifiability result and must be reported with singular
  values, visibility, selection profile, degree, and lambda.
- A successful optimizer exit is not a substitute for the independent C++
  replay threshold, and a successful replay is not a statement of physical
  adequacy or publication readiness.
