# LiouvilleCOPS localized-activity amendment: development roadmap

Prepared 2026-10-18 for the standalone `srcCG/LiouvilleCOPS` application.

This document is a development plan, not a description of features that are
already available. Every stage below is **planned** until its implementation,
tests, and acceptance evidence are recorded. The existing degree-0--4 source
model, fitting results, and data products must remain usable while this
amendment is developed.

## 1. Purpose

The amendment will test whether spatial structure finer than the current
degree-4 activity map can explain localized RG/NG pressure features during the
last ten hours of the Rosetta mission. It will do so without changing the
existing heliocentric source-rate law, surface-temperature law, CO2 activity
map, spacecraft geometry, or collisionless Liouville transport model.

The central design is to retain the current H2O activity map as a fixed
large-scale baseline and fit a separate, explicitly identified higher-degree
H2O residual. This distinction is essential: the new coefficients must not be
presented as an extension of the srcCG coefficient table because the existing
25 functions use a legacy, unnormalized convention that is defined only
through degree 4 in this application.

The development must answer three separate questions:

1. Can a higher-resolution surface-activity field reproduce the pressure
   observations better than the legacy degree-4 field?
2. Is that improvement stable under data-selection, regularization, and
   temperature/CO2 choices?
3. Is the inferred enhancement constrained by the measurements, rather than
   being an optimizer artifact or unconstrained structure on unseen facets?

A lower residual alone is not sufficient evidence of a localized physical
source.

## 2. Scope and non-goals

In scope:

- Three named, reproducible pressure-data selection profiles.
- An optional H2O activity residual containing degrees 5 through `L`, with an
  initial supported maximum of `L = 8`.
- Independent control of activity resolution, analytical-temperature fitting,
  CO2 scaling, and regularization.
- C++ forward-model support for applying a fitted residual.
- Python fitter support for estimating the residual.
- Provenance, diagnostics, regression tests, synthetic recovery tests, and
  end-to-end plotting.

Not in scope for this amendment unless separately approved:

- Changing the existing H2O/CO2 heliocentric source-rate laws.
- Adding temperature spherical harmonics.
- Fitting CO2 or CO high-degree activity residuals.
- Changing the analytical temperature formula or its solar-zenith-angle
  dependence.
- Adding collisional transport, thermal lag, or molecular-flight nucleus
  rotation.
- Replacing the SHAP5 geometry, SPICE kernels, gauge response, or visibility
  calculation.
- Claiming that angular degree maps directly to a unique physical length on a
  bilobed, irregular nucleus.

These non-goals keep the experiment identifiable: the first implementation
changes only the spatial resolution of the H2O activity multiplier.

## 3. Frozen model components and new free parameters

### 3.1 Components retained exactly

For the first controlled experiments, retain:

- The current 25-function degree-0--4 H2O activity map as the baseline.
- The current H2O heliocentric production-rate law selected by the base run.
- The current H2O solar-zenith-angle dependence.
- The current analytical temperature formula and its base-run values.
- The current CO2 and CO activity coefficients.
- The base-run CO2 scale `kappa = 1`.
- The existing mesh, facet centroids/normals/areas, shadowing, visibility,
  spacecraft attitude/position, gauge response, and cached transport kernels.
- The existing handling of nucleus orientation at each observation time.

The initial high-degree comparison must therefore use the fitter equivalents
of `--fix-temperature` and `--fix-kappa`. Later stages may release temperature
and `kappa`, but only as a separate identifiability experiment.

### 3.2 Activity parameterization

Let `A0(x)` be the *linear*, unclipped H2O field produced by all 25 legacy
coefficients at a facet centroid `x`. Define a separate high-degree residual

```text
DeltaA(x) = A_scale * sum[l=5..L] sum[m=-l..l] b[l,m] Y[l,m](x)
```

and the applied H2O activity factor

```text
A(x) = max(0, A0(x) + DeltaA(x)).
```

Important requirements:

- The legacy coefficients stay fixed when the residual-only mode is selected.
- The residual begins at degree 5 so it cannot directly refit the same basis
  functions as the baseline.
- `A_scale` is a recorded numerical scaling constant, not a fitted physical
  parameter. Its definition must make coefficient magnitudes well conditioned
  and must be identical in Python and C++.
- Clipping occurs once, after adding the baseline and residual.
- All regularization and diagnostics must state whether they act on the linear
  field, the clipped field, or the coefficients.
- A zero-length/zero-valued residual must reproduce the present forward model
  to the documented numerical tolerance.

The number of new coefficients is:

| Maximum residual degree | New coefficients, degrees 5 through `L` |
|---:|---:|
| 4 | 0 |
| 5 | 11 |
| 6 | 24 |
| 7 | 39 |
| 8 | 56 |

Degree 4 is the explicit no-enhancement control.

### 3.3 New harmonic convention

Before coding the evaluator, write a normative basis specification covering:

- Body-fixed coordinate frame used by the SHAP5 facet centroids.
- `theta = acos(z/r)` and `phi = atan2(y,x)` definitions and pole handling.
- Real spherical-harmonic definition for `m < 0`, `m = 0`, and `m > 0`.
- Associated-Legendre recurrence, normalization, and phase convention.
- Coefficient order, fixed as increasing `l`, then increasing `m` from `-l`
  through `+l`.
- Floating-point type and tolerances.

The same explicit recurrence must be implemented in C++ and Python. Do not
rely on library defaults whose normalization or Condon--Shortley phase may
differ. The convention must have its own versioned identifier in all outputs.

The new functions are intentionally independent of the 25 legacy functions.
They need not numerically reproduce the legacy low-degree functions because
no new coefficient is assigned below degree 5.

## 4. Required data-selection profiles

Implement the profiles in a machine-readable configuration directory and run
them through one common selection function. The profile names and intended
contents are:

| Profile | Intended selection | Audited expected count |
|---|---|---:|
| `legacy_base_mask` | Honor the base run's `ExcludeNG`; keep all otherwise-valid RG | 839 NG, 1,146 RG |
| `all_available` | Use every positive, outside-shape NG and RG observation, including saturated and post-transition NG | 930 NG, 1,146 RG |
| `transition_cut_100m` | Exclude NG from `2016-09-30T10:24:42` through `2016-09-30T10:37:10`, inclusive; resume at `2016-09-30T10:37:14` (audited altitude 99.48 m); keep all otherwise-valid RG | 869 NG, 1,146 RG |

The counts above are acceptance checks for the currently audited base dataset,
not universal constants. If regenerated inputs change them, the runner must
fail with an explanatory selection diff until the change is reviewed and the
profile metadata is deliberately updated.

Rules common to all profiles:

- Nonpositive/missing pressures remain unavailable, not zero observations.
- The four final SHAP5-inside points remain excluded from fitting.
- Every otherwise-valid RG observation is selected.
- The 13 audited saturated NG readings are retained unless the profile
  explicitly excludes their timestamps. Their inclusion must be recorded.
- UTC interval endpoint semantics must be explicit and unit tested.
- The selected point IDs, timestamps, pressure values, altitude, reason codes,
  and a canonical selection hash must be saved before optimization.

Do not infer a physical interpretation of the NG jump from profile membership.
The profiles are sensitivity tests for a known instrument/data interval.

## 5. Proposed user interface

Preserve the present `tools/fit_cops.py --degree 0..4` behavior. Add options
whose final spelling must be confirmed by tests and documentation:

```text
--profile NAME
--enhancement-degree {4,5,6,7,8}
--fix-legacy-activity
--enhancement-scale VALUE
--enhancement-prior FILE
```

Semantics:

- `--degree` continues to control which legacy degree-0--4 coefficients are
  free in the old fitting mode.
- `--enhancement-degree` controls only the new H2O residual.
- `--fix-legacy-activity` fixes all 25 legacy H2O coefficients; it is required
  for the first controlled residual experiments.
- Existing `--fix-temperature`, `--ordered-temperature`, and `--fix-kappa`
  remain independent controls.
- `--lambda` controls the declared residual regularization. `lambda = 0` is
  allowed as an explicit diagnostic but must produce an unregularized-fit
  warning in the manifest and report.

Add a runner such as `tools/run_fit_profile.py` to resolve a named profile,
validate prerequisites, create a new output directory, assemble the complete
fitter invocation, and save the effective command/configuration. It must never
overwrite an existing output directory.

Planned command shapes, **not valid until their implementation stage is
complete**, are:

```bash
python3 tools/run_fit_profile.py \
  --run ../../cops_output/output_fit_base \
  --profile legacy_base_mask \
  --enhancement-degree 8 \
  --fix-legacy-activity --fix-temperature --fix-kappa \
  --lambda 1e-6 3.162e-6 1e-5 3.162e-5 1e-4 3.162e-4 1e-3 \
  --jobs 4 \
  --out ../../cops_output/fit_enh_l8_legacy
```

```bash
python3 tools/run_fit_profile.py \
  --run ../../cops_output/output_fit_base \
  --profile transition_cut_100m \
  --enhancement-degree 8 \
  --fix-legacy-activity --fix-temperature --fix-kappa \
  --lambda 1e-4 \
  --out ../../cops_output/fit_enh_l8_transition
```

The actual runner must obtain `<AMPS>/cops_output` through the existing path
resolution (`cops_paths.py`/`COPS_OUTPUT_ROOT`), not by assuming that `../..`
or a particular home directory is correct.

## 6. File and cache design

### 6.1 Preserve existing formats

Do not expand `elements.bin` from 25 stored basis functions to the new maximum.
The fitter can construct the 56 degree-5--8 values at startup from the stored
facet centroids. Approximately `100,000 * 56 * 8 = 44.8 MB` is required for a
dense double-precision table, which is acceptable only after measurement on
the actual mesh.

Do not expand high-degree basis values to every molecule/facet contribution in
the transport caches. Continue to combine the point--facet sparse structure
with facet activity. This avoids multiplying the high-degree table by the much
larger contribution count.

### 6.2 New coefficient file

Keep the existing `ActivityFit.txt` 3-by-25 format unchanged. Write the new
residual to a separate, versioned file, provisionally
`ActivityEnhancementFit.txt`, containing at least:

- Format and harmonic-convention versions.
- Body-fixed frame and angle definitions.
- Minimum and maximum degree.
- `A_scale` and units/dimensional interpretation.
- One `(l, m, coefficient)` record per expected coefficient.
- Hash/identity of the fixed legacy activity file or table.
- Regularization definition and selected `lambda`.

Reject duplicates, missing modes, unexpected modes, nonfinite values, wrong
degree ranges, and incompatible convention versions.

### 6.3 C++ input options

Add an optional component rather than changing `ActivityModel=file` semantics.
Provisional input keys are:

```text
ActivityEnhancementModel = none | real_spherical_harmonics_v1
ActivityEnhancementDegree = 4..8
ActivityEnhancementFile = <path>
ActivityEnhancementScale = <value>
```

`none` must be the default and reproduce current behavior. If a file supplies
the scale and degree, duplicated input values must agree or the run must fail.

### 6.4 Provenance

Extend run and fit manifests with:

- Data-profile name, profile-file hash, selected point hash, and counts by
  gauge/reason.
- Legacy coefficient source and hash.
- Enhancement basis ID, exact degree range, ordering, scale, coefficient-file
  hash, and clipping rule.
- Fixed/free status of legacy activity, enhancement, temperature, and `kappa`.
- Regularization formula and `lambda` values.
- Optimizer termination, evaluations, bounds, Jacobian diagnostics, and rank.
- C++ verification status and prediction differences.

Any plot or table used for comparison must carry enough identity to trace it
to this manifest.

## 7. Regularization and identifiability

The first implementation should extend the existing area-weighted map-prior
idea: regularize the residual toward zero, which is equivalent to regularizing
the total linear map toward the fixed legacy map. Define the exact residual,
normalization, facet weighting, and treatment of clipped facets in equations
and tests before selecting a `lambda`.

Requirements:

- Use positive `lambda` scans for scientific comparisons.
- Keep `lambda = 0` only as a clearly labeled conditioning/overfit diagnostic.
- Report data misfit and penalty separately.
- Report singular values/effective rank for the data Jacobian and the augmented
  problem.
- Report coefficient uncertainty only when the adopted likelihood and rank
  treatment justify it.
- Report area-weighted residual power by degree and activity on facets visible
  from selected measurements versus unseen facets.
- Do not add degree-weighted roughness by assumption. It is a later candidate
  only if zero-residual map regularization proves inadequate and its formula is
  justified and tested.

The final `lambda` must not be chosen solely because it gives the smallest
training residual. Use an L-curve or another documented selection rule and
compare neighboring values for stability.

## 8. Output and diagnostic requirements

Each fit directory must contain, at minimum:

- Immutable effective profile/configuration and complete invocation.
- Selected-data table with inclusion/exclusion reason codes.
- Fit manifest and input identities.
- Legacy coefficient file/reference and enhancement coefficient file.
- Optimizer report and objective decomposition.
- Predicted NG/RG pressures and normalized/log residuals.
- A per-facet table containing centroid, area, visibility/contribution flag,
  legacy linear activity, enhancement, total linear activity, and clipped total
  activity.
- Degree-power and coefficient diagnostics.
- Pressure comparison plots, residual plots, altitude/time views, and activity
  maps that separate baseline, residual, and total.
- C++ verification results or an explicit `not verified` state.

Plotting must visually distinguish selected data, excluded data, saturated NG,
model predictions, and the instrument-transition interval. A plot may show
excluded observations for context but must not imply they contributed to the
objective.

## 9. Development stages

### Stage A0 — Baseline protection and design freeze

Goal: establish an immutable comparison before changing the parameterization.

Tasks:

1. Read `README.md`, `RUN.md`, `quick_guide.md`, `ROADMAP.md`, and current
   manifests.
2. Record Git status and preserve all unrelated user changes.
3. Regenerate or identify a valid `output_fit_base` with kernels and manifest.
4. Run the current data-only selections and a small current fitter case into
   new output directories.
5. Record pressure predictions, selected IDs/counts, objective components,
   parameter ordering, and hashes needed for zero-enhancement regression.
6. Freeze the new harmonic equations, ordering, scaling, clipping semantics,
   profile definitions, and provisional file formats in a design note or the
   code documentation.

Acceptance gate:

- The baseline is repeatable and protected from overwrite.
- The no-enhancement target and numerical comparison tolerance are recorded.
- No model or scientific default has changed.

### Stage A1 — Shared profile selection and profile runner

Implementation status (2026-10-18): **implemented and fast-tested**. The
checked-in implementation is `tools/cops_selection.py`,
`tools/run_fit_profile.py`, `fit_profiles/*.json`, and
`tests/test_a1_selection_profiles.py`. The audited base produced the required
839/1,146, 930/1,146, and 869/1,146 counts. Full-fit execution status belongs
in the generated fit manifests and must not be inferred from this roadmap.

Goal: make the three requested data choices runnable before changing physics.

Tasks:

1. Add versioned configuration files for `legacy_base_mask`, `all_available`,
   and `transition_cut_100m`.
2. Refactor selection into one importable function used by the fitter and
   runner; avoid duplicate timestamp logic.
3. Emit reason codes and canonical selection hashes.
4. Implement runner validation, dry-run/data-only behavior, unique output
   directories, and complete command recording.
5. Keep direct `fit_cops.py` invocation backward compatible.

Acceptance gate:

- Current audited data yield 839/1146, 930/1146, and 869/1146 NG/RG counts.
- Endpoint, saturation, missing-value, inside-shape, and changed-input tests
  pass.
- A profile changes selection only, not model predictions at common points.

### Stage A2 — Harmonic basis specification and dual implementation

Implementation status (2026-10-18): **implemented and validated in
isolation**. The normative definition is `HARMONIC_BASIS.md`; implementations
are `tools/activity_basis.py` and
`src/ActivityEnhancementBasis.{h,cpp}`; cross-language tests are in
`tests/test_a2_activity_basis.py`. The deterministic maximum Python/C++
absolute difference was `1.33e-15` (required <= `2e-14`) and the numerical
orthonormality error was `1.47e-14`. On the audited 100,000-face cache, a
degree-5--8 table used exactly 44,800,000 bytes and five constructions had
median time `0.0987 s`; measured maximum-RSS increase was 115,688 KiB. The
cache was not modified. Forward integration is implemented in Stage A3 and
fitter integration in Stage A4.

Goal: implement an unambiguous degree-5--8 basis independently in Python and
C++.

Tasks:

1. Add the normative mathematical specification.
2. Implement stable associated-Legendre recurrence and real harmonics in a
   small Python module.
3. Implement the identical convention in a C++ component, for example
   `src/ActivityEnhancement.{h,cpp}`.
4. Add deterministic reference coordinates including poles, axes, generic
   points, and all `(l,m)` modes through degree 8.
5. Benchmark table construction time and memory on the actual mesh.

Acceptance gate:

- Python and C++ agree for all reference points/modes within a recorded
  tolerance no looser than required by the downstream prediction test.
- Ordering, normalization, phase, pole behavior, and nonfinite inputs are
  covered by tests.
- Measured memory/time are recorded; no transport cache has been inflated by
  the number of new modes.

### Stage A3 — Optional C++ forward-model enhancement

Implementation status (2026-10-18): **implemented and validated**. The strict
file contract is `ACTIVITY_ENHANCEMENT_FILE.md`; parsing/application are in
`src/ActivityEnhancement.{h,cpp}` and `src/Source.cpp`; fast failure and
mathematical tests are in `tests/test_a3_activity_enhancement.py`. The optional
component is H2O-only, checks body frame and exact legacy-map identity, adds
the residual to the unclipped legacy field, then clips once. Enabled runs
write `ActivityEnhancementFacets.dat` and complete manifest metadata; disabled
runs do not change the table/cache schemas. An actual one-point SHAP5 run with
the checked-in all-zero degree-5 file reproduced every baseline table field
exactly, including the previously manifested Stage-A0 point; the diagnostic
contained all 100,000 facets. A single-mode test agreed with the normative
Python basis, and malformed/incomplete/incompatible files failed before
propagation. Stage A4 now generates this strict input from fitted residuals.

Goal: allow the authoritative C++ forward model to apply a fitted residual.

Tasks:

1. Add strict enhancement-file parsing and new optional input keys.
2. Compute baseline linear activity, residual, and one final clipping step.
3. Keep H2O-only scope explicit; reject attempts to assign enhancement modes
   to CO2 or CO.
4. Print a concise model summary and write all new settings to provenance.
5. Add per-facet diagnostic output without changing old output schemas when
   enhancement is disabled.

Acceptance gate:

- `ActivityEnhancementModel=none` reproduces the Stage A0 predictions within
  the frozen tolerance.
- A zero coefficient file also reproduces the baseline.
- A single nonzero mode produces the analytically expected facet field.
- Malformed/incomplete/incompatible files fail before propagation.

### Stage A4 — Python fitter residual parameters

Implementation status (2026-10-18): **implemented and validated at the unit
and mission-Jacobian level**. `tools/fit_cops.py` and
`tools/run_fit_profile.py` implement the confirmed residual-only contract;
`tests/test_a4_enhancement_fit.py` covers degrees 4--8, exact parameter
counts, derived scale, the unclipped area-weighted regularizer, analytic
derivatives, and output identity. The existing A1--A3 tests remain passing.
On the audited 100,000-facet/1,150-point base, a degree-5
`all_available` check loaded 31,672,361 contributions, retained the expected
930 NG/1,146 RG points, reproduced the base C++ pressures to `4.65e-11` NG
and `9.79e-10` RG, and compared all 11 analytic columns with finite
differences. The diagnostic reduced each step until neither perturbation
crossed the clip; the largest reported relative difference was `1.13e-6`.
Smooth, strictly positive synthetic fields agree to `rtol=3e-6`. If a clip
kink cannot be avoided after step reduction, the diagnostic explicitly uses
the one-sided direction with fewer state changes; the implemented subgradient
is zero where total linear activity is at or below zero. A complete degree-4
zero-parameter control also passed the generated C++ forward rerun and
Python/C++ validator over all 1,146 outside points (maximum relative pressure
differences `4.716e-7` NG and `4.923e-7` RG with the then-current 15-digit
legacy coefficient serialization; serialization is now 17-digit for exact
binary64 round-trip). Full synthetic-recovery and fitted nonzero Python/C++
pressure studies remain Stage A6, not an A4 claim.

Goal: estimate the new coefficients efficiently with analytic derivatives.

Tasks:

1. Construct/cached facet basis values once at startup.
2. Add residual coefficients to the parameter vector only when requested.
3. Preserve the 25 legacy coefficients as exact constants in
   `--fix-legacy-activity` mode.
4. Extend activity clipping, sparse point--facet accumulation, predictions,
   and analytic Jacobian to the residual parameters.
5. Implement residual-to-zero regularization and its Jacobian.
6. Write enhancement and facet diagnostic files without altering the legacy
   3-by-25 `ActivityFit.txt` format.
7. Add bounds/scaling deliberately; do not invent scientifically meaningful
   coefficient bounds merely to make the optimizer converge.

Acceptance gate:

- Analytic derivatives pass finite-difference checks both away from clipping
  and on documented one-sided cases near the clipping boundary.
- Degree 4 with fixed legacy activity contains zero new parameters and matches
  the baseline.
- Degrees 5--8 produce exactly the documented parameter counts.
- Existing degree-0--4 fitting tests remain unchanged and pass.

### Stage A5 — Manifest, output, and plotting integration

Implementation status (2026-10-18): **implemented and validated**. Every
completed optimizer result now has a canonical `liouville_cops.fit_report`
identity binding the selection/profile, legacy activity, fitted activity
parameters, basis/degree, regularizer, and lambda. The identity is embedded in
the prediction table, diagnostic tables, PNG metadata, `fit_report.json`, and
the corresponding optimizer record in `fit_manifest.json`. A5 adds universal
facet, degree-power, visible/unseen-region, pressure-residual, time/altitude
comparison, residual-plot, and baseline/residual/total-map products through
`tools/cops_fit_reporting.py`. `tools/plot_cops_model.py` verifies manifested
fit-table and report hashes and cross-checks profile, selection, report, and
basis identities before plotting; legacy inputs are visibly labelled. The
normal suite passes 58 tests (5 opt-in skips), and the opt-in small C++/Python
degree-4 and degree-5 fit/replay/plot integration passes. Full Stage-A6
synthetic recovery and flight-data experiments are not implied by this status.

Goal: make every result self-describing and interpretable.

Tasks:

1. Extend manifest schema and identities listed in Section 6.4.
2. Add plots for baseline/residual/total activity and selected/excluded data.
3. Update the pressure comparison plot to display the effective profile and
   enhancement degree/regularization.
4. Add tables for degree power, active/clipped surface area, and visible versus
   unseen activity.
5. Make readers reject or clearly label incompatible/legacy schemas.

Acceptance gate:

- A result can be reconstructed from its manifest and immutable inputs.
- No plot can silently mix profiles, basis conventions, or coefficient files.
- Outputs remain under `<AMPS>/cops_output` (or `COPS_OUTPUT_ROOT`), never in
  `srcCG/LiouvilleCOPS`.

### Stage A6 — Unit, synthetic, and integration validation

Implementation status (2026-10-18): **implemented and validated**. The executable driver
`tools/run_a6_validation.py` runs full Python discovery plus the C++ self-test,
and enables all external-input integration cases with `--integration`.
`tests/test_a6_validation.py` adds production-residual exact/noisy recovery,
all degree-5--6 finite-difference columns, deliberately rank-deficient unseen-
mode behavior, regularization limits, and a production Python-writer/C++-
reader round trip containing all 56 degree-5--8 coefficients. The surrounding
A1--A5 and output-path integration tests provide the basis, selection-count/
reason, zero-enhancement, report, and legacy-workflow coverage. The small
enhanced fit/replay integration enforces
and records the predeclared `1e-6` NG/RG relative-pressure tolerance.
Requirements, fixed thresholds, commands, measured results, and claim limits
are documented in `A6_VALIDATION.md`. This gate establishes implementation
consistency; it does not independently validate the physics, calibration,
full-descent fit, or source uniqueness.

Goal: prove implementation consistency before fitting flight data.

Required tests:

1. Zero-enhancement regression against Stage A0.
2. Python/C++ basis equality through degree 8.
3. Coefficient ordering and file round trip.
4. Selection-profile counts and reason codes.
5. Finite-difference Jacobian checks.
6. Synthetic recovery with one injected degree-5 mode.
7. Synthetic recovery with several modes and controlled noise.
8. Recovery behavior when the injected signal lies primarily on unseen
   facets.
9. Regularization limits: strong `lambda` approaches zero residual; zero
   `lambda` is labeled and may be rank deficient.
10. C++ verification of fitted Python pressures.

Acceptance gate:

- Python and C++ predictions agree to a predeclared relative tolerance, with
  an initial target of `1e-6`; the achieved value must be measured, not assumed.
- Synthetic tests recover only modes supported by their observation geometry
  and explicitly expose unsupported modes.
- Legacy tests and current runnable workflows remain green.

### Stage A7 — Performance characterization

**Implementation status (2026-10-18): complete for the current 100,000-facet
base and measured host.** `tools/benchmark_fit_performance.py` exercises the
production model/residual/Jacobian in degree-isolated workers, measures output
packages with deleted one-evaluation probes, benchmarks independent-job
scaling without BLAS/OpenMP oversubscription, and emits optional immutable
JSON/Markdown evidence. Deterministic policy/report tests are in
`tests/test_a7_performance.py`; the measured tables, exact command, campaign
formula, recommended 8 jobs/1 numerical thread, and interpretation limits are
recorded in `A7_PERFORMANCE.md`. The benchmark must be rerun when the host,
shape, selection, static base, libraries, or fitter changes.

Goal: define a safe compute envelope before full scans.

Tasks:

1. Measure startup, basis construction, residual/Jacobian evaluation, peak
   memory, and output size for degrees 5, 6, 7, and 8.
2. Confirm that fitter parallelism is across independent lambda fits and does
   not accidentally oversubscribe C++/BLAS/OpenMP threads.
3. Determine practical `--jobs` and thread settings on the target machine.
4. Estimate full campaign CPU time/disk from measured small runs.

Acceptance gate:

- A resource table and recommended settings exist.
- Full experiments have bounded output names, storage estimates, and no
  overwrite behavior.

### Stage A8 — Controlled flight-data experiment

Goal: isolate the effect of spatial resolution.

Run order for each profile:

1. Degree-4 no-enhancement control.
2. Degree 5 with a positive lambda scan.
3. Degree 6 with the same documented selection rule.
4. Degree 7 only if degree 6 shows stable, supported improvement.
5. Degree 8 only if conditioning, validation, and compute remain acceptable.

Keep legacy activity, temperature, and `kappa` fixed throughout this stage.
Use the same optimizer tolerances and comparable lambda grid. Save every fit in
a unique directory.

Evaluate:

- NG and RG log-residual RMS separately and together.
- Residuals versus UTC and altitude, particularly around the NG jump and
  closest approach.
- Objective reduction relative to added effective degrees of freedom.
- Jacobian spectrum/effective rank and optimizer convergence.
- Stability of enhancement location/amplitude across neighboring lambdas and
  maximum degrees.
- Activity assigned to visible/contributing versus unseen facets.
- Sensitivity to the three data-selection profiles.
- Exact C++ replay agreement.

Acceptance gate:

- Any claimed enhancement persists across a justified range of lambda/degree
  and is constrained by contributing facets.
- A result driven only by `all_available` transition/saturated NG is reported
  as selection sensitivity, not as evidence of topography.

### Stage A9 — Temperature and CO2 identifiability experiment

Goal: determine whether high-degree activity is absorbing temperature or
composition freedom.

Only after Stage A8, repeat a bounded subset of supported cases with:

1. Analytical temperature parameters free, unordered, matching literal
   current fitter behavior.
2. Analytical temperature parameters free with `--ordered-temperature`.
3. `kappa` free while temperature is fixed.
4. Both temperature and `kappa` free only if the separate cases are numerically
   identifiable.

Compare parameter correlations, Jacobian rank, bounds/floors, residual changes,
and the stability of the spatial enhancement. Do not describe fitted
temperature as a surface-temperature map; the model still uses three global
parameters in its SZA law.

Acceptance gate:

- The report states which quantities are data-constrained and which are fixed,
  bound-limited, or degenerate.
- Releasing nuisance physics does not silently become the new default.

### Stage A10 — Reproducible reporting and documentation

Goal: provide an auditable result package and truthful operating instructions.

Tasks:

1. Update `README.md` with equations, implementation, limitations, provenance,
   and validation evidence.
2. Update `RUN.md` with complete build, base-generation, profile-fit, C++
   replay, and plotting commands.
3. Update `SCRIPTS.md` with the new runner/module/files and their roles.
4. Update `quick_guide.md` with short commands and the most important controls.
5. Produce a comparison table for all profiles/degrees/lambdas actually run.
6. Archive manifests, selected data, coefficients, diagnostics, plots, logs,
   and software/input hashes under `<AMPS>/cops_output`.

Acceptance gate:

- A clean checkout plus documented external data can reproduce the selected
  result without paths tied to `/home/vtenishe`.
- Documentation distinguishes implemented behavior from planned work and
  limitations from conclusions.

### Stage A11 — Scientific interpretation decision

Goal: decide what the evidence permits the project to claim.

Possible outcomes include:

- No stable high-degree signal is supported.
- A pressure improvement exists but is not identifiable as a unique surface
  source.
- A stable body-fixed activity enhancement is supported within the current
  collisionless/source-temperature model.

Any publication statement must identify the selected profile, excluded data,
harmonic degree, regularization, fixed/free parameters, geometry visibility,
and model limitations. Local topography is only one possible explanation of a
near-surface mismatch; instrument range changes, gauge geometry/acceptance,
temperature/composition degeneracy, collisionless-model limitations, shape
resolution, and timing/attitude uncertainty must be considered before making
a causal claim.

## 10. Recommended code locations

The exact organization may change during implementation, but responsibilities
should remain separated:

| Responsibility | Existing/proposed location |
|---|---|
| Current optimizer and cached prediction | `tools/fit_cops.py` |
| Shared cache/manifest I/O | `tools/cops_io.py` |
| New harmonic math | proposed `tools/activity_basis.py`, `src/ActivityEnhancement.{h,cpp}` |
| Profile runner | proposed `tools/run_fit_profile.py` |
| Profile definitions | proposed `fit_profiles/*.json` |
| C++ source/activity model | `src/Source.{h,cpp}` |
| C++ inputs and validation | `src/Input.*`, `src/main.cpp` as applicable |
| C++ provenance | `src/Provenance.cpp` |
| Pressure/activity plotting | existing tools in `tools/`, extended rather than duplicated where practical |
| Tests | `tests/`, with small checked-in fixtures only |

Do not copy the high-order implementation into several scripts. One Python
basis module and one C++ component should be the respective sources of truth.

## 11. Risks and required safeguards

| Risk | Safeguard |
|---|---|
| More coefficients always lower training residual | Positive regularization, degree progression, rank diagnostics, synthetic and stability tests |
| Global harmonics mimic a local source with ringing | Map diagnostics, degree-power report, clipping report, neighboring-degree comparison |
| Bilobed shape makes angular modes physically ambiguous | Treat harmonics as a parameterization; report facet-space maps and visibility |
| Unseen facets acquire large activity | Area-weighted prior and visible/unseen diagnostics; do not interpret unconstrained production globally |
| NG transition drives the solution | Compare all three profiles and show excluded/selected points explicitly |
| Temperature or CO2 freedom trades against H2O structure | Keep them fixed first; release separately in Stage A9 |
| Python fit cannot be reproduced by C++ | Shared convention tests and mandatory C++ replay before scientific use |
| Legacy files/readers break | Separate enhancement file; optional default-off C++ component; regression tests |
| Memory/runtime grows unexpectedly | Compute basis per facet, not per contribution; benchmark before full campaign |
| Absolute paths make results machine-specific | Existing configurable data/output roots and manifest identities; no hardwired home path |

## 12. Stop conditions and questions that require evidence

Pause implementation or interpretation rather than guessing if any of the
following is unresolved:

- The body-fixed frame or facet-centroid coordinates used by the cached mesh
  cannot be verified.
- Python and C++ cannot reproduce the same basis or zero-enhancement pressures.
- The three profile counts differ without an understood input revision.
- The intended treatment of saturated/transition NG changes from Section 4.
- A coefficient scale or bound is proposed without a numerical rationale.
- The optimizer reports rank deficiency or bound-limited solutions that are
  hidden by aggregate residual metrics.
- C++ replay does not match Python to the declared tolerance.
- Required SPICE, COPS, source, or shape inputs lack recorded identities.

The project owner must choose any revised data policy or scientific model
assumption. The implementation may expose alternatives, but must not silently
select one.

## 13. Definition of done

This amendment is complete only when:

1. All three profiles are versioned, tested, and produce auditable selections.
2. The new harmonic convention is fully specified and matched by Python/C++.
3. Enhancement-disabled operation preserves the legacy result.
4. Degrees 5--8, coefficient files, analytic Jacobian, and regularization pass
   unit and synthetic tests.
5. C++ can replay Python fits within a measured tolerance.
6. Performance is measured and the approved flight-data campaign is complete.
7. Results include conditioning, visibility, clipping, profile sensitivity,
   and fixed/free-parameter diagnostics—not only pressure curves.
8. `README.md`, `RUN.md`, `SCRIPTS.md`, and `quick_guide.md` describe only the
   behavior actually implemented and tested.
9. All generated products live outside the source tree under the configured
   COPS output root.
10. Scientific conclusions are explicitly limited to what the validation and
    model assumptions support.
