# LiouvilleCOPS fitting: implementation and validation roadmap

Prepared 2026-10-17 from the local implementation, former fitter revision 1.9,
processed data and fit reports that were present during the documentation audit.

This is a plan for future work, not an implemented-feature list. Creating
this document did not restore/change the fitter, rerun the full descent,
launch full optimizations, or validate the instrument calibration. Reinspect
the actual working tree before executing any work package.

Subsequent path-layout update: the fitter was restored from revision 1.9 with
new output routing. New scientific products go to `<AMPS>/cops_output` (or
`COPS_OUTPUT_ROOT`), with legacy read fallback. This does not complete the
scientific validation work packages below. The formerly local scientific
results were subsequently removed; regenerate them under `<AMPS>/cops_output`.

## 1. Objective, scope, and operating rules for a future Codex session

Objective: make the standalone LiouvilleCOPS fitting workflow reproducible,
diagnostically reliable, and scientifically interpretable, then support
explicitly qualified publication claims. A low pressure residual alone is
not the completion criterion.

Scope is this application and its subdirectories only. External observation,
shape, raw COPS, and SPICE resources referenced by its inputs are data
dependencies, not authorization to modify the surrounding AMPS/srcCG code.

Operating rules:

- Read `README.md`, `RUN.md`, `quick_guide.md`, this roadmap, and applicable
  local instructions before implementing a package. Verify named functions
  and options still exist; they describe the inspected version, not a frozen API.
- Preserve all existing science outputs, inputs, and unrelated user edits.
  Use unique experiment directories. Never overwrite a historical fit to
  make a new test pass or regenerate a golden reference without explaining why.
- Use local file-editing tools and reviewable patches. Do not change source
  outside this standalone application.
- Preserve a legacy mode/baseline while changing selection, likelihood,
  parameterization, or output semantics. Record intentional numerical changes.
- Treat proposed filenames/options below as design suggestions; they do not
  exist unless explicitly listed in the current inventory. Update documentation
  only to claim features that were actually implemented and tested.
- Do not invent instrument error bars, calibration constants, collision cross
  sections, time offsets, thermal parameters, or missing data. Record unknowns
  and obtain evidence/user decisions where required.
- Do not download/install dependencies, contact collaborators, publish, commit,
  or launch large run campaigns merely because the roadmap mentions them.
  Follow the actual authorization and resource budget for the future task.
- Start with fast tests and small isolated runs. Estimate memory/disk/CPU from
  measurements before full mesh, multi-start, uncertainty, or rotation scans.
- At handoff, report changed files, tests actually run, results, remaining
  uncertainties, and the next bounded step. Do not mark a package complete
  because its code exists if its acceptance tests/evidence are missing.

Suggested kickoff request for a future session:

> Work only in LiouvilleCOPS. Read ROADMAP.md and implement W00 first (or the
> next approved incomplete package). Preserve existing outputs. Run the fast
> tests and an isolated small integration check where available. Do not change
> scientific defaults or launch a full fitting campaign without reviewing the
> relevant decision gate and resource budget. Report evidence and remaining work.

## 2. Current facts to preserve as evidence

| Item | Audited state; recheck before relying on it |
|---|---|
| Fitter | Working `tools/fit_cops.py` retained from former revision 1.9; path handling updated, fitting mathematics unchanged. No local RCS archive remains. |
| Base | The audited `output_fit_base/` had 1,150 rows and static version-1 caches, but was removed during output cleanup. Regenerate it before fitting. |
| Free parameters | Default 25 H2O coefficients + 3 analytical-temperature parameters + 1 CO2 scale. |
| Objective | Equal-gauge normalized log residuals plus optional H2O map-prior penalty; lambda defaults to zero. |
| Historical single fit | 200-evaluation limit; condition number about 1.1e16; Tmin floor makes reported temperatures effectively constant. |
| Historical selected L-curve fit | Lambda about 3.162e-3; status 0 at 1,000 evaluations; T75 at 30 K; about 64% of clipped activity integral on faces labeled unseen. Not an unseen production fraction. |
| Data selection | Saved fits use 839 NG / 1,146 RG points; current default selection can use 930 / 1,146. Do not assign current defaults retroactively. |
| Data audit | Two unmatched reference/HK stop times and one NG zero-versus-positive mismatch remain unresolved. |
| Fitted table | Total pressures replaced; species columns set to -1; production and other fields partly inherited. C++ verification needed for physical diagnostics. |
| Rotation | Observation-time orientation/illumination included; molecular-flight rotation excluded from fitting. Exact forward kernels are unsupported by fitter. |
| New boundary tool | `tools/rebuild_cops_boundary.py` already implemented; refuses existing destinations, preserves settings, regenerates/validates a static base. |
| Boundary-tool tests | `tests/test_rebuild_cops_boundary.py`: 13 fixture/unit tests and one optional integration test with three one-point C++ runs passed at implementation time. This is not full-fit validation. |

Relevant code map:

| Area | Existing locations to inspect |
|---|---|
| Observation times/pressures and selection | `src/Observations.cpp`, `src/main.cpp`; fitter `read_model_file`, `apply_data_file`, `run_intervals`, `ng_saturated`; `tools/collect_cops_sources.py`, `tools/cops_hk_table.py` |
| Cached forward model and derivatives | Fitter `read_elements`, `read_kernel`, `Model`, `Problem`; `src/Source.cpp`, `src/Liouville.cpp`, `src/Gauges.cpp` |
| Optimization and uncertainty | Fitter `fit_one`, `fit_worker`, `Problem.data_residuals`, `Problem.jacobian`, `Problem.regularization` |
| Outputs and lambda selection | Fitter `write_fit_data`, `write_fit`, `activity_diagnostics`, `max_curvature`, `main`; plotters in `tools/` |
| Geometry and numerical verification | `src/Geometry.cpp`, `src/Mesh.cpp`, `src/BVH.cpp`, `src/main.cpp` self-tests |
| Safe base regeneration | `tools/rebuild_cops_boundary.py`; reuse its safeguards, do not duplicate a destructive run launcher |

## 3. Work-package sequence and decision gates

All packages below start as **not implemented** by this roadmap. Existing
boundary tooling is a dependency, not evidence that these packages are done.

| Package | Goal | Prerequisites |
|---|---|---|
| W00 | Establish a reproducible, protected baseline | None |
| W01 | Make input/cache/output provenance and validity explicit | W00 |
| W02 | Audit observations and unify selection policy | W00; integrate manifests from W01 |
| W03 | Expand mathematical/forward-model regression tests | W00–W02 |
| W04 | Improve optimizer, parameter, and rank diagnostics | W01–W03 |
| W05 | Add controlled multi-start and convergence-aware lambda scans | W04 |
| W06 | Define/test alternative observation-error models | W02–W04; scientific decision |
| W07 | Quantify identifiability and held-out/synthetic performance | W03–W05; explicit objective from W06 or legacy mode |
| W08 | Add defensible parameter/prediction uncertainty analysis | W04–W07; likelihood decision |
| W09 | Assess numerical accuracy and physical-model sensitivity | W01–W04; feed findings into W07/W08 |
| W10 | Implement selected physics extensions only when justified | W07–W09; explicit scientific approval |
| W11 | Produce a reproducible release and claim-specific evidence package | Relevant preceding gates satisfied |

First practical milestone: W00–W03, with no silent changes to the physical
model or default objective. W09 calibration evidence gathering can proceed
before an extensive uncertainty campaign; do not condition uncertainty on
assumptions already known to be unresolved without labeling that limitation.

Decision gates requiring explicit evidence or user choice:

| Gate | Decision that must not be invented |
|---|---|
| G1: observation authority | Which reference/HK/L4 products, timestamp conventions, and missing-value policy define the scientific dataset? |
| G2: selection and likelihood | Treatment of the late NG mode change/saturation, quality flags, uncertainty weighting, and time correlation. |
| G3: parameter constraints | Whether to change ordered-temperature upper bounds, permitted temperature ranges, activity basis/prior, or default degree. |
| G4: instrument/geometry assumptions | Supported sensor-temperature mapping, species sensitivities, pressure conversion, acceptance, head offsets, and SPICE revisions. |
| G5: added physics | Whether evidence justifies thermal lag, new source laws, collisional transport, or exact-rotation fitting; define validation and compute scope. |
| G6: publication claim | Pressure reconstruction, constrained-region activity, global production, composition, or thermophysical inference require different evidence. |

These gates do not block coding neutral audit/test infrastructure. They block
silently promoting an unverified assumption into a default or physical claim.

## W00 — Protected baseline and reproducible environment

Deliverable: a small repeatable baseline experiment with an explicit source
snapshot, environment description, selected data, and numerical reference.

Implementation tasks:

1. Inspect local status and any newly generated outputs. Record the current executable,
   Makefile flags, Python/NumPy/SciPy/Matplotlib versions, and CSPICE location.
   Do not assume the parent repository's Git revision alone identifies these
   untracked application files; include source-file hashes.
2. Preserve and snapshot the working fitter and `cops_paths.py`. No local RCS
   archive remains. If either file is missing, recover the updated version from
   source control or a complete copy under an authorized implementation request.
3. Add a dependency/environment record based on a tested environment, not
   guessed version pins. Keep runtime installation separate from execution.
4. Establish a small fixture covering early descent, late descent, NG missing
   data, exclusion/saturation, and at least one geometry edge case. Keep the
   same input/prior and explicit selection flags in each comparison.
5. Capture initial pressures, residuals, objective components, parameter names,
   and Jacobian diagnostics before refactoring. Preserve historical result
   files unchanged and distinguish them from newly generated baselines.

Tests/acceptance:

- Existing boundary tests pass; self-test and a new small C++/Python comparison
  are run when dependencies exist, otherwise explicitly reported as skipped.
- Baseline can be repeated into a new directory with identical selection and
  predictions within recorded tolerances. Do not require optimizer bitwise
  identity across different libraries/hardware.
- Missing assets produce actionable errors; no tests download data silently.
- Every baseline artifact identifies the exact code and inputs that produced it.

## W01 — Provenance, cache validation, and truthful outputs

Deliverable: versioned machine-readable run/fit metadata and an output format
that cannot silently present base production as fitted production.

Primary locations: fitter readers and `write_fit`; C++ output headers;
boundary driver and plotters. A shared `tools/cops_io.py` module is a possible
new design, not an existing dependency.

Implementation tasks:

1. Define a manifest schema: run kind, input paths/hashes, effective settings,
   code/binary/dependency identities, mesh/SPICE/HK/observation identities,
   kernel format, units, selection hash, prior, bounds, objective, and parent
   run. Stream hashes; avoid loading large kernel sets into memory to hash.
2. Separate states: process completed, optimizer converged, cache verified,
   C++ prediction agreement checked, and scientific review status. A zero
   process exit status is not a scientific approval flag.
3. Validate point IDs, duplicate IDs, table/kernel correspondence, versions,
   record lengths, face-index bounds, and cross-run mesh/source identity.
   Retain explicit legacy import with warnings where old caches lack metadata.
4. Prevent new experiments from overwriting old output by default. If resume
   is later added, require matching code/input/objective/selection identities;
   do not call a partially written directory a valid checkpoint.
5. Fix fitted-table semantics: recompute fields correctly, or mark unavailable
   derived fields explicitly and version the schema. Preserve compatibility
   deliberately, updating every reader/plotter and documenting migration.
   Do not copy inherited Q values under labels suggesting fitted production.
6. Keep original observations, selected fit data, modeled pressures, and C++
   verification flags distinct. Archive generated inputs without relying only
   on mutable source inputs or stale absolute paths.

Tests/acceptance:

- Corrupt/truncated caches, duplicate rows, swapped kernels/meshes, stale inputs,
  and mismatched selection hashes fail before optimization.
- A fixture changing activity produces changed recomputed production or an
  explicit unavailable field, never an unlabeled inherited base value.
- Plotters can distinguish forward, fit, and verified products and reject
  incompatible schemas with a clear message.
- Failure/interruption retains diagnostic metadata but cannot be mistaken for
  a completed fit. Existing outputs remain untouched.

## W02 — Data audit and one explicit selection policy

Deliverable: an auditable observation table/provenance report and a common
selection implementation used by both collector and fitter.

Primary locations: `cops_hk_table.py`, `collect_cops_sources.py`, fitter input
and selection functions; `src/Observations.cpp` if schema support is changed.

Implementation tasks:

1. Audit raw BG labels/records against processed HK and actual observation
   inputs, preserving precision, units, start/stop/midpoint timestamps, and
   product identities. Use `/data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg` only if
   available; keep paths configurable and never fabricate a successful audit.
2. Investigate 09:13:04 and 09:14:04 unmatched reference times and the 01:30:04
   NG zero-versus-positive mismatch. Record evidence and unresolved status;
   do not interpolate or replace measurements silently.
3. Define stable observation IDs and explicit source types: measured,
   derived, digitized, reconstructed. Keep UTC parsing/normalization separate
   from the original timestamp text and interval-averaging convention.
4. Centralize selection with a versioned, documented precedence for
   positive/valid data, inside-body checks, FitInterval, gauge exclusions,
   saturation, and replacement-table flags. Preserve legacy behavior as an
   explicit mode until G1/G2 approves a default policy change.
5. Fix collector/fitter disagreement when ExcludeNG overlaps FitInterval;
   test boundaries inclusively and ensure re-enabling one reason cannot
   accidentally bypass a different exclusion rule in the new policy.
6. Export reasons and selection counts per gauge and interval. Detect duplicate
   IDs, empty selections, invalid flags, unsupported dates, and empty/high-range
   saturation sets without inventing a plateau threshold.
7. Keep source-comparison overlays separate from genuinely independent data.
   Do not treat digitized representations of the same observation as new
   independent samples or as automatic replacements for archive gaps.

Tests/acceptance:

- Synthetic timestamp, interval-boundary, missing-value, duplicate, and
  overlapping-exclusion fixtures exercise every decision branch.
- Collector and fitter return identical masks for a declared shared policy.
- Current audited counts are reproducible in explicit legacy mode on matching
  inputs; a corrected mode's changed counts are explained, not hidden.
- Data audit produces row-level match/mismatch evidence and never claims the
  products are identical while discrepancies remain unresolved.

## W03 — Mathematical and forward-model regression tests

Deliverable: fast automated tests plus opt-in integration tests that cover
the supported physical/parameter branches, not just the initial prior.

Primary locations: fitter `Model`/`Problem`, `Source.cpp`, `Liouville.cpp`,
`Gauges.cpp`, and tests under this application.

Implementation tasks:

1. Make fitter logic importable without launching a run; preserve CLI behavior.
   Use small fixtures for algebra and a separately tagged mission-data test.
2. Check Python/C++ predictions for flat/rough boundaries, fixed/free parameter
   combinations, supported gauge choices, multiple illumination/shadow states,
   and representative trajectory points. Include the actual fitted parameters
   from new tests, not only the original prior.
3. Verify analytic Jacobians with scale-aware finite differences at multiple
   interior parameter states. At clipping/floor boundaries test branch behavior
   and one-sided/directional limits; a central difference across a kink need
   not equal the selected branch derivative. Test the log-pressure floor too.
4. Independently test coefficient ordering/basis functions, emitted-law
   normalization, species scaling, pressure units, and the H2O prior penalty.
   Do not substitute a normalized harmonic library without an explicit basis
   transformation and regression evidence.
5. Add guards/tests for all-zero prior normalization, nonpositive temperature,
   invalid degree/bounds, empty gauge data, zero predicted pressure, singular
   Gram matrices, and nonfinite residuals. Errors should identify the problem,
   not silently generate plausible-looking results.

Tests/acceptance:

- CI/fast tests need no external mission assets; integration tests declare
  prerequisites and skip explicitly when unavailable.
- Regression tolerances reflect serialization and verified numerical error.
  Keep the existing initial-pressure consistency default 1e-6 unless evidence
  justifies changing it. Do not arbitrarily tighten the Jacobian threshold:
  the earlier sampled check had a worst relative discrepancy about 5.87e-5.
- Failures report parameter/point/branch and absolute plus relative errors,
  including safeguards for near-zero reference derivatives.

## W04 — Optimizer, temperature constraints, and rank diagnostics

Deliverable: reliable termination/constraint diagnostics and explicit handling
of weak/inactive parameter directions without changing physics silently.

Primary locations: fitter `Problem`, `fit_one`, result writer, and diagnostics.

Implementation tasks:

1. Export structured status, success, nfev/njev where available, optimality,
   active bounds, scaling, parameter coordinates, initial/final data and prior
   cost separately, and selected counts. Distinguish execution failure,
   evaluation-limit termination, and satisfied stopping criteria.
2. Report data-only and regularized Jacobian spectra/ranks with scaling and
   numerical rank thresholds. Regularization can stabilize a problem without
   creating observational information; do not label its rank as data resolution.
3. Report temperature branch occupancy and active/floored H2O regions,
   preferably weighted by their actual contribution to each selected gauge.
   Mark unidentifiable/inactive directions explicitly instead of presenting
   zero pseudoinverse sigma as exact knowledge.
4. Under G3, choose the intended ordered-temperature bounds. If all physical
   temperatures should lie in [LO, HI], design a valid coupled parameterization
   or constrained solve. Merely putting HI on each independent increment does
   not bound Tmax. Do not clamp trial temperatures after evaluating derivatives.
5. Derive/test transformed Jacobians and full covariance propagation back to
   physical Tmax/T75/Tmin, including cross-covariances. Preserve a documented
   legacy parameterization for reproducibility.
6. Consider stable SVD/QR-based uncertainty/rank calculations rather than
   relying exclusively on normal-matrix pseudoinversion. Document any changed
   rank convention and compare results on well- and ill-conditioned fixtures.

Tests/acceptance:

- A forced evaluation-limit run is never labeled converged or publishable.
- Known synthetic null directions are flagged as unconstrained, not certain.
- Ordered temperatures obey the selected physical constraints throughout the
  solve, with finite-difference agreement for the transformed derivatives.
- Inactive-temperature fixtures reproduce the diagnosed floor behavior of the
  historical fit. Legacy-mode predictions remain reproducible.

## W05 — Controlled starts, checkpointing, and lambda scans

Deliverable: repeatable search over starting conditions/regularization with
transparent ranking and failure handling.

Implementation tasks:

1. Expose initial parameters independently of the prior (the internal
   `fit_one(lam, x_start=...)` already provides a starting point for this).
   A new start must not silently change the base activity prior, fixed
   higher-degree coefficients, kernels, or observations.
2. Add seeded, bounded multi-start experiments with recorded initial vectors;
   preserve comparable objective/selection metadata. Capture every outcome,
   including failures and bound-dominated fits.
3. Define convergence-aware scan summaries. Do not select status-0/invalid
   fits as accepted L-curve candidates. If too few eligible positive-lambda
   points remain for the current curvature method, report no accepted corner.
4. Keep independent-start scans distinct from optional warm starts. Retry
   budgets and checkpoint/resume rules must be explicit, finite, and bound
   to matching manifest identities.
5. Measure peak RAM and elapsed time for one worker; cap concurrent jobs using
   the agreed resource budget. Avoid nested BLAS/OpenMP oversubscription.
   Optimize repeated cache loading only after parity tests exist.

Tests/acceptance:

- Identical seeds/configurations reproduce starting vectors and comparable
  results. Resume rejects a changed prior, dataset, objective, or base.
- A synthetic scan containing nonconverged fits cannot silently choose one
  as the accepted corner; diagnostic plots still display rejected cases.
- Reports show parameter/prediction dispersion across starts, not merely the
  best cost. Multiple similar-cost distinct solutions remain visible.

## W06 — Observation-error model and weighting experiments

Deliverable: explicitly named objective modes with documented assumptions;
retain the current equal-gauge log objective as a reproducible baseline.

Scientific gate G2: establish which uncertainty/quality information is valid
for the chosen observation product. L4 density uncertainties do not
automatically apply to raw NG/RG pressures or composition-decomposed outputs.

Implementation tasks:

1. Specify candidate residual/likelihood definitions mathematically before
   coding: legacy equal-gauge log weights, evidence-supported uncertainty
   weights, and optional duration weights as a different scientific objective.
   Document units, normalization, gauge balance, and treatment of zeros.
2. Inspect residual time correlation and instrument-mode dependence. If a
   correlated-error model is proposed, state how it is estimated without
   using held-out data and how its parameters affect likelihood normalization.
3. Treat saturation only with a justified limit/likelihood. A censored sample
   is not an ordinary measurement at the limit; excluding it and using a
   censored likelihood are distinct choices. Keep exclusion available.
4. If adding robust losses, distinguish data residuals from prior residuals;
   applying a generic robust loss to the entire stacked vector can also
   change the intended regularization. Derive the actual objective/Jacobian.
5. Export effective weights, objective definition/version, and per-interval
   influence. Do not assign fabricated uncertainty values to enable a mode.

Tests/acceptance:

- Unit tests match hand-calculated residuals/likelihoods and derivative checks.
- Missing error metadata causes a clear refusal or explicitly selected legacy
  mode, not hidden guessed values.
- Synthetic variable-cadence/correlated/saturated examples demonstrate what
  each mode estimates and expose resulting differences.
- Comparisons retain the same selection/prior or clearly identify changes;
  rescaling the objective is accompanied by reexamining lambda, not reusing
  its old numerical value as if its meaning were unchanged.

## W07 — Identifiability, synthetic recovery, and held-out validation

Deliverable: evidence of which quantities are recoverable/predictive and
which remain prior- or assumption-dependent.

Implementation tasks:

1. Design controlled recovery experiments with known activity, temperature,
   and CO2 scale, using the real sampling geometry where available. Separate
   noiseless algebra tests from realistic-error and model-mismatch experiments.
2. Test reduced degree, fixed/free temperature and kappa, alternate priors,
   lambda ranges, and multiple starts. Separate changes in a prior from changes
   in an initial guess. Preserve raw coefficient conventions.
3. Create contiguous training/validation/test time blocks and an explicit
   strategy for temporal correlation. Use only training/validation information
   to choose model complexity, weights, and lambda; reserve final test blocks.
   Do not claim random neighboring samples are independent validation.
4. Ensure held-out pressures cannot enter the fitted mask. Geometry may be
   precomputed for all times, but observed held-out pressures must not enter
   objective evaluation or tuning. Record separate selections/hashes.
5. Calculate quantitative source sensitivity/resolution using the selected
   gauge contributions/Jacobian, not only whether a face was visible at least
   once. Distinguish constrained surface regions from extrapolation.
6. Evaluate global and constrained-region production via complete C++ reruns.
   Compare predictions and derived quantities across acceptable solutions,
   not only coefficient differences in a nonnormalized basis.

Tests/acceptance:

- Injected identifiable combinations can be recovered under stated tolerances;
  intentionally nonidentifiable examples correctly report ambiguity.
- A deliberately altered held-out pressure leaves the fitted training result
  unchanged. Tuning/test separation is enforced by the runner.
- Reports show held-out errors by gauge/time/altitude and map/production
  variability across models. Same-model synthetic success is not described
  as proof of real-world calibration or collisionless validity.

## W08 — Parameter and prediction uncertainty

Deliverable: clearly labeled conditional uncertainty estimates and a separate
systematic/model-sensitivity assessment, not one unsupported universal sigma.

Prerequisites: defensible objective/likelihood and diagnostics from W04–W07.

Implementation tasks:

1. Correct local covariance reporting first: state rank, scaling, prior,
   bounds, and inactive directions; propagate full covariance through ordered
   coordinates and derived quantities where the approximation is applicable.
2. Add parameter profiles or selected physical-quantity profiles, reoptimizing
   nuisance parameters and recording convergence at each profile point.
3. If using bootstrap, justify block structure or instrument-product grouping
   using the correlation/measurement model. Do not independently resample
   strongly correlated adjacent records without acknowledging the assumption.
4. If using posterior sampling, explicitly define the observation likelihood,
   priors, coordinate-transform effects, diagnostics, and compute budget.
   The existing Tikhonov objective alone is not a complete calibrated
   probabilistic model for every physical/systematic uncertainty.
5. Propagate to pressure prediction bands, production rates, and activity
   summaries; distinguish measured-time interpolation from gap/region
   extrapolation. Keep calibration/shape/physics ensembles identifiable rather
   than folding them into an unexplained error bar.

Tests/acceptance:

- Synthetic coverage experiments and bound/null fixtures test the claimed
  meaning of intervals. Nominal coverage is reported as tested, not assumed.
- Reports explicitly state uncertainties conditioned on fixed instrument,
  geometry, prior, and physics choices, and which sources remain unquantified.
- No unconstrained direction is reported as zero uncertainty simply because
  a pseudoinverse discarded it.

## W09 — Numerical convergence and physical sensitivity

Deliverable: a structured matrix of tests linking numerical/physical choices
to pressure and inferred-quantity changes.

Implementation tasks:

1. Numerical studies: vary mesh/shape resolution, `SubdivisionAccuracy`, and
   `MaxSubdivisionLevel`; examine representative early/late/near-surface and
   inside/outside points. Define an error target relative to the scientific
   precision sought, not an arbitrary universal tolerance.
2. Run flat/rough comparisons through the existing rebuild driver with the
   same temperature/source/selection starting assumptions, then refit both
   when comparing fitted conclusions. Separate forward differences at fixed
   parameters from changes caused by reoptimization.
3. Flight-time rotation: compare static and exact forward results at the same
   fitted parameters; vary speed nodes/range and path resolution. Test frame
   consistency with nonzero spacecraft velocity and gauge offsets, since
   current analytic rotation tests do not establish that general case.
4. Instrument evidence: identify authoritative calibration documentation for
   sensitivity factors, NG acceptance/head-on factor, RG pressure conversion,
   sensor-temperature mapping, and physical gauge-head location. Cite exact
   evidence and uncertainty; do not tune these purely to reduce residuals.
5. Record the difference between SPICE-kernel revisions and validate coverage
   and frame conventions. Changes in geometry/shape require consistent bases;
   compare time-aligned predictions rather than mismatched altitude grids.
6. Assess collisionless applicability using independently justified cross
   sections, density estimates, and characteristic scales. State uncertainties
   and limit the supported altitude/time domain where evidence is insufficient.
7. Examine whether residuals correlate with local time/illumination or mode
   transitions in ways suggesting missing thermal/instrument physics; avoid
   claiming a diagnosis from correlation alone.

Tests/acceptance:

- Publish a sensitivity table with inputs, unchanged controls, outputs,
  numerical tolerances, and evidence links. Unavailable evidence is flagged.
- Chosen resolutions are supported by a convergence sequence, not a single
  larger run. Material parameter/prediction changes are propagated to claims.
- Exact-rotation forward validation is not mislabeled exact-rotation fitting;
  current version-2 kernels are still rejected by the static fitter.

## W10 — Optional physics extensions, with explicit scope approval

This is not required to make the existing approximation reproducible. Select
extensions only after W07–W09 establishes a scientific need and G5 defines
parameters, data support, independent benchmarks, and a compute budget.

Candidate tracks, each requiring its own design and acceptance plan:

- **Thermophysical source:** introduce an energy-balance/thermal-memory model
  with defined material properties, time stepping, initialization/spin-up,
  illumination history before the observed window, shadow history, and species
  production coupling. A fitted phase shift is not automatically equivalent
  to thermal inertia. Validate energy/flux units and limiting solutions.
- **Exact-rotation fitting:** implement temperature/species-dependent forward
  evaluations or a validated surrogate/cache design, with corresponding
  derivatives and frame checks. Do not relabel version-2 files as version 1
  or reuse temperature-dependent kernels while varying temperature.
- **Alternative angular/source models:** define normalization, positivity,
  parameter identifiability, and prior conventions. If changing the harmonic
  basis, transform coefficients explicitly and preserve legacy comparisons.
- **Collisional extension:** if required by applicability checks, define an
  appropriate transport treatment/interface and benchmarks. Do not append an
  arbitrary density correction to claim a validated collisional model.

Acceptance: independent analytic/numerical benchmarks, regression to the old
model in the appropriate limit, uncertainty/synthetic tests for added
parameters, and updated cache/schema/invalidation rules. Added parameters
must not be declared measured merely because an optimizer returns values.

## W11 — Reproducible scientific release and publication evidence

Deliverable: an immutable, documented experiment package and a claim/evidence
matrix. Packaging or public release requires the user's explicit instruction.

Implementation tasks:

1. Select the claim scope at G6. For each claim, identify supporting runs,
   convergence status, validation evidence, uncertainty, assumptions, and
   unsupported extrapolations. Methods/exploratory claims can be narrower than
   unique global production/composition/thermophysical claims.
2. Freeze code/environment and include manifests, exact activity-basis
   definitions/order, inputs, selection tables, physical parameters, optimizer
   settings, diagnostics, and complete C++ verification outputs.
3. Provide one documented reproducibility entry point that fails safely when
   external assets are absent. Document acquisition/provenance and redistribution
   restrictions without bundling assets that cannot be redistributed.
4. Generate figures/tables from identified verified products, not manually
   altered plots. Clearly label measurements, digitized/derived sources,
   fitted predictions, gap reconstructions, and uncertainty meaning.
5. Include rejected/nonconverged models and sensitivity results in the audit
   trail. Do not retroactively replace exploratory history with converged runs
   under the same identifiers.
6. Draft a limitations/data-code-availability section matching the actual
   release. Verify target-journal requirements when a journal is selected;
   this roadmap is not a substitute for those requirements.

Acceptance:

- A clean local replay reproduces selection and verified key results within
  documented tolerances using only declared dependencies/assets.
- No reported fitted production comes from inherited Python-table fields.
- Every physical claim has a traceable evidence/uncertainty entry; unresolved
  calibration, collisionality, uniqueness, or thermal assumptions are explicit.
- Publication-quality acceptance is a scientific/user review, not a Boolean
  assigned solely by a passing test suite or optimizer status.

## 4. Test/run tiers and resource policy

| Tier | Intended content | Permission/resource expectation |
|---|---|---|
| A: fast | Pure parsing, algebra, selection, manifest, fixtures, fake solver | Run for each implementation change; no network/mission assets. |
| B: small integration | Existing self-test; a few real geometry points; C++/Python agreement | Run with known local assets; isolated output; report executable identity and skips. |
| C: bounded fitting | Selected small datasets, derivative tests, synthetic fits, limited starts | Measure memory/time; agree point/start/evaluation limits before campaigns. |
| D: full science | Full descent, mesh/physics ensembles, multi-start scans, uncertainty runs | Explicit resource/output plan and approval for material expansion; retain failed experiments. |

Existing commands available now (not proposed future tools):

```bash
python3 -m unittest discover -s tests -p 'test_rebuild_cops_boundary.py' -v

env COPS_RUN_INTEGRATION=1 python3 -m unittest discover \
  -s tests -p 'test_rebuild_cops_boundary.py' -v

env OMP_NUM_THREADS=4 ./liouville_cops selftest.in
```

The opt-in integration command uses external resources from `fit_base.in`
and the existing executable. Record whether it was rebuilt; do not describe
testing an old binary as testing new C++ changes. For L-curve workers, separate
`--jobs` from C++ OpenMP threads and numerical-library thread pools. Inspect
RAM use before increasing either. Benchmark/reduce repeated loads only after
correctness is established.

## 5. Completion ledger for later sessions

Do not prefill test results or mark the future packages complete. After each
approved package, append an entry containing:

```text
Package ID / status: pending | in progress | implemented-unvalidated | validated | blocked
Date and source identity:
Changed files:
Preserved legacy behavior / intentional changes:
Scientific decisions and evidence (G1..G6):
Commands actually run and resource limits:
Tests passed, failed, or skipped (with reasons):
Experiment/output paths and manifests:
Acceptance criteria met / not met:
Remaining risks or required user input:
Next bounded task:
```

If blocked on a scientific choice or missing evidence, implement only
independent safe infrastructure and report the unresolved gate. Do not fill
the knowledge gap with a convenient assumption or mark the full objective
complete because all currently runnable tests pass.
