# AGENTS.md — AMPS `srcMoon` Verification and Validation

## Scope

These instructions apply to all work on the AMPS lunar application `srcMoon`, its verification/validation framework, and its use of the external validation-data repository.

Canonical validation-data root:

```text
/data/vtenishe/moon_validation_data
```

The AMPS repository root is the active checkout containing `srcMoon/`. Do not assume an absolute source path until the current checkout establishes it.

Read before working:

```text
<AMP_ROOT>/srcMoon/development_and_validation.md
/data/vtenishe/moon_validation_data/development_and_validation.md
/data/vtenishe/moon_validation_data/AGENTS.md
```

If the source-tree and data-root roadmaps differ, do not silently reconcile them. Record both SHA-256 hashes, summarize the differences, and identify which contract controls each code/test decision.

## Scientific rules

1. Never create synthetic, mock, assumed, fabricated, hand-tuned, or model-generated observations and use them as validation evidence.
2. Missing authoritative data, missing metadata, missing production capability, or unresolved observation geometry means `SKIPPED / NOT VALIDATED`, never `PASS`.
3. Never weaken an acceptance threshold because a test failed.
4. Never change calibration/training versus hold-out membership after inspecting hold-out residuals.
5. Never select observation epochs, profiles, events, energy/mass/angle cuts, or quality cuts because they agree with AMPS.
6. Everything under `moon_validation_data/**/raw/**` is immutable.
7. Every processed validation product must be reproducible by a committed script. No spreadsheet/manual editing of science values.
8. Native labels/headers/metadata control interpretation. Never guess units, scale factors, time systems, frames, fill values, quality flags, or instrument geometry.
9. Preserve real gaps. Never silently interpolate, smooth, extrapolate, gap-fill, or renormalize observations.
10. Observation processing may use only documented archive/publication corrections, never AMPS residuals.

## Status semantics

Keep these distinct:

- `PASS` — the declared criterion was evaluated and met.
- `FAIL` — capability/data existed, the test ran, and the criterion was not met.
- `ERROR` — build, runtime, parser, data-corruption, or test-infrastructure failure.
- `SKIPPED` — required capability, configuration, dependency, or authoritative data was unavailable.

Capability labels:

- `VERIFIED` — local/analytical and linked numerical gates pass.
- `VALIDATED` — verified and an independent observational comparison passes with frozen parameters.
- `EXPERIMENTAL` — implementation exists but required gates/data remain incomplete, skipped, or failed.
- `DISABLED` — not enabled/built in production.

A source/preflight PASS is never a linked-runtime PASS.

## Preserve test IDs

The roadmap freezes:

```text
U01–U26   stand-alone/local verification
I01–I46   linked AMPS integration/validation
```

Do not repurpose an existing ID.

If analysis reveals missing coverage:

1. document the gap;
2. identify whether an existing test can legitimately be extended without changing its meaning;
3. otherwise add the proposal to `srcMoon/test/TEST_GAP_PROPOSALS.md`;
4. do not assign a permanent new ID until the roadmap/test registry is updated consistently.

## Production-code rule

Tests must exercise the same production kernels, dispatch paths, geometry services, surface callbacks, movers, chemistry, sources, and samplers used by normal `srcMoon`.

Allowed:
- thin test adapters;
- deterministic fixtures;
- independent analytic/numerical reference evaluators;
- non-invasive diagnostics.

Not allowed:
- duplicate test-only force laws;
- duplicate test-only chemistry/source formulas;
- bypassing the production path;
- freezing current production output and calling it an independent physical truth.

## Analysis-before-change rule

Before changing production code or tests:

1. inventory current implementation;
2. identify production entry points;
3. identify compile-time/runtime configuration;
4. identify current tests and coverage;
5. identify data dependencies;
6. identify an independent oracle/reference;
7. write the test contract;
8. then implement.

For the initial `srcMoon` audit, do not modify production physics.

## Required component inventory

Maintain a current inventory covering at least:

- initialization/configuration;
- species table, masses, charges, particle weights, local time steps;
- mesh/AMR/refinement;
- spherical and LOLA surface geometry;
- ray/surface intersections and normals;
- illumination, terrain horizon, Earth eclipse, PSR;
- simple, Diviner-backed, and dynamic thermal state;
- lunar gravity;
- Sun/Earth differential gravity and rotating/non-inertial terms;
- Na radiation pressure/shadow gating;
- Lorentz force and ion mover;
- surface impact classification;
- sticking/adsorption/residence/re-emission;
- surface inventories and cold trapping;
- photoionization/photodissociation;
- electron-impact ionization;
- charge exchange;
- plasma/magnetic drivers;
- He source/reservoir;
- Ne source/accommodation;
- radiogenic Ar geography/transients;
- Na PSD/thermal/sputtering/impact sources;
- meteoroid forcing;
- H2O/OH migration/chemistry;
- LOS/brightness/column-density operators;
- Kaguya, LADEE, LAMP, LACE, PACE geometry/operators;
- diagnostics, sampling, budgets;
- checkpoint/restart if present;
- MPI/OpenMP/stochastic behavior;
- provenance/run manifests.

For each component record source files, functions/classes, production caller, switches, inputs, outputs, units, frames, external data, stochastic state, current U/I coverage, and uncovered failure modes.

## Test specification standard

Every test specification must contain:

```text
ID / name
purpose
production capability
production entry point(s)
configuration
inputs
independent reference/oracle
units
frame/time convention
deterministic/stochastic controls
procedure
metrics
acceptance criteria
expected PASS result
expected failure modes
required data
data hashes
predecessor gates
artifacts
status semantics
```

Distinguish:
- kernel/formula correctness;
- production wiring;
- convergence;
- conservation;
- observation-operator correctness;
- observational agreement.

## Oracle hierarchy

Prefer:

1. closed-form analytical solution;
2. independently implemented numerical reference;
3. conservation/invariant identity;
4. authoritative external benchmark;
5. reviewed regression baseline only when no first-principles oracle exists.

For stochastic tests record seed, particle counts/weights, and statistical confidence. Do not demand bitwise identity unless guaranteed.

## Build discipline

Use the normal AMPS configuration/build path discovered in the current repo.

- Rebuild generated trees when required.
- Never patch generated copies as the primary fix.
- Record compiler, MPI/OpenMP, configuration macros, source revision, and build identity.
- Preserve build logs for linked runs.
- Do not assume an old AMPS API.

## Test outputs

Use the existing source-tree convention, preferably under something like:

```text
<AMP_ROOT>/test_output/srcMoon/<test-or-campaign>/
```

Do not write test artifacts into immutable data directories.

Every linked I-series run must preserve:

```text
stdout/stderr
result.json
run_manifest.json
input/config snapshot
source revision
build/executable identity
external-data hashes
random seed
MPI/OpenMP layout
science metrics
```

Do not overwrite prior scientifically meaningful results.

## External validation data

Canonical root:

```text
/data/vtenishe/moon_validation_data
```

Before a data-dependent test:

1. identify required Dxx package(s);
2. read dataset `README.md`;
3. read `provenance.json`;
4. read `qa_report.json`;
5. verify hashes;
6. verify READY/PARTIAL/INCOMPLETE status;
7. verify selection/calibration-holdout state;
8. verify units, frames, and times;
9. verify the exact processed file required by the test.

Directory existence does not imply readiness.

## Observational validation gate

Before observational scoring:

1. linked baseline qualified;
2. numerical convergence established;
3. required surface/thermal/chemistry/source capability operational;
4. observation operator verified;
5. data package READY;
6. calibration/hold-out/event selection frozen;
7. calibration-eligible parameters declared;
8. parameters frozen before hold-out scoring.

Never tune on hold-out observations.
Never independently renormalize each observational case.

## Dataset-to-test dependency map

```text
D01  LOLA topography                    -> I14
D02  illumination/PSR                   -> I15, I20
D03  Diviner temperatures               -> I16, I17, I20
D04  PHIDRATES chemistry                -> I21, I22, I31
D05  electron-impact coefficients       -> I23
D06  ARTEMIS plasma forcing             -> I25, I26, I35, I40
D07  Kaguya TVIS Na                     -> I32
D08  LADEE NMS He/Ne/Ar/calibrated      -> I34, I36, I37, I38, I39
D09  LRO/LAMP He                        -> I35
D10  Kaguya GRS K                       -> I28, supports I37/I38
D11  IAU meteoroid streams              -> I30, supports I39
D12  Apollo 17 LACE                     -> I38
D13  Kaguya PACE pickup ions            -> I40
D14  ground-based Na tail               -> I33
D15  LADEE water events                 -> I39
```

Confirm actual local readiness and all model prerequisites before execution.

## Execution order

```text
source/code audit
→ U01–U26
→ I01–I10 baseline
→ I11–I13 convergence
→ I14–I20 surface/environment
→ I21–I24 chemistry
→ I25–I31 drivers/sources
→ I32–I40 observations
→ I41 cross-species freeze
→ I42–I46 UQ/reproducibility/release
```

Do not jump directly to observations while foundational gates remain unresolved.

### Approved development-order exception — 2026-11-12

The user explicitly authorized M2 implementation and local/supporting
production-path verification before M1 closure. This exception permits work
on I14-I20 capabilities, analytical fixtures, event diagnostics, and data QA;
it does not alter qualification order. M1 remains incomplete until I11-I13
PASS. Formal I14-I20 results still require the same-campaign I13 predecessor.
Development results must identify their limited scope, retain formal test
status as SKIPPED / NOT QUALIFIED where predecessors are unavailable, and
never substitute for convergence or observational validation. No new test ID,
acceptance relaxation, hold-out inspection, or automatic reference freeze is
authorized by this exception.

### Approved M3 development-order extension — 2026-11-17

The user approved extending early development to M3 implementation and
local/supporting production-path verification while M1/M2 remain open. The
initial electron rate model is explicitly configured Maxwellian/Voronov;
electron density and temperature are declared inputs, not inferred measured
forcing. Research and concrete proposals for charge-exchange channels and
molecular daughter kinetics precede implementation of those physics choices.
Formal I21-I24 retain their original predecessors, conservation and numerical
criteria. This extension creates no milestone PASS, permits no invented
observations/reaction data and grants no automatic review/reference approval.

## Failure handling

On failure:

1. preserve outputs;
2. classify FAIL vs ERROR correctly;
3. identify first failing invariant;
4. determine whether cause is data, infrastructure, numerics, production wiring, or physics;
5. do not weaken gates;
6. do not edit observations;
7. do not add empirical corrections silently;
8. propose the minimal production fix only after mechanism is understood.

A failed validation is evidence, not permission to change the evidence.

## Code quality

For changed code:
- detailed comments;
- explicit units;
- explicit frames/time systems;
- documented failure behavior;
- no hidden constants;
- named configuration;
- preserve baseline modes required by roadmap;
- update README;
- add/update tests.

## Completion report

At the end of every work package report:

- files inspected/changed;
- source revision;
- build/configuration;
- tests designed/run;
- exact commands;
- required data and hashes;
- PASS/FAIL/ERROR/SKIPPED counts;
- metrics;
- unresolved issues;
- missing production capability;
- incomplete datasets;
- roadmap/test-document changes required.

Never summarize a skipped or preflight-only result as validation.
