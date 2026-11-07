# AMPS core test-routine development roadmap

Status: development roadmap

Last updated: 2026-11-07

Scope: AMPS core code under `src/`; application science validation remains in
the corresponding application test directories.

## 1. Purpose

This document defines the work required to turn the portable runner in
`src/test` into a comprehensive AMPS core verification routine.  The target is
not merely to execute every source file.  Each public algorithm, generated
configuration, numerical mode, ownership transition, and persistent-data path
must be tested against an independent oracle with an explicit validity range.

The core campaign has four levels:

| Level | Purpose | Typical dependencies |
|---|---|---|
| Unit | One data structure, formula, parser, or numerical kernel | Standard library or the smallest owning library |
| Module | Several core objects through a production module interface | A configured AMPS build is allowed |
| Integration | Mesh, particles, physics, restart, or MPI ownership together | Native AMPS build and one or more MPI ranks |
| System smoke | A short application used only to test core plumbing | Native application and small retained input/reference data |

An AMPS core pass does not qualify an application-specific physical model,
observational fit, or long scientific campaign.  Those remain the
responsibility of `srcEarth/test`, `srcMoon/test`, `srcEnceladus/test`,
`srcSEP/test`, `srcSEP3D/test`, and other application-level suites.

## 2. Current baseline

As of 2026-11-07 the live global registry contains 11 executable suites and
103 registered individual cases. The three historical native family IDs and
the new T2 weight/time family ID are selection aliases of the one combined
native suite, not additional executable suites:

| Registered suite | Cases | Current scope |
|---|---:|---|
| `CORE-HARNESS-001` | 11 | Configuration, schema, registry, selection, protocol, subprocess failure, timeout/interrupt cleanup, concurrency, final report, traceability |
| `CORE-GENERAL-INPUT-001` | 4 | Input grammar, diagnostics, includes, canonicalization |
| `CORE-GENERAL-CONTAINERS-001` | 6 | Arrays, bit arrays, stacks, linked containers, strings, 4-D indexer |
| `CORE-GENERAL-NUMERICS-001` | 10 | Arrays 1-D through 5-D, managed stack, special-function identities, quadrature exactness/convergence/domains |
| `CORE-GENERAL-RNG-001` | 7 | Reproducibility, endpoints, independent streams, uniform/normal/gamma moments, isotropic directions |
| `CORE-GENERAL-MPI-001` | 5 | Two-rank scalar/vector/structured/zero-length channel traffic and rank-derived RNG streams |
| `CORE-GENERAL-LOGGER-001` | 4 | Formatting, long messages, concurrent lines, rotation/flush, fatal subprocess behavior |
| `CORE-UNITS-001` | 5 | Independently derived factors, scalar/vector round trips, dimensional identities, invalid scales |
| `CORE-SPECIES-RADIATION-001` | 4 | He/Na reference points, interpolation/scaling, acceleration units, invalid heliocentric state |
| `CORE-INTERFACE-ABI-001` | 11 | C/Fortran ABI, Geopack ordering/conventions/no-SPICE, and five typed optional-provider gates |
| `CORE-PIC-NATIVE-001` | 36 | Generated configuration/invalid decks, species tables/lookup, particle buffers, global unit initialization/output, bounded photochemistry, collisions, Maxwellian sampling, internal energy, neutral movers, sphere event, populated 3-D cells, and the host species-global weight/time API |

The component ledger contains 24 components: three are marked `covered`, ten
are `partial`, and eleven are `planned`. The harness entry remains explicitly
`partial`; registered or passing cases are not evidence of complete AMPS
coverage. T0 remains qualified by its retained 45-case evidence. The bounded
historical T1 registry is also qualified: smoke passes 51/51, and the final
98-case normal/reverse/seed-1729 shuffled campaigns each report 93 PASS, five
explicit provider SKIPs, and no FAIL/ERROR. The live registry adds five T2
cases; this does not retroactively change the T1 qualification scope.

The following limitations must remain visible until their stages close:

- `smoke` selects eight dependency-light suites (51 cases), `nightly` selects
  all 11 current suites (103 cases), `profile` selects the 36-case exclusive
  native profile, `extended` selects its 36 cases plus the 11-case interface
  coverage, and `distributed` selects only the five-case MPI channel/RNG
  family; the harness must assert these live registry sets;
- all current native physics cases are one-rank tests;
- the current Larsen--Borgnakke RT occurrence test qualifies the implemented
  hard-coded `Zrot=5`, not sensitivity to configured `RotationZnumber`;
- the production VT path declares itself unverified and must not be reported
  as covered;
- current parent photodecay does not qualify the full configured chemistry
  processor, daughter insertion, or statistical-weight bookkeeping;
- current sphere testing uses a test callback and does not qualify every
  production surface-interaction callback;
- split/merge, individual particle weights, restart, distributed ownership,
  ECSIM, gyro motion, and most model libraries remain open.

The authoritative live registry and coverage state remain
`src/test/manifest.json` and `src/test/component_inventory.json`.  This
roadmap describes how to close that ledger; it does not replace it.

## 3. Non-negotiable test contracts

### 3.1 Independent oracles

Every numerical or physical test must use one or more of the following:

1. an exact algebraic identity or conservation law;
2. an independently evaluated analytical solution;
3. a manufactured solution with a known truncation error;
4. an independently implemented quadrature or high-accuracy reference;
5. a versioned external fixture with documented provenance.

Expected values must not be generated by the production routine under test.
A snapshot produced by the same code is useful for change detection but is not
an independent physics oracle.

### 3.2 Numerical and stochastic evidence

- Conservation tests report number, physical weight, mass, charge, momentum,
  translational energy, internal energy, and total energy whenever applicable.
- Convergence tests use at least three resolutions or time steps and assert
  observed order, not a single loose-error threshold.
- Statistical tests record a deterministic case-local seed, sample count,
  analytic expectation, variance/confidence model, and false-failure budget.
- Multiple stochastic cases must not share hidden accumulators or depend on
  execution order.
- Tests cover nominal, limiting, empty, degenerate, invalid, and deliberately
  unsupported inputs.  Unsupported physics must return a typed rejection; it
  must not silently fall back to a different model.

### 3.3 Process and state isolation

The runner makes each case independently runnable by stable ID through exact
case selection.  Each case
must save and restore every process-global value it changes, including
particle-buffer pointers, species weights, time steps, random state, sampling buffers,
reaction lifetimes, callbacks, mesh pointers, and optional-data offsets.

Fatal-error paths run in child processes so they cannot terminate the suite
runner.  Timeout handling must terminate the complete child process group and
publish `ERROR`, not a partial or stale pass receipt.

### 3.4 MPI evidence

MPI tests inspect owner data and received ghost data on every participating
rank.  A rank-zero-only comparison is not distributed evidence.  Global
observables must be compared across one-, two-, and four-rank decompositions,
and particle/block migration must be counted exactly once.

### 3.5 Result protocol

Every suite must publish atomic schema-versioned JSON through the existing
runner protocol.  Stable individual IDs, descriptions, PASS/FAIL/SKIP/ERROR
states, metrics, evidence paths, and failure messages are required. Expected
case IDs now live in the global manifest and are cross-checked against runner
registries and result receipts, so a faulty adapter, missing object,
registration, or linker archive member cannot silently reduce coverage.

### 3.6 Portable machine configuration

No test may contain a user name, home directory, fixed installation prefix,
MPI launcher string, compiler path, or scratch path.  Those settings belong in
layered machine profiles such as `src/test/config/machine.example.json`.

Profiles must support:

- compiler and MPI command argument vectors;
- MPI rank flag and rank count;
- output, temporary, and reference-data roots;
- timeouts and artifact retention;
- capability declarations for optional Fortran, CUDA, SPICE, external field
  libraries, sanitizers, and schedulers.

The configuration schema now includes versioned capability declarations in
addition to `paths`, `tools`, and `execution`; harness tests reject unknown
schema versions and inconsistent required/optional capability declarations.

A missing optional capability is `SKIP`.  A machine designated for a required
nightly capability must treat that skip as a policy failure in the nightly
aggregator.

### 3.7 Native AMPS build contract

Configured native suites are exclusive and must use the established build
workflow.  Before every native rebuild the runner must:

1. confirm that the working directory is the resolved repository root (for the
   present workspace, `/home/vtenishe/Mars2/AMPS`) and reject an active build
   or test process for that checkout; the runner must discover this root and
   must not encode the present absolute path;
2. remove only the validated generated directory with `rm -rf -- build`;
3. select the application through
   `./Config.pl -application=test/core-unit-tests` and apply required options
   through `Config.pl`;
4. build through the top-level workflow with
   `env MAKEFLAGS="-j16" make amps`.

The test routine must never edit `.amps.conf` directly.  It must retain the
complete configure/build/run transcript and a copy of the executable used for
the case.  The checkout lock must span configuration, compilation, and native
execution.

## 4. Target runner organization

The completed routine should use the following tiers or campaign selectors:

| Tier | Trigger | Required content |
|---|---|---|
| Smoke | Every change and machine bring-up | Harness, parser, dependency-light general/species/units tests |
| Nightly | Every night | All deterministic one-rank unit and module tests, including native PIC physics |
| Profile matrix | Nightly or path-selected | Individual-weight, split/merge, optional-layout, field-solver, and other mutually exclusive generated configurations |
| Extended | Scheduled or pre-release | Long stochastic tests, convergence studies, restart/corruption, compiler and sanitizer matrix |
| Distributed | Scheduled/weekly and pre-release | One-/two-/four-rank mesh ownership, halos, migration, restart, and rank invariance |
| Application smoke | Scheduled after core passes | Short application plumbing runs without claiming application science validation |

The schema currently accepts `smoke`, `nightly`, `profile`, `extended`,
`distributed`, and `application-smoke`. Campaign names without registered
suites remain empty selections, not evidence that their target coverage has
been implemented.

Nonexclusive suites may run concurrently under `--jobs`.  Every configured
native profile remains behind a leading and trailing exclusive barrier.  Cases
that share one generated configuration are compiled once and selected from one
executable. The particle-buffer and physics families now belong to
`CORE-PIC-NATIVE-001`; the old suite IDs are selection aliases only. Exact
case selection must preserve those individual identities without introducing
an unsafe cross-suite build cache.

## 5. Ordered development stages

Stages are dependency ordered.  A later stage may be prototyped, but it must
not be marked qualified when an earlier required foundation remains open.

### Stage T0 -- Harden the inherited 38-case/current 45-case baseline — qualified

Objective: remove false-green and order-dependent behavior before expanding
the suite.

Required work:

- add real subprocess self-tests for timeout, interrupt, process-group cleanup,
  stale/partial receipts, malformed JSON, incompatible exit status, fail-fast,
  exclusive barriers, bounded concurrency, and machine-profile failures;
- add a source-to-test traceability ledger below the component inventory: each
  production translation unit and public header contract must map to one or
  more stable case IDs, a model-local authoritative registry, or an explicitly
  reviewed exemption; source-tree drift must reject an unclassified file;
- add a machine-readable contract/mode inventory that maps component, public
  contract, production source, compile profile, runtime mode, owning case IDs,
  status, and reviewed exemption; “all supported modes” in later gates means
  the exact entries in this inventory;

- **Implemented T0-C boundary (2026-11-07):**
  `source_traceability.json` and `configuration_contracts.json` provide those
  machine-readable records, and `validate_traceability.py` expands/classifies
  the current production source and public contracts while rejecting
  structural drift. This closes the classification mechanism only; behavioral
  gaps recorded as `partial`/`planned` and full alternative Config.pl/profile
  qualification remain assigned to their later component stages.
- globally validated expected individual-case IDs, exact `--case` selection,
  and deterministic normal/reverse/shuffle ordering are implemented and have
  retained all-order/every-case isolation evidence;
- extend the manifest/configuration schema with the selected campaign and
  required/optional capability model, including schema-version rejection and
  backward-compatibility tests;
- test `Config.pl`, `ampsConfig.pl`, generated macro headers, invalid decks,
  and mutually exclusive compile-time modes without directly editing
  `.amps.conf`;
- capture and reconcile terminal, JSON, and JUnit individual-case results;
- make every native fixture save and restore weights, time steps, sampling
  buffers, reaction lifetimes, random state, mesh pointers, and callbacks;
- after particle deletion, verify both `Prev` and `Next` links with the
  production particle-list integrity checker;
- run current native cases individually and in several registry orders to
  prove order independence;
- replace contradictory “current” evidence in status documents with one
  chronological evidence table;
- combine compatible native families into one configured build to eliminate
  the duplicate rebuild without merging their result identities;
- copy native executables, generated-configuration identity and cryptographic
  hashes into `artifact_dir`; successful work directories may be deleted and
  therefore cannot serve as retained evidence;
- define a `CoreGlobalStateSnapshot`/guard for mutable native-fixture state.
  It must exactly compare before/after pointers, counters, weights, time steps,
  offsets, callbacks, mesh, reaction and sampling state.  State that cannot be
  restored safely requires process-per-case isolation.

Acceptance gate:

- all 45 currently registered cases (including the inherited 38-case baseline)
  pass in normal, reverse and deterministic
  shuffled order, and every case passes when individually selected;
- injected timeout, crash, missing result, stale result, and hidden case
  failure are all reported with the correct non-PASS state;
- smoke, nightly, and extended selections are distinct and documented;
- no fixture leaves changed process-global state after its scope exits;
- every current production translation unit has an initial tested/partial/
  planned/exempt classification without treating compilation as behavioral
  coverage;
- exact commands exist for case selection and deterministic normal/reverse/
  shuffle execution;
- the `src/test` inventory state agrees with the narrower harness contracts
  actually demonstrated by retained evidence.

T0 acceptance result (2026-11-07): **PASS**. The final normal, reverse, and
seed-1729 shuffled campaigns at `/tmp/amps-t0-all-final-normal`,
`/tmp/amps-t0-all-final-reverse`, and
`/tmp/amps-t0-all-final-complete-shuffle-1729` each report 45/45 individual
cases and 4/4 suites PASS with zero FAIL/SKIP/ERROR. All 21 dependency-light
cases also pass as one-case global selections under
`/tmp/amps-t0-individual-<CASE-ID>`.
Each of the 22 executable native cases is invoked with one exact `--case` in a
fresh MPI process from the campaign's single clean build; the two generated-
configuration cases additionally pass as exact global one-case selections at
`/tmp/amps-t0-individual-CONFIG-GENERATOR-001` and
`/tmp/amps-t0-individual-CONFIG-INVALID-DECK-001`. The three all-order
retained native executables have the same SHA-256
`446e920d95fd2f174c7cfce8fadede469ed2af7905116a37fb608c3b9529484a`.
This gate qualifies the T0 harness baseline only; it does not change any
`partial` or `planned` AMPS component into complete behavioral coverage.

### Stage T1 -- General numerics, species, units, and interfaces — bounded scope qualified

Implemented suites and aliases:

- `CORE-GENERAL-NUMERICS-001`;
- `CORE-GENERAL-RNG-001`;
- `CORE-GENERAL-MPI-001`;
- `CORE-GENERAL-LOGGER-001`;
- `CORE-SPECIES-RADIATION-001`;
- `CORE-SPECIES-001` (selection alias within `CORE-PIC-NATIVE-001`);
- `CORE-UNITS-001`;
- `CORE-INTERFACE-ABI-001`.

Required cases:

- arrays 1-D through 5-D: indexing, bounds policy, resize, deep copy, move,
  alias/ownership transfer, empty arrays, and asymmetric dimensions;
- managed stack and container exhaustion/reuse paths;
- random streams: fixed-seed reproducibility, endpoint exclusion, uniform,
  Gaussian, gamma and other used distribution moments, and independent
  cell/rank streams;
- special functions at reference points, asymptotic limits, invalid domains,
  and independently evaluated identities;
- quadrature polynomial exactness and observed convergence for smooth and
  endpoint-sensitive integrands;
- logger formatting, concurrent writers, long messages, rotation/flush, and
  fatal diagnostics in a subprocess;
- MPI channel scalar/vector/structured payloads and zero-length messages;
- atomic/planetary/species constants, generated mass and charge tables,
  parser-to-compiled index mapping, duplicates and unknown names;
- base and derived unit conversions, SI/normalized round trips, dimensional
  consistency, and invalid units;
- C/C++/Fortran ABI round trips, interpolation argument ordering, units and
  coordinate conventions, and typed missing-library behavior;
- ARMS, BATSRUS, CCMC, Kameleon and SWMF adapter fixtures: component ordering,
  unit/coordinate conversion, time interpolation, malformed inputs and typed
  behavior when the optional provider is unavailable.

Independent oracles: standard mathematical identities, hand-computed small
fixtures, high-precision constants with provenance, exact polynomial
integrals, and ABI echo functions that do not call the wrapper under test.

Acceptance gate: `src/general`, `src/species`, and `src/pic/units` have no
unlisted public algorithms; required interfaces pass on at least one supported
Fortran toolchain and cleanly skip on profiles that declare Fortran absent.

Implementation status (2026-11-07): all listed suites and 53 new T1 case IDs
are registered. `src/test/t1_public_algorithm_inventory.json` explicitly
classifies tested, partial, unavailable, and planned public families. Focused
receipts pass 10 numerics, seven RNG, five MPI, four logger, five units, four
species-radiation, five generated/configured species cases, and two production-
linked global-unit cases. The interface
receipt reports six executable ABI/Geopack passes and five explicit provider
SKIPs; those SKIPs preserve the absence of independent ARMS, BATSRUS, CCMC,
Kameleon, and SWMF fixtures rather than claiming adapter qualification.

Qualification status: **QUALIFIED FOR THE DECLARED BOUNDED T1 CONTRACTS**.
`/tmp/amps-t1-smoke-qualified` reports 51/51 PASS. The complete retained
normal, reverse, and seed-1729 shuffled campaigns under
`/tmp/amps-t1-all-qualified-*` each report 93 PASS, five typed optional-provider
SKIPs, and no FAIL/ERROR across all 98 registered cases. All three clean native
builds retained the identical executable SHA-256
`17b98e5ecd22915fc181cce51b8d71bd6abef517802a456ff2efba659c581919`.
The legacy `cLogger` daemon lifecycle, CRC/debug streams, thread/spinlock/CUDA
helpers, optional provider fixtures, and T01/T05/T96/TA15/TA16/KMAG numerical
references remain explicitly outside the bounded implemented T1 subset.

### Stage T2 -- Particle records, weights, time steps, and split/merge

T2 family IDs (the weight/time ID is implemented as a native-suite alias;
the other two remain proposed):

- `CORE-PIC-PBUFFER-LAYOUT-001`;
- `CORE-PIC-WEIGHT-TIME-001`;
- `CORE-PIC-SPLIT-MERGE-001`.

Implementation status (2026-11-07): the first bounded weight/time slice is
registered as five cases selected through the `CORE-PIC-WEIGHT-TIME-001` alias
of the combined native suite.  The current host profile is species-global for
both statistical weight and time step, with individual weight correction off.
It directly tests per-species/block-visible round trips and minimum reduction,
scaled species-distribution copies, and typed rejection of zero, negative,
nonfinite, and out-of-range inputs.  Production validation was added at the
public setters/copy helpers rather than relying on downstream memory access.

Evidence is test-led and retains the failure that found the fixture defect.
`/tmp/amps-t2-weight-time-initial` passes 5/5 with executable SHA-256
`46678c6054daad2ff9b7e3a1aabce9792d12139515e02f93bfaa0296210f5ffe`.
The first full run at `/tmp/amps-t2-all-normal` reports 96 PASS, two FAIL, and
five SKIP because the test fixture restored AMPS's `-1` uninitialized sentinel
through a production setter that now correctly rejects nonpositive values.
After correcting the fixture to restore the raw saved lifecycle state without
invoking physical-value setter/reducer semantics,
`/tmp/amps-t2-weight-time-collision-recheck` passes the five new cases plus the
two affected populated-cell collision cases (7/7), and
`/tmp/amps-t2-all-qualified-normal` reports 98 PASS, five typed provider SKIPs,
and zero FAIL/ERROR across the live 103-case registry.  The entire 36-case
native owner also passes reverse and deterministic seed-1729 shuffled order in
`/tmp/amps-t2-native-qualified-reverse` and
`/tmp/amps-t2-native-qualified-shuffle-1729`; all corrected-order executables
have SHA-256
`6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`.
This qualifies the first host-profile slice and its native order isolation;
full-registry reverse/shuffle, alternate
local/single-global/individual-weight profiles, particle-layout coverage, and
split/merge remain open. `split_merge_t2.cpp` is not registered and carries no
qualification evidence.

Required generated profiles:

1. fixed species weight with the standard optional layout;
2. individual particle weight enabled;
3. anonymous optional fields with mixed sizes and alignments;
4. every optional field needed by IDF, reactions, and sampling together.

Required cases:

- clone and initialize every enabled field, verify offsets/alignment, free-list
  reuse, stale-state clearing, list insertion/removal, and empty/full buffer;
- verify buffer exhaustion and invalid record access in a subprocess;
- global/species/local weights and time steps, round trips, rescaling
  invariance, and zero/negative/nonfinite rejection;
- port all fixed- and variable-weight cases from
  `test/srcUnitTests/split_merge.cpp`;
- after split/merge, assert particle count, total physical weight, species,
  weighted position, momentum, translational/internal energy, optional-data
  policy, and full list integrity;
- cover empty, one-particle, two-particle, target-unreachable, and unrelated
  species cases.

Acceptance gate: all supported record layouts compile and pass; physical
moments are invariant under an allowed split/merge cycle; fatal paths cannot
terminate the global runner.

### Stage T3 -- One-rank AMR topology, geometry, and interpolation

Proposed suites:

- `CORE-MESH-TOPOLOGY-001`;
- `CORE-MESH-ADAPTATION-001`;
- `CORE-MESH-CUTCELL-001`;
- `CORE-PIC-INTERPOLATION-001`;
- adapter for `test/srcInterpolateAMRTest`.

Required cases:

- manufactured 1-D, 2-D, and 3-D uniform/nonuniform trees;
- exact parent/child relations and face/edge/corner neighbors;
- point search at interiors, faces, edges, corners and domain limits;
- refine/coarsen cycles, associated-data transfer, total-volume conservation,
  modification counters, load-list uniqueness, and maximum-level failure;
- analytic plane and sphere cut-cell volume/area fractions;
- surface normals, orientation, triangle intersections, near-tangent and
  degenerate geometry, and NASTRAN read/write fixtures;
- interpolation partition of unity, nonnegative weights where required,
  exact constant/affine reproduction, refinement-interface behavior,
  duplicate-cell merging and permutation invariance.

Acceptance gate: deterministic mesh checksums and conserved quantities agree
after refine/coarsen cycles; geometric error converges at its documented order;
no assertion relies on MPI rank zero or application-specific mesh code.

### Stage T4 -- Movers, boundaries, injection, and sampling

Proposed suites:

- `CORE-PIC-MOVER-CHARGED-001`;
- `CORE-PIC-MOVER-RELATIVISTIC-001`;
- `CORE-PIC-GYRO-001`;
- `CORE-PIC-BOUNDARY-001`;
- `CORE-PIC-INJECTION-001`;
- `CORE-PIC-SAMPLING-001`;
- adapter for `test/srcPeriodicBCTest`.

Mover cases:

- uniform electric field against
  `v=v0+qEt/m` and `x=x0+v0*t+qE*t^2/(2m)`;
- uniform magnetic field against analytic gyrofrequency, radius, phase,
  parallel motion and constant kinetic energy;
- crossed uniform fields against analytic `E x B` drift;
- relativistic Boris against constant gamma and relativistic gyrofrequency;
- guiding-center/gyro limits: parallel motion, magnetic moment, analytic
  `E x B`, grad-B and curvature drifts, plus typed rejection outside the
  method's validity domain;
- three or more time-step levels for observed convergence and reversals of
  charge sign, field sign and mass scaling.

Boundary cases:

- same-level and coarse/fine block crossings, periodic wrapping and outer loss;
- production sphere/circle/rotation-body/Nastran intersections;
- specular, diffuse, absorbing, tangent/no-hit and start-inside behavior;
- moving plane/sphere solutions in a stated inertial frame;
- exact event time, hit point, remaining substep, callback count, final cell,
  and particle-list integrity.

Injection and sampling cases:

- half-space Maxwellian flux with zero, inward and outward bulk flow;
- Poisson event counts and uniform face position;
- uniform-volume sources and the spherical radial CDF proportional to `r^3`;
- cell-volume and source-shape normalization, ownership and species weight;
- known weighted populations yielding exact density, velocity, pressure,
  temperature, flux and energy moments;
- distribution histogram normalization and under/overflow policy;
- sampling interval boundaries, collect/swap/flush and no stale carry-over.

Acceptance gate: every production mover selected by a supported configuration
has a convergence/invariant case; every production boundary mode has exact
accounting; sampled and injected physical rates are invariant under legal
weight/time-step rescaling.

### Stage T5 -- Collisions and internal degrees of freedom

Proposed suites:

- `CORE-PIC-COLLISION-ALGORITHMS-001`;
- `CORE-PIC-COLLISION-VHS-VSS-001`;
- `CORE-PIC-IDF-001`.

Required cases:

- apply the same independent `n_i*n_j*<sigma*g>` fixture to `ntc`,
  `ntc_simplified`, `mf_simplified`, `mf`, and `mf_Yinsi`;
- verify the documented same-species factor and per-particle/directed-event
  convention independently for each selector;
- compare VHS/VSS cross sections and angular moments with their analytical
  laws and independent quadrature of `<sigma*g>`;
- cover unequal masses, weights and time steps; update-flag permutations;
  zero relative speed; empty, one-particle and cold cells; majorant and
  collision-limiting branches;
- require pairwise or weighted momentum and total-energy conservation;
- heat/cold translational and rotational populations in both directions,
  predict final equilibrium from conserved energy, and compare the relaxation
  curve with the configured transition probability;
- build `RotationZnumber` profiles such as 2, 5 and 10 and require measured
  rate sensitivity.  Do not change the production law merely to satisfy the
  test; review the physical cross-species combination first;
- qualify mixed monatomic/linear/nonlinear species and multi-mode vibrational
  low/high-temperature limits;
- keep VT/VV as a typed unsupported case until implementation, independent
  conservation, detailed balance, and equilibrium tests exist.

Acceptance gate: every compiled collision selector has rate, conservation,
edge and statistical tests; configured relaxation parameters demonstrably
control the production path; no selector is inferred from NTC evidence.

### Stage T6 -- Reactions and chemistry models

Proposed suites:

- `CORE-PHOTOCHEM-EXECUTION-001`;
- `CORE-CHARGE-EXCHANGE-001`;
- `CORE-ELECTRON-IMPACT-001`;
- `CORE-CHEMICAL-SOLVER-001`.

Photochemistry cases:

- run the full configured chemistry processor on populated cells;
- verify parent survival, both particle-list links and no-daughter mode;
- add a daughter-producing individual-weight profile and assert number, mass,
  charge, momentum, energy and physical-weight budgets;
- test multiple heliocentric distances, branching invariance, three-body
  angular laws, unavailable daughters, disabled reactions, multiple parents,
  invalid lifetime and reaction-network cycles.

Charge-exchange cases:

- Lindsay--Stebbings and Maher--Tinsley values at table points and boundaries;
- interpolation behavior and Maxwellian rate coefficients from independent
  quadrature;
- reaction identity, product mapping and all conservation/weight budgets.

Electron-impact cases:

- exact zero below threshold, table knots/midpoints, rate integration,
  branching, product mapping and conservation.

Chemical-solver cases:

- exact `A -> B` exponential decay;
- closed-form reversible `A <-> B` evolution;
- a documented stiff Robertson benchmark;
- positivity, elemental/charge conservation and time-step convergence;
- invalid, disconnected and singular network rejection.

Acceptance gate: each reaction model closes its physical budgets and agrees
with an independent rate/solution; daughter insertion and statistical-weight
handling are exercised through production hooks, not only isolated kernels.

### Stage T7 -- Surfaces, exosphere, sputtering, and dust

Proposed suites:

- `CORE-SURFACE-001`;
- `CORE-EXOSPHERE-SOURCES-001`;
- `CORE-SPUTTERING-001`;
- `CORE-DUST-001`;
- adapters for `test/srcModelDustTest` and relevant model-local tests.

Required cases:

- surface-element indexing, positive areas and total analytic area;
- adsorption/desorption state and flux accounting;
- source-shape integration and total-rate normalization;
- uniform and cosine-law surface distributions;
- source, impact and sampled-flux closure;
- sputtering yield tables, threshold/energy/angle laws and product weights;
- exosphere source-process limiting cases and surface-temperature/orbital
  inputs where those are core rather than application contracts;
- dust size-bin mapping, charging equilibrium/transients, force limits,
  injection, Mie-table interpolation, sampled moments and split/merge budgets.

Acceptance gate: every enabled surface/source mode integrates to its requested
physical rate and closes its accounting; tabulated models have boundary and
interpolation fixtures with provenance.

### Stage T8 -- ECSIM, deposition, field solvers, and waves

Proposed suites:

- `CORE-ECSIM-OPERATORS-001`;
- `CORE-ECSIM-BOUNDARIES-001`;
- `CORE-PIC-DEPOSITION-001`;
- `CORE-PIC-FIELD-SOLVER-001`;
- `CORE-PIC-AUXILIARY-001`;
- adapters for `test/srcStencilTest`, `test/srcLinearSystemSolver`,
  `test/srcGMRESPeriodicBCTest`, `test/srcFastWave`, and
  `test/srcLightWaveTest`.

Required cases:

- curl, grad-div, Laplacian and curl-curl on constant, affine, polynomial and
  sinusoidal manufactured fields for orders 2, 4, 6 and 8;
- discrete `div(curl)=0` and `curl(grad)=0` identities within their stencil
  domains;
- periodic, Dirichlet, Neumann and conducting boundary stencils where
  supported;
- one-particle analytic charge/current deposition, total deposited charge,
  and discrete continuity during a cell crossing;
- manufactured Poisson/Maxwell solutions, residual histories and convergence;
- vacuum wave dispersion and a documented particle-plus-field energy budget;
- straight and manufactured-curved field-line interpolation, vertex/edge
  exchange, particle splitting on field lines, and conservation during a
  field-line repartition;
- ray tracing and column integration through constant and affine analytic
  fields, including exact path length and boundary clipping;
- stopping-power tables at knots, interpolated points and validity limits;
- background/fluid-field interpolation and update-mode dispatch against small
  committed interface fixtures, with typed behavior when an optional coupling
  library is unavailable;
- runtime-state and user-particle callback dispatch, invocation order and
  exactly-once accounting; deterministic trajectory-record formatting;
- remote halo stencil tests are deferred to T10 but must reuse these same
  manufactured fields.

Acceptance gate: each advertised spatial order achieves its expected observed
order; solver residual and conservation gates are independent of production
diagnostic output.

### Stage T9 -- Restart, data files, and output integrity

Proposed suites:

- `CORE-PIC-RESTART-001`;
- `CORE-PIC-DATAFILE-001`;
- `CORE-PIC-OUTPUT-001`;
- adapter for `test/srcDataReaderTest`.

Required cases:

- inventory every persistent format and record whether it currently carries
  magic, format version, byte order, integer widths and generated-layout
  identity; where absent, add reviewed metadata and a compatibility policy
  before claiming cross-build portability;
- round-trip mesh topology, owner metadata, particles, optional fields,
  weights, time steps, internal energy, chemistry state, boundary state,
  sampling counters and RNG state;
- compare one step after restart with an uninterrupted run;
- empty mesh/blocks/cells and maximum populated records;
- truncated, corrupt, unknown-version and incompatible-layout failures;
- endianness and integer-width handling, atomic publication, interrupted-write
  behavior, and explicitly supported cross-build compatibility;
- background-data interpolation against small committed fixtures;
- Tecplot/data output variable order, units, dimensions and finite values;
- where supported, write with one rank and read with a different decomposition.

Acceptance gate: restart continuation is equivalent within stated exact or
floating-point tolerances; corruption cannot be accepted as a valid state;
every persistent format has an explicit version/compatibility test.

### Stage T10 -- MPI ownership, halos, block and particle migration

Proposed suites:

- `CORE-MESH-MPI-OWNERSHIP-001`;
- `CORE-PIC-MPI-HALO-001`;
- `CORE-PIC-MPI-MIGRATION-001`;
- `CORE-PIC-MPI-RANK-INVARIANCE-001`.

Required matrix: one, two and four ranks for the same manufactured domain.

Required cases:

- deterministic owner maps and unique block load lists;
- Replace/Add halo operations at faces, edges and corners;
- same-level and coarse/fine send masks;
- explicit owner values and received-ghost values on every rank;
- packing/unpacking blocks with particle lists and optional data;
- block redistribution and particle migration exactly once;
- dynamic load balancing and affinity/repartition decisions on deliberately
  imbalanced block/particle costs, including conservation across repartition;
- periodic migration and outer loss at rank boundaries;
- distributed cut-cell consistency;
- restart under the same and a changed decomposition;
- rank-count-invariant global number, weight, mass, momentum, energy, sampled
  rates and checksums where deterministic ordering permits exact comparison.

Acceptance gate: every rank publishes its local assertions; no test relies on
rank-zero summaries alone; the complete rank matrix passes without duplicate
or lost blocks/particles.

### Stage T11 -- Integrate existing model-local registries

Required adapters:

- `src/models/sep_common` tests;
- `src/models/swcme/test` routine and validation profiles;
- `src/models/sep_corona_swcme/test` including the reduced-front tests;
- the authoritative staged `src/models/sep_coronal_cme` registry;
- remaining reusable legacy applications under `test/`.

The adapter must invoke the authoritative local runner or exact local case
registry and import its result; it must not duplicate expected values or
silently narrow the local qualification scope.  Full-volume, reduced-front,
background-only, particle, and application claims remain distinct.

Acceptance gate: every imported case has a stable global identity and evidence
path; the global report preserves local typed failures and limitations; no
local pass is renamed into a broader physical qualification.

### Stage T12 -- Portability, sanitizers, performance, and coverage closure

Required matrix, limited to configurations actually supported by AMPS:

- supported GNU/Clang compiler profiles and required Fortran compiler;
- debug and optimized builds;
- AddressSanitizer and UndefinedBehaviorSanitizer; ThreadSanitizer for selected
  thread-safe dependency-light modules where compatible with MPI;
- leak checking for long-lived native fixtures;
- CPU and supported accelerator backends without equating one backend's pass
  with another's;
- serial and supported hybrid/OpenMP profiles, including thread-count
  invariance, concurrent particle-list operations and deterministic global
  invariants where ordering permits;
- at least two materially different MPI implementations or scheduler profiles
  when available to the project;
- deterministic performance guardrails for mesh construction, particle push,
  collision processing, halo exchange and restart I/O, with performance
  failures reported separately from physics correctness.

Coverage closure requires:

1. every `src/test/component_inventory.json` entry is `covered` or has an
   explicitly reviewed exemption with rationale;
2. source-tree drift detection rejects new core modules without an inventory
   decision;
3. `make -C src/test strict-coverage` passes;
4. smoke, nightly, profile, extended and distributed receipts are retained;
5. required tiers have zero unapproved SKIP, FAIL and ERROR;
6. all test commands, machine profiles, references, known limitations and
   result locations are documented and reproducible.

Only after these conditions are met should strict coverage become a required
nightly gate.

## 6. Legacy-test migration map

The following existing applications are starting points, not automatically
trusted coverage.  Each must be reviewed for independent oracles, failure exit
propagation, hard-coded paths, deterministic state, and exact case reporting.

| Existing application | Target global area |
|---|---|
| `test/srcUnitTests` | Particle buffer, distributions, IDF, collisions, split/merge provenance |
| `test/srcInterpolateAMRTest` | AMR interpolation |
| `test/srcPeriodicBCTest` | Periodic particle/field boundaries |
| `test/srcStencilTest` | ECSIM/discrete operators |
| `test/srcLinearSystemSolver` | Linear solver residual/convergence |
| `test/srcGMRESPeriodicBCTest` | Periodic GMRES behavior |
| `test/srcFastWave` | Fast-wave dispersion/convergence |
| `test/srcLightWaveTest` | Electromagnetic-wave dispersion/convergence |
| `test/srcLogger` | Logger and diagnostics |
| `test/srcDataReaderTest` | External/background data readers |
| `test/srcModelChemistryTest` | Reaction/chemistry integration |
| `test/srcModelCollisionTest` | Collision-model integration |
| `test/srcModelDustTest` | Dust-model integration |

Migration is complete only when the global runner requires the exact migrated
case set and the legacy executable can no longer hide a failed assertion
behind process success.

## 7. Test data and reference policy

Every committed reference must record:

- stable case ID and owning component;
- physical quantity, SI units, coordinate frame and sign convention;
- source publication, analytical derivation, external dataset, or generating
  script version;
- raw-input checksum where applicable;
- interpolation and preprocessing procedure;
- tolerance rationale and expected convergence order;
- validity limits and known exclusions.

Small immutable references belong beside the owning suite.  Large or licensed
data use `data_root` from the machine profile and a committed manifest with
checksums.  Missing optional data produces a visible `SKIP`; required nightly
data must be provisioned by the site profile.

Reference-regeneration scripts may create candidate files in a temporary or
artifact directory.  They must never overwrite an accepted reference during a
normal test run, and generated candidates require review before promotion.

## 8. Evidence retained by every run

The global run directory must contain:

- resolved machine configuration and immutable selection;
- one complete log per suite;
- atomic aggregate JSON and JUnit reports;
- normalized individual-case results and suite summaries;
- generated input/configuration files;
- configure/build transcript, retained executable copy and hash, and generated
  configuration identity for native suites;
- case metrics, seeds, tolerances and conservation/convergence residuals;
- reference/data checksums;
- MPI rank count and per-rank evidence for distributed cases;
- failure work directories when retention policy requests them.

The final terminal report must list every selected case that executed.  With
fail-fast, selected but unexecuted suites must remain visible from the
selection/result difference and must not be mistaken for passing tests.

## 9. Definition of done

### Individual case

- Stable ID and clear description.
- Independent oracle and documented tolerance.
- Nominal and relevant limiting behavior.
- Deterministic isolation or documented stochastic confidence.
- Typed result and retained metrics/evidence.
- Runnable alone and in arbitrary suite order.

### Suite

- Exact expected case registry.
- Portable prerequisites and machine configuration.
- No hidden application path or global-state dependency.
- Full runner protocol, timeout and cleanup behavior.
- README with build/run commands, oracles, limitations and evidence.

### Component

- All public algorithms and supported compile/runtime modes are mapped.
- Unit, module and required integration behavior are covered.
- Invalid and unsupported modes have typed tests.
- MPI/restart/configuration variants are included where the component owns
  those behaviors.
- The component ledger is changed to `covered` only after retained evidence.

### Comprehensive AMPS core routine

- All 24 current inventory components satisfy their component definition of
  done, and new components cannot bypass inventory review.
- One-rank deterministic nightly, generated-profile matrix, extended
  numerical checks, and one-/two-/four-rank distributed campaigns all pass.
- Required sanitizer/compiler profiles pass.
- There are no unreviewed skips or known false-green paths.
- `strict-coverage` is enabled and required.

Until that point, reports must state the exact passing case and component
counts and must not describe the current routine as complete AMPS-core
qualification.
