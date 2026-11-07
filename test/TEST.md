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

As of 2026-11-07 the live global registry contains 17 executable suites and
210 registered individual cases. The historical particle-buffer, physics,
species, packed-fixed-layout, and species-global weight/time family IDs are
selection aliases of the combined baseline native suite, not additional
executable suites. Five additional native owners select mutually exclusive T2
generated profiles, while the T3 mesh owner selects three compile-time
dimensional profiles; every such owner performs its required clean rebuild:

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
| `CORE-PIC-NATIVE-001` | 43 | Baseline packed/fixed layout, generated configuration, species/global units, bounded physics, and species-global weight/time |
| `CORE-PIC-SPLIT-MERGE-001` | 13 | Packed individual-weight layout plus conservative split/merge |
| `CORE-PIC-PBUFFER-LAYOUT-ALIGNED-FIXED-001` | 7 | Aligned-16 fixed-weight record layout and lifecycle |
| `CORE-PIC-PBUFFER-LAYOUT-ALIGNED-INDIVIDUAL-001` | 13 | Aligned-16 individual-weight layout plus split/merge |
| `CORE-PIC-WEIGHT-TIME-LOCAL-001` | 5 | Species-local block weight/time storage and rescaling |
| `CORE-PIC-WEIGHT-TIME-SINGLE-GLOBAL-001` | 4 | One weight/time slot shared across species |
| `CORE-MESH-NATIVE-001` | 58 | One-rank 1-D/2-D/3-D AMR topology/adaptation, 3-D interpolation, cut-cell geometry, and NASTRAN I/O |

The component ledger contains 24 components. The harness entry remains explicitly
`partial`; registered or passing cases are not evidence of complete AMPS
coverage. T0 remains qualified by its retained 45-case evidence. The bounded
historical T1 registry is also qualified: smoke passes 51/51, and the final
98-case normal/reverse/seed-1729 shuffled campaigns each report 93 PASS, five
explicit provider SKIPs, and no FAIL/ERROR. T2 is qualified for its declared
one-rank host profile matrix by the 43-case baseline receipt and the 85-case
six-profile reverse and shuffled campaigns. T3 is qualified for its declared
one-rank Cartesian scope by corrected-tree normal, reverse, and seed-1729
shuffled 58-case campaigns. This does not retroactively change an earlier
stage's scope.

The following limitations must remain visible until their stages close:

- `smoke` selects eight dependency-light suites (51 cases), `nightly` and
  `--all` select all 17 current suites (210 cases), `profile` selects seven
  exclusive native owners (143 cases), `extended` selects the baseline,
  packed-individual, interface, and mesh owners (125 cases), and `distributed`
  selects only the five-case MPI channel/RNG family; the harness must assert
  these live registry sets;
- all current native physics cases are one-rank tests;
- the current Larsen--Borgnakke RT occurrence test qualifies the implemented
  hard-coded `Zrot=5`, not sensitivity to configured `RotationZnumber`;
- the production VT path declares itself unverified and must not be reported
  as covered;
- current parent photodecay does not qualify the full configured chemistry
  processor, daughter insertion, or statistical-weight bookkeeping;
- current sphere testing uses a test callback and does not qualify every
  production surface-interaction callback;
- CUDA/device invalid-state diagnostics, semantic reduction of opaque
  application extension bytes, restart, distributed ownership, ECSIM, gyro
  motion, and most model libraries remain open.

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

Status (2026-11-07): **QUALIFIED FOR THE DECLARED ONE-RANK HOST CONTRACTS**.
The six executable owners and 85 profile cases are:

| Suite/profile | Cases | Qualified contract |
|---|---:|---|
| `CORE-PIC-NATIVE-001` | 43 | Packed fixed-weight production layout and species-global weight/time |
| `CORE-PIC-SPLIT-MERGE-001` | 13 | Packed individual-weight layout and split/merge |
| `CORE-PIC-PBUFFER-LAYOUT-ALIGNED-FIXED-001` | 7 | Aligned-16 fixed layout with anonymous 1/2/4/8/16-byte fields |
| `CORE-PIC-PBUFFER-LAYOUT-ALIGNED-INDIVIDUAL-001` | 13 | Aligned-16 individual layout and split/merge |
| `CORE-PIC-WEIGHT-TIME-LOCAL-001` | 5 | Species-local block storage and physical-distribution rescaling |
| `CORE-PIC-WEIGHT-TIME-SINGLE-GLOBAL-001` | 4 | Single-global slots and cross-species aliasing |

The record tests check typed aligned allocation, offsets and stride, clone,
free-list reuse with stale-payload clearing, intrusive-list integrity, exact
capacity, invalid host access, and real buffer exhaustion in an isolated
nonzero subprocess that must report `The particle buffer is full`. Individual
weights and every public host weight/time setter reject nonpositive or
nonfinite values. Split/merge independently checks model count, total physical
weight, species, weighted position, momentum, translational energy, rotational
energy, each vibrational mode, boundary/target-unreachable behavior, and list
integrity. Anonymous extension bytes follow the documented clone/survivor
policy because their semantics are unavailable to the generic kernel.

The retained qualification commands and receipts are:

```sh
python3 src/test/run_tests.py --test CORE-PIC-NATIVE-001 \
  --output-dir /tmp/amps-t2-native-full-normal
python3 src/test/run_tests.py --campaign profile --order reverse \
  --output-dir /tmp/amps-t2-profile-reverse
python3 src/test/run_tests.py --campaign profile --order shuffle \
  --shuffle-seed 1729 --output-dir /tmp/amps-t2-profile-shuffle-1729
```

The baseline passes 43/43. Both full profile campaigns pass 85/85 across all
six clean native rebuilds. Focused passing receipts are
`/tmp/amps-t2-layout-packed-fixed-recheck` (7/7),
`/tmp/amps-t2-split-merge-packed-recheck` (13/13),
`/tmp/amps-t2-layout-aligned-fixed-initial` (7/7),
`/tmp/amps-t2-layout-aligned-individual-recheck` (13/13),
`/tmp/amps-t2-weight-local-recheck` (5/5), and
`/tmp/amps-t2-weight-single-global-initial` (4/4).

Initial failures remain part of the evidence. The packed-fixed build initially
used `int` test offsets where the checked API requires `long int&`; the packed
split/merge fixture had an unreachable two-particle split bin and a merge
velocity above the debugger guard; the aligned-individual oracle incorrectly
required one particular source tuple to survive; and the local profile omitted
the target-free exosphere block needed to generate internal-boundary base
types. The corrections changed types or independent fixtures/profile input,
not the asserted production invariants. Earlier species-global sentinel
restoration evidence remains documented in `src/test/STATUS.md`.

Residual boundaries are explicit: device-side/CUDA invalid-state diagnostics,
semantic averaging of opaque application fields, unequal-weight collision
physics, injection/sampling normalization, restart, and distributed ownership
are not qualified by T2. T3 qualification below is independent evidence and
does not broaden this T2 contract.

### Stage T3 -- One-rank AMR topology, geometry, and interpolation

Status (2026-11-07): **QUALIFIED FOR THE DECLARED ONE-RANK CARTESIAN
CONTRACTS**. `CORE-MESH-NATIVE-001` is the single executable owner; the
roadmap family IDs `CORE-MESH-TOPOLOGY-001`, `CORE-MESH-ADAPTATION-001`,
`CORE-MESH-CUTCELL-001`, and `CORE-PIC-INTERPOLATION-001` remain selection
aliases rather than duplicate build owners. Its 58 cases comprise 14 topology
and adaptation cases in each of `DIM=1`, `DIM=2`, and `DIM=3`, plus eight 3-D
production-interpolation cases and eight 3-D cut-cell/NASTRAN cases.

The topology and adaptation fixtures use manufactured dyadic trees and
independent integer, volume, and data-transfer oracles. They check
exact-boundary search, parent/child relations, same-level and 2:1
face/edge/corner neighbors, canonical checksums, refine/coarsen callbacks,
associated-data conservation, modification counters, unique leaf load lists,
maximum-level rejection, and repeated split/coarsen cycles. Lower-dimensional
compilation is real evidence: shared AMR/PIC sources were hardened to avoid
3-D-only loop bounds, neighbor access, node-ID types, surface declarations,
mover locals, and field-line scratch sizes in the 1-D and 2-D profiles. The
post-campaign audit additionally bounded global-domain initialization and
`buildMesh` reads/writes to active dimensions while retaining zeroed 3-vector
diagnostics. Each dimensional uniform case now runs production `buildMesh` on
a private topology-only mesh and checks the mandatory root split, generation,
maximum level, and unique `2^D` leaf list.

The interpolation cases execute production row/stencil code on uniform and
2:1 AMR meshes, including all 256 legacy octant-refinement masks. They require
unit-sum, nonnegative leaf support, affine-coordinate reproduction, duplicate
row merging, insertion-order invariance, cache invalidation after adjacent
refinement, and immutable boundary coordinates. The fixes select the `<=`
level path when a coarse cell is surrounded by fine support, inspect all
face/edge/corner neighbor slots before choosing the coarsest candidate, and
correct the edge-1/edge-3 geometric predicates. These are production selection
repairs, not relaxed test tolerances.

Cut-cell checks use independent triangle, ray, prism, and sphere equations.
Production now distinguishes contact with the open block interior from mere
boundary contact, preventing an aligned interface from being counted as a
volume cut. Recursive integration depth follows the interface-dominated
`O(A*h)` error bound rather than treating the error as one isolated subcell
volume. Sphere recursion preserves unsplit coordinate bounds, ray/edge metrics
are finite and translation invariant, and sparse out-of-order NASTRAN IDs are
resolved as keys before deterministic read/write/read comparison.

The final corrected-tree qualification commands and receipts are:

```sh
python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order normal \
  --output-dir /tmp/amps-t3-qualified-normal
python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order reverse \
  --output-dir /tmp/amps-t3-qualified-reverse
python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order shuffle \
  --shuffle-seed 1729 --output-dir /tmp/amps-t3-qualified-shuffle-1729
```

The earlier non-final paths `/tmp/amps-t3-native-normal`,
`/tmp/amps-t3-native-reverse`, and `/tmp/amps-t3-native-shuffle-1729` passed
58/58 before an independent audit applied an additional lower-dimensional
bounds correction. They are retained bring-up evidence, not the final
qualification authority. The three corrected-tree commands above each pass
58/58. Their retained executable SHA-256 values are identical across normal,
reverse, and shuffled order: 1-D
`54e79ddf1c34e2910a059a23e38232bea6322ab7a1ee6d31e35204acd8dbf4e0`, 2-D
`8d5311b99bfba3812599336b2cdbb4fbcd53ab16a3e625122b6db44a30b4aad9`, and
3-D `96f13a4be39ed221b7a252bb8313db3371ed49b2fd464ea1b5d09a545bf42e34`.
Each command performs the required clean
`Config.pl`/top-level build separately for `test/core-mesh-1d`,
`test/core-mesh-2d`, and `test/core-mesh-3d`, then launches every selected case
as a fresh one-rank MPI process. `/tmp/amps-t3-final-fixes-fourth` additionally
passes the four focused interpolation/cut-cell regressions that closed the
preceding 3-D failures.

The final-tree cross-stage regression

```sh
python3 src/test/run_tests.py --campaign profile --order normal \
  --output-dir /tmp/amps-t2-t3-profile-qualified
```

passes 143/143 cases in all seven configured suites with no FAIL, SKIP, or
ERROR. It revalidates the six T2 owners (85 cases) together with the complete
T3 owner (58 cases); it does not expand either stage's declared scope.

The bounded T3 result does not qualify MPI owner/ghost exchange, restart,
non-Cartesian coordinates, arbitrary application-associated data, deep
multi-level imbalance, open/nonmanifold/self-intersecting or large production
surfaces, or multi-rank surface ownership. The legacy `OutputLayerNodeMap`
path remains hard-coded for 3-D and is neither exercised nor qualified in the
1-D/2-D profiles. The three input decks are ignored by
the repository's existing `*.input` rule; Codex did not alter Git control files,
so a human preparing a commit must explicitly force-add
`input/test/core-mesh-{1d,2d,3d}.input`.

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

Current implementation checkpoint (2026-11-07): the registered baseline
profile supplies 19 T4 cases and the separate periodic profile supplies one.
The covered contracts are nonrelativistic uniform-E, uniform-B, and crossed-
field Boris convergence; relativistic uniform-B convergence and gamma
invariance; the public first-order guiding-center RHS for parallel,
`E x B`, grad-B, and curvature terms with typed invalid domains; static sphere
reflection/absorption/tangent/start-inside behavior; exact inertial moving-
plane and moving-sphere events; outer loss; periodic wrapping and associated-
data synchronization; Maxwell half-space flux; exponential waiting times;
uniform cell, triangle, and spherical-volume positions; full raw density,
flux, pressure-tensor, kinetic-energy, speed, temperature moments; legal
weight rescaling; and collect/swap/flush lifecycle. Every new equation,
unit/frame convention, stochastic tolerance, and unsupported domain is
documented next to the implementation or fixture.

This is a substantial bounded T4 qualification, not the universal acceptance
gate above. Still open are explicit coarse/fine mover crossing; full
particle-level gyro/adiabatic convergence; production circle, rotation-body,
Nastran, and diffuse-surface dispatch; multiple-boundary competition;
injection ownership/species-rate normalization; histogram under/overflow;
and MPI migration. The final retained baseline receipts pass 51/51 combined
T4/T5 cases in normal, reverse, and seed-1729 shuffled order at
`/tmp/amps-t45-baseline-final-{normal,reverse,shuffle-1729}`. The periodic
case passes in `/tmp/amps-t45-profile-variants-final`.

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

Current implementation checkpoint (2026-11-07): all five compiled selectors
(`ntc`, `ntc_simplified`, `mf`, `mf_simplified`, and `mf_Yinsi`) use the same
independent unlike/same directed-rate fixture and have zero-relative-speed,
empty/singleton, aggressive collision-limiting, and conservation evidence.
The advanced selectors are also checked with unequal global weights/time
steps; NTC and MF have a separate individual-record-weight profile. Typed
capabilities state that simplified selectors require fixed record weight and
Yinsi individual weights are unimplemented. Collision limiting now caps
trials per model particle and applies an unbiased sample multiplier without
mutating its global threshold from the Yinsi cell kernel.

Separate clean VHS and VSS profiles validate generated diameter, reference-
speed, power-index, and VSS-alpha tables, independent Maxwell quadrature,
the Koura--Matsumoto angular first/second moments, and elastic invariants.
LB tests cover conserved hot-to-cold and cold-to-hot RT partitions, all update-
flag permutations, monatomic/linear/nonlinear rotation, multi-mode low/high-
temperature vibration, and generated `Zrot=2,5,10` sensitivity. Unequal
rotating cross-species `Zrot` remains a typed ambiguity until a physical
mixing law is selected; VT is compiled but explicitly unverified/terminating,
and VV is absent. Those are unsupported capabilities, not test passes.

The baseline T5 portion is included in each 51/51 order-isolation receipt
above. `/tmp/amps-t45-profile-variants-final` additionally passes all seven
selected variant cases across six clean native builds: periodic, VHS, VSS,
Zrot 2, Zrot 10, and two individual-weight collision cases. Within the
declared supported capability domain, the T5 acceptance gate is satisfied;
the typed unsupported domains remain future implementation work.

The final cross-stage command `python3 src/test/run_tests.py --campaign
profile --order normal --output-dir /tmp/amps-t2-t5-profile-qualified` passes
194/194 cases in 12/12 T2-T5 profile suites. Its predecessor retained at
`/tmp/amps-t2-t5-profile-final` found a stale seven-species test oracle after
the charged T4 deck became nine-species; the independent table was extended
to H+ and electron mass/charge contracts, passed its focused 5/5 rerun, and
then passed the complete profile campaign. The diagnostic 193/194 receipt is
not qualification evidence.

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
