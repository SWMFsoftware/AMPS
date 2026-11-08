# AMPS core testing design

## Scope and boundaries

The core campaign verifies code under `src`. Its purpose is to localize defects
before an application-level regression fails. It preserves four distinct
levels:

| Level | Meaning | Expected dependencies | Typical duration |
|---|---|---|---|
| unit | One algorithm, data structure, parser, or numerical kernel | Standard library or the smallest owning library | milliseconds–seconds |
| module | Several core objects through a public module contract | Configured compile is allowed; no physical application campaign | seconds–minutes |
| integration | Two or more core modules, restart, MPI ownership, or generated configuration | Native AMPS build and 1/2/4 ranks | minutes |
| system | A short application calculation used only to verify core plumbing | Native application, retained reference data | minutes–hours |

Application validation remains in `srcEarth/test`, `srcMoon/test`,
`srcEnceladus/test`, `srcSEP/test`, and `srcSEP3D/test`. A core pass does not
qualify a physical model, observational comparison, or scientific campaign.
The global runner may orchestrate an explicitly allow-listed application root
and translate its exact local outcomes. This supplies common scheduling and
evidence without moving the scientific oracle or promoting a local SKIP/PASS
beyond the application suite's documented scope.

## Frozen T1 receipt and live T2 registry

The bounded T1 registry is frozen at 11 suites and 98 stable individual cases,
including a 31-case configured native suite.  It is qualified by retained
normal, reverse, and seed-1729 shuffled campaigns: each has 93 PASS, five
typed provider SKIPs, and no FAIL/ERROR.  Its retained native executable has
SHA-256
`17b98e5ecd22915fc181cce51b8d71bd6abef517802a456ff2efba659c581919`.
Five external provider contracts in the interface suite are
deliberately `SKIP` unless their actual provider libraries,
configured coupling state, and independently versioned data fixtures are
available. Those skips remain visible limitations, not passes.

The live T2 manifest contains the same 11 suites and 103 cases.  The
dependency-light `smoke` campaign remains eight suites and 51 cases;
`nightly` and `--all` select all 103 cases; `distributed` owns the current
five-case two-rank channel/RNG evidence; `profile` selects the 36-case
configured native suite; and `extended` selects that native suite together
with the 11-case interface suite.  The first T2 slice adds five host
species-global particle weight/time-step cases.  Its retained normal evidence
at `/tmp/amps-t2-all-qualified-normal` has 98 PASS, five provider SKIPs, and
no FAIL/ERROR; focused and collision-regression receipts are in
`/tmp/amps-t2-weight-time-initial` and
`/tmp/amps-t2-weight-time-collision-recheck`.  The complete 36-case native
owner also passes reverse and deterministic seed-1729 shuffled order at
`/tmp/amps-t2-native-qualified-reverse` and
`/tmp/amps-t2-native-qualified-shuffle-1729`; its corrected normal, reverse,
and shuffled executable hash is
`6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`.
This evidence does not qualify species-local or single-global profiles,
individual correction and split/merge, unequal-weight/time collision physics,
or CUDA device diagnostics.

## Test oracle rules

- Exact algebraic identities and invariants are preferred over snapshot files.
- Numerical kernels use independently evaluated analytical solutions,
  manufactured solutions, or externally versioned fixtures with provenance.
- Stochastic tests fix and report seeds, test moments/confidence intervals, and
  choose sample sizes from a stated false-failure probability.
- Conservation tests separately report number, mass, charge, momentum, energy,
  and statistical-weight residuals as applicable.
- Convergence tests use at least three resolutions and test observed order,
  rather than comparing one discretization with a loose tolerance.
- MPI tests validate owner values and received ghosts on every participating
  rank. A rank-zero-only output comparison is insufficient.
- Invalid configurations must fail with the expected typed diagnostic; a crash
  or silent fallback is not a pass.
- Every reference records units, frame/sign convention, generating source,
  version, and tolerance rationale.

## Registry and execution contracts

`manifest.json` is the only execution registry. `component_inventory.json` is
the coverage ledger. The global runner validates both before executing a test,
so a renamed/deleted module, missing runner, duplicate ID, or stale coverage
reference cannot silently reduce a nightly.

`source_traceability.json` is the T0 rule-based source/public-contract ledger,
and `configuration_contracts.json` is the T0 machine/generated-mode ledger.
`validate_traceability.py` expands the former across the current source tree
and cross-checks both ledgers against the manifest, component inventory, stable
case registries, and model-local registries. Longest-path ownership makes, for
example, `src/pic/ecsim` more specific than `src/pic`. Tests/build trees and
declared local test registries are excluded from production-file discovery.
The validator proves classification integrity only; it cannot turn a compiled
or named source into behavioral coverage.

`t1_public_algorithm_inventory.json` is the complementary T1 family ledger.
It enumerates public arrays, numerical helpers, RNG distributions, MPI
channels, logging, constants, unit normalization, species radiation,
generated species, ABI/Geopack interfaces, and optional providers. A family
may be `tested`, `partial`, `unavailable`, or `planned`; the latter three
states must survive even when every executable case in an aggregate run is
green.

An individual runner receives a schema-1 context containing resolved
repository/work/artifact paths, tool argv lists, execution settings, and its
manifest entry. It must atomically write a schema-1 result with matching ID,
status, message, cases, evidence, and optional metrics. Exit codes are fixed:
`0` for PASS/SKIP, `1` for FAIL, and `2` for ERROR. The global runner enforces
that correspondence.

`--full-validation` is a non-narrowable execution mode over the complete
manifest. It differs from `--all` by propagating qualification intent in the
schema-1 context, allowing reviewed application adapters to select their long
or qualification paths. Selection and final reports retain the flag and full
campaign union; missing site resources keep their typed status.

Before invoking an individual case, its executor emits and flushes
`[CASE] RUN  <stable-id> <description>`. Suite adapters stream child output,
and the global runner relays these records while preserving the complete log.
Case-start output is informational rather than evidence: the exact result
registry and atomic JSON receipt remain authoritative if a process terminates
between its start record and final result.

The scheduler derives one core stream and one stream per registered
application. `--streams N` bounds concurrently active streams; `--jobs N`
bounds consecutive nonexclusive suites inside each stream. Configured AMPS
builds remain globally exclusive through a writer-priority checkout gate: an
exclusive suite waits for every nonexclusive reader in every stream and blocks
new readers while it is waiting. Results are restored to selected order even
when stream completion order differs. Each stream retains its own text/JSON
receipt and the global JSON/JUnit reports remain the combined authority.

### Result and failure semantics

The four states answer different questions:

| State | Meaning | Runner exit | Nightly aggregate |
|---|---|---:|---|
| PASS | The test executed and every asserted behavior met its oracle | 0 | success |
| FAIL | The test executed, but a software/scientific assertion or required case failed | 1 | failure |
| SKIP | A declared machine capability or data prerequisite was unavailable | 0 | visible, policy-dependent |
| ERROR | The harness could not produce trustworthy test evidence | 2 | failure |

A timeout, malformed receipt, missing receipt, protocol mismatch, duplicate
case ID, hidden failed case, and incompatible exit code are `ERROR`. A compile
failure in a suite that is explicitly testing whether production code builds
is `FAIL`. A compiler that was declared but is not installed is `SKIP`. This
separation prevents infrastructure incidents from being misreported as model
regressions and prevents missing prerequisites from being silently treated as
passes.

The global report records the selected manifest IDs independently of the
results. With `--fail-fast`, unexecuted selected IDs remain evident from that
difference. JUnit contains one test case per individual case when available;
the JSON report retains suite paths, metrics, artifacts, component ownership,
and the full coverage snapshot.

### Isolation and portability contract

- Tests never infer the checkout location from a username, home directory, or
  fixed filesystem prefix. The resolved repository root arrives in context.
- Machine tools are argv vectors. Runners call them without a shell and never
  reconstruct scheduler/compiler commands by string concatenation.
- Relative profile paths are anchored at the repository root; environment and
  home expansion occurs once during configuration loading.
- Dependency-light suites must not modify generated repository configuration.
  Generated headers needed only for compilation are mirrored into `work_dir`.
- Native configured suites must hold the repository lock, reject known active
  builds/runs, remove only the validated `repository/build` path, use
  `Config.pl`, and retain the full configure/build/run transcript.
- A suite may read shared references, but generated expectations must not be
  written back into the source or reference tree.

### C++ harness contract

`native/test_harness.h` uses a function-local registry, stable string IDs, and
exception-based assertion boundaries. It intentionally does not catch fatal
signals or memory corruption: abnormal termination makes the suite adapter
report failed/untrustworthy execution rather than converting corruption into
a normal assertion. Duplicate IDs fail during registration.

Standalone tests link their translation units directly. Configured AMPS
applications place objects in `main.a`; each test translation unit therefore
needs an externally referenced link anchor or an equivalent explicit link
mechanism. Suite adapters require the exact ID multiset so link omission is
detected without relying on unspecified cross-translation-unit initialization
order.

## Component plan

The inventory is normative for status. The complete target for each production
area is:

| Component | Required test families before `covered` |
|---|---|
| general | arrays/containers; random streams and reproducibility; distribution moments; special functions; quadrature order; logger and fatal diagnostics; runtime input |
| interface | ABI/header consumers; C/Fortran round trips; interpolation adapters; unavailable-library typed behavior |
| meshAMR | tree construction; neighbor topology; refine/coarsen; load ownership; ghost exchange at 1/2/4 ranks; cut-cell volume/face geometry; interpolation; restart |
| species | compiled tables; mass/charge constants; parser mapping; duplicate/unknown species failures |
| PIC core | particle-buffer/list lifecycle; weights and time steps; movers and convergence/invariants; injection; sampling; boundary conditions; collisions; splitting/merging; restart; chemistry hooks |
| PIC ECSIM | discrete curl/grad/div/laplacian identities; boundary stencils; solver residuals; halo exchange; wave dispersion/convergence |
| PIC gyro | uniform-field orbit; drifts; adiabatic invariants; time-step convergence; invalid state |
| PIC mesh/parallel | block collectors; owner/ghost synchronization; send masks; migration; rank-count invariance |
| PIC units | base/derived conversion tables; round trips; dimensional consistency; unknown unit failures |
| charge exchange | rate lookup/interpolation; event probability; branching; momentum/energy/weight bookkeeping |
| chemical solver | network assembly; exact linear decay; stiff benchmark; positivity; elemental/charge conservation; rejected invalid network |
| dust | size-bin mapping; charging; force limits; source normalization; moments and split/merge |
| electron impact | thresholds; cross-section interpolation; rate integration; products and weights |
| exosphere/surface/sputtering | surface indexing; source shape normalization; flux accounting; adsorption/desorption; reflection; yield/energy/angle laws |
| photolytic reactions | lifetimes; branching; parent removal; daughter kinematics; statistical weights; no-daughter mode |
| SEP/CME libraries | adapt existing authoritative registries without duplicating expected values or overstating physical qualification |

T1 closes bounded portions of `general`, `PIC units`, `species`, and
`interface`.  The first T2 slice additionally closes the host species-global
weight/time-step API subset of `PIC core`; its other compile-time profiles and
particle-physics consumers remain open. It does not make the broad `general`,
`interface`, or `PIC core` inventory entries complete. The legacy
fork/shared-memory logger, CRC/debug streams,
thread/spinlock and CUDA helpers, older grids/distributions, optional provider
numerics, and T01/T05/T96/TA15/TA16/KMAG reference comparisons remain explicit
gaps. `PIC units` and `species` are covered only for the exact public contracts
recorded in the component and T1 inventories.

## Campaigns at the T12 checkpoint

- `smoke`: eight suites and 51 fast dependency-light cases; expected on every
  change and machine bring-up.
- `nightly`: the established 42 suites and 375 deterministic core cases,
  including typed optional-provider states and exclusive configured owners;
  no observational, imported model, or portability-closure campaign runs.
- `distributed`: two suites and 21 cases: the five-case two-rank channel/RNG
  module plus T10's 16-case collective 1/2/4-rank AMR/PIC matrix.
- `profile`: 19 configured/native suites and 233 cases.
- `extended`: 32 suites and 910 cases, adding interface/optional-provider,
  all 637 T11 imported identities, and the 16-case bounded T12 matrix to the
  configured profiles.
- `--all`: all 49 suites and 1,028 cases. It is not expected to be green while
  imported local FAIL/ERROR evidence and the strict-coverage T12 case remain
  red; aggregate execution must preserve those states rather than hide them.

The current scheduler runs nonexclusive suites concurrently within explicit
exclusive barriers. Future distributed sharding may place those batches on
different workers; each shard must publish its own JSON/JUnit receipt, and a
final aggregator must require the exact manifest selection count before
reporting success.

T11's out-of-line `case_registries.json` is identity data, not an expected-
physics database. The manifest loader requires a unique exact case list and
each adapter separately reconciles its complete ordered local/global registry
before selection. This two-sided check prevents either local drift or central
narrowing. Imported SCCM software verification, reduced-front/background,
full-volume deferred BG3D-4, observational, and legacy disposition claims are
kept distinct in IDs, messages and evidence.

T12 is designed as a bounded, machine-profile-aware matrix rather than one
host's tool autodetection script. Required tools use the normal configuration;
optional compiler/MPI/CUDA/leak commands accept argv vectors from a machine
overlay before conventional fallbacks. Every optional backend retains its own
case ID and cannot borrow another backend's PASS. Performance records require
an independent checksum before applying a broad elapsed-time ceiling, keeping
correctness and capacity failures distinct. `T12-COVERAGE-STRICT-001` invokes
the real strict inventory gate and remains FAIL until all coverage conditions
are satisfied; an anticipated failure is never converted into PASS or SKIP.

## Generated configuration contracts

Generated configurations are test inputs and must be identified as precisely
as runtime fixtures. Each entry in `configuration_contracts.json` records the
authority files, profile, selectors, generated outputs, owning suite/cases,
reviewed status, rationale, and open gaps. `tested` requires direct assertions
and no known gap; `partial` identifies a bounded demonstrated profile;
`planned` has no accepted direct case; `exempt` requires a reviewed reason.

The current configured core application is one bounded partial profile. Its
single exclusive suite must select the application with `Config.pl`, never by
editing `.amps.conf`, and retain the generated configuration and executable
identity. The bounded macro checks and invalid-deck probe do not qualify all
generator options. Complete Config.pl/ampsConfig.pl qualification still needs
isolated positive and negative fixtures, independent generated macro/species
oracles beyond the now-registered duplicate-`SpeciesList` and compiled-table
cases, reproducibility checks, and one contract per mutually exclusive profile.
A failure to configure is `FAIL` only when the tool executed and
violated an asserted generator contract; missing tools or an untrustworthy
receipt keep their normal `SKIP`/`ERROR` meanings.

## Migration sequence

1. Framework, portable profiles, inventory audit, input parser, core containers,
   and a configured particle-buffer application.
2. Import the existing `srcUnitTests` cases into dependency-free/native suites;
   remove hard-coded GoogleTest paths and fix failure exit propagation. T1 now
   includes the bounded general, units, species, ABI, and native species
   families recorded in `t1_public_algorithm_inventory.json`; remaining public
   families retain their explicit status there.
3. Adapt existing stencil, linear solver, interpolation, collision, chemistry,
   dust, SWCME, and SEP/CME runners without reimplementing their oracles.
4. Add missing general/species/unit/model tests and invalid-input coverage.
5. Add mesh/PIC integration with 1/2/4-rank owner/ghost evidence, restart, and
   rank-count invariance.
6. Add sanitizers, compiler matrix, deterministic stochastic budgets,
   performance guardrails, sharding, and nightly aggregation.
7. Change every inventory entry to `covered` only with retained evidence, then
   enable `strict-coverage` as a required nightly gate.

## Current native physics boundary

`CORE-PIC-NATIVE-001` implements a bounded one-rank portion of the PIC,
species-generation, and photolytic-reaction plan. It rejects duplicate
top-level species declarations and independently checks compiled constants,
generated ordering/mass/charge/type tables, and typed unknown-name lookup. Its
five T2 weight/time-step cases exercise host species-global round trips,
minimum reduction, scaled species-distribution copies, and typed rejection of
nonpositive, nonfinite, or out-of-range inputs. Its
physics cases cover table/radial rates, sampled branching, configured missing
daughters, exponential scheduling and 3-D parent survival;
two/three-body conservation and two-body angular isotropy; drifting Maxwellian
moments; hard-sphere scattering invariants; populated-cell NTC rate,
conservation, and RT equilibration; classical LB equilibrium/partition/RT
moments and occurrence rate; neutral ballistic Simple/Boris movement; and an
analytic specular spherical-surface event. These tests use independent
algebraic/statistical oracles and are not application snapshots.

This does not make the relevant inventory components complete. In particular, a
zero-force neutral Boris limit does not qualify charged orbit accuracy; the
NTC cases do not qualify the other collision selectors or unequal statistical
weights/time steps; and parent-only photochemistry does not qualify daughter
insertion or weight bookkeeping.  The weight/time-step cases do not qualify
species-local or single-global compile profiles, individual particle-weight
correction and split/merge, unequal-weight/time collision physics, or CUDA
device diagnostics.  The exact remaining items stay enumerated in
`component_inventory.json` and `test/srcCoreUnitTests/README.md`.
