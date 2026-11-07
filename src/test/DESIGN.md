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

An individual runner receives a schema-1 context containing resolved
repository/work/artifact paths, tool argv lists, execution settings, and its
manifest entry. It must atomically write a schema-1 result with matching ID,
status, message, cases, evidence, and optional metrics. Exit codes are fixed:
`0` for PASS/SKIP, `1` for FAIL, and `2` for ERROR. The global runner enforces
that correspondence.

Before invoking an individual case, its executor emits and flushes
`[CASE] RUN  <stable-id> <description>`. Suite adapters stream child output,
and the global runner relays these records while preserving the complete log.
Case-start output is informational rather than evidence: the exact result
registry and atomic JSON receipt remain authoritative if a process terminates
between its start record and final result.

Configured AMPS builds are exclusive. The implemented `--jobs N` scheduler
executes consecutive nonexclusive suites in a bounded worker pool (default
`N=1`) and places a leading/trailing barrier around every exclusive suite.
Results are restored to manifest order even when completion order differs.
This prevents checkout-wide application configuration from racing while
allowing isolated dependency-light suites to use otherwise idle cores.

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

## Nightly tiers

- `smoke`: harness integrity and fast dependency-light tests; expected on every
  change and machine bring-up.
- `nightly`: deterministic unit/module coverage suitable for routine CI; no
  observational or long campaign runs.
- `extended`: currently the exclusive configured native profile; future MPI
  rank matrices, restart, sanitizer, convergence, and long statistical suites
  belong here as they are implemented.

The current scheduler runs nonexclusive suites concurrently within explicit
exclusive barriers. Future distributed sharding may place those batches on
different workers; each shard must publish its own JSON/JUnit receipt, and a
final aggregator must require the exact manifest selection count before
reporting success.

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
oracles, reproducibility checks, and one contract per mutually exclusive
profile. A failure to configure is `FAIL` only when the tool executed and
violated an asserted generator contract; missing tools or an untrustworthy
receipt keep their normal `SKIP`/`ERROR` meanings.

## Migration sequence

1. Framework, portable profiles, inventory audit, input parser, core containers,
   and a configured particle-buffer application.
2. Import the existing `srcUnitTests` cases into dependency-free/native suites;
   remove hard-coded GoogleTest paths and fix failure exit propagation.
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

`CORE-PIC-NATIVE-001` implements a bounded one-rank portion of the PIC and
photolytic-reaction plan: table/radial rates, sampled branching, configured
missing daughters, exponential scheduling and 3-D parent survival;
two/three-body conservation and two-body angular isotropy; drifting Maxwellian
moments; hard-sphere scattering invariants; populated-cell NTC rate,
conservation, and RT equilibration; classical LB equilibrium/partition/RT
moments and occurrence rate; neutral ballistic Simple/Boris movement; and an
analytic specular spherical-surface event. These tests use independent
algebraic/statistical oracles and are not application snapshots.

This does not make either inventory component complete. In particular, a
zero-force neutral Boris limit does not qualify charged orbit accuracy; the
NTC cases do not qualify the other collision selectors or unequal statistical
weights/time steps; and parent-only photochemistry does not qualify daughter
insertion or weight bookkeeping. The exact remaining items stay enumerated in
`component_inventory.json` and `test/srcCoreUnitTests/README.md`.
