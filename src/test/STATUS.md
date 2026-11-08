# AMPS core-test status

Last updated: 2026-11-07

## T7 checkpoint

T7 is complete for the bounded one-rank spherical scope selected for this
stage. The live registry contains 31 suites and 311 cases; structural
traceability passes with 24 component areas, 209 production translation
units, 186 public contracts, and 23 configuration contracts. The final T7
campaign reports 29 PASS, no FAIL/SKIP/ERROR, and six passing suites. Its 26
dependency-light cases also pass in normal, reverse, and deterministic
seed-707 shuffled order.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| Spherical surface | Flat/index round trips, positive exact curved areas, total `4*pi*R^2`, unit normals, default inventory/sample-buffer semantics, deterministic accommodation and adsorption-minus-desorption exchange | PASS for spherical one-rank contracts; nonspherical surfaces, target chemistry/temperature laws and MPI reduction remain open |
| Exosphere sources | Area-weighted shapes, uniform and default Knudsen-cosine laws, stationary Poisson event times, cumulative process selection, physical/sampled ledgers, and production particle insertion | PASS for the five target-independent compiled spherical modes through native AMR ownership and the mover; target-specific, external-domain, background-plasma, SPICE and MPI modes remain open |
| Sputtering | Rubin/Teolis boundaries and interpolation; product catalogue; H2O/O2 speed/weight laws; selectable angles with Knudsen-cosine default; incident-rate, daughter-weight and number/mass/momentum/energy ledgers | PASS for the shared source layer; native sputtered-particle insertion and exact page-level provenance of inherited numeric tables remain open |
| Dust kernels | Log bins/power laws, SI potential/Lorentz limits, spherical direction/mass, Fink Mie interpolation, weighted moments, and dust-owned individual-weight split/merge within explicit exact-intrinsic size groups | PASS for the deterministic physical and conservation APIs; the two-output merge cannot generally preserve the pressure tensor, higher moments or raw speed moment |
| Dust charging and injection | Independent Horanyi-current fixtures at three potentials and one equilibrium root; transient convergence for four radii; translated spherical production injection, DUST records, native AMR owner/mover, sampled mass rate and signed reservoir closure | PASS for bounded one-rank profiles; multi-rank arbitration/migration/reduction, restart and long stochastic convergence remain open |

Production corrections include bounds-safe Rubin/Teolis interpolation, the
Teolis interval-offset correction, a single shared sputtering product
catalogue with idempotent initialization, dimensionally correct Mie bounds in
log size parameter, dimension-safe dust injection normalization, translated
sphere geometry, initialized particle records and checked AMR ownership, the
removable `p=1` dust-size limit, and shared checked dust bin/force/moment
kernels used by production call sites. The legacy exosphere injector now uses
the shared exact angular inverse CDFs, stationary Poisson scheduler and
cumulative compiled-process selector. New headers and tests document
equations, SI units, coordinate conventions, numerical bounds, ownership and
unsupported domains.

Final retained T7 evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t7-final` | PASS 29/29 cases and 6/6 suites; no FAIL/SKIP/ERROR |
| `/tmp/amps-t7-kernels-current` | PASS 26/26 dependency-light cases in normal order |
| `/tmp/amps-t7-kernels-reverse` | PASS 26/26 dependency-light cases in reverse order |
| `/tmp/amps-t7-kernels-shuffle-707` | PASS 26/26 dependency-light cases in deterministic shuffled order |
| `/tmp/amps-t7-exosphere-native-final` | PASS native five-mode spherical exosphere source, ownership and ledger case after final relative-rate closure hardening |
| `/tmp/amps-t7-dust-charging-native-current` | PASS native signed-current/equilibrium/transient charging case |
| `/tmp/amps-t7-dust-injection-native-current3` | PASS native translated-sphere dust injection and mass-ledger case |
| `/tmp/amps-t7-harness-final` | PASS 11/11 harness schema, campaign, selection, scheduler, protocol and traceability checks |
| `python3 src/test/validate_traceability.py` | PASS: 24 areas, 209 translation units, 186 public contracts, 23 configuration contracts, 31 suites, 311 cases |

The final native exosphere transcript emitted 618 numerical particles in four
source intervals: 122 impact-vaporization, 141 PSD, 112 thermal-desorption,
122 sputtering and 121 vertical-injection records. The accumulated calculated,
surface-sampled and global production ledgers each report
`6.1799999999999988e11 s^-1`; source occurrence counts are observability checks,
not rate estimators. The dust injection receipt emitted one represented grain
in the first half-grain interval and none in the second, closing the signed
reservoir from `-4.1887902047863901e-17 kg` to exactly zero while reporting
`8.3775804095727796e-16 kg s^-1` sampled mass rate.

Earlier dust receipts that reported two capability SKIPs are retained only as
diagnostics. They motivated the dedicated DUST-species native injection
profile and the dust-owned size-group split/merge API; the all-green final
campaign supersedes them. Likewise, the first native injection attempt exposed
radius-record initialization order and the second exposed a floating endpoint
at `exp(log(a_max))`; both defects were corrected before the retained receipt.

Every configured T7 run performed the mandatory root/process preflight,
removed only `build`, selected its application through `Config.pl`, and
compiled through top-level `make amps` with `MAKEFLAGS=-j16`. Focused
warnings-as-errors and ASan/UBSan runs passed the dependency-light dust and
sputtering/exosphere kernels; LeakSanitizer is unavailable under this
environment's ptrace restrictions. No Git state operation was performed.

## T6 checkpoint

T6 is qualified for the implemented, bounded reaction domains. The live
registry at this checkpoint contains 25 suites and 282 cases, and structural
traceability passes with 24 component areas, 209 production translation units,
181 public contracts, and 20 configuration contracts. The gate covers the
dependency-light chemical integrator, configured H2O photolysis execution,
one individual-weight daughter-production profile, Burger electron-impact
rates/products, and the two compiled charge-exchange rate laws.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| Chemical integrator | `A -> B`, reversible `A <-> B`, stiff decay, coupled chain, positivity/conservation, invalid networks | PASS for the dependency-light finite network API; runtime deck assembly and cell coupling remain open |
| Photolysis execution | Populated-cell parent survival/decay, intrusive-list integrity, disabled/invalid lifetime paths, idempotent initialization | PASS for configured H2O quiet-Sun execution; application multi-parent networks remain open |
| Photolysis daughters | Individual-weight H2O channel, daughter identity/count, weight, charge, mass, momentum and kinetic-energy accounting | PASS for the selected one-channel native profile; MPI/threaded insertion and general runtime networks remain open |
| Charge exchange | Lindsay--Stebbings and Maher--Tinsley formulas, finite-domain values and invalid-domain sentinel | PASS for the exposed H/H+ lifetime functions; event selection, products and particle weights are not exposed by this module |
| Electron impact | Burger total rates, thresholds, eligible-channel selection and H2O product statistics | PASS for compiled parent-rate dispatch and H2O products; runtime electron coupling and other-parent product insertion remain open |

The populated-cell fixture now calls
`PIC::DomainBlockDecomposition::UpdateBlockTable()` after allocating owner
blocks. This is a substantive test correction: without rebuilding that table,
whole-mesh chemistry processors observe `nLocalBlocks == 0` and a nominally
passing parent-decay case is a no-op. The final receipts therefore supersede
the earlier diagnostic runs.

Final retained T6 evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t6-chemical-integrated` | PASS 7/7 chemical-solver cases |
| `/tmp/amps-t6-photochem-baseline-blocktable` | PASS 6/6 populated-cell photochemistry cases |
| `/tmp/amps-t6-photochem-daughter-blocktable` | PASS 1/1 individual-weight daughter case |
| `/tmp/amps-t6-reactions-build2` | PASS 7/7 charge-exchange/electron-impact cases |
| `/tmp/amps-t6-gate-reverse` | PASS 21/21 cases and 4/4 suites in reverse order |
| `/tmp/amps-t6-gate-shuffle-706` | PASS 21/21 cases and 4/4 suites in deterministic shuffled order |
| `/tmp/amps-t6-harness-prebuild-2` | PASS 11/11 harness self-tests |
| `python3 src/test/validate_traceability.py` | PASS: 24 areas, 209 translation units, 181 public contracts, 20 configuration contracts, 25 suites, 282 cases |

Each native T6 owner used the mandatory root/process preflight, removed only
`build`, selected its deck through `Config.pl`, and compiled through top-level
`make amps` with `MAKEFLAGS=-j16`. No Git state operation was performed.

## T4/T5 checkpoint

The live registry contains 22 suites and 261 cases. Structural traceability
passes with 24 component areas, 209 production translation units, 181 public
contracts, and 19 configuration contracts. T4 is implemented and qualified
for the bounded one-rank profiles below but remains partial against the full
roadmap acceptance gate; T5 satisfies its supported-capability gate while
retaining typed unsupported domains.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| Charged movers | Uniform-E velocity/position convergence; uniform-B phase/energy convergence; crossed-field drift; relativistic gamma and gyrofrequency | PASS for constant-field Boris profiles; coarse/fine crossing and additional production movers remain open |
| Guiding center | Parallel electric/mirror force, `E x B`, grad-B and curvature drift; neutral, zero-B, negative-mu and non-finite statuses | PASS for the field-local RHS; full orbit integration, adiabatic validity and convergence remain open |
| Boundaries | Static sphere reflection/absorption/tangent/start-inside, outer loss, exact moving plane/sphere, periodic real/ghost wrap and data sync | PASS for registered modes; circle/rotation-body/Nastran mover dispatch, diffuse law, multiple boundaries and MPI migration remain open |
| Injection | Maxwell half-space flux; Poisson waits; uniform cell, triangle, and ball positions | PASS for distribution kernels; ownership, species-rate/weight normalization and distributed injection remain open |
| Sampling | Exact density, flux, full pressure tensor, energy, speed, temperature, rescaling invariance, collect/swap/flush | PASS for raw cell moments; histogram under/overflow and distributed reduction remain open |
| Collision selectors | Same independent rate/conservation fixture for NTC, simplified NTC, MF, simplified MF, and Yinsi; zero/edge/limiting paths; unequal global W/dt on advanced selectors | PASS in declared domains; simplified selectors remain fixed-weight-only and Yinsi record weights are typed unimplemented |
| Molecular laws | Generated VHS/VSS power tables; Maxwell quadrature; VSS alpha moments and elastic invariants | PASS for the two compile profiles; broader transport calibration and alternate mixture laws remain open |
| LB internal energy | RT hot/cold equilibrium, update flags, mode families, vibrational limits, generated Zrot 2/5/10 | PASS for RT; unequal rotating-pair Zrot is typed ambiguous, VT is unverified/terminating, VV is absent |

Production corrections made under these tests include the previously omitted
VHS/VSS power-index emission, VSS angular dispatch for every selector,
unbiased per-particle collision limiting without Yinsi global-threshold
mutation, generated rather than literal `Zrot`, a typed guiding-center RHS,
and shared checked Poisson/uniform-ball source samplers. New code and fixtures
carry equations, SI units, frames, stochastic tolerances, validity limits, and
typed-failure comments.

Final retained evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t45-baseline-final-normal` | PASS 51/51, one clean baseline build |
| `/tmp/amps-t45-baseline-final-reverse` | PASS 51/51, reverse case order |
| `/tmp/amps-t45-baseline-final-shuffle-1729` | PASS 51/51, deterministic shuffled order |
| `/tmp/amps-t45-profile-variants-final` | PASS 7/7 cases and 6/6 suites: periodic, VHS, VSS, Zrot 2, Zrot 10, individual-weight NTC/MF |
| `/tmp/amps-t2-t5-profile-qualified` | PASS 194/194 cases and 12/12 suites on the corrected final tree, covering every T2-T5 configured profile |
| `python3 src/test/validate_traceability.py --output /tmp/amps-t45-traceability-pretest.json` | PASS: 24 areas, 209 translation units, 181 public contracts, 19 configuration contracts, 22 suites, 261 cases |

The native adapter performed the mandatory root/process preflight, removed
only `build`, selected each application with `Config.pl`, and compiled with
top-level `make amps` and `MAKEFLAGS=-j16` for every retained build. No Git
state operation was performed.

The first full-profile diagnostic at `/tmp/amps-t2-t5-profile-final`
reported 193 PASS and one FAIL: the pre-T4 species oracle still asserted a
seven-neutral deck after H+ and electrons were appended for the charged-mover
profile. The oracle now documents and checks all nine symbols, indices,
masses, charges, and gas types. `/tmp/amps-t45-species-nine-profile-fix`
passes the focused 5/5 correction, and the complete
`/tmp/amps-t2-t5-profile-qualified` rerun supersedes the diagnostic with
194/194 PASS. The failed receipt is retained as discovery evidence and is not
reported as qualification.

## Implemented

- Portable global runner with validated registry, isolated work/artifact
  directories, timeouts, typed PASS/FAIL/SKIP/ERROR states, JSON, and JUnit.
- The global runner now ends every execution with a stable
  `FINAL INDIVIDUAL TEST RESULTS` table.  It prints one PASS/FAIL/SKIP/ERROR
  row for each executed case, then separate individual-case and suite
  summaries.  The JSON receipt preserves the same normalized rows in
  `individual_results` and their totals in `individual_counts`; suites which
  cannot publish case evidence are explicitly labeled as suite-level rows.
- Layered machine configuration with external JSON profiles, environment
  overrides, relative-path resolution, and no hard-coded installation paths.
- Explicit component inventory with structural drift audit and a separate
  strict-completeness gate. `source_traceability.json` now classifies every
  discovered production translation unit and public compile contract, while
  `configuration_contracts.json` records the machine, Config.pl, generated-
  file, compile-profile, and runtime-mode boundaries.
- The dependency-free `validate_traceability.py` expands those rules, rejects
  unclassified/overlapping files and stale suite/case/local-registry owners,
  and cross-checks the 24-area component inventory. This is structural
  governance evidence, not behavioral physics coverage.
- Eleven registered harness self-tests cover the original configuration,
  registry, selection, protocol, scheduler, and report contracts plus schema,
  hostile-subprocess failures, real timeout and SIGINT process-tree cleanup,
  bounded concurrency/fail-fast, and the full source/configuration
  traceability contract.
- Generic runtime-input parser adapter (four existing independent fixtures).
- General utility tests for arrays, bit arrays, stack allocation/reuse, linked
  containers, string helpers, and sparse 4-D indexing.
- T1 adds 10 general-numerics cases for asymmetric arrays 1-D through 5-D,
  empty/copy/move/alias/resize/bounds behavior, managed-stack reuse,
  special-function identities, and independent quadrature exactness,
  convergence, and invalid-domain checks.
- T1 adds seven deterministic RNG cases for recurrence/reseeding, endpoint
  exclusion, independent interleaved streams, uniform/normal/gamma moments,
  and isotropic directions, plus a separate two-rank RNG-stream case.
- T1 adds four two-rank `CMPI_channel` payload cases, four diagnostic-logger
  cases, five SI/CGS/normalized-units cases, and four He/Na radiation-table
  cases. The associated production hardening covers safe multidimensional
  array ownership/extent arithmetic, complete hybrid RNG reseeding, checked
  MPI framing/alignment, finite and representable unit scales, Gamma/bulk-
  normal domains, GCD `INT_MIN`, finite radiation-table clamping, and invalid
  or unrepresentable heliocentric inputs.
- T1 adds `CORE-INTERFACE-ABI-001`: six executable C/Fortran/Geopack cases and
  five explicit optional-provider gates. ARMS, BATSRUS, CCMC, Kameleon, and
  SWMF remain `SKIP` until independent libraries and versioned fixtures are
  configured; these skips are not behavioral passes.
- New dependency-free C++ harness shared by standalone and configured tests.
- One combined `CORE-PIC-NATIVE-001` registered owner for the configured
  `test/srcCoreUnitTests` application. It replaces two executable native suite
  entries while preserving their old IDs as selection aliases, so compatible
  particle-buffer and physics cases require one configured build. T1 extends
  it from 24 to 31 cases and adds `CORE-SPECIES-001` as a third alias for
  duplicate-species rejection, compiled constants, generated mapping/tables,
  and unknown-name diagnostics. Two additional production-linked cases cover
  transactional/idempotent global AMU/charge initialization and conversion-
  table output.
- T2 extends the baseline exclusive native owner from 31 to 43 cases and adds
  five profile owners, for six suites and 85 native profile cases. The matrix
  covers packed/aligned-16 and fixed/individual-weight records,
  species-global/species-local/single-global weight and time-step modes, and
  packed/aligned split/merge. Production record requests now enforce typed
  power-of-two alignment and checked offset/stride arithmetic; record reuse
  clears stale payload; invalid host pointers are contained; deletion marks
  the addressed record; and host weight setters require finite positive
  values before accessing storage. Species-local distribution copies rescale
  physical density when weight or time step changes.
- T2 split/merge cases assert model count, physical weight, species, weighted
  position, momentum, translational energy, rotational energy, and each
  vibrational mode. The central-variance implementation avoids catastrophic
  cancellation. Opaque application bytes follow a clone/survivor policy
  because generic PIC code has no field semantics. Buffer exhaustion is
  isolated per MPI subprocess and must exit nonzero with the expected
  diagnostic.
- T3 adds `CORE-MESH-NATIVE-001`, a 58-case owner that performs clean 1-D,
  2-D, and 3-D configured builds. Fourteen topology/adaptation cases execute in
  each dimension; eight production interpolation and eight cut-cell/NASTRAN
  cases execute in 3-D. The former roadmap suite IDs remain selection aliases.
  Production hardening covers lower-dimensional AMR/PIC compilation, safe
  split/coarsen data transfer, full neighbor-slot interpolation selection,
  affine reproduction across all 256 legacy refinement masks, open-interior
  cut-surface contact, interface-error-driven integration depth, translated
  sphere recursion, ray/edge metrics, and sparse NASTRAN identifiers.
- Expanded 19-case native physics suite:
  - four photochemistry kernel cases for Huebner H2O rates/yield/radial
    scaling, exponential decay probability/event time, product momentum/
    energy and two-body angular isotropy, and the seven-channel branching
    distribution plus neutral-only compiled-product availability;
  - two hard-sphere collision-kernel cases for unequal-mass momentum/energy/
    relative-speed conservation and spherical isotropy;
  - one drifting-Maxwellian distribution case for component means, analytic
    mean thermal speed, and recovered kinetic temperature;
  - six Larsen-Borgnakke cases for optional particle storage, rotational
    equipartition, vibrational harmonic-oscillator means, beta energy
    partition, two-stage RT redistribution/conservation, and the analytic
    `1/Z_rot=1/5` relaxation occurrence rate;
  - three mover cases for neutral ballistic Simple/Boris limits, destination-
    cell insertion, and an analytic specular reflection at a registered sphere;
  - three populated 3-D-cell cases: NTC/LB H2O rotational equilibration with
    momentum/total-energy conservation, the independent Maxwellian H2O/O
    collision rate `n0*n1*sigma*<g>`, and 2,400-particle photolytic parent
    survival with no daughter production.
- Native physics statistics use fixed case-local seeds, documented sample
  counts, analytical expectations, and explicit confidence/tolerance bounds.
- The photoproduct angular test exposed and fixed uniform-polar-angle sampling
  in `GetProductVelocity_2`; the implementation now samples `cos(theta)`
  uniformly for the correct solid-angle measure.
- The neutral-only branching test exposed and fixed a second production defect:
  product IDs were compacted around unavailable charged daughters while their
  velocity vectors retained original channel slots. Both ordinary and
  given-product generators now compact IDs and velocities together, with a
  NaN-sentinel regression for a surviving slot-1 daughter.
- Native-run safety includes a repository-specific advisory lock, an active
  build/run preflight, exact validation of the generated `build` target,
  application selection through `Config.pl`, fixed `MAKEFLAGS=-j16`, and a
  retained configure/build/MPI transcript.
- The native lock now uses the shared checkout's read-only `Config.pl` inode,
  so cooperating jobs on different nodes or with different `TMPDIR` values
  cannot reconfigure the same checkout concurrently. Explicit output reuse
  invalidates the old per-suite receipt before launch; an adapter must publish
  fresh evidence. PASS receipts must contain at least one executed case.
- Exact order-independent case registries detect a missing, renamed,
  duplicated, or unlinked C++ test even when the remaining executable returns
  success, without relying on unspecified cross-translation-unit static
  initialization order.
- Exact `--case` selection and deterministic `--order normal|reverse|shuffle`
  with a recorded shuffle seed are implemented. The harness checks selection,
  aliasing, and deterministic order construction. The retained T0
  normal/reverse/shuffle campaigns each pass their historical 45-case
  registry. All 21 T0
  dependency-light cases also pass through one-case global selections; the 22
  executable native cases are launched as exact one-case MPI processes from
  one clean build, while the two configuration cases execute as individually
  identified adapter cases.
- The global runner accepts `--jobs N` (default `1`) and executes only
  consecutive nonexclusive suites concurrently. Every `exclusive: true`
  native suite has a leading and trailing barrier, so repository-wide
  configuration/build state cannot overlap another suite. Results are restored
  to manifest order; concurrent case-start lines carry their suite ID.
- The runner now streams `[CASE] RUN <ID> <description>` immediately before
  every Python, dependency-light C++, parser, and configured native case. The
  terminal shows case-level progress while full compiler/MPI output remains in
  the suite log; streaming retains timeout enforcement and kills the complete
  POSIX child process group on timeout or user interrupt. The harness exercises
  both paths with a real descendant process.
- The legacy `test/srcUnitTests` entry point now returns GoogleTest's result
  instead of reporting process success after failed assertions.
- `SRC_UNIT_TESTS_MIGRATION.md` maps every legacy source/case to its registered
  replacement or an explicit open item. It records the legacy out-of-bounds RT
  assertion, allocator-mismatched fixture teardown, non-reset statistical
  accumulators, and still-unported individual-weight split/merge profile.
- `array_1d` now has explicit deep-copy/resizing-copy and ownership-transferring
  move semantics. The container tests exposed that the previous implicit copy
  constructor aliased owned storage and could double-delete when copy elision
  was unavailable.
- CRC diagnostic formatting no longer uses overlapping source/destination
  `sprintf`; bounded append formatting preserves termination on long paths.
- The manifest distinguishes smoke, nightly, profile, extended, distributed,
  and application-smoke campaigns. The retained T1 qualification registry
  assigned eight suites and 51 cases to smoke and 11 suites/98 cases to
  nightly. At the T3 checkpoint the registry had 17 executable suites and 210
  cases: smoke contained 51, nightly/`--all` contained 210, `profile`
  contained seven exclusive native suites and 143 cases, and `extended`
  contained 125. The current T7 registry instead has 31/311 total and 17/211
  in `profile`; `distributed` still retains its five-case MPI selection. The
  harness pins the current distinct sets. Empty campaign selections do not
  imply implemented coverage; every native profile remains a scheduler
  barrier wherever selected.
- Framework, suite-adapter, native-fixture, and core-change comments document
  protocol/failure semantics, path containment, portability, concurrency,
  memory ownership, MPI reduction, test oracles, and validity limits. The
  framework/design and configured-application READMEs carry the matching
  maintainer and machine-configuration contracts.

## Not yet complete

Only `src/general/input` and `src/pic/units` are currently marked `covered`.
`src/species` has complete computational path coverage for this T1 scope but
remains `partial` until a separately versioned paper/digitization fixture
independently revalidates the retained He/Na tables. The test harness is
`partial`, and all other `partial` and
`planned` entries in
`component_inventory.json` remain open. Existing model/application runners
outside the adapted T7 spherical-exosphere and dust applications still remain
to be integrated;
distributed mesh ownership, restart, and rank-invariance coverage are open.
T2 closes the declared host particle-layout, full-buffer, weight/time-mode,
individual-weight, and split/merge contracts, but not every use of those
facilities by later physics.

The PIC and photolytic-reaction entries remain `partial`. Production NTC cell
selection, a cross-species analytic rate, and H2O RT equilibration now have
one-rank populated-cell evidence. Species-global, species-local, and
single-global weight/time modes, individual particle weights, local rescaling,
and conservative split/merge now have T2 host evidence. T5 covers all five
collision selectors in their declared weight domains, collision limiting, and
the configured VHS/VSS laws; Yinsi individual-record weights and simplified-
selector unequal-weight extensions remain open. T4 covers constant-field
charged movers plus bounded static/moving/periodic surfaces, but coarse/fine
crossings, additional production surface types and rank migration remain open.
T6 covers three-body photolysis angular statistics and one individual-weight
H2O daughter-insertion profile; quiet-Sun/runtime networks and MPI/thread-safe
daughter processing remain open. VT/VV relaxation is not claimed because the
current VT implementation explicitly exits as unverified.

The RT occurrence case deliberately qualifies the current LB closure's
hard-coded `Zrot=5`. The input generator emits and initialization validates
`RotationZnumber`, but `LB::RedistributeEnergy` does not currently use that
array. Selecting a cross-species combination rule is a production-physics
decision, so this test campaign records the discrepancy instead of inventing
one. A later correction must set a non-5 input and verify the configured rate.

`HARNESS-TRACEABILITY-001` imports and runs the same exhaustive validation used
by the standalone `validate_traceability.py` command. It expands and audits
every current production translation unit/public contract and validates the
configuration ledger, so structural drift is part of smoke/nightly execution.
The explicit standalone command remains useful for retaining the expanded
per-file JSON receipt.

T1 is qualified for its declared bounded contracts. The stale campaign oracle
from the first smoke run was corrected; `/tmp/amps-t1-smoke-qualified` passes
51/51. The final 98-case normal, reverse, and seed-1729 shuffled campaigns each
report 93 PASS, five typed optional-provider SKIPs, and no FAIL/ERROR. This
qualification does not promote the explicit `partial`, `planned`, or
`unavailable` contracts described above.

T2 is qualified for the declared one-rank host profile matrix. The baseline
owner passes 43/43, and the six-suite profile campaign passes 85/85 in reverse
and deterministic seed-1729 shuffled orders. The matrix covers four record
layouts, three supported weight/time organizations, individual finite-positive
weight validation, local rescaling, and packed/aligned conservative
split/merge. Earlier species-global sentinel restoration and the initial T2
layout/profile fixture failures remain retained diagnostics; each was
corrected without weakening the production invariant. CUDA/device diagnostics,
semantic reduction of opaque fields, unequal-weight collision physics,
injection/sampling normalization, restart, and distributed ownership remain
outside T2. The separately implemented T3 scope is described below.

T3 is qualified for its declared one-rank Cartesian scope. The first normal,
reverse, and seed-1729 shuffled 58/58 receipts were followed by an independent
audit that found and corrected one additional lower-dimensional array-bound
defect, so those receipts are retained only as superseded bring-up evidence.
The corrected tree then passed all 58 cases in normal, reverse, and seed-1729
shuffled order with identical per-profile executable hashes. The scope does not
include MPI owner/ghost exchange, restart, non-Cartesian coordinates, arbitrary
application-associated data, deep imbalance, or general nonmanifold/open/
self-intersecting and multi-rank surfaces. The legacy 3-D-hard-coded
`OutputLayerNodeMap` is not exercised or qualified in the lower-dimensional
profiles.

The lower-dimensional bring-up was substantive rather than a source-only
compile: dimension-bounded mesh-signature and associated-data loops, compact
neighbor access, the 2-D neighbor call signature, dimension-generic mesh and
node-ID types, mover local scope, and always-three-component field-line vectors
were corrected before 1-D/2-D cases could execute. The final audit also bounded
global-domain initialization and `buildMesh` reads/writes to active dimensions,
with explicitly zeroed three-component diagnostics. The per-dimension uniform
case now executes production `buildMesh` on a private topology-only mesh and
asserts its root split, generation/maximum level, and unique `2^D` bottom list.
The interpolation audit then corrected the coarse-surrounded-by-fine `<=`
selector, searched every neighbor
slot before selecting the coarsest support, and repaired the edge-1/edge-3
predicates. Cut-cell work distinguishes open-interior intersection from face
contact and selects recursive depth from the interface error `O(A*h)`; sphere
subdivision preserves unsplit bounds, and NASTRAN lookup supports sparse,
out-of-order IDs.

`make -C src/test strict-coverage` is expected to fail until those items are
implemented and qualified. A passing smoke run must not be described as
complete AMPS coverage.

## Current validation evidence

The live source tree registers 31 executable suites and 311 individual cases.
Its `profile` campaign selects 17 configured owners and 211 cases, including
the expanded 93-case baseline owner, historical T2 variants, the 58-case T3
dimensional mesh owner, T4/T5 compile variants, T6 reaction profiles, and
three T7 applications. The earlier seven-suite/143-case combined T2+T3
profile remains retained historical evidence rather than the current campaign
membership. The historical bounded T1
qualification remains 11 suites/98 cases with 31 native cases; its normal,
reverse, and shuffled native rebuilds retain the identical executable SHA-256
`17b98e5ecd22915fc181cce51b8d71bd6abef517802a456ff2efba659c581919`.
T2 has a 43/43 normal baseline receipt plus 85/85 reverse and deterministic
shuffled profile receipts. T3 has final corrected-tree 58/58 normal, reverse,
and seed-1729 shuffled receipts; its first three 58/58 receipts remain
superseded bring-up evidence. The final-tree combined profile regression passes
all 143 T2+T3 configured cases. T0's earlier 45-case receipts remain
historical evidence, and the retained T1 receipts remain the authority for T1
itself.

| Date | Command/evidence | Outcome | Scope |
|---|---|---|---|
| 2026-11-07 | `python3 src/test/run_tests.py --campaign profile --order normal --output-dir /tmp/amps-t2-t3-profile-qualified` | PASS 143/143, 7/7 suites, 0 FAIL/SKIP/ERROR (426.294 s) | Historical T2+T3 checkpoint of all six T2 configured owners (85 cases) plus the T3 mesh owner (58 cases); the live profile campaign has since expanded. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order normal --output-dir /tmp/amps-t3-qualified-normal` | PASS 58/58 | Final corrected-tree normal order. Executable SHA-256: 1-D `54e79ddf1c34e2910a059a23e38232bea6322ab7a1ee6d31e35204acd8dbf4e0`; 2-D `8d5311b99bfba3812599336b2cdbb4fbcd53ab16a3e625122b6db44a30b4aad9`; 3-D `96f13a4be39ed221b7a252bb8313db3371ed49b2fd464ea1b5d09a545bf42e34`. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order reverse --output-dir /tmp/amps-t3-qualified-reverse` | PASS 58/58 | Final corrected-tree reverse order; all three executable hashes match the normal receipt exactly. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t3-qualified-shuffle-1729` | PASS 58/58 | Final corrected-tree deterministic shuffled order; all three executable hashes match normal and reverse exactly. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order normal --output-dir /tmp/amps-t3-native-normal` | Superseded bring-up PASS 58/58 | All three dimensional profiles passed, but an independent audit subsequently corrected a lower-dimensional array-bound defect; superseded by the qualified rerun above. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order reverse --output-dir /tmp/amps-t3-native-reverse` | Superseded bring-up PASS 58/58 | Reverse case order across the same three clean profile builds; not final authority after the post-campaign fix. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t3-native-shuffle-1729` | Superseded bring-up PASS 58/58 | Deterministic shuffled order before the audit correction; final authority is the qualified receipt above. |
| 2026-11-07 | `/tmp/amps-t3-final-fixes-fourth` | PASS 4/4 focused | Final 3-D interpolation-mask and planar cut-cell regressions before the independent lower-D audit. |
| 2026-11-07 | `/tmp/amps-t3-1d-escalated-initial`, `/tmp/amps-t3-2d-second` | PASS 6/6 + 6/6 focused | Search, neighbors, coarsen, data transfer, and cycles after lower-dimensional compilation/neighbor hardening. These are focused bring-up receipts, not the final matrix. |
| 2026-11-07 | `/tmp/amps-t3-3d-initial`, `/tmp/amps-t3-3d-second` | Historical FAIL 16/21, then 3/5 | Exposed missing interpolation configuration, sphere-recursion bounds, stale/incorrect AMR interpolation selection, and planar cut-cell contact/depth errors; later focused and 58-case bring-up runs passed. |
| 2026-11-07 | `python3 src/test/validate_traceability.py --output /tmp/amps-t3-doc-audit.json` | PASS structural audit | 24 areas, 209 production translation units, 181 public contracts, 16 configuration contracts, 17 suites, and 210 cases. This structural receipt is distinct from the final behavioral campaigns. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-PIC-NATIVE-001 --output-dir /tmp/amps-t2-native-full-normal` | PASS 43/43 | Complete baseline packed/fixed, species-global native owner in normal order. |
| 2026-11-07 | `python3 src/test/run_tests.py --campaign profile --order reverse --output-dir /tmp/amps-t2-profile-reverse` | PASS 85/85, 6/6 suites | All T2 generated profiles and cases in reverse suite/case order; each owner performed its required clean rebuild. |
| 2026-11-07 | `python3 src/test/run_tests.py --campaign profile --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t2-profile-shuffle-1729` | PASS 85/85, 6/6 suites | All T2 generated profiles and cases in deterministic shuffled order. |
| 2026-11-07 | `/tmp/amps-t2-layout-packed-fixed-recheck` | PASS 7/7 | Packed fixed record offsets, clone/reuse, list, capacity, invalid access, and exhaustion after correcting test offset types to the checked `long int&` API. |
| 2026-11-07 | `/tmp/amps-t2-split-merge-packed-recheck` | PASS 13/13 | Packed individual layout and split/merge after making the two-particle bin reachable and keeping the high-bulk fixture below the debugger speed guard. |
| 2026-11-07 | `/tmp/amps-t2-layout-aligned-fixed-initial` | PASS 7/7 | Aligned-16 fixed record profile with mixed-size anonymous extensions. |
| 2026-11-07 | `/tmp/amps-t2-layout-aligned-individual-recheck` | PASS 13/13 | Aligned-16 individual layout and split/merge after the opaque-state oracle accepted an exact clone of either valid source tuple. |
| 2026-11-07 | `/tmp/amps-t2-weight-local-recheck` | PASS 5/5 | Species-local block-associated weight/time and rescaling after adding the target-free exosphere generator block required by mesh headers. |
| 2026-11-07 | `/tmp/amps-t2-weight-single-global-initial` | PASS 4/4 | Single-global round trips, cross-species slot aliasing, and invalid host inputs. |
| 2026-11-07 | `/tmp/amps-t2-harness-final-prebuild` | PASS 11/11 | Harness registry, campaign-count, schema, protocol, subprocess, and traceability regression before the final native profile campaign. |
| 2026-11-07 | `/tmp/amps-t2-layout-packed-fixed-initial` | Historical build failure | The fixture used `int` momentum offsets with a production API requiring checked `long int&`; the type mismatch was corrected and the seven cases passed in the recheck. |
| 2026-11-07 | `/tmp/amps-t2-split-merge-packed-initial` | Historical FAIL: 11 PASS, 2 FAIL | Independent fixture exposed an unreachable two-particle split bin and a merge bulk speed above the debugger guard; corrected with identical-velocity input and a `1e9` bulk offset. |
| 2026-11-07 | `/tmp/amps-t2-layout-aligned-individual-initial` | Historical FAIL: 12 PASS, 1 FAIL | Oracle over-specified which opaque tuple member was cloned; corrected to accept either exact existing record under the documented clone/survivor policy. |
| 2026-11-07 | `/tmp/amps-t2-weight-local-initial` | Historical build failure | Minimal local deck omitted the exosphere-generated internal-boundary base classes; corrected by adding the target-free exosphere block. |
| 2026-11-07 | `/tmp/amps-t2-weight-time-initial` | PASS 5/5 | Focused species-global weight/time API cases from one clean native build; executable SHA-256 `46678c6054daad2ff9b7e3a1aabce9792d12139515e02f93bfaa0296210f5ffe`. |
| 2026-11-07 | `/tmp/amps-t2-all-normal` | Historical FAIL: 96 PASS, 2 FAIL, 5 SKIP | First 103-case normal campaign. `PHYS-3D-COLLISION-RELAXATION-001` and `PHYS-3D-COLLISION-RATE-001` exposed a test-fixture defect: restoring the `-1` uninitialized sentinel conflicted with the new finite-positive production setter contract. |
| 2026-11-07 | `/tmp/amps-t2-weight-time-collision-recheck` | PASS 7/7 | Corrective focused run: the two populated-cell collision cases and all five T2 weight/time cases pass together after raw saved lifecycle-state restoration; executable SHA-256 `6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`. |
| 2026-11-07 | `/tmp/amps-t2-all-qualified-normal` | 98 PASS, 5 SKIP, 0 FAIL/ERROR | Qualified normal-order receipt for the live 103-case registry; all 36 native cases pass and the five optional providers remain typed SKIPs. |
| 2026-11-07 | `/tmp/amps-t2-native-qualified-reverse` | PASS 36/36 | Complete native owner in reverse case order after a clean rebuild; executable SHA-256 `6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`. |
| 2026-11-07 | `/tmp/amps-t2-native-qualified-shuffle-1729` | PASS 36/36 | Complete native owner in deterministic shuffled case order after a clean rebuild; executable SHA-256 `6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`. |
| 2026-11-07 | `python3 src/test/validate_traceability.py` | Historical T2 PASS | Then-current T2 structural classification: 24 areas, 209 production translation units, 181 public contracts, 15 configuration contracts, 16 suites, and 152 registered cases. No behavioral case execution. |
| 2026-11-07 | `python3 src/test/run_tests.py --campaign smoke --jobs 4 --output-dir /tmp/amps-t1-smoke-qualified` | PASS 51/51, 8/8 suites | Corrected campaign-selection oracle and complete dependency-light T1 smoke gate. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order normal --jobs 4 --output-dir /tmp/amps-t1-all-qualified-normal` | 93 PASS, 5 SKIP, 0 FAIL/ERROR | Complete 98-case normal order; native 31/31 PASS; provider SKIPs remain typed and explicit. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order reverse --jobs 4 --output-dir /tmp/amps-t1-all-qualified-reverse` | 93 PASS, 5 SKIP, 0 FAIL/ERROR | Complete reverse suite/case order with a mandatory clean native rebuild. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order shuffle --shuffle-seed 1729 --jobs 4 --output-dir /tmp/amps-t1-all-qualified-shuffle-1729` | 93 PASS, 5 SKIP, 0 FAIL/ERROR | Complete deterministic shuffled order with a mandatory clean native rebuild. All three final native executables have SHA-256 `17b98e5ecd22915fc181cce51b8d71bd6abef517802a456ff2efba659c581919`. |
| 2026-11-07 | `/tmp/amps-t1-rng-hybrid-final`, `/tmp/amps-t1-mpi-integrated-final`, `/tmp/amps-t1-numerics-integrated-final`, `/tmp/amps-t1-units-integrated-final`, `/tmp/amps-t1-species-radiation-final-2` | PASS: 7 + 5 + 10 + 5 + 4 focused cases | Focused post-audit evidence for RNG, MPI channel/rank streams, arrays/numerics, units, and radiation tables. The complete campaigns above are the qualification authority. |
| 2026-11-07 | `/tmp/amps-t1-smoke-final` | Historical FAIL: 50 PASS, 1 FAIL | Preserved first live smoke receipt that exposed the stale `HARNESS-SELECTION-001` oracle; resolved by `/tmp/amps-t1-smoke-qualified`. |
| 2026-11-07 | `/tmp/amps-t1-interface-global` | 6 PASS, 5 SKIP | Executable C/Fortran/Geopack contracts pass; five external provider capabilities remain explicitly unavailable. |
| 2026-11-07 | `python3 src/test/validate_traceability.py --output /tmp/amps-core-traceability.json` | PASS | Structural classification: 24 areas, 209 production translation units, 181 public contracts, 12 configuration contracts, 4 suites, 45 registered cases. No behavioral case execution. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order normal --output-dir /tmp/amps-t0-all-final-normal` | PASS 45/45, 4/4 suites | Final-tree normal order; JSON/JUnit reconcile; 22 native executable cases each ran in a fresh exact-case MPI process. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order reverse --output-dir /tmp/amps-t0-all-final-reverse` | PASS 45/45, 4/4 suites | Final-tree reverse suite/case order; exact reverse receipt and identical native executable hash. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t0-all-final-complete-shuffle-1729` | PASS 45/45, 4/4 suites | Final-tree deterministic shuffled suite/case order; includes real timeout/SIGINT cleanup, terminal/JSON/JUnit reconciliation, distinct campaign assertions, effective generated-mode checks, and an explicit retained native preflight transcript. |
| 2026-11-07 | `/tmp/amps-t0-harness-final-complete` | PASS 11/11, 1/1 suite | Final harness-only receipt after campaign-set, interrupt, traceability, and JSON/JUnit reconciliation hardening. |
| 2026-11-07 | `/tmp/amps-t0-individual-<CASE-ID>` for all 21 dependency-light IDs | PASS 21/21 | Each dependency-light case selected and executed alone through the global `--case` interface. Native executable cases use exact one-case MPI subprocesses in every all-order campaign. |
| 2026-11-07 | `/tmp/amps-t0-individual-CONFIG-GENERATOR-001`, `/tmp/amps-t0-individual-CONFIG-INVALID-DECK-001` | PASS 1/1 each | Both adapter-owned configuration cases pass as exact one-case global selections after independent required clean rebuilds. Together with the 21 dependency-light probes and 22 one-case native MPI processes, every registered case has individual-selection evidence. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-normal` | PASS 21/21, 3/3 suites | Final dependency-light registry in normal order, including exact campaign-set assertions. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-reverse` | PASS 21/21, 3/3 suites | Same dependency-light registry in reverse suite/case order. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-shuffle-1729` | PASS 21/21, 3/3 suites | Same dependency-light registry with deterministic shuffle seed 1729. |
| 2026-11-07 | `/tmp/amps-core-native-t0-b` | PASS 24/24, 1/1 suite | Combined native normal order from one clean configured build; includes both Config.pl/generated-deck cases, three particle-buffer cases, and 19 physics cases. |
| 2026-11-07 | `/tmp/amps-core-native-t0-b-reverse` | PASS 24/24, 1/1 suite | Combined native reverse order from one clean configured build; JSON case selection/results reconcile exactly. |
| 2026-11-07 | `test_output/core-tests/run-20261107T072523Z-3412060/` | PASS 38/38, 5/5 suites | Last complete pre-consolidation execution; historical registry only. |
| 2026-11-07 | `/tmp/amps-core-report-final` | PASS 16/16, 3/3 suites | Last retained dependency-light receipt before the 11-case harness expansion. |
| 2026-11-07 | `/tmp/amps-core-physics-audit-final` | PASS 19/19 | Last retained native-physics family receipt before suite consolidation. |
| 2026-11-07 | `/tmp/amps-core-pbuffer-expanded` | PASS 3/3 | Last retained particle-buffer family receipt before suite consolidation. |

The following records are retained as historical bring-up/progression. Their
different case counts are expected for their then-current registries; none
supersedes the consolidated final-tree 45-case evidence above.

## Historical validation evidence

- `/tmp/amps-t0-final-report-reconciliation` correctly failed
  `HARNESS-REPORT-001` while adding direct JUnit coverage because the synthetic
  self-test receipt omitted the real writer's required `seconds` field. The
  fixture was corrected without changing the writer; the immediate rerun at
  `/tmp/amps-t0-final-report-reconciliation-2` and the final 11-case harness
  receipt both pass.

- `python3 src/test/run_tests.py --all` — PASS on 2026-11-07 after both
  required clean configured native builds: 38/38 individual tests and 5/5
  suites passed with 0 FAIL/SKIP/ERROR.  The command printed the final 38-row
  individual result table and the two summaries.  JSON/JUnit receipts and all
  suite transcripts are retained in
  `test_output/core-tests/run-20261107T072523Z-3412060/`; the JSON report
  independently records `individual_counts={TOTAL:38,PASS:38,FAIL:0,SKIP:0,ERROR:0}`
  and suite `counts={TOTAL:5,PASS:5,FAIL:0,SKIP:0,ERROR:0}`.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-report-final` — PASS after strengthening the formatter
  regression: `HARNESS-REPORT-001` directly captures and checks the heading,
  per-case states, and typed suite-level fallback; all 16 dependency-light
  cases were printed and recorded, with 3/3 suites and 0 FAIL/SKIP/ERROR.

- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-smoke` — PASS: the three nonexclusive smoke suites
  executed concurrently, all 15 individual cases passed, suite-qualified case
  starts were streamed, and aggregate results remained in manifest order.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-expanded-4` — PASS after the required clean configured
  build and one-rank MPI run: the exact 19-case registry passed with 0
  FAIL/SKIP/ERROR in 54.629 seconds. JSON/JUnit reports and the full native
  transcript are retained in that directory. This was the expanded
  native-physics evidence for that registry state.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-pbuffer-expanded` — PASS after the required
  clean rebuild with RT and collision-frequency storage enabled: exact 3-case
  registry, 0 FAIL/SKIP/ERROR in 37.360 seconds.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-final` — PASS after interruption/scratch/fail-fast
  hardening: 3 suites and 15 cases, 0 FAIL/SKIP/ERROR in 3.007 seconds.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-expanded-final` — PASS after the independent audit
  fixes and required clean rebuild: exact 19-case registry, 0
  FAIL/SKIP/ERROR in 54.818 seconds. This supersedes
  `/tmp/amps-core-physics-expanded-4` as the then-current native evidence.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-final-2` — PASS after the final scheduler timeout,
  interrupt-cleanup, fail-fast, and scratch-identity changes: 3 suites, 15
  cases, 0 FAIL/SKIP/ERROR in 3.015 seconds. This supersedes the earlier
  concurrent smoke receipts.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-audit-final` — PASS after all final independent-audit
  fixes and the required clean rebuild: exact 19-case registry, 0
  FAIL/SKIP/ERROR in 55.407 seconds. This was the final native-family evidence
  before consolidation.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-smoke-audit-final` — PASS: final scheduler/protocol tree, 3
  suites and 15 cases, 0 FAIL/SKIP/ERROR in 2.987 seconds. This was a prior
  dependency-light receipt.
- `mpiexec -n 1 ./amps --filter PIC-PBUFFER-` — PASS on the executable from
  the final clean native build: 3/3 particle-buffer cases after hybrid-thread
  fixture restoration was added.

- `python3 src/test/run_tests.py --tier smoke --output-dir /tmp/amps-core-test-final-2`
  — PASS: 3 suites, 14 individual cases, 0 FAIL/SKIP/ERROR. JSON and JUnit
  receipts are in that output directory.
- After the detailed code/documentation audit,
  `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-test-final-3` — PASS: the same 3 suites and 14 individual
  cases, 0 FAIL/SKIP/ERROR. This was the final dependency-light receipt for
  that earlier registry state.
- `AMPS_TEST_OUTPUT_ROOT=/tmp/amps-core-top-level make core_test_nightly`
  — PASS: validated the top-level Make/nightly entry point with the same 14
  cases.
- `SCRATCH=/tmp python3 src/test/run_tests.py --config
  src/test/config/machine.example.json --test CORE-HARNESS-001 --output-dir
  /tmp/amps-core-machine-profile-2` — PASS and verified that transient work is
  routed through the configured scratch root.
- `python3 src/test/run_tests.py --coverage --require-complete-coverage` —
  expected exit 1: 2 complete and 22 incomplete inventory components.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-native-final-3` — PASS after a clean configured
  native build and one-rank MPI run: 1 suite, 3 cases, 0 FAIL/SKIP/ERROR. The
  JSON/JUnit reports, full build log, and particle-buffer transcript are in
  `/tmp/amps-core-native-final-3`.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-5` — PASS after a clean configured build and one-rank
  MPI run: 1 suite, 11 cases, 0 FAIL/SKIP/ERROR. JSON/JUnit reports and the
  complete native transcript are retained in that directory.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-pbuffer-physics-final` — PASS after a second
  clean configured build: the pre-existing exact 3-case registry remains green
  with the expanded species and optional particle-data layout.
- `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-smoke-physics-final` — PASS: 3 dependency-light suites, 14
  cases, 0 FAIL/SKIP/ERROR after registering the physics suite.
- `python3 src/test/run_tests.py --all --list` and `python3
  src/test/run_tests.py --coverage` — PASS registry/inventory audit: 5 suites,
  24 components, 2 covered, and 22 honestly partial/planned.
- `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-case-progress-final` — PASS: all 14 cases were announced by
  stable ID and description before execution, with unchanged JSON/JUnit
  results.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-case-progress-native` — PASS after the required
  clean native build: all 3 MPI-launched cases were announced live and the
  exact registry remained green.
- `git diff --check` — PASS after the implementation and documentation audit.
- `PYTHONPYCACHEPREFIX=/tmp/amps-core-pycache python3 -m compileall -q
  src/test` — PASS for every Python runner/framework module without writing
  bytecode into the source tree.

## Native bring-up record

The first clean native attempt (`/tmp/amps-core-native-final`) failed while
compiling `meshAMR`: the new input deck had not invoked `exosphere.pl`, which
is the established generator for AMPS' default internal-boundary user-data
classes. The fix adds a target-free, SPICE-free parser block. It supplies type
definitions only; exosphere physics and source modules remain excluded.

The second clean attempt (`/tmp/amps-core-native-final-2`) compiled and linked
successfully but reported no cases. AMPS stores application objects in the
static archive `main.a`, and the linker did not load `particle_buffer.o` when
it contained only static registration constructors. The fixture now exports
`LinkParticleBufferTests()`, called by `main.cpp`, to make archive membership
deterministic. The third clean attempt passed the exact three-case registry.
Both failed attempts remain useful diagnostic evidence and were not described
as test passes.

## Native physics bring-up record

- `/tmp/amps-core-physics-1` failed during input generation because the
  collision table used a single backslash instead of AMPS' double-backslash
  continuation token.
- `/tmp/amps-core-physics-2` reached the IDF parser and rejected unsupported
  `RotZnum all` shorthand. The unused rotation-collision and pair-index values
  were removed because every relaxation switch is off in this kernel fixture.
- `/tmp/amps-core-physics-3` built and passed all nine non-mover cases. The
  mover fixture incorrectly selected an unallocated AMR internal/root node;
  the second case inherited it and triggered the expected signal handler.
- `/tmp/amps-core-physics-4` resolved the owner leaf and produced correct
  ballistic states, but inspected the starting leaf's temporary list for an
  endpoint exactly on a leaf face. The fixture now avoids an ambiguous face
  endpoint and resolves the destination owner explicitly.
- `/tmp/amps-core-physics-5` passed the exact 11-case registry. None of the
  failed bring-up runs is counted as validation evidence.
- `/tmp/amps-core-physics-expanded-1` executed all 19 expanded cases and
  reported 16 PASS / 3 FAIL. The failures were test-fixture/oracle defects:
  the neutral-only build intentionally omits unavailable ionic daughters, the
  RT oracle initially used one molecule's rotational exponent rather than the
  pair's six rotational degrees, and the Boris free-flight trajectory crossed
  the newly registered sphere. All three definitions were corrected without
  relaxing a production tolerance.
- `/tmp/amps-core-physics-expanded-2` was manually interrupted because the
  240,000-draw branching check also generated full two/three-body kinematics
  on every draw and was unsuitable for routine execution. Its orphaned adapter
  and MPI processes were explicitly cleaned up. The final test separates the
  high-statistics production channel selector from a smaller product-contract
  sample; product kinematics retain their dedicated conservation case.
- `/tmp/amps-core-physics-expanded-3` produced a typed `ERROR` because the
  interrupted adapter still held the repository-specific native lock. No test
  result from that invocation is counted as scientific validation.
- `/tmp/amps-core-physics-expanded-4` passed the final exact 19-case registry
  in 54.629 seconds. It supersedes the three diagnostic expanded bring-up runs;
  those failures/interruption remain recorded and are not described as passes.
