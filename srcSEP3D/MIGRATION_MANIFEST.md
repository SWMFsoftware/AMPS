# srcSEP3D Migration Manifest

This manifest records what Phases R0–R2, M, B, T, P, A, O, and V removed, retained, or
replaced. It is an
evidence document: completion claims must correspond to an executable test or
an explicitly identified external build gate.

## Phase R0 — Rebaseline the production tree

### Deleted production files

| File | Disposition | Reason |
|---|---|---|
| `SEP3D.cpp` | deleted | contained the axisymmetric y=0 mover and legacy global energy sampler; neither belongs in a three-dimensional transport application |
| `amps/amps_sampling.h` | deleted | declaration-only placeholder with no production implementation; retaining it implied sampling existed when it did not |
| `core/sep3d_energy_distribution.{h,cpp}` | deleted | unintegrated replacement kernel had no registered acceptance test or production consumer; Phase O now provides one complete sampling implementation under `output/` |
| `test/stage1` | removed from deliverable | generated executable from the uploaded archive, not source; rebuilt locally by the runner |

### Retained and rewritten files

| File | R0 state |
|---|---|
| `main_lib.cpp` | retains only AMPS/Exosphere link hooks and explicit execution guards; wedge mesh, fixed resolution, Maxwellian prepopulation, placeholder sampling, and placeholder outputs are removed |
| `main.cpp` | clean AMPS driver source retained for the `main.a` archive; strict-warning clean by construction |
| `SEP3D.h` | declares only current application hooks and includes the AMPS mover-status adapter; no retired mover or sampler declarations |
| `core/sep3d_types.h` | replaces integer-valued mover return codes with semantic `ParticleMotionOutcome` |
| `makefile` | truthful production manifest, enclosing-AMPS production build, archive symbol audit, AMPS-independent test targets, and absolute path discovery valid from both `srcSEP3D` and copied `build/main` locations |
| `test/run_tests.py` | rebuilt around the selector vocabulary used by `srcSEP` and structured JSON/JUnit evidence |

### Mover return-code correction

The prototype core assigned local integers to mover results. They did not match
the AMPS definitions in `pic.h`: AMPS uses deleted-on-face `0`, left-domain
`2`, and motion-finished `3`. A model extension also occupied a value used by
AMPS. This is corrected as follows:

1. L0 returns `Core::ParticleMotionOutcome`, whose enumerators are semantic and
   are never used as AMPS integers.
2. `amps/amps_mover_status.h` is the only conversion boundary.
3. The adapter maps `Advanced` to `_PARTICLE_MOTION_FINISHED_` and
   `LeftDomain` to `_PARTICLE_LEFT_THE_DOMAIN_`.
4. `static_assert` binds the adapter to the actual AMPS ABI at production
   compile time.
5. `BLDL3D03` also reads the supplied AMPS `pic.h` and records the observed
   values in the test evidence.

Invalid background, underflow, and internal error statuses are not mapped to a
valid AMPS particle result. They must be handled before the adapter is called.

### Production execution policy

The host installs an immutable configuration through `ConfigureApplication`.
Phases M/B/T build the AMPS mesh, reserve/fill the frozen cell layout, publish
an analytic or coupled ambient snapshot, and validate prescribed or coupled
scattering input. Phase P supplies the numerical transport; Phase A supplies
the single AMPS dispatcher and common SWCME source conversion; Phase O supplies
read-only products and complete restart state. `amps_time_step()` now executes
the typed Runtime/PIC step. It does not silently substitute a local-state
resolver, shock schedule, or MPI observation gather that the host has not
installed.

### R0 evidence

| Gate | Implementation | Source-only behavior |
|---|---|---|
| `BLDL3D01` | delegates through `make strict-production` to the enclosing top-level `make amps`, then audits `build/main` archives | SKIP when a real `Makefile.conf` is unavailable |
| `BLDL3D02` | production-tree and makefile-manifest inspection | executable without AMPS |
| `BLDL3D03` | adapter contract plus real `pic.h` macro inspection | executable with an AMPS source tree; SKIP if `pic.h` is unavailable |
| `BLDL3D04` | compile-only AMPS `Pi`-macro collision regression | executable without a configured AMPS build |
| `BLDL3D05` | source/copy makefile path-resolution fixture | executable without AMPS libraries; invokes both makefiles from an unrelated working directory |
| `LAY01`, `LAY02`, `BLD01` | three-layer AMPS-dependency guard | executable without AMPS/MPI |

## Earlier Step 1–4 material: corrected status

The uploaded tree described Steps 1–4 as complete. R0 corrects that claim:

- the directory layout and standalone registry are useful and retained;
- the original production source had not actually retired the axisymmetric
  mover until this R0 change;
- `sep_common` now owns the canonical shared sources and archive beside both
  applications; the former `srcSEP/util` copies are absent;
- SWCME now has one canonical `src/models/swcme` implementation, which
  srcSEP3D consumes directly without inspecting or requiring srcSEP;
- a configured enclosing AMPS executable link is still required before the
  production-build portion of R1 can be declared complete;
- R2 owns validated snapshot metadata and lifecycle transitions; Phases M/B/T
  now provide the physical mesh/storage, ambient field, and scattering input
  that those transitions describe.

No later physics phase should use the old Step 1–4 completion labels as release
evidence. The current runner reports the gates that are actually executable.

## AMPS macro compatibility correction

`SEP3D::Core::Const::Pi` was renamed to `kPi`. AMPS exposes `Pi` as a global
macro in `general/constants.h`; after `pic.h` includes the application header,
the preprocessor rewrote the namespaced declaration and caused
`cell_centered_linear_interpolation_cpp.cpp` to fail with “expected
unqualified-id before numeric constant.” `BLDL3D04` now reproduces this include
order in a compile-only regression.

## Copied-build makefile path correction

AMPS copies `srcSEP3D` into `AMPS/build` and renames that application directory
to `AMPS/build/main`. A literal `../Makefile.conf` is therefore wrong in the
copied tree: it resolves to `AMPS/build/Makefile.conf` instead of the canonical
`AMPS/Makefile.conf`.

The makefile now obtains the absolute directory containing the active makefile
from `MAKEFILE_LIST`, identifies whether it is the source or copied location,
and resolves `AMPS_ROOT` once. `AMPS_CONFIG`, `SEP_COMMON_DIR`, `SWCME_DIR`,
both shared archives, and all shared-object paths are absolute descendants of
that root.
This logic does not depend on the current working directory. `BLDL3D05`
constructs both layouts, invokes each makefile from a third directory, and
requires identical canonical paths.

## Enclosing production-build correction

The application makefile is designed to compile production objects after AMPS
has copied the application to `AMPS/build/main`. A direct production submake in
`AMPS/srcSEP3D` lacks the include flags for generated headers such as
`AMPS/build/pic/pic.h`; such a compile is not equivalent to the production
build. `strict-production` and `BLDL3D01` therefore delegate to the top-level
`make amps` target and audit the resulting absolute `build/main` archive paths.

## Forbidden production identifiers and operations

`BLDL3D02` rejects the following in production sources:

- `Mover_Axisymmetric_SecondOrder`
- `TotalParticleAcceleration`
- `GlobalEnergyDistribution`
- `inject_particle_onto_field_line`
- `CMPI_channel`
- the former wedge bounds `8.760e+08` and `9.445e+08`
- `PrepopulateDomain`, placeholder mesh output, or mesh-file creation in the
  R0 application boundary

Documentation and tests may name retired identifiers when explaining or
detecting them; production headers and translation units may not.

## Phase R1 — shared kernel and SWCME integration

| Contract | Implementation | Evidence |
|---|---|---|
| one dependency-free SEP kernel archive | `src/models/sep_common/makefile` produces exactly seven members | `ARCH3D02`, `UTIL02`, archive `verify` |
| one compiled SWCME archive | `src/models/swcme/makefile` produces exactly `swcme3d.o` | `ARCH3D02`, `R1SW01` |
| srcSEP3D consumes canonical objects | its makefile rebuilds `mainlib.a` from local objects plus the canonical object lists | `ARCH3D02` |
| bounded common 1-D/3-D coupling runner | `src/models/swcme/test/run_tests.py` exposes srcSEP-style selectors and JSON/JUnit | `SWCME3D01`, runner unit tests |

The uploaded SWCME subset omits translation units referenced by its historical
extended Makefile. R1 does not report those missing campaigns as passes. The
new runner exposes only the four executable common-library gates, while the
larger catalog remains pending source restoration.

## Phase R2 — immutable configuration and Runtime lifecycle

| File | Ownership and invariant |
|---|---|
| `runtime/run_configuration.{h,cpp}` | validates host options, derives the complete pre-mesh byte layout, and publishes immutable physics/output manifests |
| `runtime/runtime.{h,cpp}` | owns lifecycle state, snapshot metadata, pinned generation, and restartable step/output/checkpoint counters |
| `runtime/runtime_adapters.{h,cpp}` | maps analytic Parker and SWMF authorities onto identical acquisition/publication calls |
| `main_lib.cpp`, `SEP3D.h` | own and expose the process Runtime without parsing process arguments or parameter files |

The lifecycle is `Created → Configured → MeshReady → WaitingForSnapshot →
SnapshotReady`, with `Running` and `Checkpointing` as transactional excursions
back to `SnapshotReady`, and `Finalized` as the terminal state. Every mutating
method validates all inputs before committing state. `LIFE3D02` exercises all
80 operation/state pairs and proves rejected calls preserve state, counters,
and snapshot generation.

The physics fingerprint includes authorities, transport/domain selection,
physical radii and timestep, seed, background cadence, and storage-layout
identity. Output directory, prefix, and cadence remain in the resolved manifest
but are intentionally excluded from physics identity. `LIFE3D03` protects that
separation and the exact storage offsets.

The corrected species boundary treats AMPS `SpeciesList` as immutable and
complete. Runtime no longer carries an asserted index, species name, mass, or
charge and no longer calls molecular-data setters. Immediately after AMPS base
initialization, the application enumerates every generated index through the
chemical-symbol, mass, and signed-charge accessors and validates the complete
table before mesh or after-parser initialization. The configured timestep and
base weight are then installed for every species and every owner-local block.
Injection iterates the same table, assigns the exact per-species sample count,
and reconstructs momentum bounds from each AMPS mass. `CFG3D07` supplies
positive mixed ion/electron coverage and independent negative controls for
count, index continuity, symbol uniqueness, mass, charge, and observer range;
`R3D05` covers the species-dependent energy-to-momentum conversion.

## Phase M — mesh and storage

| File | Ownership and invariant |
|---|---|
| `mesh/mesh_model.{h,cpp}` | owns domain presets, radial/tube resolution, standalone balanced octree, deterministic identities, memory estimate, owner-only storage, and mixed-resolution gradients |
| `runtime/run_configuration.{h,cpp}` | fingerprints Phase-M options and freezes the complete background/turbulence cell layout before allocation |
| `main_lib.cpp` | registers exact AMPS static/sampling requests, drives `localResolution`, partitions/allocates the production mesh, and binds the frozen layout |
| `MESH_STORAGE.md` | records formulas, byte order, AMPS initialization order, and evidence |

The standalone octree is a verification oracle, not a competing production
mesh. Production allocation stays in AMPS. Both call the same resolution
function, and production cell population uses only owner-local decomposition
blocks. `MSH3D01–09` are the release evidence.

## Phase B — background providers and snapshots

| File | Ownership and invariant |
|---|---|
| `background/bg_provider.{h,cpp}` | complete SI sample and transactional per-point batch interface |
| `background/bg_parker.{h,cpp}` | analytic Cartesian Parker field/plasma state and derivatives with finite polar limits |
| `background/bg_swmf.{h,cpp}` | read-only, unit-aware, frame/epoch-checked SWMF/AWSoM ambient import |
| `background/background_snapshot.{h,cpp}` | validate-before-publication immutable snapshot and compatible current/next interpolation |
| `runtime/runtime_adapters.{h,cpp}` | publishes physical Parker/SWMF snapshots through the existing neutral Runtime path |
| `BACKGROUND_FIELD.md` | complete-field, coupling-unit, ownership, and atomicity contract |

No provider writes AMPS memory. The L3 boundary maps a validated immutable
snapshot to the owner-local cell list and then publishes its metadata. Failed
candidates and failed batch indices do not modify active output.
`BGP3D01–06` and `SNAP3D01–08` are the release evidence.

## Phase T — turbulence and scattering inputs

| File | Ownership and invariant |
|---|---|
| `turbulence/turbulence_provider.h` | separate turbulence authority, typed missing/ballistic state, directional magnetic variance |
| `turbulence/turbulence_models.{h,cpp}` | sep_common-independent provider API, normalized prescribed spectrum, explicit AWSoM w+/w− conversion, and resonance policy |
| `turbulence/coefficient_bridge.h` | opt-in adapter declarations for canonical sep_common coefficient types; intentionally excluded from the AMPS-facing provider header |
| `runtime/run_configuration.{h,cpp}` | validates/fingerprints amplitude, spectral band/index, correlation scale, missing policy, and resonance policy |
| `TURBULENCE_SCATTERING.md` | units, sign convention, normalization, policies, and shared-kernel boundary |

The AWSoM convention is field aligned: plus travels along `+B`, minus against
it, and `deltaB²=mu0*w`. Heliocentric outward/inward aliases are determined
from `B·r`, not assumed from a variable name. The coefficient bridge calls the
canonical `sep_common.a`; `COEF3D02` checks bitwise identity with a direct
shared-kernel call. Self-consistent 3-D turbulence remains a configuration-time
reserved feature. `TUR3D01–04` and `COEF3D01–02` are the release evidence.
`BLDL3D06` additionally proves that `SEP3D.h` does not transitively expose
source/restart headers to generic AMPS translation units, that the provider
header compiles without a sep_common include path, and that only the opt-in
bridge requires the canonical one. Its production-recipe probe deliberately
ignores `CPPFLAGS`, `CXXFLAGS`, and `INCLUDE`, then verifies that the
target-scoped `CPLUS_INCLUDE_PATH` still supplies `src/models/sep_common` and
`src/models/swcme` while compiling the root-level `main_lib.o` target. This is
the exact build-system boundary implicated by installed AMPS generic recipes.

## Phase P — transport cores

| File | Ownership and invariant |
|---|---|
| `transport/parker_transport.{h,cpp}` | complete tensor-Parker Itô step and exact frozen-divergence cooling |
| `transport/focused_transport.{h,cpp}` | bounded split focused SDE with declared pitch scheme |
| `transport/time_step.{h,cpp}` | named stability/accuracy limits with typed, unclamped underflow |
| `transport/keyed_random.{h,cpp}` | particle/purpose-keyed restartable stochastic streams |
| `TRANSPORT_CORES.md` | governing equations, numerical split, limits, and reproducibility contract |

The core contains no AMPS list manipulation or source injection. Perpendicular
diffusion and drift inputs are explicit zero-only guards. `COEF3D03–05`,
`PRK3D01–08`, `FTE3D01–07`, and `RNG3D01–03` are the Phase-P release evidence.

The production AMPS resolver now fills the previously missing
`dKappaParallelDsMPerS` input by sampling the same canonical coefficient path
at field-aligned neighboring positions. Centered differences are preferred;
one-sided differences are explicit at boundaries, and the resolver fails
closed when neither neighbor is usable. The host-neutral arithmetic and its
boundary policy are registered as `COEF3D06` and run in `phase-t`.

## Phase A — AMPS mover and source adapters

| File | Ownership and invariant |
|---|---|
| `adapters/transport_adapter.{h,cpp}` | complete-record validation, exact two-core dispatch, named substep, boundaries, and expanding-shock crossing |
| `adapters/particle_ledger.{h,cpp}` | exact integer particle conservation per step/species |
| `adapters/swcme_source_adapter.{h,cpp}` | canonical SWCME source to shared DSA sampler; dimension-independent semantic keys |
| `amps/amps_particle_adapter.{h,cpp}` | packed AMPS particle state, ABI return mapping, deterministic velocity reconstruction, and destination-list insertion |
| `AMPS_ADAPTERS.md` | physics mapping, algorithms, host responsibilities, and evidence |

The AMPS allocator slot is never a particle identity. New particles must be
initialized with the stable ID and gyrotropic state returned by the common
source adapter. The generated mover configuration must route to the one Phase-A
dispatcher and make its declaration visible to the AMPS mover translation
unit. `ADP3D01`, `NAT3D04–05/08`, and `SHK3D01–04` are the standalone evidence;
`BLDL3D01/03` remain the linked ABI gate.

## Phase O — sampling, output, and restart

| File | Ownership and invariant |
|---|---|
| `output/sampling.{h,cpp}` | stable-ID ordering, compensated cell/spacecraft/field-line reductions, and closed-ledger shock diagnostics |
| `output/publication.{h,cpp}` | unit-bearing CSV schemas, identity manifest, artifact hashes, and atomic directory publication |
| `output/restart.{h,cpp}` | versioned canonical little-endian codec for all runtime/stochastic/particle/generation/sampling/ledger state |
| `output/output_coordinator.{h,cpp}` | Runtime-owned cadence, checkpoint commit/rollback, and pre-mesh restore |
| `SAMPLING_OUTPUT_RESTART.md` | products, equations, schemas, transactional rules, and evidence |

No sampler changes a particle or consumes a random stream. No restart parser
mutates caller state before all checks succeed. A missing coupled snapshot is
either rejected or awaited with an explicit bounded callback; another
generation is never substituted. `NAT3D06–07` and `RST3D01–03` are the Phase-O
evidence.

## Phase V — integration and scientific validation

| File | Ownership and invariant |
|---|---|
| `validation/validation_metrics.{h,cpp}` | canonical stable-ID rank merge, exact global ledger reduction, resource budgets, and positive-series scientific metrics |
| `validation/case_registry.json` | immutable IDs, evidence classes, release/diagnostic roles, and quantitative thresholds |
| `validation/run_validation.py` | linked invocation, checksum/provenance validation, series/convergence evaluation, and atomic JSON/JUnit reports |
| `validation/templates/` | non-passing instructions for series and convergence evidence bundles |
| `test/individual-test/test_validation.cpp` | controlled integration plus analytical Parker/focused/SWCME validation |
| `test/test_validation_runner.py` | CLI, SKIP, checksum-failure, and convergence-runner regression tests |
| `INTEGRATION_SCIENTIFIC_VALIDATION.md` | algorithms, equations, evidence boundary, case roles, and physical interpretation |

Phase V preserves application independence: cross-model evidence is exported
by its owner and consumed as immutable bytes; neither runner constructs a path
to the sibling `srcSEP` application. Controlled `INT3D`/`VFY3D` results cannot
satisfy linked `NAT3D`/`MPI3D`, cross-model `XM3D`, or observational `OV3D`
gates. Missing external evidence is reported as `SKIP`, while malformed or
checksum-invalid evidence is `ERROR`.

The controlled release prerequisites are `INT3D01–03`, `VFY3D01–05`, and
`VALRUN3D01`. Linked and scientific release closure requires the real target
executable and reviewed evidence bundles; the source package intentionally
contains templates rather than manufactured passing event records.

## Remaining production boundary

The remaining work is release evidence rather than a placeholder numerical
phase: configure the AMPS mover declaration/macro, install a pinned-snapshot
coefficient resolver, connect SWCME event scheduling, feed deterministic global
observations into the Phase-V audit, and run the registered linked MPI and
scientific campaigns.
