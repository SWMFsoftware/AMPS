# Standalone SEP integration roadmap: existing core shock/CME models

## Purpose and scope

Integrate the existing core shock/CME implementations into `srcSEP`, starting
with the shared reduced shock-front provider selected by the user. Add
adapters, application orchestration and tests while retaining the core models
as the physics authorities. Developing new physical shock models is outside
this integration roadmap.

The first deliverable is a qualified Parker-spiral transport application
supporting all three production movers with that provider:

- `parker`: pitch-angle-averaged Parker transport.
- `fte-dmumu`: focused transport with pitch-angle diffusion.
- `fte-mfp`: focused transport with discrete scattering events.

Other existing core backends follow through their maintained interfaces and
declared capabilities. The MHD-driven application remains a separate
integration path; existing coupling regressions remain required checks for
unintended shared-code changes.

This is a staged work plan, not authorization to implement every stage at
once. The user authorized stages 0–4 and subsequently stages 5–6; their implemented scope, executed tests,
source identities and outstanding validation boundaries are recorded in
[STANDALONE_INTEGRATION_STATUS.md](validation/STANDALONE_INTEGRATION_STATUS.md).
The approved stages 5–6 source/boundary decisions are recorded below. Stages
7–9 remain separate work; neither a short native campaign nor a passing
unsupported-mode rejection qualifies those later capabilities.

## Requirements that apply to every stage

- AMPS' physical timestep remains constant. Geometry integration, tracing and
  particle numerical subdivisions do not change the physical timestep or
  create extra sampling opportunities.
- Particle motion relative to plasma-comoving field lines excludes added
  solar-wind displacement. Physical flow and its derivatives remain available
  to momentum, scattering, source and shock-relative calculations.
- Move every existing vertex, including both endpoints. Extend inward along
  B; new segments start empty and existing source rules populate them. Remove
  trailing segments only when completely outside the maximum radius, deleting
  attached particles once and recording represented losses.
- Sampling exposure belongs to each segment. Preserve the global output
  schedule and distinguish an unmeasured segment from measured vacuum.
- Do not relax test tolerances, replace missing physics with assumptions, or
  describe unsupported downstream transport as implemented. Request physical
  clarification before dependent implementation when a contract is missing.
- Core models own their equations, configuration and physical state.
  Application adapters must not introduce competing shock or ambient laws.
- Use SI and explicit frames, composition, epochs, identities and validity.
  Preserve immutable particle-read realizations and transactional publication.
- Add detailed comments explaining physics, numerical methods, ownership,
  units, invariants, validity limits and failure paths to implementation code.
- Preserve unrelated local changes and regression references. Reference
  regeneration requires separate physical verification and reviewed provenance.

## Stage 0 — Freeze the integration scope and establish the baseline

### Work

Create a capability inventory for the exact source revision being integrated:

| Existing component | Integration role | Admission condition |
| --- | --- | --- |
| Shared reduced shock-front provider | First production target: ambient, analytical front and local shock states | Required core and native-background gates pass |
| SWCME1D | Additional analytical backend through its existing adapter | Selected profile supplies every primitive and derivative required by the requested transport mode |
| Shared coronal-CME implementation | Existing capabilities through its maintained neutral field-line interface | Exported state, identities and temporal coverage satisfy the consumer contract |
| Full corona–SWCME regional model | Future consumer of completed core capabilities | Required sheath/ejecta/material-state qualifications are complete |

Capture the commit, current patch, required untracked source files, compiler
configuration, input assets and reference checksums. Build in isolated
workspaces so an active nightly campaign and unrelated local changes remain
untouched.

### Tests

Establish the current individual results of the registered suites:

- `CORE-MODEL-SEP-CORONA-SWCME-REDUCED-001`
- `CORE-MODEL-SWCME-001`
- `CORE-MODEL-SEP-CORONAL-CME-001`
- `CORE-MODEL-SEP-CORONA-SWCME-FULL-001`
- `CORE-PIC-FIELD-LINE-001`

Use the executable registries rather than documentation totals as the source
of selected case identities and observed results.

**Exit gate:** a reproducible baseline and a capability matrix identifying
supported, incomplete and unqualified modes. A failed prerequisite blocks its
dependent integration work; unrelated failures remain visible.

## Stage 1 — Establish the core-to-SEP state contract

### Work

Extend existing thin adapters rather than reproduce core equations in
`srcSEP`. Consumer-facing records need:

- Provider/configuration identity, epoch, generation and validity interval.
- Cartesian position, B, U and their declared coordinate basis.
- Mass density and explicitly identified species densities.
- Pressure, temperature and composition definitions.
- Velocity-gradient tensor and magnetic derivatives.
- Shock geometry, normal speed, classification and available one-sided states.
- Explicit capability and failure status.

The reduced adapter currently identifies its number density as electron
density. Proton injection must use the core composition to obtain proton
density; electron and proton density must not be treated as interchangeable.

The shock calculation and particle transport must use the same upstream state.

### Tests

- Compare native adapter values with direct core queries at identical
  positions and times.
- Independently check SI conversions, density/composition mapping and vector
  transformations.
- Detect mismatched event identity, stale generation and invalid coverage.
- Inject a late query failure and verify that no partial state is published.
- Verify explicit rejection of missing downstream or thermal capabilities.

**Exit gate:** one coherent, validated state reaches the native field-line
interface without alternate analytical formulas or hidden fallbacks.

## Stage 2 — Connect material field-line motion to the core ambient

### Work

Supply the implemented geometry container with application callbacks for:

- Vertex trajectories under the selected core velocity.
- Target-epoch magnetic tracing.
- Physical state at moved and newly inserted vertices.
- Stable material identities and source identities.

Use the core's actual velocity, including any radial variation or azimuthal
component. The exact constant-radial-speed map is only a specialization for a
profile that declares those assumptions.

Preserve the approved vertex-motion, empty-insertion and outer-removal rules.
Validate compatibility between B, flow and moving geometry. Do not resample a
stationary field that is inconsistent with the advected line. Provider coverage
must include the actual inserted endpoints, including any below the nominal
minimum radius after a complete tracing increment.

### Tests

- Exact uniform-translation and constant-radial-flow references.
- Movement of both endpoints.
- Repeated insertion, whole-segment trimming and large-displacement cases.
- Magnetic-field/tangent agreement under spatial refinement.
- Flux-tube area and volume consistency.
- Failed preparation leaves live geometry, particles and sampling unchanged.
- No solar-wind term enters relative particle displacement.
- New segments have empty particle lists and zero sampling exposure.

**Exit gate:** the actual native container advances correctly under the selected
analytical background, with unchanged AMPS timestep.

## Stage 3 — Integrate analytical shock connectivity in space and time

### Work

Add native field-line connectivity for the selected finite SSE front. Find
candidate intersections, reject unsupported roots and obtain local physics
through the core front evaluator. An origin-centred spherical routine must
not substitute for finite-front geometry.

Connectivity must support:

- Multiple physical intersections.
- Segment-endpoint hits.
- Tangency and loss of connection.
- Entry into and exit from the supported front.
- Explicit queries within the complete AMPS interval.

A source born at time `t_b` must be attached to the correct material coordinate
at `t_b`. It must not combine a future shock position with departure-time
geometry. Provide a read-only representation of geometry evolution over the
interval, even if live native topology is committed only at joined boundaries.

### Tests

- Independent plane, sphere and SSE intersection references.
- Unsupported rear-surface rejection.
- Endpoint ownership consistent with particle coordinate conventions.
- Connection changes during a step.
- Surface residual, normal and normal-speed verification.
- Queries do not mutate the committed provider epoch.
- Source location is consistent with its physical birth time.

**Exit gate:** every accepted connection has a valid core shock record and a
consistent material coordinate and epoch.

## Stage 4 — Qualify the three transport solvers on the same moving background

### Work

Audit the complete native formulations before comparing results.

For Parker transport, establish where the flux-tube metric, coefficient
gradient, statistical weighting and momentum measure enter. Removing bulk
coordinate advection does not remove geometric diffusion terms.

For focused transport, verify streaming, focusing, strain, cooling and
scattering using the same immutable background. Preserve the distinction
between physical solar-wind quantities used in momentum and shock physics and
relative coordinate advection, which vanishes on plasma-comoving lines.

### Tests

| Solver | Required checks |
| --- | --- |
| `parker` | Constant and spatially varying diffusion; tube-area metric; expansion cooling; independent deterministic transport reference |
| `fte-dmumu` | Ballistic streaming; focusing; pitch-angle equilibrium; bounded pitch-angle treatment; expansion and strain |
| `fte-mfp` | Event-rate and waiting-time distributions; residual-hazard persistence; ballistic limit; declared scattering-frame conservation |
| Cross-solver | Parker/focused agreement in the common diffusion limit using equivalent coefficients and source normalization |

Discrete-event and continuous-diffusion movers need not reproduce identical
early anisotropic distributions merely because they share a mean free path.

Hold AMPS' timestep fixed. Refine internal numerical controls and compare
against high-accuracy references at the same physical endpoints. If fixed-step
coupling errors exceed the accuracy requirement, improve the integration
composition rather than shorten the physical timestep.

**Exit gate:** each solver reproduces its governing equation, and comparisons
use physically equivalent inputs.

## Stage 5 — Connect existing core source and acceleration capabilities

### Work

Inventory existing core source/acceleration interfaces and determine which
are compatible with the selected provider. For each admitted source, declare:

- Species and upstream population.
- Energy/momentum distribution.
- Pitch-angle distribution and frame.
- Efficiency and physical normalization.
- Release side and timing.
- Whether acceleration is represented by the source or resolved by transport.

Reuse an existing core implementation only after validating that it consumes
the selected provider's actual state. Do not import another backend's ambient
density or compression through its source helper.

Account for the physical intersection area between a flux tube and the shock.
Tangency must receive a valid finite-footprint treatment or an explicit
unsupported status; an obliquity floor is not a solution.

### Tests

- Independent incident and injected particle counts.
- Weight normalization and zero-efficiency limit.
- Species-density mapping.
- Supported versus sub-fast/unsupported front regions.
- Birth time and remaining movement duration.
- Source identity unchanged by geometry insertion or MPI ownership.
- No double counting between prescribed acceleration and resolved compression.
- Manufactured verification sources remain distinguished from production
  source physics.

**Exit gate:** a selected existing core prescription produces the correct
represented population through the native source path.

**Prerequisite:** if no existing source prescription supplies the required
physics, record that gap and request the physical selection. This integration
task must not invent one.

## Stage 6 — Deliver the first end-to-end reduced-front transport application

### Work

Add an explicit particle-enabled consumer mode while preserving the current
background-only mode and its zero-particle tests. Integrate the full sequence:

1. Prepare the immutable core realization for the physical interval.
2. Establish material geometry and shock connectivity.
3. Evaluate source events at their declared times.
4. Advance particles through the complete AMPS interval.
5. Complete sampling and joined particle work.
6. Commit geometry/background updates.
7. Refresh connectivity, ownership and observer caches.

The reduced model has no downstream volume. Its particle mode therefore needs
an explicit policy at that boundary. Immediate jump states must not be extended
throughout the interior.

### Tests

Run the same event and source with all three movers:

- Particle-free background.
- Manufactured monoenergetic source.
- Selected production source.
- Connected and disconnected observer cases.
- Front crossing a field-line vertex.
- Inner insertion and outer removal during active transport.
- Invalid provider query during preparation.
- Number accounting: injected, retained, escaped and boundary-removed.
- Energy accounting with separately identified injection, background work and
  escape.

For Parker diffusion, test moving-boundary first passage; endpoint-only
crossing detection is insufficient.

**Exit gate:** a complete one-rank native campaign passes independent source,
transport, sampling and accounting checks.

**Completed 2026-11-22 for the declared prescribed reduced-SSE/Parker case.**
Parker, FTE-Dmumu and FTE-MFP each completed 400 unchanged 300-s physical steps
and passed the complete-driver gate. See the
[qualification report](validation/STANDALONE_INTEGRATION_STATUS.md) for source
identities, independent checks, numerical/MPI results and remaining scope.

## Stage 5–6 source and boundary decisions

The user selected an **already accelerated population**, using the existing
core local-compression test-particle DSA spectrum. The transport calculation
must not add another prescribed shock-acceleration operator to this population.
Ambient adiabatic momentum evolution remains part of the transport equations.
For the core's accepted local compression `r`, retain the distinct measures

```text
q = 3 r/(r-1),    f(p) proportional to p^(-q),
dN/dp proportional to p^(2-q).
```

This is the selected spectral prescription, with the usual steady,
nonrelativistic-shock, test-particle diffusion assumptions; it does not solve
time-dependent acceleration, determine injection efficiency or establish an
energy cutoff. Failed/sub-fast/non-forward local shock records do not acquire
an invented DSA source. Preserve the core's numerical rejection and do not
floor compression or cap the spectral index to force a sample.

Two independent application selectors are implemented in the maintained
parser. Their presence does not close the stage-6 complete-application gate:

| Selector | Requested alternatives | Default/status |
| --- | --- | --- |
| `shock-return-policy` | `delete`, `pass` | **`delete`**, explicitly selected by the user |
| `shock-source-normalization` | `gross-emission`, `net-upstream-escape` | **`net-upstream-escape`**, accepted for the first production path; the alternative retains its distinct physical meaning and qualification requirements |

`delete` removes a returning particle at its first resolved upstream-to-shock
encounter and records represented number, birth energy and return-time energy
separately. A Parker Brownian path can hit the moving shock and return upstream
within the same numerical interval: endpoint signs alone cannot implement this
boundary. The physical AMPS step and sampling schedule remain unchanged.

The user resolved the `pass` alternative: it is recognized but **not
implemented**. Selecting it must terminate at the native configuration
preflight, before injection or particle stepping, through the AMPS host trap:

```cpp
exit(__LINE__,__FILE__,
    "Shock return policy 'pass' is not implemented yet for reduced SSE/Parker transport. "
    "The provider has no downstream transport volume. Select 'delete'.");
```

Neutral physics/configuration helpers retain typed unsupported status; the
native application boundary owns this host-specific `exit` call. It must not
silently fall back to deletion, continue in ambient fields, wait until the
first particle crosses, or reinterpret local RH limits as a downstream volume.
The requested default remains `delete`.

The normalization alternatives describe different physical populations, not a
rescaling identity. Gross emission and net upstream first passage need their
own declared dimensional rate/efficiency and release interpretation. Preserve
separate candidate, committed, returned and surviving-number/energy ledgers;
never renormalize away losses to obtain the requested net source. In particular,
the existing coronal-CME schema-5 production source accepts net first passage
at a finite reference surface. A new reduced-application gross mode must not
silently relax that existing core gate or claim its source qualification.

The question of **shock-frame isotropy** was explicitly raised for verification.
Standard DSA assumes near-isotropy in the local scattering-centre frame;
neglecting scattering-centre drift gives the local plasma frame. It does not
imply isotropy in the shock frame. The [Drury derivation, sections 2.1 and 2.3](https://www.oa.uj.edu.pl/~mio/Ast-Wys-En/Literatura/drury.pdf)
uses plasma-frame momentum with shock-frame spatial coordinates; an explicit
boost changes angular dependence. [Summerlin and Baring, figure 5 and appendix](https://arxiv.org/pdf/1110.5968)
also distinguish fluid-isotropic populations from their shock-frame angular
and crossing distributions. For an oblique relative boost, representability
by the plasma-frame gyrotropic mover must be established; transforming
individual four-vectors alone does not justify discarding gyrophase dependence.
The maintained [core frame contract](../src/models/sep_coronal_cme/model/physics.md)
requires that check. The user has now selected **upstream-plasma-frame
isotropy** for the parent population; shock-frame isotropy is not the selected
source prescription. Together with the already accepted DSA spectral choice,
this settles the parent population's momentum/angular frame.

For particles crossing from an isotropic reservoir into the upstream transport
domain, the crossing distribution is positive-normal-flux weighted, not uniform
in pitch cosine. Reuse the existing core focused-escape machinery with the
actual signed front-relative guiding-centre speed and magnetic polarity.
Physical U enters this relative source flux; it is still absent from relative
particle displacement on a plasma-comoving field line.

Finally, the default absorbing boundary makes release location a physical
input. A Parker diffusion process born exactly on an absorbing shock is killed
immediately. A finite upstream reference distance or qualified source/boundary
operator is therefore required; an arbitrary mesh-dependent epsilon is not a
physical closure. The user accepted a finite upstream reference surface. Its
positive physical distance is an explicit run input and must stay fixed during
numerical refinement.

The approved first production path reuses the existing **finite upstream
reference first-passage** source, with its dimensional
rate defined at that release surface and later returns absorbed at the shock.
For a 1-D consumer, release points must be intersections of the same material
field line with the upstream offset surface at the actual birth epoch.
A Cartesian displacement along the shock normal generally leaves an oblique
magnetic field line and must not be used as its particle placement rule.
The offset's signed normal distance `L_ref` is a required positive SI input;
its numerical value is not inferred from a cell size, a fraction of the AMPS
timestep, or a fixed diffusion-depth Péclet number. Geometrical admissibility
and diffusion/kinetic applicability are checked independently. A calibrated
physical range should be propagated as source uncertainty while holding each
chosen `L_ref` fixed under numerical refinement.

**The source/release interpretation and net-first-passage default are now
accepted; no numerical production default for `L_ref` is assigned.** Particle-
enabled runs must supply a finite, strictly positive SI release distance and
pass geometrical/domain checks before allocation. Missing or invalid release
input must produce a clear configuration error rather than a guessed offset.
An uncalibrated test value can serve only a named verification/sensitivity
case. This required runtime input does not prevent implementing and testing
the agreed baseline with explicitly declared manufactured test distances.

The two normalization definitions remain separate choices with separate
physical meanings; agreement on the net-first-passage baseline does not supply
a gross-to-net escape conversion or qualify the alternative automatically.
The switches and particle-enabled source are implemented for stages 5–6.
The bounded numerical/native checks and complete three-mover stage-6 campaign
pass. See the maintained
[qualification status](validation/STANDALONE_INTEGRATION_STATUS.md) for actual
executions and limitations. This decision record is not itself qualification.

## Stage 7 — Qualify parallel execution, sampling and restart

**COMPLETE (2026-11-22)** for the prescribed reduced SSE/Parker material
profile: all 44 checkpoint/thread, 58 material restart/repartition and 48 native
real-MPI jobs pass. All three movers qualify same-decomposition continuation,
1↔8 ownership changes and two-/four-worker restart; Dmumu additionally passes
on eight MPI ranks with four workers per rank. Late rank-local preparation
failure preserves complete live state and permits retry. Position refresh,
worker timestep ownership, sampling addition order and complete material-record
serialization were repaired without changing physical dt, acceptance bounds or
reference solutions. Exact final images, actual process receipts and the frozen
build closure are recorded in
[the maintained status report](validation/STANDALONE_INTEGRATION_STATUS.md).
The subsequent stages 8–9 have their own scoped qualification record below.

### Work

Complete the production state lifecycle for changing geometry, including
ownership tables, migration, source counters, stochastic histories and output
continuation.

Moving-geometry checkpoints must retain topology, material labels, background
identity and segment exposure. Existing fixed-geometry checkpoint
qualification is insufficient.

### Tests

- One, two, four and the production-required MPI rank count.
- Multiple supported thread counts.
- Migration during geometry evolution.
- Rank-local preparation failure without partial publication.
- Repartitioned ownership without source duplication.
- New and old segments sampling the same density with different ages.
- Measured vacuum distinguished from no measurement.
- Continuous execution versus separate-process checkpoint/resume.
- Output cadence and exposure denominators unchanged.
- Complete represented-number and signed energy-ledger comparisons.

Require exact identity where the algorithm promises determinism; use
independently justified statistical criteria for stochastic observables.

**Exit gate:** real native MPI and restart campaigns pass. A serial MPI
dependency adapter cannot satisfy this gate.

## Stage 8 — Admit the other existing core backends

**COMPLETE (2026-11-22) for the approved admission scope.** SBACK01–05
pass strict and ASan/UBSan checks; SBACK06–08 pass actual two-rank preflight
for all three movers. The maintained diagnostic/exchange interfaces remain
available. Rejected backends gain no propagation qualification from these
passing rejection tests.

**Authorized scope clarification, 2026-11-22:** the user selected “Reject
unsupported modes; qualify existing supported transport.” The executable
SWCME1D, neutral-bundle and full-regional interfaces do not currently establish
all required material contracts. Preserve their diagnostic capabilities, add
explicit admission/rejection checks and publish the exact
[support matrix](validation/STANDALONE_BACKEND_SUPPORT.md). This scope does
not authorize an invented core plasma evolution law or qualify rejected
transport. Implementation and current regression receipts remain distinct.

### Work

Reuse the completed consumer machinery through each backend's maintained route:

- **SWCME1D:** extend its existing adapter for the selected supported profile.
- **Coronal-CME:** use its versioned neutral field-line exchange interface.
- **Full corona–SWCME:** consume only capabilities whose core qualification
  is complete.

A single-epoch field-line bundle does not establish moving-background support.
Required material identities, temporal coverage and state reconstruction must
be available explicitly. Backend selection must not change mover physics,
normalization, units or boundary policy implicitly.

### Tests

For every admitted mode:

- Direct-core versus native-consumer equivalence.
- Configuration, frame, composition and epoch checks.
- Field-line identity and temporal-coverage checks.
- Shared controlled-input comparisons.
- Native transport and restart campaign.
- Rejection of unimplemented regional or thermal capabilities.

**Exit gate:** a published support matrix identifies the exact
model/profile/mover combinations qualified for production.

A passing rejection test qualifies rejection behavior; it does not qualify
transport with that unsupported model.

## Stage 9 — Register and run the complete regression campaign

**COMPLETE (2026-11-22) for the supported reduced SSE/Parker scope.** The
fresh 150-job native/restart/thread campaign, three 400-step applications,
eight admission/compatibility gates, four coupling builds and nine coupling
comparisons pass. Physical dt stays 300 s; the 1e-10 budget bound, original
references and byte-exact continuation gates are unchanged. Eleven original
application-verifier controls, eleven maintained launcher controls and 45
orchestrator checks pass. Exact commands and source/build identities are in
[the qualification report](validation/STANDALONE_INTEGRATION_STATUS.md).
The broad core/geometry scope is 341 PASS / 181 SKIP / one ERROR: the
full-regional refinement's unchanged 180-s timeout remains unqualified and
blocks its dependent backend. This scoped completion does not certify the
whole shared registry or any rejected moving-background model.

The new maintained `test/run_standalone_application.py` executes the original
SAPP01–03 complete-window gate; previously only its external-evidence verifier
was registered as a Make audit. Shared suites now include the full application
and SBACK01–08 admission cases, with matching application-orchestrator readers.
The unchanged dt/window, source, reference and budget contracts and all other
required campaign roles are recorded in
[STAGES8_9_TEST_CONTRACT.md](validation/STAGES8_9_TEST_CONTRACT.md).

Every new test must be reachable through the existing validation procedure.
New `srcSEP` runners must be included in the shared orchestrator with
individual result readers, following [test/README.md](test/README.md).

Use:

- The shared registry in [src/test/manifest.json](../src/test/manifest.json).
- The complete application procedure documented in
  [srcSEP/test/README.md](test/README.md).
- Existing core-model suites and geometry tests.
- Existing AMPS–BATS-R-US coupling regressions to detect shared-code damage:
  `AMPS-BATSRUS-Europa`, `PS-BATSRUS-OH`, `MPS-BATSRUS-PIC-coupling`, and
  `AMPS-BATSRUS_dynamic_fluxrope`.

Each stage must produce JSON/JUnit results, per-test logs, exact commands,
source/configuration hashes, seeds, metrics, acceptance criteria and reference
provenance. Report `PASS`, `FAIL`, `SKIP` and `ERROR` distinctly. Missing output,
crashes and MPI startup failures cannot be reported as physics passes.

Add isolated negative controls demonstrating detection of wrong density
mapping, stale epochs, extra bulk displacement, unsupported shock roots and
global sampling denominators. Independent references must not be generated
from the production routine being tested.

**Final acceptance:** all required supported-mode gates pass; every limitation
remains explicit; no reference is regenerated merely to accommodate changed
output. Preserve the original physical and numerical acceptance criteria.

## Execution order and remaining physical decisions

Stages 0–4 have the bounded reduced-core implementation and numerical/native
evidence described in the stage report. The subsequent stages 5–6 authorization
adds the approved finite-reference net source and absorbing return boundary to
that same core authority. Their numerical, native MPI and complete one-rank
application gates now pass for the declared reduced case; implementation alone
was not used to close them.
Additional existing backends enter after that path is qualified.
Full CME-volume transport remains conditional on the corresponding core
capabilities; it is not hidden inside adapter work.

Before dependent implementation, resolve any missing selections for the core
ambient/profile, compatible source prescription, release side and particle
boundary rule at the unmodelled downstream region. Use the selected core
profile's own thermal/composition definitions wherever available; do not
introduce a competing closure. Any additional acceleration or downstream
capability must come from an existing qualified implementation or remain an
explicit prerequisite outside this integration task.

## Maintained interfaces and supporting documentation

- [Reduced shock-front provider](../src/models/sep_corona_swcme/shock_front/provider.h)
- [Reduced shock-front physics contract](../src/models/sep_corona_swcme/docs/shock-front/model.md)
- [Corona–SWCME implementation and capability status](../src/models/sep_corona_swcme/README.md)
- [SWCME model interfaces](../src/models/swcme/README.md)
- [Coronal-CME maintained architecture](../src/models/sep_coronal_cme/README.md)
- [Reduced-background native adapter](adapters/reduced_shock_background_adapter.h)
- [SWCME1D adapter](adapters/swcme1d_adapter.h)
- [Neutral field-line bundle adapter](adapters/field_line_bundle_adapter.h)
- [Material geometry implementation contract](../src/pic/README.md)
- [Material geometry design](../src/pic/MOVING_FIELD_LINE_GEOMETRY_DESIGN.md)
- [Transport numerics](TRANSPORT_NUMERICS.md)
- [Application development requirements](AGENTS.md)


## Extension: complete and integrate the remaining shock/CME backends

**Added 2026-11-22, after commit 80466e6. Status: PLANNED.** This extension
defines the remaining development and integration work for SWCME1D, SWCME3D,
the Coronal-CME ellipsoidal front/piston and the full corona–SWCME regional
model. It does not change the completed, bounded stages 0–9 or mark any
additional backend as implemented or qualified. A later implementation request
must name its stage or bounded work package.

The intended deliverable is a separately identified, physically complete
background profile for each requested backend, connected to the three existing
field-line movers and their source, sampling, ownership and restart lifecycle.
A diagnostic API, linked archive, geometrical shock surface or passing rejection
test is insufficient. Completion is declared per profile and supported region,
rather than for every option in a model family at once.

### Existing ownership and mandatory dependencies

| Authority | Reuse and remaining work |
| --- | --- |
| [SWCME core](../src/models/swcme/README.md) | Preserve canonical 1-D/3-D kinematics, geometry, MHD jump and checked evaluators. Complete a conservative material background and its thermodynamics; retain the historical empirical/phenomenological modes as explicitly identified diagnostics. |
| [Coronal-CME core](../src/models/sep_coronal_cme/README.md) | Reuse fixed-orientation ellipsoid/piston geometry, smooth/tabulated histories, PFSS/open/closed plasma, local shocks and source kernels. Complete event assembly and the selected live consumer contract. |
| [Corona–SWCME composite](../src/models/sep_corona_swcme/README.md) | Own the common regional sheath/ejecta/contact state, compatible material evolution, continuous handoff and event lineage. A surface-local RH solution cannot replace a volume closure. |
| [Neutral SEP contracts](../src/models/sep_common/README.sep_coronal_cme.md) | Own provider-independent exchange, status, identity and units. Extend records only for actual required capabilities, with versioned compatibility and transactional import. |
| Native consumers | Keep PIC/MPI storage, clocks, geometry publication, particles and sampling in thin application adapters. Three-dimensional shock geometry does not require three-dimensional particle motion in srcSEP. |

The existing
[background-first launch roadmap](../src/models/sep_coronal_cme/docs/CME_SURFACE_LAUNCH_ROADMAP.md)
controls the coronal/composite order:

1. **G1:** complete spatial background in srcSEP3D; actual zero-particle
   native launch, handoff and 1-AU evolution.
2. **G2:** the same qualified background in srcSEP; actual evolving
   field-line publication, still zero particles.
3. **G3:** nonzero-particle source/transport qualification in srcSEP3D.
4. **G4:** nonzero-particle source/transport qualification in srcSEP, with
   comparisons restricted to common physical limits.

G1 precedes G2, G2 precedes G3, and G3 precedes G4 for these routes. Existing
qualified reduced-SSE tests remain active throughout. Independent SWCME-only
profiles may progress through their own background and srcSEP gates without
claiming completion of coronal G1–G4. Work on srcSEP3D here is an explicit
cross-application prerequisite, not permission to add Cartesian particle movers
to srcSEP.

The old bundle-only architecture and the later launch roadmap's live-provider
target must be reconciled in the maintained architecture contract and its
tests before changing linkage. Preserve snapshot replay as a separate,
explicit capability. Do not bypass an architecture assertion or relabel an old
single-epoch bundle to make live evolution appear available.

### Contracts preserved throughout stages 10–21

- AMPS physical dt is fixed at initialization and remains unchanged throughout
  stepping and sampling. Retain the existing 300-s/400-step campaign and its
  output cadence. Background ODE refinement, event roots and internal numerical
  subdivisions must not shorten an AMPS step or add sampling opportunities.
  Any proposed different fixed-dt particle campaign requires an explicit
  decision before use; it is not an automatic response to a failed test.
- Every particle step consumes an immutable owning realization of the complete
  physical interval. Dense-time geometry, background and source queries share
  its configuration, frame, identities and coverage; they do not mutate a live
  provider or mix endpoint fields with departure-time geometry.
- Particle displacement along every srcSEP field line is relative to the
  plasma-comoving geometry and contains no extra solar-wind advection term.
  Physical U and all nine Cartesian components of grad(U) remain inputs to
  line motion, shock inflow, compression, strain and work.
- Move both endpoint vertices. Trace inserted inner vertices along the actual
  magnetic field, not minus U. New segments start empty with zero exposure;
  existing source rules populate them. Remove an outer segment only when
  completely outside the declared radial domain, account for attached
  particles once, and retain segment-specific sampling exposure.
- Use one authoritative event/background and one declared shock authority per
  profile. Composition, species density basis, EOS, thermal law, signed B,
  physical frames, source interpretation and region classification are shared
  by geometry, shock, source and movers.
- Retain the approved accelerated-population source contract where applicable:
  local-compression DSA momentum spectrum, parent isotropy in the upstream
  plasma frame, equation-appropriate angular/flux projection, an explicit
  finite upstream release reference, net upstream escape and absorbing returns.
  Shock-frame flux calculations require the declared transformation; they do
  not make the parent distribution shock-frame isotropic.
- The default shock-return policy remains deletion. The pass option continues
  to terminate through the existing native error handler until separately
  specified downstream particle physics and its tests are implemented.
  Completing a downstream background volume does not by itself enable passage
  or authorize an additional acceleration mechanism.
- Prescribed read-only turbulence remains the first integration target.
  Local wave evolution, feedback, gross-emission conversion, additional shock
  families and observational qualification require their own completed
  contracts; none is supplied by a new backend adapter.
- Keep existing references, physical bounds, acceptance criteria and test IDs.
  Distinguish a physical no-shock/no-connection result from numerical failure.
  No density/pressure/B floor, hidden ambient fallback, lost support area or
  changed runtime ceiling may turn an incomplete model into an admitted one.

### Stage 10 — Rebaseline every backend and freeze the capability matrix

**Recorded status:** baseline audit complete; the exact baseline and open decisions are
in [STAGE10_BACKEND_CAPABILITY_AUDIT.md](validation/STAGE10_BACKEND_CAPABILITY_AUDIT.md).
The retained full-regional timeout is not a passing closure gate.

**Dependencies:** completed stage 9; no new model admission.

**Work.** Inventory actual production APIs, parsers, compiled objects, model
tests, native entry points and current failures at the new starting revision.
Classify each item as implemented, tested in isolation, wired natively or
qualified in production. Preserve failed fixtures and source/build/input
identities. Reconcile historical documentation with executable code, including
the neutral bundle schema-1/schema-3 discrepancy and the revised regional
closure selection. Historical Mars1/Mars2 paths and old test totals are
provenance, not assumptions about the active checkout.

Define individual rows for SWCME1D SHOCK_ONLY/FULL_ICME, SWCME3D
sphere/ellipsoid/finite-SSE, the coronal ellipsoid/piston, and each intended
composite region/profile. Record geometry, evolution law, regions, density
basis, thermodynamics, derivatives, source, return boundary and restart support.
A diagnostic FULL_ICME mode must not acquire a production label by aliasing it.

**Tests and exit gate.** Run the exact existing core catalogs and bounded
application baselines; record PASS/FAIL/SKIP/ERROR separately. Retain
CSWC-FULL-CSWC0628-MASS-TIME and the original CSWC0628 physical refinement
failures until their actual causes and independent gates are resolved. Freeze
new quantitative acceptance and performance criteria before examining new
candidate residuals. Publish the matrix and a concrete decision record for
every unresolved closure. Stage 10 cannot close by counting dormant tests as
executions.

### Stage 11 — Complete the common material-background interface

**Recorded status:** complete for the supported reduced SSE/Parker profile.
The interface and reduced native binding are implemented,
including explicit initial geometry-generation presence and qualified-region
validation at every new trajectory API evaluation. The positive legacy
trajectory, physical timestep and RK controls are unchanged. Final matching
strict/sanitizer, shared-registry, 48 native MPI and repeated two-rank CLI
checks pass, as do all 58 material restart, 44 checkpoint/threaded and three
400-step application gates. The legacy heterogeneous native startup failed its
material-owner prerequisite. The integrated routine driver selects each case's
reviewed fixture and retains every original gate. The maintained `make test`
passes all 222 SEP Routine cases, all source prerequisites and all four embedded
SWCME Routine cases. The controlled analytical target passes all 18 comparisons.
Focused results and current execution evidence are recorded in
[STAGE11_MATERIAL_BACKGROUND_CONTRACT.md](validation/STAGE11_MATERIAL_BACKGROUND_CONTRACT.md).
No later model profile is admitted by this work.

**Dependencies:** stage 10 and the relevant approved core physics contracts.

**Work.** Define the smallest shared capability record required by the
consumers: profile/revision identity, immutable physical interval and dense-time
coverage, material map/trajectory and inverse where needed, region/side,
rho, species number densities, pressure/temperatures, U, B, grad(U), required
magnetic derivatives, EOS, shock/contact records and status. State SI units
and the coordinate/velocity frames in each public field.

Separate a snapshot-exchange capability from material evolution and from
physical downstream-volume support. A stored divergence alone is insufficient
for field-aligned strain; density history is not a substitute for grad(U).
Smooth derivatives must stay on their physical branch. Discontinuity crossings
need one-sided states and an explicit interface contract.

Extend admission from the current coarse backend inventory to exact qualified
profile capabilities. Preserve existing rejected empirical/diagnostic profiles
and unknown-profile rejection. Do not change RequireMovingTransport for an
entire backend family to unconditional success. Reuse the existing CLI,
configuration fingerprint and publication transaction; do not add a competing
parser, clock or mutable model cache.

**Tests and exit gate.** Independent species/EOS/frame/derivative checks,
canonical identical-point comparisons, missing capability/coverage/side
negatives, foreign fingerprints and stale generations, failed candidate
rollback, public-header/ABI checks and versioned serialization rejection.
Every required quantity must have an identified owner and tested availability;
an unimplemented field cannot be initialized to a plausible default.

### Stage 12 — Complete SWCME1D material evolution and thermal export

**Recorded decision:** the first `SHOCK_ONLY` material candidate uses a prescribed
constant radial wind, inverse-square mass density, ideal-induction Parker winding
and an explicit sustaining force and work term. Its controlled verification case
uses provider coverage of 0.30–1.20 AU, adiabatic expansion with `gamma = 5/3`,
and explicit electron/proton reference values at 1 AU. The approved parameters
and analytical reference values are in [D10-1](validation/STAGE10_BACKEND_CAPABILITY_AUDIT.md#decision-records-and-frozen-acceptance).

**Upstream implementation.** The canonical `swcme::material1d::Event` supplies
the prescribed adiabatic profile, explicit sustaining force and work, full
Cartesian gradients, fixed-step material map and local fast-shock limits.
Its detached SEP adapter exposes the canonical fields and species. A distinct
profile identity preserves both original diagnostic modes. Native source,
geometry and restart ownership are still required before admission. The
ten independent `MAT1D` cases and the neutral parity checks exercise this
upstream scope. See [the SWCME API and equations](../src/models/swcme/README.md#prescribed-swcme1d-material-upstream).

**Dependencies:** stages 10–11 and the complete upstream input and thermal
contract for the selected profile. `FULL_ICME` additionally requires its own
regional material and contact closure.

**Work.** Preserve the legacy diagnostic law
n = C2/r² + C4/r⁴ + C6/r⁶ with constant radial U. Its nonzero div(rho U)
must remain visible. Use the approved alternative in a new versioned profile.
The canonical core owns the coupled flow, density, induction, thermal state
and sustaining force. Adapters export that same state and derivatives.

Implement the approved versioned profile in the SWCME core. Complete the
material trajectory, its temporal support, radial derivatives and composition
mapping. For FULL_ICME, qualify the intended sheath/ejecta/transition state
rather than exporting an upstream temperature throughout the disturbed volume.
Extend the 1-D adapter to export the actual regional state and derivatives.
Maintain the existing ballistic/DBM/data-driven kinematics and checked status
semantics; each admitted combination needs explicit coverage.

**Tests and exit gate.** Independent radial mass inventory/flux and material
continuity, magnetic flux/induction, thermal/work and regional boundary tests;
direct-core versus adapter agreement for density basis, signed B, U,
pressure/temperatures and grad(U); ordered-time and no-extrapolation checks.
Exercise SHOCK_ONLY and FULL_ICME separately. A corrected upstream profile
does not qualify FULL_ICME. Admit only the exact profile whose complete
background contract passes.

### Stage 13 — Complete SWCME3D geometry and conservative material state

**Dependencies:** stages 10–11 and stage 12's shared laws wherever reused;
the 3-D case retains its own independent qualification.

**Work.** Reuse canonical sphere, self-similar ellipsoid and finite-SSE geometry,
apex histories and local ideal-MHD shock evaluation. Complete three-component
material motion, all required spatial derivatives, composition and thermal
export for each selected profile. Shared 1-D/3-D laws must share their core
implementation rather than evolve independently.

A radial density/velocity correction does not establish 3-D induction,
latitudinal consistency or oblique downstream evolution. Qualify those
explicitly. Remove a consumer's sphere-only restriction only after its actual
nonspherical support, interpolation, side handling and derivatives pass.
Keep a finite shock cap distinct from a closed ejecta body.

Build a thin native field-line consumer of the 3-D core state. Keep native
geometry, tube metrics and particle coordinates in the existing field-line
machinery. A linked swcme3d object or a Cartesian mesh provider is useful reuse,
not evidence that this material-line consumer already exists.

**Tests and exit gate.** Exact spherical and matched 1-D limits; independently
reconstructed ellipsoid/SSE normals, local normal speeds and surface areas;
tilted axes, both magnetic polarities, finite-support misses, multiple roots,
tangencies and weak/sub-fast classification. Check Cartesian derivative
convergence and mass/flux/thermal budgets on independently deformed material
volumes. Record shape/profile-specific support; a sphere PASS cannot admit all
three geometries.

### Stage 14 — Complete Coronal-CME ellipsoid/piston event assembly

**Dependencies:** stages 10–11 and the maintained coronal/composite contracts.
This is shared-core background work; native srcSEP integration still awaits G1.

**Core implementation.** The coronal library materializes independent front
and piston histories and publishes immutable background epochs only after
geometry, support, ambient, and shock checks succeed. Solar cut-cell
quadrature integrates exposed area even when the original cell midpoint is
buried. The composite library provides the retained PFSS/open/closed/Parker
ambient owner, including signed sectors, physical source-surface side
selection, composition, and branch-preserving derivatives. `S14BG01`–`10`
exercise the assembly and its transactional failure paths. G1 regional and
native SEP qualification remain separate. See [the coronal core](../src/models/sep_coronal_cme/README.md).

**Work.** Assemble the existing fixed-orientation triaxial front and independently
configured piston from the canonical input histories. Preserve center/apex
equivalence, independent semiaxis rates, analytic normal velocities, stable
patch identities and solar clipping. Qualify attachment/detachment, field/plasma
support and source-eligible radial support as separate concepts.

Bind the actual PFSS/open/closed/Parker ambient, composition/EOS, local oblique
MHD jump and complete histories to one event. Preserve delayed shock formation:
a geometric sub-fast patch carries no invented compression or particle source.
Keep dynamic attitude/deflection and research-only histories unavailable unless
their own maintained production contracts and tests are completed.

Deliver coherent background-only core epochs without requiring source plans,
particle records, RNG state or wave feedback. Complete source-surface and
one-sided topology handling using canonical physics. Integrate regional state
only through stages 15–16; an ambient-plus-front result remains intermediate.

**Tests and exit gate.** Independent triaxial normal/speed/area and clipping
references, smooth/Hermite history extrema between knots, launch support,
patch identity, delayed activation and weak-shock RH checks. Check signed
PFSS sectors, closed/open boundaries, EOS and one-sided derivatives.
Transactional failure must preserve the prior complete epoch. This stage
qualifies the assembled components, not the missing regional volume or G1.

### Stage 15 — Complete the full regional sheath and material closure

**Dependencies:** stages 10–11 and 14; the actual selected regional closure.

**Work.** Follow the later
[conservative cross-ray scope](../src/models/sep_corona_swcme/docs/CONSERVATIVE_CROSS_RAY_TRANSPORT.md)
and
[analytical material-map design](../src/models/sep_corona_swcme/docs/ANALYTICAL_LATERAL_TRANSPORT_DESIGN.md).
Retain the three-dimensional ambient and analytical/semi-empirical evolution
requirement; do not revive fixed-ray radial dynamics as a complete 3-D closure
or introduce a global MHD runtime. Existing HLL, static CT, Cauchy and piston
kernels remain reusable references, not an already selected production solver.

First qualify the actual ambient reference flow, inverse and variational
derivatives on supported branches. Then test whether a proposed deformation
can represent the required oblique RH downstream velocity, field and contact
motion. A scalar deformation that fails representability must reject; do not
drop transverse components, discard required rays or fit a mass source.
A richer vector map/force/heating law requires a concrete physics decision.

Complete one common curved geometry, compatible material labels, shock
admission times, initial inventory and inverse map. For ideal source-free
cohorts, enforce rho = rho0/J and B = F B0/J using the same F = dX/da,
positive J and U = dX/dt. Thermal evolution must follow the selected EOS and
heat/work law. These identities alone do not prove momentum balance, correct
shock entropy or a qualified analytical plasma solution.

Record one physical shock authority per profile and retain the computed-shock
requirement wherever the selected contact-driven closure specifies it. A
prescribed-front alternative requires a distinct approved profile; it cannot
close a computed-shock gate. Do not independently impose incompatible front,
contact and downstream states. Preserve a specified zero-added-heat material
law where applicable, while accounting separately for any sustaining heat/work
required by an isothermal ambient. Account for all physical influx, contact/lateral/
solar/outer losses and sustaining force/heating/work. Use local absolute
contact leakage as well as signed totals so cancellation cannot hide penetration.

**Tests and exit gate.** Independent planar and spherical limits, oblique and
quasi-perpendicular shocks, startup inventory and admission accounting,
positive/invertible maps, signed magnetic flux, induction, entropy and
force/work residuals. Require mass/time/angular refinement and the original
minimum-J, leakage and coverage gates. Fix the existing refinement timeout
through a measured algorithmic correction while preserving its physical
assertions and original ceiling. No local kernel PASS closes BG3D-4.

### Stage 16 — Complete ejecta, regional joins and continuous handoff

**Dependencies:** stages 14–15 and the existing BG3D-5/6 contracts.

**Work.** Supply the separate ejecta reference state, material map, composition,
thermal state and magnetic topology required by the selected full profile.
Use the actual checksummed reference/vector-potential assets where specified.
Do not infer a spheromak, an empty initial body or quiet ambient behind a shock.
Reconcile any older contact-fraction geometry with the selected common material
contact; independent interface definitions cannot coexist as one physical body.

Qualify normal velocity, normal magnetic field, normal traction and permitted
tangential discontinuities at the sheath/ejecta contact. Declare every covered
wake, transition and post-event state. A smoothing layer needs a compatible
mass/flux/energy construction; component-wise B blending is not sufficient.

Transfer the complete surface, velocity and regional material state into the
outer SWCME continuation. Match the full geometry and required derivatives,
not just apex radius/speed. Differentiate a single transition construction,
including time-dependent weight derivatives. Preserve event identity,
inventory, magnetic flux and work across the handoff without resetting clocks,
source history or material labels.

**Tests and exit gate.** Independent ejecta mass/flux/EOS checks, contact
traction/slip and topology tests, matched direct-handoff limits and deliberately
incompatible equal-apex histories. Refine the transition and check positive
maps and bounded temporal/spatial derivatives throughout it. Missing regions,
foreign assets, unsupported topology and out-of-coverage queries must reject.
BG3D-5/6 close only with complete qualified regional coverage.

### Stage 17 — Qualify native three-dimensional background evolution: G1

**Dependencies:** completed shared coronal/full-regional stages 14–16.
SWCME-only backgrounds retain their separately declared scope.

**Work.** Install the complete shared provider in an independently configured
srcSEP3D background-only application. Stage actual current/previous native
fields, derivative/category data, owner values and received ghost values.
Publish one common epoch only after collective readiness. Add background
restart for material state, handoff and integrated evolution as required.
Retain zero physical particle allocation and zero source counts.

Deliver a short pre/inside/post-handoff test and a separate low-corona-to-1-AU
case covering the actual event duration and every declared region. A sampled
trajectory reaching 1 AU is not a native propagation receipt. Do not reduce the
new full background gate to the existing outer-only 400-step particle fixture.

**Tests and exit gate.** Actual one-/four-rank evolving native runs, independent
spatial probes in ambient/sheath/ejecta/transition regions, owner/received-ghost
agreement, stale-slot and rank-local failure controls, background restart and
mesh/time/surface/material convergence. Quantify physical residuals separately
from discretization error. G1 passes only when all required regions and the
actual launch/handoff/1-AU sequence pass without particles.

### Stage 18 — Integrate and qualify native moving field-line backgrounds: G2

**Dependencies:** stages 11–13 for a selected SWCME-only profile; G1 for the
coronal/composite routes. No source or particle milestone is inferred.

**Work.** Add thin profile-specific adapters using the already-qualified core
state. For the coronal/composite live route, update the authoritative linkage
and architecture tests together and retain versioned bundle replay explicitly.
A replay capable of material motion needs actual temporal primitives,
trajectories/lineage and derivative coverage; old bundles remain snapshots.

Connect native vertex trajectories, B tracing, flux-tube metrics, inner insertion
and outer trimming to the same immutable interval. Preserve one-sided region/
sector information, both polarities and all physical derivatives. Certify
coverage for inserted vertices and all numerical stencils, including a complete
tracing increment below the nominal inner sampling bound. Topology-changing
events require their own physical mapping; otherwise reject rather than join
unrelated field lines.

**Tests and exit gate.** Independently built srcSEP, real one-/four-rank
multi-epoch zero-particle runs, identical-point canonical comparisons with G1,
B/tangent and flux-tube-area convergence, material continuity, empty insertion,
segment exposure and retired-segment accounting. Exercise temporal coverage,
late-rank rollback/retry, restart and an actual background-only 1-AU case for
the coronal/composite profiles. The moving-line representation itself must
pass before any corresponding particle motion is enabled.

### Stage 19 — Close the required three-dimensional particle prerequisite: G3

**Dependencies:** G1 and G2 for coronal/composite profiles.

**Work.** Complete the existing srcSEP3D source/allocation/transport milestone
through that application's maintained interfaces. Freeze the source's physical
meaning, units, finite release support, spectrum, scattering prescription,
supported regions and boundaries against the qualified background. Preserve
the single acceleration mechanism: an already accelerated injected population
must not silently acquire a second resolved shock-acceleration model.

**Tests and exit gate.** Actual nonzero births and transport, independent
number/energy/source/return oracles, domain and shock contact checks,
time/space/statistical convergence under the approved fixed-step campaign
contract, real MPI/repartition and fresh restart. A wave-feedback claim requires
its separate physical gate. Only a complete G3 receipt unlocks new srcSEP
particle work for these routes. This stage introduces no Cartesian mover into
srcSEP and makes no claim that general 3-D transport equals field-line transport.

### Stage 20 — Complete standalone sources, all three movers and restart: G4

**Dependencies:** stage 18 for each SWCME-only profile; G3 additionally for
coronal/composite profiles. Physical source mapping must be explicit.

**Work.** Reuse the approved DSA/finite-release/net-escape/absorbing-return
contract through a neutral source interface. Determine physical source rates
and finite tube/front overlap from the selected backend's actual units and
geometry. SWCME's relative-only source weight is dimensionless; it cannot be
treated as an absolute injection rate or silently converted to the finite-
reference net source. Freeze that conversion and its physical input before
native allocation.

Find every supported moving front/line contact at its actual event time.
Retain partial footprints, multiple roots and tangencies with explicit
ownership/support decisions. New particles attach to the birth-time material
coordinate and use only their physical birth remainder within the unchanged
host step. Do not apply additional wind displacement or additional sampling.
Use consistent shock/source histories across Parker, FTE-D_mumu and FTE-MFP;
the equation-specific angular projection must preserve physical normalization.

Complete geometry/particle migration, source clocks, immutable birth history,
keyed RNG, exposure and sampling-buffer phase in a versioned checkpoint.
Preserve old valid reduced-profile restart behavior and reject older images
which cannot represent a new profile's required state.

**Tests and exit gate.** Independent dimensional source and finite-tube flux
checks, delayed activation, first-return absorption, causal births and correct
frames. Exercise ballistic/focusing/compression and the common diffusion limit
with independent oracles for all three movers. Run actual MPI 1/2/4/8,
two-/four-worker continuation, combined eight-rank/four-worker cases, 1↔8
repartition and late-rank failure/retry. Compare complete promised-deterministic
images exactly and stochastic observables using frozen independent criteria.

Each admitted backend/profile must also execute the complete 400-step,
300-s application output gate and the event-specific launch/handoff duration
required by its declared scope. Independently reconstruct every surviving
segment's exposure and signed number/birth/current-energy ledgers. Keep the
existing 1e-10 application budget limit. A short run, an empty population or a
matching source-only fixture does not close G4.

### Stage 21 — Make every admitted profile an automatic regression gate

**Dependencies:** the exact profile's completed core and native gates above.

**Work.** Register new stable positive tests in the core manifest, catalogs,
traceability, source/Make entry points and the complete SEP orchestrator with
individual readers. Preserve existing SBACK rejection IDs for the legacy
unsupported configurations; add distinct positive cases for new qualified
profiles. New stage numbers in this document are work packages, not registered
test IDs or evidence of implementation.

Supply profile-specific physical decks and independent output verifiers while
retaining the original reduced SAPP gate unchanged. A new verifier must check
its own backend/region identities and physical source/sampling expectations;
changing a label on reduced-model output cannot qualify another backend.

The scheduled nightly SEP runner must select every required new profile without
manual --only flags. Provision and audit the native profile and two-/four-worker
build prerequisites through the routine; a configured test with a missing build
is ERROR, not PASS. Keep long campaigns explicit and automatically selected by
the full nightly path; --quick or a core smoke/extended selection alone cannot
claim complete nightly coverage. Require a selection/reader control which
fails if a profile, mover, role, process result or output window disappears.

Measure aggregate build/run cost and provision sufficient scheduled capacity.
If the expanded campaign cannot fit the current nightly runner, add an explicit
scheduled full-campaign runner and a coverage guard. Do not shorten the physical
window, hide an omitted profile or relax a per-test runtime ceiling to fit it.

**Tests and exit gate.** Strict warnings, ASan/UBSan and relevant thread checks;
all retained reduced-SSE numerical/native/application gates; each new backend's
complete positive and negative catalogs; CLI PCLI01–12 before/after runtime
changes; catalog/traceability and source-distribution audits. Rebuild and run
AMPS-BATSRUS-Europa, PS-BATSRUS-OH, MPS-BATSRUS-PIC-coupling and
AMPS-BATSRUS_dynamic_fluxrope with their original inputs/references/limits.

Publish per-profile JSON/JUnit, commands, seeds, MPI/workers, physical dt,
compiler/dependency/object/archive/executable hashes, physical input provenance,
independent metrics, runtime and retained failures. Demonstrate one complete
scheduled-orchestrator execution on the final committed source. Passing
individual manual campaigns alone does not prove automatic nightly coverage.
Update the support matrix and release notes only for exact combinations whose
required gates pass; unrelated FAIL/ERROR/SKIP results remain explicit.

### Decisions and completion checklist for this extension

Already selected: three srcSEP field-line movers, fixed physical stepping,
material-relative displacement, the approved accelerated-population source
contract and default absorbing shock return, prescribed initial turbulence,
three-dimensional ambient, analytical/semi-empirical regional evolution without
an MHD runtime, and coronal G1 → G2 → G3 → G4 order. Do not reopen those
decisions merely because an adapter is inconvenient.

| Decision to complete before dependent coding | Required concrete record |
| --- | --- |
| SWCME conservative production profile | Compatible rho/U/B/thermal equations, region support, source terms if any, SI/frame/composition definitions and independent conservation/induction/thermal gates. Preserve the old diagnostic profile. |
| Full regional lateral/map closure | Actual ambient-flow consistency and RH/contact representability evidence; then the selected vector/map evolution, initialization, inverse, force/heating/work and interface laws. The current candidate investigation is not a selected physical closure. |
| Source normalization for each new profile | Physical meaning and dimensions of the backend rate and release-reference conversion, tube/front measure, explicit inputs and zero/inapplicable support rules; no dimensionless-to-absolute assumption. |
| Profile-dependent missing initial/boundary data | Actual ejecta assets, initial sheath inventory, support/topology and thermal/flux boundary definitions wherever the maintained contract does not already supply them. No invented prehistory or default empty region. |
| Optional downstream particle transport | Separate specification of crossing/scattering/energy exchange and acceleration accounting, followed by native qualification. Until then preserve deletion by default and rejection of pass. |

A missing decision blocks only its dependent work. Continue independent audits,
contract design and tests whose expected results are established; ask a focused
physics question with the measured conflict and candidate consequences before
implementing an unresolved closure.

For each backend, completion requires: **approved core closure → independent
physics/numerical gates → actual native background → actual nonzero-particle
source and all three movers → MPI/thread/restart/sampling → automatic regression
evidence**. The coronal/composite path additionally retains G1–G4 ordering.
A fully covered background may be complete while an optional downstream
particle capability remains explicitly unavailable. Observation-based
validation is a separate later claim requiring independent data, frozen
calibration/hold-out roles and response metrics.
