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
Stages 8–9 remain pending.

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
