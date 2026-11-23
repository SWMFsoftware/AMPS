# Immutable background snapshots and the simulation clock

Step 2 establishes one explicit state boundary between background providers and
particle transport. It does not change the governing transport equations or
select a different mover.

Current implementation and qualification are tracked in [STATUS.md](STATUS.md).
The later CP03 correction supersedes the original density-history derivative
reconstruction described in earlier reports.

## Authoritative clock

`PIC::SimulationTime` exposes application time. srcSEP reads it
through `SEP::Background::SimulationTimeSeconds()`, which returns seconds from
the PIC simulation-time origin and exposes no setter or independent advance
operation.

For SWMF single-global/counter-OFF execution the coupler owns the completed
physical-time ledger and supplies one explicit full nominal interval. SEP
prepares/freezes at its departure; the coupler presents the declared end during
PIC work, and commits its ledger only after successful application completion.
The initialized timestep remains constant. A host-limit acknowledgement on a
skipped call is not physical evolution and does not advance sampling exposure.

The former standalone driver kept a global elapsed time and a second static
SWCME launch-time counter. Those clocks were incremented in a different order
from PIC and caused the prepared SWCME state to lag the particle step. They have
been removed. Shock motion, Tecplot titles, spectral output, SWCME preparation,
and background validity now use the PIC clock adapter.

Local integration counters inside a single mover are not application clocks.
They measure substeps within the mover's supplied `dtTotal` and are intentionally
unchanged.

## Snapshot contents

`SEP::Background::BackgroundSnapshot` is immutable after construction and
records:

| Field | Meaning |
|---|---|
| provider | `analytic`, `swcme`, `reduced-shock`, `swmf`, or `local-evolution` authority |
| ownership | model-owned, imported read-only, or an explicit handoff copy |
| epoch | provider-state epoch in simulation seconds |
| previous/current physical epochs | provider realization times used for temporal derivatives; never a particle substep |
| validity | closed interval of simulation seconds for which particles may consume the state |
| field-line generation | nonzero identity of the imported/constructed field-line state |
| configuration fingerprint | stable 64-bit FNV-1a identity of canonical background-affecting configuration |
| provenance | human-readable origin or handoff record |

FNV-1a is used only as a deterministic identity/checking hash; it is not a
cryptographic signature. The canonical configuration deliberately excludes the
particle mover. Changing from Parker transport to either focused-transport
parameterization must not change the background provider, time origin, shock,
IMF, or turbulence representation.

The application fingerprint's content schema is `srcsep-background-v3`. The
existing full run/checkpoint fingerprint remains unchanged; background identity
uses a private projection that removes only the run mover value and its source
provenance. Its fixed canonical mover sentinel is never installed or executed.
All other existing run identity inputs are retained. For analytic backgrounds,
the same declared native origin/speed enters as full binary64 scientific text;
distinct representable winds therefore receive distinct identities. For SWMF,
the standard BATSRUS wave contract identity comes from its single typed adapter
authority. A two-wave import supplies no circular polarization merely because
its normalization is now explicit.

The snapshot does not deep-copy the large AMPS field-line arrays for every
particle. Instead, `SnapshotStore` freezes their provider/generation identity.
Publication is prohibited during the global particle read phase, and every
entry through the common `SEP::ParticleMover` wrapper acquires a const snapshot
reference that remains valid for the full mover call. The enclosing read-phase
guard owns the snapshot, and an atomic phase marker avoids a mutex/reference-count
operation in every particle call. This gives all scheduler threads one read-only
realization without a per-particle field-array copy.

`EvaluateLocalBackgroundAt` obtains an explicit Cartesian velocity gradient
`G_ij=partial U_i/partial x_j` [s^-1] from the declared provider. Divergence is
`trace(G)` and field-aligned strain is `b_i G_ij b_j`, using the physical B
direction. A derivative of `U dot b` along a curved line also contains a
curvature term and is not interchangeable with that strain. Temporal density
history alone supplies neither quantity without a material-coordinate,
continuity/source and composition contract.

The analytic radial and supported canonical SWCME/reduced branches have explicit
derivative contracts. The prescribed SWMF path reconstructs the full Cartesian
tensor internally from received position/velocity samples at one epoch. It
gathers across the original PT communicator, independently of PIC's disconnected
communicators, and uses a distance-weighted affine fit with column-pivoted
Householder QR. The deterministic stencil expands from 32 to at most 512 samples,
requires four distinct producer origin-grid lines, and rejects a rank-deficient
or ill-conditioned fit (`cond_1(R)>1e6`). Affine exactness, deficient-stencil
rejection, and convergence on a symmetric cubic field are covered by CP03G05–07.
Those manufactured tests do not establish accuracy for every received mesh.
Endpoint tensors are prepared before publication and only interpolated during
the immutable particle phase. Missing coverage rejects; no radial law, density
history, zero transverse derivative, or pseudoinverse fills a missing dimension.
The canonical SWCME1D
FULL_ICME interface also lacks regional thermodynamics needed for a coherent
native publication and is rejected before mutation.

Every mover composes snapshot validity with its departure-relative substep
time and verifies the field-line generation. A scheduled birth additionally
requires one snapshot covering its whole explicit parent window. Schema-4
native birth integration is in progress; its next native gate is pending.
Reaching a validity endpoint does not authorize silently changing the provider
or using a second background realization inside a particle read phase.

For a fixed analytic generated domain, shared `amps_init()` now validates the
whole native realization collectively before installing the existing inverse-square density
and initial histories. Every current magnetic vector must be finite/nonzero;
every current velocity must match its explicit radial declaration; and all
current/previous B/U/n slots must exist. Only then are current n and previous
B/U/n written, with each previous value equal to its current fixed-background
value. Every rank finishes preparation before any rank writes; the first failing
rank's status/context is returned consistently to all ranks. A bad later endpoint
on one rank leaves every rank's earlier primitives/history untouched. No
temperature or pressure profile is created.

`util/sep_analytic_radial_background.*` evaluates the explicit constant-speed
field `U=v (x-origin)/r` and its Cartesian tensor
`G=(v/r) (I-rhat*rhat)`, in SI. Its native declaration reads the existing active
initialization origin and wind speed; without runtime input it preserves the
documented legacy origin zero and 400 km/s. It has no independent mutable
configuration, fingerprint or snapshot publisher. Generated Parker/straight
lists now carry declared U themselves; mesh registration copies it exactly and
also preserves arbitrary imported list velocities. An imported line does not
thereby receive an analytic derivative or fixed initial history contract.
Selecting an imported domain with analytic history ownership fails explicitly
instead of replacing its primitives.

`ARB01`–`ARB06` provide independent translated-origin/nondefault-speed values,
Cartesian tensor and velocity-difference checks, invalid-input atomic rejection,
the geometry-only Parker builder contract, and origin/speed identity sensitivity.
`ARBN01`–`ARBN04` inspect actual
native startup histories, configured generator/list registration and atomic
whole-history rejection and frozen-mover/provenance identity exclusion. Their executed build/results and broader limitations
remain in [STATUS.md](STATUS.md). These fixed-background repairs do not qualify
time-dependent material geometry, source-frame mapping or production closure.

## Provider ownership and handoff

The store applies these rules before accepting a publication:

- analytic, SWCME, and reduced-shock state must be `ModelOwned`;
- SWMF state must be `ImportedReadOnly`;
- normal updates cannot change provider or configuration fingerprint;
- epochs and field-line generations cannot move backwards;
- a provider cannot replace another provider through the ordinary `Publish`
  operation;
- SWMF state can become locally evolved only after the caller copies the
  imported arrays and calls `PublishLocalEvolutionHandoff`, which records a new
  generation with `HandoffCopy` ownership;
- no update or handoff is permitted while particles are reading a snapshot.

Step 11 applies the same rule to wave ownership. `Turbulence::State` records
the imported epoch, source checksum, integrated/spectral representation, and a
one-time handoff marker. `swmf-read-only` is immutable;
`swmf-initial-then-local` can evolve only after that explicit handoff. A later
field-line remap increments its generation while conserving branch and
spectral-bin energy on the covered physical arc length.

The private SWCME adapter's `Configure` and `PrepareState` entry points call
`AssertProviderMayWrite` before changing their backing cache. This prevents a
provider update during a particle phase or while a different provider owns the
background, without exposing the provider's state type through `sep.h`.

### D01 fail-closed SWCME queries

`SW1DAdapter::QueryAtRadius` returns a `QueryResult`, not a Boolean plus three
write-through references. The result contains one immutable SI sample only when
the operation succeeded. Its status distinguishes an unprepared model,
non-finite input/output, an out-of-domain radius, non-positive density, invalid
speed, invalid divergence, and a canonical provider failure. A rejected query
therefore cannot leave density updated while speed or divergence still contains
an older value.

The default `strict` policy rejects every invalid query. Two recovery policies
exist for controlled experiments and are never selected implicitly:

- `clamp-radius` maps only a radius below the canonical `1.05 R_sun` boundary
  to that boundary; it does not repair provider output;
- `diagnostic-fallback` returns a user-specified sample whose density and speed
  must be positive and whose three fields must be finite.

The selected policy and fallback SI values are frozen in `Run::Configuration`,
included in its restart fingerprint, and copied into the background metadata
fingerprint. Process-local atomic counters record successful queries, rejected
queries, radius clamps, and diagnostic fallbacks. Rank zero prints these
counters in the final run summary. A failure message records MPI rank,
simulation epoch, requested/evaluated radius, failed field, immutable source
state ID, and canonical detail.

Preparation is transactional as well. `PrepareState` constructs a candidate
cache and publishes neither it nor new snapshot metadata until canonical
validation succeeds. If preparation fails, the previous valid cache and state
ID remain unchanged. The caller receives an explicit status and terminates the
step, so particles can never consume old fields under a new epoch label.

The handoff API records ownership but does not itself copy AMPS arrays. A future
caller enabling local evolution is responsible for making that one-time private
copy before publishing the handoff. No existing Step 2 path requests local
evolution.

## Update order

For a standalone SWCME step:

1. Read the upcoming epoch from `PIC::SimulationTime`.
2. Ask the private SWCME adapter to prepare its state at that exact epoch.
3. If preparation failed, report rank/epoch/previous state ID and stop without
   changing the current snapshot. Otherwise publish model-owned metadata valid
   through the upcoming global time step.
4. Enter a `ParticleReadPhase` and call `PIC::TimeStep()`.
5. Each particle mover acquires the same const snapshot in
   `SEP::ParticleMover`.
6. End the read phase before shock/turbulence/background writers run.

For SWMF coupling, `PrepareSnapshotForParticleStep()` observes
`AMPS2SWMF::MagneticFieldLineUpdate::LastCouplingTime`. A changed coupling epoch
is published between particle phases as a new imported, read-only field-line
generation. The latest import remains valid until another coupling generation
arrives. A missing, future, stale, or cross-provider state is fatal; srcSEP does
not silently continue with unrelated arrays.

The coupled clock binding does not relax these checks. In particular, a newer
import may lie after the actual AMPS departure when coupling intervals are not
divisible by the nominal timestep. A temporal reconstruction/retained-background
policy must be explicitly established before qualifying that evolving case;
neither a shortened AMPS step nor altered metadata supplies the missing state.

### Reduced shock-front ambient publication

With `--reduced-shock-event`, the application resolves the strict HCI event
and its checksummed assets before native allocation. At each cadence-aligned
epoch it prepares the shared front transactionally, stages ambient samples for
every native field-line vertex, preserves the outgoing vertex state in the
previous datums, installs the candidate, and publishes one `ReducedShock`
snapshot using the provider's exact generation. A failed query changes neither
the arrays nor snapshot metadata.

This provider's spatial contract is ambient-only. Accepted shock records and
their immediate upstream/downstream limits remain attached to the separate
shared surface epoch; they are not a downstream volume and are never copied
behind the front. Every rank evaluates the replicated line and exact MPI
min/max agreement is required for an endpoint plasma/IMF record. The native
post-step reduction independently requires zero actual AMPS particles.

The implementation and runnable case are
`adapters/reduced_shock_background_adapter.*` and
[`examples/reduced-shock/README.md`](examples/reduced-shock/README.md).

## Public implementation surface

- `util/sep_background_snapshot.*`: dependency-light immutable value, store,
  read-phase guard, ownership checks, and deterministic fingerprint.
- `util/sep_background_runtime.*`: PIC clock adapter, runtime fingerprint,
  standalone publication, SWMF refresh, and explicit local handoff.
- `SEP::ParticleMover`: common per-mover const acquisition point.
- `amps_time_step`: one global RAII read phase around `PIC::TimeStep()`.

## Verification

Run the dependency-light focused tests with:

```sh
make test-state-unit
make test-swcme-fail-closed-unit
```

The test compiles the production snapshot implementation using C++11,
`-Wall -Wextra -Werror`, and pthread support. It checks construction domains,
compile-time non-assignability, interval enforcement, provider isolation,
explicit SWMF handoff, generation/epoch monotonicity, publication exclusion,
concurrent mover views, deterministic fingerprints, and the single-clock source
invariant.

The D01 test links the production adapter to the canonical SWCME model and
checks unprepared access, the exact inner boundary and its adjacent invalid
point, NaN/infinity, negative density, zero speed, failed-preparation atomicity,
explicit clamp/fallback counters, complete diagnostics, and unchanged fast and
slow preset trajectories.

Those assertions are implemented once in `util/sep_swcme_validation.cpp` and
registered as extended native ID `D01`. `make test-swcme-fail-closed-unit`
compiles a dependency-light launcher around that descriptor; a linked complete
run discovers the same callback with `amps --list-tests` and executes it through
`test/run_tests.py --amps ../amps --all`. The callback restores a prepared fast,
strict state before returning, while the complete runner additionally gives it
process isolation.

The complete native gate remains:

```sh
make -j test SEP_EXECUTABLE=/path/to/amps
```

It additionally requires the enclosing AMPS configuration, PIC/field-line
headers, MPI toolchain, and linked executable.

## Current limits

Step 2 freezes identity and publication lifetime; it does not redesign AMPS
vertex storage as intrinsically const C++ memory. Background writers must remain
outside `PIC::TimeStep()`, as enforced by the store API. Later steps still need
to move particle-wave source accumulation to thread-local storage and perform a
deterministic post-mover reduction. That work is intentionally not folded into
this state/clock change.

## Explicit analytic shock departure and completed host geometry

The fixed analytic host prepares the existing shock trajectory explicitly at
the physical departure before capturing the previous radius or opening the
particle read phase. `Shock::EvaluateSourceAtEpoch` returns an immutable
configuration/epoch/radius/speed view; native injection requests it at its
frozen step-start birth epoch. Count, compression, intersection and swept speed
consume that same view. No mutable PIC counter query or lazy trajectory write
is permitted in this path. Known injection S is passed to the existing core
allocation helper, so initial S and Cartesian X have one coordinate authority.

After actual PIC work, migration and the read guard have joined, the host closes
the birth window and collectively advances the analytic owner to that window's
completed epoch. Only then does it derive the wave-work shock-radius interval.
This prevents a cached departure radius from masquerading as completed geometry
and prevents the old uninitialized sentinel-to-launch displacement from being
counted as physical shock work. A restored previous radius inconsistent with
the declared departure rejects instead of bridging a stale interval. Candidate
validation and MPI agreement precede native cache publication.

This repair preserves the existing radius-dependent speed/compression/density
laws. That legacy source remains global-zero-centered; a translated analytic
wind has no newly inferred source count/sphere contract and is rejected when
this source host path is requested. SWCME retains its prepared provider-owned
trajectory, and imported/SWMF clocks are never replaced by this analytic API.
The fixed-background accessor tests may use translated origins without running
a shock-source host step. No thermal or guiding-center source closure is added.

The source-specific native query and consumers require an active frozen driver
window and exactly matching prepared owner/departure/full trajectory config.
They reject unprepared or foreign views without consulting PIC's mutable clock
or modifying the output. General arbitrary-time trajectory mathematics remains
available separately through `EvaluateSourceAtEpoch`; it cannot by itself
authorize native source consumption.

ASE01–05 check explicit-time values and owner binding only. HOST01–03 are separate, exclusive
fresh-process diagnostics of real PIC/host steps; ordinary catalog SKIP is not
native evidence. The complete input protocol, immutable per-rank traces,
counter-policy distinction and exact qualification limits are documented in
[`test/README.md`](test/README.md). Native completion is tied to the staged
binary/status record, rather than inferred from the source-only trajectory gate.

All srcSEP particle coordinates use `plasma-comoving-v1`, regardless of whether
the data provider is analytic, SWCME, reduced-front, SWMF or local evolution.
The provider identifies data authority, not a choice of particle reference
frame. Imported SWMF metadata must agree with this contract. Unknown providers
and incompatible imported frames reject before transport. Physical U and
grad(U) remain distinct from the zero coordinate-advection value; the material
particle equation never zeros thermodynamic compression. The legacy analytical native
geometry/shock update ordering remains unqualified; see TRANSPORT_NUMERICS.md
for the identified conditional-update and endpoint limitations.


## Reduced-core material owner

`ReducedShock::IntervalView` retains the actual core provider and explicit
closed physical interval. Configuration changes cannot rebind an existing
view. `CaptureInterval` preserves the strict prepared-cache validity API;
`CapturePhysicalStep` queries the core's continuous analytical functions over
its declared event support. It neither interpolates RH states nor advances the
committed front as a side effect. Cache cadence, geometry generation and
physical epoch are distinct metadata; none changes the AMPS physical timestep.

`Native::InstallBackground` stages a copy of each native stored-data prefix,
all current/previous B,U,np,p,T and optional ideal E=-U cross B, and the complete
interval connectivity. Initial previous state is the same initial condition,
not an invented earlier epoch. A repeated physical epoch preserves genuine
previous state. `CadencePublication` stages the core's ordinary Prepare path
on a value copy sharing immutable configuration/ambient; only its Current
pointer is independent. No live core cache changes before readiness.
`SnapshotStore::Publication` likewise preallocates/validates metadata; commit
uses a single-use expected-owner token. Every native physical state is checked
for exact replicated MPI agreement before publication. Local sampling
numerators are not physical clocks and are never summed as extra exposure.

`Native::Advance` stages every moved vertex and target-epoch primitive; a new
vertex gets previous state from the actual inverse flow map. Old diagnostics
are copied; no wave remap or source population is invented. Only the Parker
exterior is admitted as a material native background. The core's PFSS query
capability alone does not prove that branch satisfies the required material
induction/continuity contract. Stages 0–4 preserve a prescribed, read-only
ambient rather than extending one-sided RH limits through a downstream volume.
# Reduced source ownership and PIC phase ordering

The reduced core remains a prescribed ambient/front provider. Its
`particle_mode=disabled` capability is unchanged. A separate application-owned
consumer is enabled explicitly and bound into the run configuration; no
downstream volume or core particle engine is implied.

The host freezes the full physical source window before PIC changes its
counter. PIC samples departure particles, then invokes the reduced source
callback inside the immutable read phase. Future births are therefore absent
from the earlier sample. Their own material-line positions/front states are
queried at explicit birth epochs, and their remaining first movement is stored
as a duration. Neither absolute timestamp rounding nor the provider's cache
cadence changes AMPS dt, output cadence or exposure. After all workers join,
geometry retirement and number/birth-energy/event-energy ledgers are audited.
