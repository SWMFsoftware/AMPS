# Phase A: AMPS Mover and SWCME Source Adapters

Phase A connects the AMPS particle buffer and SWCME shock source to the
AMPS-independent Phase-P transport cores. It deliberately does not introduce
another Parker equation, focused equation, scattering closure, or shock model.
The adapter boundary performs representation conversion, validation, list
bookkeeping, and conservation accounting only.

## Production mover path

`AMPS::Movers::MoveParticle` is the sole AMPS mover entry point. The immutable
run configuration selects exactly one of three registered cores:

| Canonical name | Configuration | Core |
|---|---|---|
| `parker` | `TransportModel::Parker3D` | `AdvanceParker` |
| `focused-diffusion` | `TransportModel::FocusedDiffusion3D` | `AdvanceFocused` |
| `focused-scattering` | `TransportModel::FocusedScattering3D` | `AdvanceFocusedScattering` |

The path for one particle is:

1. Read the AMPS Cartesian position, species, statistical-weight correction,
   and the packed srcSEP3D persistent extension.
2. Reject an absent schema tag or zero stable particle ID. The allocator slot
   is never used as a stochastic identity because slots change after deletion,
   restart, and MPI migration.
3. Require `Runtime::Running`, then enter the requested-time loop. Before each
   accepted substep the bridge finds the current AMR node and calls the
   installed `LocalRecordResolver`. The resolver uses only the Runtime-pinned
   background/turbulence generations and returns a complete local record.
4. `AdvanceParticleRequestedTime` selects a named substep, dispatches one
   Phase-P core, advances the moving-shock radius by consumed time, and repeats
   until the complete AMPS interval is consumed or a terminal state occurs.
5. Classify the final radius as active, inner-boundary absorbed, outer-boundary
   escaped, or failed. Evaluate an expanding-shock intersection on the accepted
   segment and suppress duplicate crossings of the same shock generation.
6. For an active particle, write position, gyrotropic persistent state, and a
   deterministic Cartesian velocity; then insert the particle into AMPS's
   temporary destination-cell list and return `_PARTICLE_MOTION_FINISHED_`.
   Terminal records are deleted and return `_PARTICLE_LEFT_THE_DOMAIN_`.

The local resolver is explicit because a coupled run must never read mutable
SWMF arrays while mover workers are active. For a Parker run it must provide
`kappaParallelM2PerS` and the full field-aligned derivative
`dKappaParallelDsMPerS`; setting the derivative to zero is valid only for a
physically constant coefficient. For focused diffusion it supplies both
`D_mumu` and `dD_mumu/dmu`; for focused scattering it supplies
`lambda_parallel` and the directional turbulence fractions. All are produced
by the shared coefficient bridge under the immutable input selectors.

### Persistent particle schema

`RequestParticleStorage()` reserves one packed AMPS extension before the
particle buffer is frozen. The extension contains:

- a 64-bit schema tag;
- stable particle ID;
- completed global step and transport substep;
- last crossed shock generation;
- total relativistic momentum, pitch cosine, and gyrophase;
- residual exponential scattering optical depth and next scattering-event
  index.

All access uses `memcpy`; AMPS extension offsets are not assumed to satisfy C++
structure alignment. AMPS's own particle restart therefore carries the same
stochastic tuple. `InitializeParticle()` is the only supported way to attach
this state to a newly allocated source particle.

### Velocity reconstruction

The transport state is gyrotropic. AMPS still expects a Cartesian velocity for
generic diagnostics. The adapter computes

\[
\mathbf v=v\left[\mu\mathbf b+\sqrt{1-\mu^2}
(\cos\phi\,\mathbf e_1+\sin\phi\,\mathbf e_2)\right],
\]

where `e1` is formed by crossing `b` with the Cartesian axis least aligned with
it and `e2=b×e1`. This deterministic basis is nonsingular at magnetic poles and
does not consume random numbers.

## Expanding-shock intersection

For a straight accepted substep
`x(u)=x0+u(x1-x0)`, `0<=u<=1`, and a spherical shock
`R(u)=R0+u Vsh dt`, `FirstShockIntersection` solves

\[
|\mathbf x_0-\mathbf c+u\Delta\mathbf x|^2
=(R_0+uV_{sh}\Delta t)^2.
\]

The smallest root in `[0,1]` is the physical first crossing. A linear fallback
handles a numerically vanishing quadratic coefficient. The mover also adds a
conservative pre-step limiter based on radial gap divided by particle, plasma,
and shock closing speeds; the post-step quadratic remains authoritative.

## SWCME injection path

`MakeShockSourceRecord` consumes the canonical
`swcme::sep::SEPSourceState`, already validated by SWCME shock acceleration.
For isotropic diffusive-shock acceleration,

\[
f(p)\propto p^{-q},\qquad \frac{dN}{dp}\propto p^2f(p)
\propto p^{-(q-2)}.
\]

Therefore the shared `sep_common` linear-momentum spectrum receives
`powerIndex=q-2`. In `local-compression-dsa` mode, q is the patch-local value
published by SWCME. In `fixed-phase-space-power-law` mode,
`ConfigureSpeciesSpectrum` replaces only that shape with the reviewed positive
q from `[source].phase_space_power_index`; shock position, normal, generation,
patch weight, and compression diagnostics remain provider-owned. This is how a
published \(f(p)\propto p^{-5}\) boundary is represented without pretending
that its slope came from a compression ratio.

The minimum and maximum kinetic energies are converted to relativistic SI
momentum separately with every compiled species mass.
`SampleInjectedParticle` then uses independent keyed streams for momentum,
pitch angle, gyrophase, and stable ID. Isotropic launch uses `mu=2u-1` and
`phi=2 pi u`. The spectrum selector does not derive an absolute birth rate
from a boundary phase-space density; the source-number contract remains
explicit and auditable.

Each macroparticle represents

\[
w_i=\frac{w_{patch}\,\epsilon_{inj}}{N_{macro}}.
\]

The adapter is dimension independent: a 1-D or 3-D application supplied with
the same common source record, campaign, generation, patch ID, species, and
macro index obtains identical spectrum fingerprints, random keys, momenta,
and weights. No `srcSEP` source file is inspected or linked.

### Schema-3 standalone provider and exact per-species count

`CreateStandaloneSwcmeShockProvider` re-resolves the frozen raw `[swcme]`
assignments and requires the resulting canonical manifest/fingerprint to match
the immutable application configuration byte-for-byte. It then constructs the
canonical `swcme::sep::Interface3D`, evaluates the complete surface at the
declared first valid epoch, and rejects an invalid MHD solve, invalid mesh,
empty active source, or insufficient computational sample count before AMPS
mesh allocation. Before `event.valid_from` it publishes a valid inactive state;
it does not invent a shock.

For schemas 3 and 4, `source.samples_per_step` is an exact integer for each compiled
AMPS species, not an independent expectation for every patch. For each species,
`AllocateExactPatchMacroparticles` reserves one representative for each
positive-weight active patch, apportions the remaining integer samples in
proportion to canonical physical patch weight, and assigns largest remainders
with stable source-ID/index tie-breaking. The sum is exactly the input count for
that species on every active step. A count below the active patch cardinality
fails closed because silently omitting a nonzero source patch would not be a
conservative representation.

Before allocation, `ConfigureSpeciesSpectrum` converts the declared total
kinetic-energy bounds to momentum using the current compiled AMPS mass. This
prevents SWCME's reference-particle momentum interval from being reused across
unlike species. `SourceRequest::prescribedMacroparticles` carries each exact patch allocation
through `BuildInjectionPlan`. That path never stochastically rounds or caps the
count. Instead, every particle receives

\[
W_{patch}=\frac{\dot N_{seed}\,\epsilon_{inj}\,
w_{patch}\,\Delta t}{N_{patch}},
\]

through AMPS' individual statistical-weight correction. Schemas 1–2 retain the
keyed stochastic-rounding/cap behavior for compatibility. Both paths record the
represented physical population, energy, momentum, count, rejections, and caps
in the source ledger.

## Exact particle ledger

Each globally reduced `(step,species)` row must satisfy the integer identity

\[
N_{start}+N_{injected}=N_{end}+N_{escaped}+N_{absorbed}+N_{failed}.
\]

`advanced` and `shockCrossings` are diagnostics and are not additional sinks.
A mismatched close returns an error and leaves the row open, preserving the
evidence needed to locate a list or return-code bug.

Production opens rank-local rows immediately before `PIC::TimeStep`. Movers
record exactly one final disposition after consuming the complete requested
time. After AMPS finishes list exchange, start/outcome/end counters are summed
with `MPI_Allreduce` and imported only if the global row closes exactly. This
ordering makes ordinary rank migration invisible to conservation while still
detecting an invalid final-list insertion or an unaccounted deletion.

## Particle splitting and merging

### AMPS-core audit

AMPS supplies two population-control layers in
`src/pic/pic_particle_spliting.cpp`:

- the legacy global `ParticleSplitting::Mode`, called automatically near the
  end of `PIC::TimeStep`; and
- the newer `MergeParticleList(spec, head, target)` and
  `SplitParticleList(spec, head, target)` linked-list routines.

The generic routines are useful allocation/list mechanisms, but they are not a
safe complete SEP resampler without application callbacks. The audit found:

1. historical empty-bin/list paths could dereference `particles[0]` or call
   `std::next(end)`; guards now return for non-positive targets, empty merge
   lists, and split lists with fewer than two particles;
2. `CloneParticle` correctly preserves AMPS `next`/`prev` links, but it copies
   every application extension byte, including srcSEP3D stable identity and
   residual scattering state;
3. the generic 3-to-2 merge conserves a nonrelativistic weighted `v^2`
   moment. SEP particles may be relativistic, so that is not conservation of
   their kinetic energy;
4. the generic merge uses process-order `rnd()` and does not update
   application momentum, pitch, gyrophase, shock-generation, or semantic RNG
   fields after changing Cartesian velocity.

The recommended core-level completion is a registered application resampling
policy with callbacks for (a) energy/momentum reconstruction, (b) extension
state after split/merge, and (c) semantic random direction. Until such a policy
exists, srcSEP3D disables the automatic legacy mode immediately after
`PIC::Init_BeforeParser()` and does not call the generic high-level merge.

### SEP-aware boundary controller

`ApplyPopulationControl` still uses the tested AMPS primitives
`GetNewParticle(head)`, `CloneParticle`, and `DeleteParticle(ptr, head)` so list
ownership remains in AMPS. The application layer supplies the missing physics:

- split the current heaviest representative into an original and clone with
  exactly half the input individual weight;
- merge three deterministic low-weight representatives into two equal-weight
  outputs at their weighted position centroid;
- solve the output momentum displacement by bisection against exact
  relativistic energy
  `sqrt(p^2 c^2 + m^2 c^4)-m c^2`, while conserving vector momentum and total
  statistical weight;
- resolve the magnetic direction at the actual output centroid and recompute
  pitch, gyrophase, and Cartesian velocity from each new momentum;
- assign semantic stable IDs to merge products, a distinct ID to the split
  child, reset event optical depth under the new histories, and retain the
  maximum last-shock generation;
- restart the semantic operation ordinal independently within each physical
  cell/species population, so changing MPI block ownership cannot change a
  merge direction or post-resampling stable ID; and
- report maximum relative weight, momentum, and relativistic-energy residuals.

Hysteresis is per active cell and compiled species. An empty cell is never
populated. If a nonempty count is below `minimum` or above `maximum`, repeated
operations drive it to `target`; a count inside the band is untouched. The
controller runs only at a joined iteration boundary, after shock injection and
before observer publication/checkpointing. Thus no mover, source allocator, or
MPI list exchange can run concurrently, and the controlled representation is
the one that enters the next `PIC::TimeStep`.

The limits are deliberately not one global particle count. Moving weight
between spatial cells merely to satisfy a global target changes density and
the phase-space distribution. Per-cell/per-species maxima provide a concrete
upper bound proportional to the occupied active mesh, while local minima
protect statistics only where a represented population already exists.

## Production configuration requirement

Run the hook after configuring AMPS and before compiling `pic_mover.cpp`:

```bash
make -C srcSEP3D prepare-production
```

`amps/install_mover_hook.py` inserts the exact declaration and maps the
generated `_PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_` macro to
`SEP3D::AMPS::Movers::MoveParticle`. It is idempotent and refuses to overwrite
an unrelated mover. `strict-production` depends on this target and audits the
result. During `amps_init()`, srcSEP3D installs the resolver, substep cap,
ledger, and initial shock state as one immutable `AMPS::Movers::Context`; a
coupled host supplies providers, not a second mover context.

Repeated source cadences under one physical shock generation receive distinct
`injectionSequence` keys derived from `(shock generation, Runtime tick)`. The
physical generation remains unchanged in `lastShockGeneration`, so stochastic
identity cannot collide and crossing de-duplication retains its intended
meaning. A source row is retained only on the AMPS rank that owns the patch
position, avoiding replicated physical totals at checkpoint gather.

## Evidence

- `ADP3D01`: exact three-entry registry and validating dispatch.
- `POP3D01`: relativistic 3-to-2 weight, momentum, energy, and position
  centroid conservation.
- `NAT3D04`: inner, outer, and invalid-background dispositions.
- `NAT3D05`: exact ledger closure and transactional mismatch.
- `NAT3D08`: first moving-shock root and generation de-duplication.
- `SHK3D01–04`: dimensional source identity, analytic shock geometry, source
  ownership guards, and source-weight normalization.
- `BLDL3D01/03`: configured AMPS compilation and actual mover-return ABI.
- `R3D01–02`: installed generated hook and complete re-resolved requested-time
  advancement.
- `R3D05`: species-dependent energy-to-momentum conversion, physical source
  normalization, cap/disconnection policy, and unique cadence identity.
- `R3D08`: canonical provider preflight, delayed activation, deterministic
  largest-remainder allocation, exact per-species count, and no downstream cap.
