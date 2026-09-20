# Phase A: AMPS Mover and SWCME Source Adapters

Phase A connects the AMPS particle buffer and SWCME shock source to the
AMPS-independent Phase-P transport cores. It deliberately does not introduce
another Parker equation, focused equation, scattering closure, or shock model.
The adapter boundary performs representation conversion, validation, list
bookkeeping, and conservation accounting only.

## Production mover path

`AMPS::Movers::MoveParticle` is the sole AMPS mover entry point. The immutable
run configuration selects exactly one of two registered cores:

| Canonical name | Configuration | Core |
|---|---|---|
| `parker3d-tensor` | `TransportModel::Parker3D` | `AdvanceParker` |
| `focused3d-split` | `TransportModel::Focused3D` | `AdvanceFocused` |

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
physically constant coefficient. For focused transport it supplies both
`D_mumu` and `dD_mumu/dmu` from the shared coefficient bridge.

### Persistent particle schema

`RequestParticleStorage()` reserves one packed AMPS extension before the
particle buffer is frozen. The extension contains:

- a 64-bit schema tag;
- stable particle ID;
- completed global step and transport substep;
- last crossed shock generation;
- total relativistic momentum, pitch cosine, and gyrophase.

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
`powerIndex=q-2`. The minimum and maximum kinetic energies are converted to
relativistic SI momentum by the common SWCME helper. `SampleInjectedParticle`
then uses independent keyed streams for momentum, pitch angle, gyrophase, and
stable ID. Isotropic launch uses `mu=2u-1` and `phi=2 pi u`.

Each macroparticle represents

\[
w_i=\frac{w_{patch}\,\epsilon_{inj}}{N_{macro}}.
\]

The adapter is dimension independent: a 1-D or 3-D application supplied with
the same common source record, campaign, generation, patch ID, species, and
macro index obtains identical spectrum fingerprints, random keys, momenta,
and weights. No `srcSEP` source file is inspected or linked.

### Schema-3 standalone provider and exact global count

`CreateStandaloneSwcmeShockProvider` re-resolves the frozen raw `[swcme]`
assignments and requires the resulting canonical manifest/fingerprint to match
the immutable application configuration byte-for-byte. It then constructs the
canonical `swcme::sep::Interface3D`, evaluates the complete surface at the
declared first valid epoch, and rejects an invalid MHD solve, invalid mesh,
empty active source, or insufficient computational sample count before AMPS
mesh allocation. Before `event.valid_from` it publishes a valid inactive state;
it does not invent a shock.

For schema 3, `source.samples_per_step` is a global integer, not an independent
expectation for every patch. `AllocateExactPatchMacroparticles` reserves one
representative for each positive-weight active patch, apportions the remaining
integer samples in proportion to canonical physical patch weight, and assigns
largest remainders with stable source-ID/index tie-breaking. The sum is exactly
the input count on every active step. A count below the active patch cardinality
fails closed because silently omitting a nonzero source patch would not be a
conservative representation.

`SourceRequest::prescribedMacroparticles` carries each exact patch allocation
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

- `ADP3D01`: exact two-entry registry and validating dispatch.
- `NAT3D04`: inner, outer, and invalid-background dispositions.
- `NAT3D05`: exact ledger closure and transactional mismatch.
- `NAT3D08`: first moving-shock root and generation de-duplication.
- `SHK3D01–04`: dimensional source identity, analytic shock geometry, source
  ownership guards, and source-weight normalization.
- `BLDL3D01/03`: configured AMPS compilation and actual mover-return ABI.
- `R3D01–02`: installed generated hook and complete re-resolved requested-time
  advancement.
- `R3D05`: physical source normalization, cap/disconnection policy, and unique
  cadence identity.
- `R3D08`: canonical provider preflight, delayed activation, deterministic
  largest-remainder allocation, exact global count, and no downstream cap.
