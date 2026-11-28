# Field-line transport and turbulence numerical contract

The stage-8 support matrix rejects backgrounds without a qualified material
evolution contract independently of mover selection. Stage-9 full application
execution uses the original fixed dt=300 s, 400-step output window and signed
number/birth-energy/current-energy ledgers. The maintained launcher verifies
matching build/input identity and actual receipts before the unchanged
independent checker reconstructs each segment's age/exposure. No test runner
may shorten the physical step, insert extra samples or use a global denominator
for younger material segments. See the
[stage-8/9 contract](validation/STAGES8_9_TEST_CONTRACT.md).

Steps 6–10 define one numerical vocabulary for every public srcSEP mover and one
canonical implementation each for the Parker and coefficient-driven focused
transport equations. The kernels under `util/` intentionally depend on neither
PIC nor MPI, while `transport_common.*`, `coefficient_providers.*`, and the
three canonical mover adapters connect those kernels to the host application.

The 2026-11-16 user decisions require shared polarized-spectrum scattering for
both focused movers, with the finite-band resonance gap preserved. Its new
joint momentum/pitch operator is under development. The retained plasma-frame
Dmumu kick and isotropic MFP reset described below remain implementation
history; they do not qualify that selected physical coupling. The current
requirements and study design are recorded in
`validation/development_requirements.json` and
`validation/studies/weak_alfven_domain.json`.

## State, frames, and units

The common particle state contains species, field-line identifier, host line
coordinate, particle mass [kg], and local parallel/normal velocities [m/s]. The
host coordinate is not assumed to be arc length. Every core returns a physical
displacement [m], and only `PICAdapter::AdvanceAlongFieldLine` passes that
distance to `FieldLine::move`; this is the single metric conversion.

Momentum is relativistic SI momentum [kg m/s]. `MomentumFromSpeed` requires
`0 <= v < c`; `SpeedFromMomentum` accepts finite nonnegative momentum. Neither
routine clips an invalid state. In the Parker/isotropic limit, cooling is
evaluated in the local plasma frame:

`dp/dt = -(p/3) div(U)`, therefore `p(t+dt)=p(t) exp[-div(U) dt/3]`.

The geometric coordinate and stored parallel velocity follow the line tangent
`t`. All srcSEP field-line particle coordinates are plasma-comoving, independent
of the provider of magnetic field, plasma state, or shock data.
The native adapter preserves physical `U dot t` and the full physical grad(U),
but supplies a separate coordinate-advection value to all three mover cores:
exactly zero for every srcSEP provider, including analytic, SWCME, reduced-front,
SWMF and local-evolution data. A provider cannot select an Eulerian particle
equation. The shared fixed-grid helper remains available for other applications
and manufactured controls; it is not a selectable srcSEP provider mode.
Both focused predictor/corrector paths and geometry/shock-crossing limits use
that same coordinate-relative value. Parker retains its separate Ito diffusion
drift. Physical compression, field-aligned strain and wave-frame transforms
continue to consume physical plasma inputs. The public core background field
`plasmaAdvectionMPerS` retains its historical name for fixture compatibility;
its argument is explicitly advection relative to the declared geometry.

This particle equation requires a material geometry owner as well. The current
analytical path is not yet qualified as a coherent moving-grid host: the generic
`PIC::FieldLine::Update()` is conditional on DATAFILE/global-clock build modes,
executes inside `PIC::TimeStep()` after departure shock distances/snapshots are
prepared, and skips the last vertex. A shock intersection with current native
vertices alone does not prove that the vertices represent the required epoch.
The approved material boundary rule moves every existing node, including both
endpoints, extends the inner end along B until it covers the minimum radius,
and removes wholly exterior trailing segments. New segments start empty and
use existing source rules. The callback API and required geometry, sampling,
conservation and restart work are specified in
`../src/pic/MOVING_FIELD_LINE_GEOMETRY_DESIGN.md`. The generic native transaction
and sampling implementation are available; the analytical owner still needs
its declared prescribed-versus-material physical state callback and qualified
geometry/shock publication ordering. Preflight rejects an unregistered owner.
The fixed-geometry checkpoint codec rejects material topology/exposure rather
than synthesizing restart history. This is not yet a complete analytical host
or a full coupled/native-MPI qualification.

The PT/CON_bline wrapper declares its plasma-advected label convention on the
wire; this declaration enters the immutable source namespace and background
fingerprint. It cannot change on a live population. Prefix insertion/removal
maps S by oldFirstLabel-newFirstLabel, retaining the material segment fraction.
Surviving particle X is refreshed from the newly published geometry before
sampling; weights, momentum, identity and RNG state are unchanged. Removed
segments follow the existing absorbing endpoint policy and emit represented
particle-loss counts/weights. Noncontiguous producer labels reject rather than
silently invent a remeshing map. New footpoint labels are boundary construction,
not permission to move surviving particles to unrelated labels. This remains
the existing departure-snapshot/operator-split geometry policy; it does not
introduce temporal field extrapolation or qualify arbitrary relabelling/remeshing.
The AMPS global timestep and sampling exposure remain unchanged.

For the registered reduced-material owner, checkpoint schema 5 retains topology,
labels, the actual background/cache epochs, all particle history, diagnostic
buffer phase and each segment's exposure. Native MPI continuation is qualified
separately from the retained fixed-line codec. After a head insertion the stored
binary64 coordinate is authoritative: Cartesian X is interpolated from its
translated fractional part, because integer addition can round that fraction.
An unrepresentable segment address rejects the private geometry transaction;
no particle is clipped to another segment.

Segment-attached workers select the configured global physical timestep without
looking up an Eulerian AMR block. Field-line ownership and block ownership can
belong to different MPI ranks. Material sampling reads a canonical index ordered
by species and the exact bits of S, parallel/normal velocity and individual
weight. These are the per-particle sampling inputs; equal keys yield equal
contributions. The index changes neither particle/list state nor any sampling
formula. It removes thread-arrival order from binary64 diagnostic addition,
while retaining one exposure tick at the same fixed dt. These contracts are
graded by SINT08's reversed-order native sample and the TCR01–03 fresh-process
threaded restart cases. The isolated wave/codec CP05O/OF cases do not qualify
particle motion on a frozen nonzero-U line or self-consistent wave feedback.

The keyed reduced-material checkpoint serializes each complete normalized
particle record in immutable-ID order. Thread arrival and pool allocation are
not physical ordering variables. Sampling, inventory and conservative
population control already use deterministic orders, so this changes neither
transport nor a ledger. The fixed-line shared-stream operator codec retains
its semantic list order. Native material continuation additionally reverses
actual lists and requires the complete image to remain exact.

## Status, boundaries, and substeps

Numerical boundaries return `Status`/`StatusCode`. Invalid particle state,
coefficient, argument, domain exit, and timestep underflow remain distinct.
Production adapters translate an unexpected non-OK result to the host fatal
error path; focused tests inspect it directly.

Physical line-coordinate boundaries declare either absorption or reflection.
Pitch angle always uses repeated reflection into `[-1,1]`; the modulo-folding
implementation handles a stochastic increment that crosses multiple boundaries
and keeps `mu=0` regular.

Each mover assembles named limits for its active operators and selects their
minimum with `SelectSubstep`. The selector records the limiting name and count.
If that minimum falls below the declared threshold it reports `StepUnderflow`;
it never increases the interval with an undocumented floor.

This complete composable-limit contract is used by all three canonical movers.
The event-driven mover additionally splits at the current shock and immutable
snapshot-validity boundaries.

## Deterministic stochastic streams

`KeyedRandomStream` is constructed from campaign seed, particle identity,
operator identity, and event/substep index. MPI rank and worker number are not
keys. The open-interval uniform generator supports Box–Muller normals without a
singular logarithm, and an identical key reproduces an identical stream
independently of scheduling.

Transport schema 4 persists the complete three operator streams, including
key identity, algorithm version, cached normal, and accepted-interval count.
Stream schema 2 uses rejection of zero mantissas to return strictly open
uniforms; the previous endpoint-rounding algorithm is rejected on restore.
MFP residual optical depth/event count and exact-split serial are also stored.
Fresh source identities come from persisted source-kind/physical-line/species/
event/macroparticle tuples and exclude sampled momentum. Exact splitting keeps
the parent's history and gives only the actual child an independent stream.
Migration byte copies preserve the entire record. Complete field-line physics
checkpoint qualification remains separate from these storage/kernel checks.

### Native birth windows and first movement (schema 4)

The common native adapter now stores a padding-free 64-byte `Birth::Record`
with every particle. It distinguishes an explicitly created initial condition,
a source applied at the beginning of a timestep, and a uniform-in-time source
birth. A scheduled birth records the frozen parent start/end epochs [s], its
physical parent duration [s], its positive remaining first-move duration [s],
and whether that movement remains pending. Default/unknown records are errors;
loading an existing particle never supplies a fresh epoch or event zero.

The shared host driver prepares one background snapshot, collectively freezes
the source window, opens the particle read phase, calls `PIC::TimeStep()`,
releases the read guard, and closes the window before wave evolution. Freeze
uses the authoritative PIC clock and the active global timestep. Audited
standalone single-global builds advance the clock at start (counter OFF) or end
(counter ON). Species-global builds currently advance it only for one compiled
species with counter ON; equal multi-species durations alone do not implement
core clock synchronization and are rejected. Unequal species durations and
unbound externally owned clocks return `UnsupportedConfiguration`. SWMF
single-global/counter-OFF execution instead requires the coupler's explicit
full-step binding. Its duration must equal the initialized, unchanged global
timestep. The application freezes at the physical departure, and the coupler
then presents the end counter during PIC execution without another PIC clock
increment. Source-window closure verifies that endpoint and all pending births
before the wrapper can commit completed physical time. Skipped calls add no
particle movement or sampling exposure; deferred time survives normal session
transitions. Invalid background coverage still rejects rather than changing
dt, relabelling an import or inventing temporal reconstruction.

For this PT origin-grid producer, source channels carry its explicitly supplied
global line label rather than a local array index or moving geometry hash.
The namespace records the declared origin grid/frame. Missing or changed
labels still reject. A declared fresh run initializes the existing scheduler
after the first actual import and before particle prepopulation; continuation
must restore source history instead. Source event counters, identity hashes
and keyed random algorithms are unchanged by the adapter.
The same snapshot must cover the complete parent interval; replacement after
freeze or closing while particle readers/workers remain active is rejected.
Closure also requires the PIC clock at its declared end and no pending birth.
Failed closure preserves the active window and stored histories for diagnosis
or a valid retry.

All three movers resolve their effective duration before coefficient or random
draws. A pending uniform birth uses `(1-phase)*parent_duration`, once, instead
of receiving both helper pre-drift and a full host step. The remaining duration
is authoritative even when its absolute departure epoch rounds to the parent
end at a large PIC epoch. Snapshot time remaining is evaluated as
`(snapshot_until-parent_end)+(first_duration-elapsed)`, retaining that small
interval independently of the PIC counter's start/end update order. A declared
positive terminal first interval below the adaptive refinement floor is
completed whole; a physical/stability limit that requires a smaller interval
still returns `StepUnderflow`. No phase is clipped, redrawn, or discarded.
Every later substep of a move follows the same rule
(`PICAdapter::SelectParticleSubstep`): the movers sum their substeps in
floating point, so after a substep set by another limit (for example the
snapshot boundary, computed from absolute epochs) the time left in the move
can be a round-off residue of about one ulp of the step (7e-15 s of a 60 s
step). That remainder is completed whole instead of being reported as a
stability underflow; a limit smaller than the remainder still returns
`StepUnderflow`, and a move that leaves no such remainder selects the same
substeps as before.

Successful phase-space, complete RNG/hazard history, and first-move consumption
are committed together. A statistical split copies the parent's physical birth
provenance to both representatives; the child receives its separately audited
lineage identity and future stream without inventing a new physical age.
Weight-only merging preserves the birth records of surviving phase points.
Native migration copies pending/consumed history exactly. Completed physics
checkpoints reject active driver windows and any pending first move; schema-3
records and incompatible layouts are rejected without migration or guessed age.
The retained checkpoint codec also rejects a live four-polarization authority
until its complete serializer and collective restore transaction are qualified.

Current nonthermal field-line sources declare their existing no-pre-drift,
step-start timing explicitly before reserving source counters. Two-argument
fresh initialization is initial-condition-only and cannot run inside an active
source window. The new overload accepts an explicit birth record and existing
local-plasma-frame phase space; it does not infer or transform a velocity frame.
The fresh coupled SWMF volume initializer is a separate path from boundary
injection. It consumes the received proton density, temperature and bulk velocity,
requires a single compiled H+ species, and samples a Maxwell–Jüttner **volume**
distribution in the local plasma frame. For proper density `n0`, simulation-volume
density is `gamma_U*n0`; the rest-momentum draw includes the Lorentz phase-measure
factor `1+U dot v'/c^2`. Spatial weights use the existing magnetic-flux area
`A=Phi/|B|` and the full Jacobian `n0*gamma_U*A*L`. Simpson quadrature controls
proposal counts only; each sampled point receives its actual importance weight,
so count rounding and the 200-macro proposal ceiling do not renormalize density.
Purpose-keyed source streams and persisted line-event identities initialize the
complete transport schema as an explicit initial condition. Existing records
cannot acquire invented source/RNG history through this path. MJV01–03 verify
the momentum kernel against independent quadrature and Lorentz moments; native
population normalization and the full coupled regression are separate gates.

Maxwellian beginning and optional thermal solar-wind boundary helpers remain gated before
RNG/counter/allocation mutation. The user selected drifting relativistic
Maxwell–Jüttner flux injection with verification of its Maxwellian limit.
The separate full-particle material-plane kernel is now implemented and its
`MJ01`–`MJ06` strict-warning/ASan/UBSan callbacks pass **6/6**, including
the actual cold flux-weighted moments and expanded cold/high-Gamma rejection.
Independent review's cold-boost counterexample failed the old code; its new
precision preflight/error budget now rejects unresolved encoded crossings
before RNG/output commit. The expanded MJ01–06 also pass natively in build8
on two MPI ranks.
Its explicit proper density/temperature,
bulk velocity and boundary velocity/normal inputs do not establish their native
provider mapping. The source count, material face geometry, time integration,
native frame/storage and guiding-center boundary bridge remain separate from
temporal bookkeeping. The retained classical/truncated helpers
are not the selected law and are not enabled by new timing metadata.

`util/sep_maxwell_juttner_flux.*` samples the exact material-plane radial
gamma-two/gamma-three mixture and cosine-weighted proper-frame hemisphere,
retaining full rest and boosted momenta. It stages six keyed uniforms per
event and commits RNG/output together; an unrepresentable event returns an
explicit status without a tail redraw or clipping. The SI rate uses proper
density and scaled K2 normalization. Independent high-precision references,
original positive integrals, boosted moving-plane quadrature and campaign/
repeat/order controls pass. A tilted cap's full-particle gyromotion crossing
flux cannot simply be projected into the gyroaveraged mover state. This API
enables no native source until a matching reduced boundary law or explicit
proper-frame cap/B alignment domain is established. Numerical methods, tested
range, physical limits and reproducible commands are in
[MAXWELL_JUTTNER_FLUX_IMPLEMENTATION.md](validation/physics/MAXWELL_JUTTNER_FLUX_IMPLEMENTATION.md).

The retained analytic moving-shock query still advances its trajectory using
the mutable PIC epoch. Its source entry now compares that query epoch to the
declared frozen departure and rejects a mismatch before even reserving a source
event. This affects the normal counter-OFF/end-counter source ordering; fixed
boundary sources remain available. No end-epoch shock is labeled as a
step-start birth. A pure trajectory query at departure plus a separate end-step
geometry/energy-ledger update is required before enabling the affected path;
no new temporal quadrature is inferred by this guard.

`BIRTH01`–`BIRTH11` exercise the shared value core under ASan/UBSan, including
the late-epoch counterexample and atomic rejection. The new native descriptors
`BIRTHN01`–`BIRTHN06` compare actual Parker/Dmumu/MFP remaining-time movement
against explicit-duration references, retain a full-host negative control, and
exercise read-phase close rejection, a tiny positive MFP interval, actual
clone/split/MPI history preservation, and pending-checkpoint rejection. These
native schema-4 cases pass **6/6 on one and two MPI ranks in build6b** and the
expanded Cartesian-view/atomic rejection gates pass **6/6 on two ranks in
build8**. Initial
build6 evidence is retained: `BIRTHN04` passed while five fixtures correctly hit
duplicate-publication rejection at startup's epoch/generation. Static fixtures
now explicitly run one second after the saved startup epoch; the immutable
publication invariant is unchanged. These cases execute actual mover managers,
core clone/split and MPI migration, with an explicitly controlled PIC clock.
An actual `PIC::TimeStep`/host-driver counter-advancement gate remains separate,
as do real moving-shock and relativistic thermal-source qualification. See
[STATUS.md](STATUS.md) for current gates.
The schema-4 CP05 cases and separate fresh-process pairs also pass in build8;
those software diagnostics do not qualify actual thermal-source physics.

### Conserved-flux checkpoint owner and process boundary

The corrected native codec's inner schema 4 captures actual flux-owner
state: configured reference seed area [m²], per-line conserved magnetic flux
[Wb], each record's generation/provenance, and exact lazy/uninitialized table
slots. Capture is observational; it does not call the lazy flux getter. Incoming
records are validated against the fixed physical line count and frozen seed
area before staging. Decoded B/flux must also give finite positive production
Simpson areas/volumes. An opaque area callback rejects Save/Load before any
publication or native replacement; its address is neither serialized nor
silently cleared. Staged vectors own all allocations, a nonphysical mutation
serial detects stale plans, and the collective commit swaps ready values.
Replicated epoch/configuration/source/flux state must agree across ranks on
both save and load. This schema also persists fieldwise
vertex XYZ triples outside their native associated buffers. They must be finite
and exactly match the declared fixed geometry; they are the decoded oracle for
finite/exact particle Cartesian views. A drawn residual optical depth must be
strictly positive, matching `LoadParticle`; zero and negative zero reject before
replacement. Unknown earlier corrected schemas are rejected.

`FLCP01`–`FLCP04` pass the dependency-light fieldwise-codec/validation tests with
ASan/UBSan. New native `CP05C05` exercises the actual value owner, stale commit
and opaque callback. Paired `CP05F01`–`CP05F03` use separate MPI processes:
the writer installs an explicit flux/provenance different from startup, saves
after one step, and records uninterrupted final state; the reader configures
operators without warming caches or creating sources, loads, advances, and
compares final checkpoint bytes. Hardened schema-4 `CP05C01`–`CP05C05` pass
**5/5 on two ranks in build8**, including CRC-valid Cartesian/hazard negatives.
All three separate-process `CP05F01`–`CP05F03` writer/reader pairs pass with
byte-identical final artifacts and actual nondefault flux restoration. This
qualifies the bounded fixed-geometry diagnostic API, not a general production
resume interface or real-source integration.

The mover commit now preflights canonical `GetCartesian(coordinate)` and
publishes PB::X together with phase, RNG and birth state. The coordinate remains
the sole authority. Build7a's independent tests first reproduced stale views
after all three movers and partial publication on nonfinite geometry; build8
passes these independent assertions without tolerance changes. Future material geometry
publication must also refresh derived views for records whose coordinate does
not move; the fixed-geometry diagnostic does not establish that transaction.

Run the paired API fixture through the existing test-input interface:

```sh
python3 -B srcSEP/test/run_fresh_process_checkpoint_tests.py \
  --amps ./amps --mpi-ranks 2 --output /tmp/sep-fresh-checkpoint-new
```

The output directory must be new/empty; every role owns separate retained
scratch and structured reports. This qualifies only the explicitly bounded
fixed analytic diagnostic API when it passes. It is separate from a production
resume CLI, material-grid geometry continuation, four-polarization persistence,
and output/sample continuation. Current execution gates remain in
[STATUS.md](STATUS.md).

The build8 evidence point is binary SHA256
`049a1959918dd154e5a81d75ed7ebb18a57764f29c09c3917c81b4068f9c58eb`,
with its captured 162 native/23 canonical source manifest. The combined two-rank
campaign passed 116/116 without SKIP; the fresh checkpoint cases used separate
writer/reader processes with that same binary. The source hashes/reports are
tracked by STATUS. The later independently exposed extreme-shear material-flux
conditioning case remains outside this baseline pass and does not acquire
qualification from its six original MTG callbacks.

WP24 additionally makes key composition ordered and purpose tagged. The source
tuple is `(campaign,event,field-line,species,macroparticle,purpose)`; swapping
two fields changes the stream, and count/spectrum/pitch/gyrophase draws cannot
perturb one another. The primary field-line shock source no longer consumes the
process-global RNG.

## Shock/source and diagnostic numerical contracts

WP21 integrates a piecewise-linear `v_shock(r)` exactly from a serialized launch
datum, and WP22 intersects its spherical surface with the complete 3-D
field-line polyline using oriented/tolerance-aware typed results. WP25 evaluates
the swept shell volume and upstream normal-relative kinetic energy, queues the
result, and lets the authoritative turbulence core apply and ledger it.

WP27 diagnostics reject invalid or luminal particle state rather than clipping
it. Relativistic Larmor radius is `p_perp/(|q|B_local)`. WP28 uses half-open bins
with separate underflow/overflow/invalid counters, preserves `sum(w^2)` for
effective sample size/error, and constructs number density, omnidirectional
differential intensity, directional crossing flux, and normalized shape as
separate unit-bearing products.

`CampaignRandomSeed` is the production adapter boundary. A campaign that needs
a nonzero seed must set it before particle motion; leaving it at zero is still
deterministic and explicit.

## Parker equation

Fixed analytic generated-domain startup now uses one explicit radial origin and
speed for list U and the Cartesian tensor
`G=(v/r)(I-rhat*rhat)`. Parker/straight generators declare U; native registration
preserves declared list/imported U. Shared startup preflights the whole analytic
realization before applying the existing inverse-square density law and
identical initial previous B/U/n. The runtime initializer is the configuration
authority; imported/time-dependent providers and thermal profiles are not
synthesized. `ARB01`–`ARB05` and `ARBN01`–`ARBN04` track pure/native scope;
production closure and material-grid qualification remain separate. The
coefficient provider now consumes the explicit radius/frame in the same
background view: analytic `|x-origin|`, declared SWCME/reduced global radius,
and unavailable SWMF/handoff origin. This is established independently of
controlled derivative fixtures, with no new endpoint-U validation. Its versioned
frame contract enters the existing background fingerprint and coefficient
provenance; matching wire layout does not authorize an older physical
interpretation. The optional parallel-diffusion bridge consumes the view's
solar-relative Cartesian vector in the unchanged B/U basis and rejects an
unknown origin before its numerical stencil wrapper. `RADIUSN01`/`RADIUSN02`
use explicitly tagged schema6 `--test-input` in separate actual immutable
startup processes, through the same LoadFile/Install owner. Their independent
QLT1 and speed/radius power-law references have strong global-zero negative
controls. Production `--input` remains incompatible with test mode. Both
old-producer callbacks fail as intended on two MPI ranks in isolated BUILD10a:
QLT1 and the optional library match their independently computed old-global
negative controls. Post-fix native execution remains pending; see
[PROVIDER_RADIUS_CONTRACT.md](validation/physics/PROVIDER_RADIUS_CONTRACT.md).

`SpatialDiffusionProvider` returns `kappa_parallel` [m²/s],
`d(kappa_parallel)/ds` [m/s], status, and provenance together. The Itô update is

`ds = [U_parallel + d(kappa_parallel)/ds] dt
      + sqrt(2 kappa_parallel dt) dW`.

The production provider delegates to the configured srcSEP `GetDxx` contract;
a missing configured diffusion callback is the explicit zero-diffusion case.
Convection, variable-coefficient drift, Gaussian noise, and exact cooling are
composed in `AdvanceParker`. Segment-scale drift and diffusion constraints are
applied outside the pure core by the PIC adapter.

## Coefficient-driven focused transport

`PitchAngleDiffusionProvider` returns `Dmumu` [s^-1], `dDmumu/dmu` [s^-1],
status, provenance, and the immutable turbulence-state identity used to obtain
the coefficient. A configured analytical derivative is used directly; the
numerical mode evaluates a symmetric stencil except for one-sided endpoints.

One substep uses a symmetric deterministic/stochastic/deterministic split:

1. half of the focusing, parallel-flow-gradient, and cooling operators;
2. one Itô kick `dmu = (dDmumu/dmu)dt + sqrt(2 Dmumu dt)dW`;
3. the remaining deterministic half; and
4. streaming with the midpoint pitch angle and local plasma advection.

The deterministic pitch-angle half-steps use explicit midpoint integration.
The production adapter limits streaming, focusing, diffusion amplitude, and
diffusion drift separately.

## Event-driven mean-free-path focused transport

`MeanFreePathProvider` returns `lambda_parallel` [m], branch rates `nu+` and
`nu-` [s^-1], status, provenance, optional momentum validity bounds, and the
turbulence-state identity. Positive infinite lambda is the explicit zero-rate
ballistic limit. A particle stores the remaining unit-exponential optical depth
`tau=-log(xi)`. Each deterministic interval integrates the nonhomogeneous
hazard `integral(nu+ + nu-) dt`; a monotone root solve locates an event inside
the interval, and only a completed event draws the next optical depth. Segment,
shock, snapshot, and outer-timestep partitions therefore do not resample the
physical process.

Between events, the mover applies midpoint focusing and parallel-flow-gradient
drift, exact plasma-frame cooling, and midpoint streaming. At an event the
branch is sampled with probabilities `nu+/(nu+ + nu-)` and
`nu-/(nu+ + nu-)`; an empty branch is never selected. The particle direction
is redistributed isotropically in that signed Alfvén wave frame. Exact Lorentz
velocity transforms preserve wave-frame speed and keep the returned
plasma-frame state subluminal.

## Unified coefficient registry

`util/sep_coefficient_registry.*` publishes canonical names, SI units, and
parameter schemas for spatial, pitch-angle, and mean-free-path coefficients.
The PIC-facing providers add background ownership validation and provenance.
The only supported conversions are centralized as `kappa=v*lambda/3` and the
explicit isotropic `Dmumu` closure. A spatial-from-MFP/MFP-from-spatial cycle is
rejected. Self-consistent and SWMF MFP transport must use `from-spatial`; an
SWMF source is valid only with a read-only SWMF snapshot.

Invalid lambda output follows an explicit policy. `fail` returns a structured
error; `ballistic` substitutes positive infinity and increments both invalid-
sample and ballistic-substitution counters. The canonical mover has no silent
minimum/maximum lambda clamp.

## QLT coefficient normalization

### Standard BATSRUS imported wave view

The user selected standard BATSRUS exports: raw Wave1/Wave2 contain **total**
Alfvén-wave energy [J/m³], Wave1 propagates along actual B and Wave2 opposite B.
`util/sep_standard_batsrus_waves.*` converts to a private local magnetic view
with `W_B+=Wave1/2`, `W_B-=Wave2/2`; the canonical magnetic conversion remains
`deltaB_branch²=2 mu0 W_B,branch`. Thus the independent raw reference is
`deltaB²=mu0(Wave1+Wave2)`. Local integrated/spectral magnetic energy is
unchanged. This explicit producer contract supersedes the earlier claim that
raw SWMF columns themselves were magnetic-only.

`ReadStandardBatsrusWaveView` is the actual coefficient-path native accessor.
It requires an active read-only SWMF snapshot matching the context generation,
valid coordinate/segment and allocated two-field vertex storage. Both endpoints
are validated before interpolation, including a zero-weight endpoint; invalid,
nonfinite and unrepresentable positive energy/variance reject before output
replacement. Imported PIC arrays remain total energy and are never halved in
place. The standard contract, raw endpoints, fraction and generation bind the
view checksum; representation/checksum feed coefficient provenance, and the
SWMF background fingerprint uses the same contract authority.

Wave labels are relative to B, so neither field polarity nor geometric line
ordering swaps Wave1 and Wave2. No circular polarization or spectral shape is
invented; two integrated exports cannot initialize four independent imported
spectra. A local evolution handoff needs its own explicit exactly-once
conversion/owner transaction. This adapter does not enable that handoff or
resolve SWMF velocity gradients, composition/density/temperature basis, frames,
material geometry or the full focused-transport provider.

`SBW01`–`SBW05` pass the strict-warning/ASan/UBSan pure contract suite;
`SBWN01`–`SBWN02` use actual native vertex arrays and this same accessor, with
controlled read-only metadata, polarity/endpoint-order controls, raw-byte and
RNG/PB-count preservation. All five pure and both native cases pass on two
MPI ranks in build8; the pure sanitizer gate also passes. `WN02` remains a
synthetic already-magnetic, unknown-direction view; it makes no raw SWMF claim
and retains the generic unknown-mapping controls. Details and commands are in
[STANDARD_BATSRUS_WAVE_ADAPTER.md](validation/physics/STANDARD_BATSRUS_WAVE_ADAPTER.md).

For `P(k)=C k^(-5/3)` on `[kmin,kmax]`, the normalization is

`C=(2/3) deltaB^2/(kmin^(-2/3)-kmax^(-2/3))`,

so `integral P(k) dk = deltaB^2`. With
`k_res=Omega/(v |mu|)`, the slab expression is

`Dmumu=(pi/2) Omega^2 (1-mu^2) P(k_res)/(B^2 v |mu|)` [s^-1].

An invalid input or a resonance outside the represented band returns zero from
the historical QLT API. Provider adapters subsequently reject any non-finite or
negative coefficient before advancing a particle.

## Full focused equation and authoritative background view

Production `fte-dmumu` and `fte-mfp` use the local plasma-frame gyrotropic SDE

`dmu/dt = (1-mu^2)/2[-v dln|B|/ds + mu(divU - 3 bb:gradU)] + dDmumu/dmu`,

`dln(p)/dt = -1/2[(1-mu^2)divU + (3mu^2-1)bb:gradU]`.

The stochastic pitch increment is `sqrt(2 Dmumu dt) dW`. The explicit
`ReducedFieldAligned1D` mode preserves the earlier one-dimensional closure for
documented comparisons, but production adapters select `FullGyrotropic`.
`LocalBackgroundView` is the single SI record supplying magnetic field/vector,
unit field direction, line tangent, plasma velocity, number/mass density,
`divU`, `bb:gradU`, Alfvén speed, provider provenance, generation, and validity.
Providers evaluate requested physical arc-length displacements and the Dmumu
core samples its predicted midpoint rather than a hard-coded location.

Density-derived divergence uses the interval between physical background
epochs stored by `BackgroundSnapshot`, never a particle timestep. All three
movers compose the same snapshot validity/generation limit before every shell;
a stale generation or an attempt beyond the validity endpoint is an error.

## Wave feedback and concurrency

The numerical core deposits only a signed streaming record labelled with the
same turbulence identity carried by its coefficient. PIC-facing movers never
write shared G+/G- wave arrays from worker threads. They append a typed,
self-contained record after every completed deterministic interval/event.
Species and statistical weight are copied while the particle is alive; the
post-step queue never retains or dereferences a PIC handle. If an absorbing
boundary is crossed, the attempted trajectory is clipped to the exact endpoint,
interval time is scaled to the in-domain fraction, the boundary record is
enqueued, and only then is the particle deleted.

Step 12 formalizes the next reduction boundary in
`util/sep_reproducible_reduction.*`. Thread-local records carry a versioned
schema plus source, field line, segment, branch, spectral bin, species, stable
particle, step, event, interval, and purpose keys. Duplicate complete keys are
errors. After workers join (and after rank-local records are gathered), a
complete-key sort establishes one fixed
accumulation order with long-double intermediates. Non-additive authoritative
wave state uses owner/gather synchronization, never an `MPI_Allreduce`
overwrite. Canonical evidence hashes are bitwise invariant to the tested worker
and synthetic rank partitions.

The retained turbulence driver, ownership, boundary, ledger, remap, and restart
contract is documented in [TURBULENCE_MODEL.md](TURBULENCE_MODEL.md).

## WP33 population control and WP36 system accounting

`RunConfiguration::populationControl` owns the actual spatial/momentum/pitch
bin counts, minimum/maximum population, deterministic-selection requirement,
invariant tolerance, and lineage schema passed to PIC merging and splitting.
`sep_population_control.*` provides the independent SI moment observer used by
controlled and native fixtures. Split children inherit phase-space state and
partition weight with a final-child residual; stable child IDs hash parent,
generation, ordinal, and schema without rank or thread.

`sep_system_ledger.*` closes number, charge, energy, and parallel momentum over
named signed seam transactions. Its checksum-protected checkpoint carries
field-line generation, stable IDs, lineage, RNG counters, residual scattering
hazards, turbulence-state identity, and output-manifest identity. This source
contract is ready for native adapters; complete injection/escape/coupling/remap/
restart equality still requires a linked AMPS campaign.

## Verification boundary

The dependency-light suites are `make test-transport-common-unit`,
`make test-parker-unit`, `make test-fte-dmumu-unit`,
`make test-fte-mfp-unit`, `make test-coefficients-unit`,
`make test-turbulence-core-unit`, and `make test-reproducibility-unit`. They use strict C++11
warnings, AddressSanitizer, and UndefinedBehaviorSanitizer and cover `CORE01`–
`CORE07`, `PARK01`–`PARK07`, `FTED01`–`FTED08`, `FTEM01`–`FTEM08`, and
`COEF01`–`COEF05`, `TURB02`–`TURB23`, `TURBOWN01`, and `PAR01`–`PAR05`.
The mover and turbulence cases are descriptors in the same registry consumed
by `--test` and `--test-group`; their source-only runners no longer carry
parallel copies of the callbacks. `FTED08` adds the analytical Legendre-mode
decay reference, while `TURB21` and `TURB22` add translated-profile advection
and time-dependent wave-energy-source references. `TURB22` verifies how a
known source is integrated and accounted; it does not replace the native
particle-distribution-to-growth-rate test. Source checks additionally
verify public dispatch and the common adapter boundary.

Linked-application cases CV02-CV05 exercise these same pure cores through the
production `amps --test` registry. CV02 checks the constant-kappa Green
function over independent seed, particle-count, and timestep refinements. CV03
checks the conservative Ito `d(kappa)/ds` drift, uniform equilibrium, analytic
versus numerical gradients, and a reversed-drift control against an independent
finite-volume solver. CV04 verifies frozen-divergence cooling and second-order
spherical midpoint convergence across species and energy. CV05 verifies
constant-gradient focusing, pitch-angle bounds, and `(1-mu^2)/B` conservation.
Run them with `make test-cv02-cv05-unit SEP_EXECUTABLE=/path/to/amps`; absence
of that linked binary is SKIP and is never replaced by source-only evidence.

Refinement evidence is evaluated uniformly as
`p=log(e_coarse/e_fine)/log(h_coarse/h_fine)`. The common helper rejects zero,
negative, nonfinite, or reversed-resolution inputs so an empty/error-free
special case cannot masquerade as convergence. `PARK07` adds exact absorbing
first-passage statistics, and `FTEM08` adds the finite-time persistent-flight
MSD and its `kappa_parallel=v*lambda_parallel/3` diffusion asymptote.

PIC-facing adapter, shell, MPI, restart, SWMF, long-run conservation and
observational gates require their actual configured application runs. Source
checks and dependency-light tests do not replace them. The current native build
and campaign dispositions are recorded in [STATUS.md](STATUS.md); each passed
software diagnostic remains distinct from production/scientific qualification.


## Reduced-core Parker material integration

For physical-number representatives, the forward Parker measure is
`N ds dp = 4 pi p^2 A f ds dp`, with conserved tube flux `A |B|=Phi`.
Transforming `(1/A) partial_s(A kappa partial_s f)` to the forward
number-density equation gives the Itô drift `partial_s kappa + kappa
partial_s ln A`. Consequently the material-line step is

```text
ds = [partial_s kappa - kappa partial_s ln|B|] dt + sqrt(2 kappa) dW
p_new = p_old exp[-div(U) dt/3].
```

The segment-plus-fraction conversion supplies units/arc metric only; it cannot
replace this tube-area term. No extra particle weight renormalization or U
coordinate drift is applied. The native shell includes this same drift when
selecting its internal geometry limit. The core helper defaults to constant
area for existing manufactured constant-area fixtures; the native flux tube
passes its actual log-area gradient explicitly. Invalid/nonfinite metric and
overflow reject before a random draw.

Both focused solvers retain the full plasma-frame gyrotropic coefficients:
`dmu/dt=(1-mu^2)/2[-v partial_s ln|B|+mu(div(U)-3 bb:grad(U))]` and
`dln(p)/dt=-[(1-mu^2)div(U)+(3mu^2-1)bb:grad(U)]/2`. The latter averages to
Parker cooling for isotropic pitch. Physical U, its tensor, Alfvén speed and
shock inflow remain available even though relative coordinate advection is
zero. Internal subdivision refines these operators at the same physical
endpoint; it does not change physical dt or the sampling interval.

The reduced native owner prepares geometry over the completed full driver
window. Its read-only interval view maps every departure vertex to an explicit
query time before finding finite SSE intersections. The returned segment
fraction labels that same material edge at the requested time. Prefix
insertion changes live indices only at commit; stage 5 sources must map these
interval labels to the committed topology rather than mix future fronts with
old positions. Derived vertex shock distances are diagnostics. The reduced
MFP limiter uses the complete supported root list, with cut speed
`(Vn-U1 dot n)/(tangent dot n)` relative to material, avoiding a false zero
where two nearest-root diagnostic branches switch. Tangency returns explicit
unsupported proximity estimation rather than an obliquity floor. This limiter
is not a stochastic first-passage solver or an acceleration prescription.

SINT06 verifies core material continuity and induction by an independent
worldline derivative; SINT07R checks curved-field tangent convergence after a
fixed 600-s motion. SINT14 grades the tube drift against exact conditional and
Gaussian moments. SINT15 grades focusing/cooling against analytical
coefficients, angular quadrature and constant-gradient exact characteristics.
SINT16 uses equivalent `D_mumu=nu(1-mu^2)/2`, `lambda=v/nu` and
`kappa=v^2/(3 nu)`. It compares independent finite-time persistent-flight and
Parker second moments with fixed seeds/confidence criteria; different early
angular distributions of the two scattering processes are preserved.

Actual native SINT14N–16N checks use the unchanged 300-s example step and
manufactured slow ions on the real core ambient. The MFP path is explicitly
conditioned to survive without a scattering event; existing FTEM tests grade
waiting times, residual hazards and wave-frame conservation. These diagnostics
qualify the native adapters with the pure solver tests, not a physical shock
source, turbulent feedback or downstream transport.

The native SINT15N comparison exposed a separate time-centering error in
`fte-dmumu`: using the final cooled speed throughout streaming created a
first-order displacement bias. Each deterministic half now advances log(p)
and mu together by midpoint, predicts coefficient position using half-step
speed, and integrates parallel velocity by endpoint trapezoid. SINT15 adds an
exact relativistic exponentially cooling streaming reference and a second-order
internal-refinement check. The same full AMPS timestep and original native
trajectory bounds are retained; no output/sample opportunity is inserted.
# Finite-reference reduced source and absorbing return

The core Parker exterior skips winding-angle quadrature only for an exact
Cartesian rotation-axis position: rotating that direction about the axis is
the identity for any magnetogram. It still evaluates the same PFSS footpoint,
full plasma primitive, magnetic-null condition and declared coverage. Nonzero
transverse coordinates retain the ordinary winding integration, including
finite-difference queries near the axis. SINT26 checks independent axial-dipole
and radial-wind relations; this optimization introduces no constant-wind or
radial-field approximation to the general provider.

The core MHD jump solver similarly uses the exact gas-dynamic compression
identity only to locate two of its existing scan-grid candidates when both
transverse B and incident velocity are componentwise zero. Both candidates
must still be physical and bracket the energy residual. The unchanged
bisection and all entropy, characteristic-ordering and RH-residual checks
then determine admission. Nonzero transverse components, weak roots near the
identity branch or an unsuccessful proposed bracket use the original scan.
SINT27 grades independent parallel-shock identities; the core RH3D catalog
grades oblique, perpendicular, weak and inadmissible cases. No shock model,
compression cutoff or residual tolerance changes.

Repeated native boundary queries memoize successful roots per worker using the
retained immutable material-view owner and the exact floating-point epoch.
There is no timestamp quantization, trajectory interpolation or shared mutable
cache. A geometry publication clears the worker's entries, including a common
physical endpoint whose arc label changed after prefix insertion. Retaining
the shared owner prevents address reuse from aliasing an old view. Invalid
epochs still take the provider's strict error path; the 4096-entry memory bound
only causes a full recomputation on eviction. Native SINT25N–28N compare a
repeated query bitwise and retain the independent full-intersection 1-mm bound;
their second physical step also exercises owner replacement after motion.

An application-owned net first-passage source reuses the core DSA number law
`dN/dp ∝ p^(2-3r/(r-1))`, with local core RH compression `r`. Its physical
normal offset `L_ref > 0` is fixed under numerical refinement. Placement is an
intersection of the same material line with the finite leading sphere's
upstream offset at the actual birth epoch. A normal Cartesian displacement of
the connected shock point would generally leave an oblique magnetic line.

Birth proposals stratify the complete physical window: `t_j=t0+dt*(j+u_j)/N`.
Their physical weights are `R(t_j)*dt/N`. This is unbiased Monte Carlo time
integration of a prescribed rate, not an assumption that `R` is constant or a
Poisson physical counting model. Disconnected/sub-fast proposals create no
particles. Unsupported or failed preparations return errors before allocation
and source-counter publication. Zero efficiency creates no source particles.

The positive-normal-flux angular law is proportional to
`max(U_ref·n - V_ref,n + v*polarity*mu*(b·n),0)/2`. The core inserts its exact
zero-flux threshold into the integration nodes and inverts the integrated
linear density, including partial loss cones. It does not interpolate a
quadratic CDF linearly or redraw unsupported angles. Parker retains the same
source's angular marginal and isotropic diagnostic projection; the focused
equations retain the escaping angular distribution.

The admitted **net** spectrum requires positive outward angular support at its
minimum momentum at the actual reference point (or at the explicit verification
monoenergy). The maximum normal particle velocity is monotone in speed, so this
qualifies the finite band. An empty cone rejects the source transaction before
allocation and counter commit. Dropping slow proposals would reduce the net
rate; redrawing or renormalizing would change its prescribed DSA marginal.
Neither operation supplies the missing gross-to-net escape law. SINT32N uses
an incompatible 1–2 eV band and requires atomic rejection. The retained
`no_focused_escape` ledger field stays zero for an admitted net source.

MFP boundary callbacks receive explicit beginning/end offsets in each numerical
interval. The final end is the caller's exact closed endpoint. After a scatter,
`elapsed+(duration-elapsed)` can exceed that endpoint by an ULP. The core now
carries the endpoint through its callback and commits it directly, without
clamping queries or increasing the boundary tolerance. Absorbed clocks
interpolate from the nearer endpoint and preserve fractions zero/one exactly.
SINT28 uses 4096 fixed-key non-dyadic internal durations and early event times:
the retained old kernel misses 66 endpoints; the correction misses none.
Physical AMPS dt and the sampling schedule remain unchanged.

Material trajectory RK4 now bisects the complete requested interval for both
maximum-step and accuracy refinement, retaining adjacent exact endpoints.
Repeated maximum-size steps formerly left an arbitrarily short final tail;
at heliocentric coordinates its proportional error budget could be smaller
than coordinate roundoff, and further halving made that budget smaller still.
The captured iteration-84 case exhausted 100000 attempts with about 0.00049 s
remaining. The endpoint-preserving partition removes this artificial tail.
The local estimate still receives positionTolerance*abs(h)/totalDuration;
there is no roundoff acceptance floor or unused-budget redistribution.
Maximum-step subdivision also consumes the unchanged work bound. SINT30
reproduces the old failure and grades the correction against a separately
converged midpoint trajectory. Internal RK intervals do not change physical
AMPS dt or segment sampling exposure.

The polar shock-root search retains the actual first/last advected endpoints
and all provider validity checks. It uses the departure mesh only to suggest
a narrower search bracket, validates both suggested vertices with the same
core trajectory map at the exact query epoch, and otherwise continues binary
search in a proved range. Native tests require exact equality with the retained
full-range binary search as well as the existing full-polyline 1-mm check.

The immutable core ambient memoizes canonical PFSS results only at exactly
equal axial source points. It retains both hemispheres and all signed-zero
azimuthal limits used by the existing Cartesian pole evaluator. Five nearby
representable radii retain separate canonical results; unmatched feet use
the ordinary PFSS query. Nonaxisymmetric SINT31 controls
compare with uncached canonical evaluations and retain magnetic-null rejection.
No time, position or field value is rounded to obtain a cache hit.

Shock diagnostics resolve the same reduced authority as transport. The
`reduced-SSE/Parker` title's `R_sh` is the heliocentric leading apex radius,
evaluated at the actual output epoch, rather than the SSE generating-sphere
radius or an inactive legacy shock's radius. The shared diagnostic query also
serves spectral output. Native tests query a non-cadence epoch so an old
prepared-front cache cannot satisfy the metadata check accidentally.

Field-line vertex primitives store B in tesla. The established text columns
advertise nanotesla, so serialization multiplies current and previous B by
1e9. This is an output conversion only; coefficients, transport, geometry and
live vertex records retain SI values. The native fixture compares the actual
written column against live B with the writer's declared text precision.

On material coordinates, relative particle displacement contains no added U.
The moving shock's coordinate follows its HCI trajectory relative to advected
vertices. The background's actual U remains in front-relative flux, strain,
cooling and scattering-frame physics. Existing departure-snapshot transport
composition remains in force; explicit-time source/front geometry does not
change the background snapshot's ownership or turn RH limits into a volume.

For a frozen-coefficient Parker numerical interval with signed upstream
endpoint distances d0,d1 and K=kappa*dt, a bridge with two positive endpoints
hits with probability `exp(-d0*d1/K)`. The conditional first-hit time is sampled
without replacing the Gaussian endpoint. Endpoint signs alone miss returned
paths. The linear moving-boundary approximation is refined before stochastic
draws using quarter/midpoint material-root checks and an explicit SI tolerance;
the host timestep and sample count remain unchanged. Internal interval/curve
refinement and coefficient convergence are numerical qualifications, not an
observational or universal diffusion-approximation qualification.

Each particle move retains its last checked boundary interval as a numerical
proposal hint. Every new proposal still evaluates the actual start/end and
quarter/midpoint material roots with the same spatial tolerance, topology
checks, 64-attempt budget and minimum-step/progress guards, before any particle
noise/event draw. After eight checked intervals a factor-two growth probe
allows recovery from a short interval near a mesh vertex. Reusing a hint never
authorizes an unchecked future interval or interpolated background.

The history is local to one particle move and resets at its next physical
call. It changes no AMPS timestep, sampling exposure, birth endpoint or
persistent RNG record. The internal numerical mesh can change individual
fixed-seed trajectories; the conditional first-passage laws and numerical error
bounds remain fixed. SINT32 grades the production refiner with s(t)=t^2:
128 identical accepted intervals at a 1-mm bound require 149 checks with the
hint versus 1024 cold checks. Failed preparation must leave history unchanged.

Focused intervals resolve the first crossing of their accepted streaming
interpolant. The MFP kernel invokes the boundary for every interval between
events, before the endpoint scattering operation. An absorbed particle draws
no later scattering angle, branch or optical depth. Its terminal momentum uses
the same interval's logarithmic-momentum interpolant. Existing field-line
endpoint exit rules are retained and recorded separately from shock return.

The joined ledger checks candidate minus no-focused-escape exclusions against
committed source number and birth energy, then checks committed equals retained
plus shock-return, domain-exit and geometry-retirement measures. Birth energy
is immutable per physical representative; event energy is evaluated at the
actual numerical loss. Signed deterministic flow work and prescribed wave-frame
scattering work are separate. Source energy is an explicitly supplied already
accelerated population, not a newly resolved accelerator energy budget.
The ledger prints retained birth energy separately from retained current
kinetic energy, so external checks can recompute both balances independently.

The first particle consumer admits one transverse supported leading front.
Tangency, multiple simultaneous absorbing fronts, occupied connectivity/
criticality changes and unsupported downstream continuation return explicit
errors. They do not receive obliquity floors, guessed branches or ambient
continuation. The existing conservative merger lacks a birth-energy invariant;
if its actual trigger is reached for this source, it rejects before changing
weights. Cohort-preserving resampling, remapping and moving restart belong to
stage 7. Exact splitting already copies the complete cohort record.


## Common material admission and the fixed physical timestep

The stage-11 neutral interval requires `end == begin + configured_dt` and full
physical coverage before native material advancement. It never shortens a
physical step to a cache cadence or event endpoint. Internal trajectory
quadrature remains the established endpoint-preserving adaptive RK method;
its controls do not change the sampling timestep. Relative particle transport
and moving field-line geometry retain their existing frame convention.

Capabilities distinguish primitives, species/EOS, full grad(U)/grad(B), dense
time, trajectory, inverse/deformation, local shock/contact limits and physical
downstream volume. An unsupported quantity, stale realization, unpreserved
branch or out-of-coverage query returns explicit status. The native snapshot
store remains the sole publication/read-phase authority. No velocity, density,
heating, mapping or interpolation law is invented by the neutral adapter.

MBC15–19 compare the reduced export with the separately constructed canonical
provider, charge/mass/ideal-gas identities, independent radial flight-time
quadrature/bisection in the actual accelerating Parker wind and actual local
shock states. These are synthetic/canonical verification, not observations.
The original physical residual, sampling, conservation and continuation gates
remain controlling; see `validation/STAGE11_MATERIAL_BACKGROUND_CONTRACT.md`.
