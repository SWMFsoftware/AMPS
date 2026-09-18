# Immutable background snapshots and the simulation clock

Step 2 establishes one explicit state boundary between background providers and
particle transport. It does not change the governing transport equations or
select a different mover.

## Authoritative clock

`PIC::SimulationTime` is the only owner of application time. srcSEP reads it
through `SEP::Background::SimulationTimeSeconds()`, which returns seconds from
the PIC simulation-time origin and exposes no setter or independent advance
operation.

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
| provider | `analytic`, `swcme`, `swmf`, or `local-evolution` authority |
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

The snapshot does not deep-copy the large AMPS field-line arrays for every
particle. Instead, `SnapshotStore` freezes their provider/generation identity.
Publication is prohibited during the global particle read phase, and every
entry through the common `SEP::ParticleMover` wrapper acquires a const snapshot
reference that remains valid for the full mover call. The enclosing read-phase
guard owns the snapshot, and an atomic phase marker avoids a mutex/reference-count
operation in every particle call. This gives all scheduler threads one read-only
realization without a per-particle field-array copy.

`EvaluateLocalBackgroundAt` derives density-based `div(U)` from the previous
and current physical epochs. A static zero-length epoch pair is valid only when
the two density values agree. Every mover composes the same snapshot-validity
limit and verifies the field-line generation before each substep shell; reaching
the closed validity endpoint requires the driver to end the read phase before
continuing with a new snapshot.

## Provider ownership and handoff

The store applies these rules before accepting a publication:

- analytic and SWCME state must be `ModelOwned`;
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
