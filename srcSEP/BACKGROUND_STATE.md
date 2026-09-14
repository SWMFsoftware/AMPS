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

SWCME's legacy `SetModelAndState` entry point also calls
`AssertProviderMayWrite` before changing its backing cache. This prevents code
that bypasses the normal publisher from mutating SWCME state during a particle
phase or while a different provider owns the background.

The handoff API records ownership but does not itself copy AMPS arrays. A future
caller enabling local evolution is responsible for making that one-time private
copy before publishing the handoff. No existing Step 2 path requests local
evolution.

## Update order

For a standalone SWCME step:

1. Read the upcoming epoch from `PIC::SimulationTime`.
2. Prepare `swcme1d::StepState` at that exact epoch.
3. Publish model-owned metadata valid through the upcoming global time step.
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
```

The test compiles the production snapshot implementation using C++11,
`-Wall -Wextra -Werror`, and pthread support. It checks construction domains,
compile-time non-assignability, interval enforcement, provider isolation,
explicit SWMF handoff, generation/epoch monotonicity, publication exclusion,
concurrent mover views, deterministic fingerprints, and the single-clock source
invariant.

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
