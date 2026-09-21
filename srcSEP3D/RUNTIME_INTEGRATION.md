# R01–R07 Production Runtime Integration

This document is the implementation contract for the seven runtime
improvements that connect the portable srcSEP3D physics to a configured AMPS
application.  The word *tick* always means the integer index owned by
`Runtime`; physical time is derived as `tick * run.time_step_s`.  Coupled data
are immutable while particles move, and all publication, source, observer, and
checkpoint operations occur only after AMPS has joined at a step boundary.

## R01: one AMPS mover hook and one immutable context

`amps/install_mover_hook.py` edits the configured
`AMPS/build/pic/picGlobal.dfn` after AMPS configuration and before
`pic_mover.cpp` is compiled.  It inserts the exact declaration of
`SEP3D::AMPS::Movers::MoveParticle` and maps
`_PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_` to that function.  The script
is idempotent and refuses to replace an unrelated existing mover mapping.
`make prepare-production` invokes it, while `strict-production` makes the
preparation an explicit prerequisite and audits the generated mapping.

`amps_init()` installs one copied `AMPS::Movers::Context`.  The context contains
the local-state resolver, the hard substep limit, the current immutable shock
record, and a host-owned ledger pointer.  Replacement is forbidden while
`Runtime` is `Running`.  The resolver samples only the currently published
background/turbulence generation; it never reads a mutable SWMF work array.
AMPS return integers are isolated in `amps/amps_mover_status.h` and checked
against the configured `pic.h` ABI.

Before each AMPS particle phase, every rank opens one rank-local ledger row per
species.  After AMPS has completed particle exchange, the application sums
start, disposition, crossing, and end counters collectively.  Only a globally
closed row is retained.  Thus rank migration is neither a source nor a sink,
and a particle deleted after being classified active is detected immediately.

## R02: consume the complete requested time

One AMPS mover call requests a physical interval `dt_requested`; it does not
request one numerical substep.  `AdvanceParticleRequestedTime` repeats:

1. resolve the cell, background, turbulence, and coefficients at the current
   particle position;
2. compute all named stability limits and select an accepted substep no larger
   than the remaining time;
3. advance the selected Parker or focused core with a semantic random key;
4. evaluate inner/outer boundaries and the first moving-shock intersection;
5. accumulate consumed time and continue from the accepted state.

The loop terminates only on exact requested-time completion (apart from a
machine-roundoff remainder), a semantic terminal disposition, a resolver/core
error, or the configured hard substep count.  Local state is re-resolved before
every accepted substep, and an expanding shock radius is advanced by the
already consumed time.  Active particles store `completedStep = startTick + 1`
so the next global step cannot reuse a keyed random stream.

## R03: transactional background and turbulence generations

Snapshot updates follow an explicit side state:

`Idle -> Requested -> Filling -> Staged -> Idle`

A provider failure moves to `Failed`; acknowledging the failure returns to
`Idle` without modifying the active generation.  Background and turbulence
are prepared and evaluated into candidate storage.  An `MPI_Allreduce` proves
that every rank accepted the complete candidate before `Runtime` publishes its
descriptor and the application swaps the immutable shared pointers.  Movers
therefore see the old complete pair or the new complete pair, never a mixture.

The associated AMPS cell bytes are a cache and diagnostic copy.  The immutable
snapshot pointer is the mover authority, so a failed collective preparation
cannot expose partially filled cache bytes to particle workers.  Generation,
provider identity, frame, configuration fingerprint, epoch, and validity
interval are part of the publication and restart provenance.

## R04: authoritative tick clock and integer event schedule

`RuntimeCounters.currentTick` is authoritative. `CurrentTimeS()` and
`NextStepEndTimeS()` are derived from the frozen base time step; repeated
floating additions are not used. Species identity is likewise single-source:
AMPS `SpeciesList` generates the complete count, index order, chemical table,
mass table, and signed-charge table. Between `PIC::Init_BeforeParser()` and
mesh/after-parser initialization, `BindCompiledSpeciesTable()` enumerates every
index, copies those immutable values, and calls
`ValidateCompiledSpeciesBinding`. The validator rejects count disagreement,
non-contiguous indices, duplicate/empty symbols, non-positive mass, neutral or
non-finite charge, and observer indices outside the generated table. It never
writes AMPS molecular data and does not depend on a chemical-species macro.

`amps_init()` applies the common explicit timestep and base statistical weight
to every generated global slot and every owner-local block. Each due source
event then allocates the exact requested count independently for every compiled
species. Its kinetic-energy interval is converted to momentum with that
species' AMPS mass before sampling, and the species index remains part of the
random key and source/particle ledgers. Before every particle phase,
`VerifyClockAgreement` checks Runtime time, PIC time, PIC time step, snapshot
time, and shock time.

Background refresh, injection, sampling, and checkpointing each have an
integer cadence and a persisted next-event tick.  `CompleteStep` advances the
tick once and advances due events by their cadence.  Restart validates both the
counters and event schedule before mesh mutation.  Startup prints the restored
tick, physical time, base step, and all next-event ticks.

## R05: shock providers, physical source rates, and source ledgers

`ShockProvider` separates analytic or host-published shock reconstruction from
particle allocation.  It returns one immutable `ShockState` with generation,
coverage interval, geometry, compression, provenance, and zero or more source
patches.  Source creation occurs at the configured injection tick after
transport and snapshot publication.

For each connected active patch, the expected macroparticle count is

`physical_rate_per_s * cadence_interval_s / macroparticle_weight`.

The fractional remainder is stochastically rounded with a semantic key.  A
configured cap is explicit; accepted particle weights are renormalized so the
represented physical number remains exact.  Momentum, pitch, gyrophase, and
identity use separate streams.  The injection opportunity combines physical
shock generation with the authoritative tick, preventing stable-ID reuse when
one shock generation remains valid across multiple source cadences.  The
physical shock generation remains separate for crossing de-duplication.

AMPS allocation uses `InitiateParticle`, then writes the packed srcSEP3D state.
Only the rank owning the patch position records the physical source row.
Ledgers retain represented number, energy, momentum, macro count, rejection,
cap, inactive, and disconnected counters; a failed allocation aborts the run
instead of silently changing source normalization.

## R06: observers and commit-only sampling windows

Observers are immutable configuration objects.  Supported geometry includes
fixed Cartesian, fixed heliographic, moving Cartesian, spherical shell, and a
host-updated field-connected location.  Each observer declares collection
radius, accepted species and pitch range, logarithmic energy bounds/bin count,
cadence, products, and normalization.

At a due joined boundary, every rank exports read-only particle and cell
records.  Root gathers complete records, and the sampler sorts particles by
stable ID before compensated reduction.  Spacecraft spectra include the sum
of statistical weights, sum of squared weights, macro count, and standard
uncertainty.  Sampling is ordered after motion and source injection, so the
timestamp names an unambiguous post-step state.

`ObserverRuntime` separates `Capture`, `PreparePublication`, and
`CommitPublication`.  A failed file transaction leaves the pending window
intact; counters and accumulators reset only after the staging directory and
manifest are atomically published.

## R07: complete, versioned restart and deterministic continuation

Restart schema 2 uses magic `SEP3DR02`, an explicitly serialized
little-endian payload, payload length, and checksum.  It stores:

- physics fingerprint, resolved manifest, code identity, storage-layout
  fingerprint, and saved MPI rank count;
- integer clock/counters, next-event ticks, base time step, active snapshot,
  and background/turbulence/source generations;
- shock state, campaign seed, next stable ID, every active particle's complete
  stochastic tuple, global closed particle rows, and physical source rows;
- observer sampling counters including a pending-window summary.

Loading builds and validates a candidate before modifying Runtime or AMPS.
Identity, layout, checksum, record limits, finite state, unique stable IDs,
ledger closure, snapshot availability, and rank policy must all pass.  Runtime
counters and events are restored before mesh creation.  Particles are then
reinserted on their spatial owner, immutable providers are re-based to the
saved generation, and PIC simulation time is set from the integer clock.

Checkpointing begins only at a joined due tick.  Particle/source tables are
gathered, root writes a sibling staging file and renames it atomically, and all
ranks either complete or abort the Runtime checkpoint transition together.
The sequence counter advances only after successful publication.  With the
deterministic repartition policy, a different rank count may restore by spatial
ownership without changing particle random histories.

## Acceptance evidence

`R3D01` through `R3D07` are registered native C++ tests and are included by
`test/run_tests.py --all` and `--suite improvements-r`.  They cover the build
hook, complete subcycling, transactional publication, integer clocks/events,
source conservation and cadence identity, observer commit behavior, and the
full restart round trip.  `BLDL3D01`, `BLDL3D03`, and `BLDL3D05` additionally
exercise the configured AMPS boundary and copied `build/main` layout.
