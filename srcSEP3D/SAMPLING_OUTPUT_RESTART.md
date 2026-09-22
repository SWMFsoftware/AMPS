# Phase O: Sampling, Output, and Restart

Phase O produces diagnostics without modifying particle trajectories and makes
restart state sufficient to reproduce the next stochastic draw and output
sequence. The implementation is divided into a pure sampler, a transactional
publisher, a canonical restart codec, and a Runtime coordinator.

## Read-only deterministic sampling

`SamplingRequest` contains immutable particle observations, cell definitions,
virtual-spacecraft volumes, field-line projection bins, closed ledger rows, and
the previous sampling counters. `Sample()` copies and sorts observations by
stable particle ID before any reduction. Duplicate or zero IDs are errors.
Consequently AMPS linked-list order, OpenMP scheduling, and MPI migration do
not define the arithmetic order.

Floating sums use Neumaier compensation. The fixed stable-ID order plus a
specified compensated algorithm makes repeated sampling bitwise repeatable for
one global observation set and reduces loss when statistical weights span many
orders of magnitude. Production MPI collection must form that global set in
stable-ID order before calling `Sample`; rank-local partial sums are not a
substitute because their tree order would change roundoff.

### Cell moments

For cell volume `V` and particle weights `w_i`, the products are

\[
n=\frac{\sum_iw_i}{V},\qquad
\mathbf F=\frac{\sum_i w_i\mu_i v_i\hat{\mathbf r}_i}{V},\qquad
E=\frac{\sum_iw_iK_i}{V},\qquad
\langle\mu\rangle=\frac{\sum_iw_i\mu_i}{\sum_iw_i}.
\]

The energy and speed use the exact relativistic momentum relations. The radial
flux direction is an explicitly named diagnostic approximation because the
minimal observation does not duplicate the background `b` vector; a future
field-aligned flux product should add `b` to the observation schema rather than
silently changing this column.

### Virtual spacecraft and field-line projections

A virtual spacecraft selects particles inside a declared spherical collection
radius, bins exact kinetic energy, reports represented particles per joule,
and reports the dipole anisotropy `3<mu>`. A field-line projection computes
`s=(x-origin)·direction`, bins in physical metres, and reports represented
particles per metre. Every bin edge is strict, finite, and increasing; the
rightmost edge belongs to the final bin.

R06 resolves each configured observer at the authoritative post-step time.
Fixed Cartesian/heliographic, moving Cartesian, spherical-shell, and
host-updated field-connected locations share the same immutable sampling
definition. Species and pitch-angle acceptance are applied before logarithmic
energy binning. Each spectrum retains represented weight, squared-weight sum,
accepted macro count, declared normalization, and the standard uncertainty
derived from the squared weights.

Production gathers complete particle/cell records to root and lets `Sample()`
perform the stable-ID ordering; it does not reduce rank-local floating sums.
`ObserverRuntime` uses capture/prepare/commit phases. A failed staging write
leaves the pending observation window intact and only a successful atomic
publication clears its counters.

Shock diagnostics are copied only from closed ledger rows. Thus injected,
escaped, absorbed, failed, and crossing counts have the same conservation
authority as the mover rather than being inferred later from floating output.

## Transactional output publication

`Publish()` writes four CSV artifacts into a new `.staging` directory:

- `cells.csv`;
- `spacecraft.csv`;
- `field_lines.csv`;
- `shocks.csv`.

Column names carry SI units (`_m`, `_J`, `_m-3`, `_m-2_s-1`, and so on). After
each stream is closed, the writer computes a 64-bit FNV-1a artifact digest. It
then writes `manifest.txt` with schema version, output sequence, simulation
time, configuration fingerprint, code identity, snapshot generation and
fingerprint, sampling counters, and artifact hashes. One directory rename is
the publication commit. An existing sequence is never silently reused.

FNV-1a is an integrity and reproducibility checksum, not a cryptographic
signature. `ParseAndVerifyPublication()` is independent of the writer path and
checks required manifest keys, exact unit-bearing headers, artifact count, and
every digest before replacing caller output.

`PublishIfDue()` derives cadence solely from `RuntimeCounters`: a completed
step sets `stepsSinceOutput=0` exactly when it increments `outputSequence`.
There is no second static modulo counter that can drift after restart.

## Canonical restart image

The restart file is not a memory dump. It contains:

1. eight-byte magic `SEP3DR02`;
2. little-endian payload length;
3. versioned payload;
4. FNV-1a checksum of the exact payload bytes.

Schema 2 serializes:

- physics fingerprint, resolved configuration manifest, code identity,
  storage-layout fingerprint, snapshot fingerprint, and saved MPI rank count;
- authoritative current tick, completed steps, output/checkpoint counters, the
  base time step, and every persisted next-event tick;
- complete active snapshot descriptor plus background, turbulence, and source
  generations and the complete shock state;
- campaign seed and next stable particle ID;
- committed and pending sampling counters;
- every active particle's stable ID, species, Cartesian position, momentum,
  pitch cosine, gyrophase, weight, completed step, substep, and last shock
  generation;
- every globally closed particle-ledger row and physical source-ledger row.

Integers and IEEE-754 bit patterns are written field by field in little-endian
order. Structure padding, native enum width, and host ABI therefore cannot
alter the file. Particles and ledger rows are sorted by physical identity
before writing. The keyed random generator needs no opaque engine dump: the
campaign/particle/step/substep/purpose tuple reconstructs its future exactly.

Writing uses a sibling `.staging` file and a final rename. Reading builds a
candidate and validates magic, length, schema, checksum, fingerprints, record
limits, finite particle state, unique sorted IDs, and exact ledger closure.
Caller output is changed only after all checks pass. Optional load expectations
also verify the resolved manifest and storage layout. Rank-count mismatch is
either rejected or admitted only under the explicit deterministic repartition
policy; particles are then reinserted on their spatial AMPS owner without
changing their stochastic tuple.

## Snapshot policy and lifecycle

A restart may name a background generation not yet published by a coupled
host. `MissingSnapshotPolicy::Reject` fails immediately.
`MissingSnapshotPolicy::Wait` requires an explicit read-only availability
callback and bounded timeout; it never substitutes the newest or closest
generation.

`WriteRestartAtBoundary()` enters Runtime's `Checkpointing` state, writes the
post-commit checkpoint sequence, then commits the matching Runtime increment.
A failed write calls `AbortCheckpoint()` and returns to `SnapshotReady` without
changing the sequence. `RestoreRestartBeforeMesh()` is legal only from
`Configured`, so counters are restored before mesh allocation and background
publication.

The production coordinator gathers complete particles and rank-owned source
rows at a joined due tick. Root writes a sibling staging file and renames it;
the success flag is broadcast before every rank completes or aborts the
checkpoint transition. Restore validates the whole candidate before mesh
mutation, restores clock/events first, re-bases provider generations, sets PIC
simulation time, then inserts particles. This ordering makes a split run use
the same next `(campaign,particle,step,substep,purpose)` random key as the
uninterrupted run.

## Evidence

- `NAT3D06`: repeated/reversed input produces identical products and caller
  observations are untouched.
- `NAT3D07`: manifest, SI schemas, atomic bundle, independent parsing, and
  corruption detection.
- `RST3D01`: complete round trip and identical future keyed normal draws.
- `RST3D02`: fingerprint/checksum failures are transactional.
- `RST3D03`: explicit missing-snapshot reject/wait behavior.
- `R3D06`: resolved observers, uncertainty, and commit-only window reset.
- `R3D07`: complete schema-2 clock/event/provider/shock/RNG/ledger/sampling
  round trip.
- `R3D09`: finite native Tecplot records for invalid background padding,
  not-yet-sampled cells, completed empty particle cells, and occupied cells.

## Native Tecplot empty-data semantics

AMPS particle moments and the srcSEP3D background are sampled independently.
A completed particle sampling window containing zero particles is therefore a
valid measurement: particle number, density, velocity moments, energy, and
temperature are written as finite zeros. `particle_sample_present=0` records
that condition, while `particle_sampling_window_valid=0` means that the file
was written before any sampling interval completed (as in initialization-only
mode).

The Cartesian AMR tree also contains padding cells outside the configured
heliocentric shell. Background physics is undefined there. Those columns use
finite zero placeholders with `background_valid=0`; post-processing must mask
them with that flag. This avoids `NaN` propagation in Tecplot without assigning
physical meaning to the placeholder values.
