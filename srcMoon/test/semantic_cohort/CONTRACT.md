# Supporting normal-startup and semantic-particle composition contract

This extends the existing U15/U16 chemistry execution support; it creates no
new permanent test ID. It is a software and mechanical fixture, not sodium
kinetic qualification, an observation or a replacement for I21–I24.

The production entries are `NativeRegistration::SetBeforeParticleBuffer` and
the actual `amps_init()` pre-buffer boundary, `SemanticRandom`'s opaque PIC
storage, `ParticleTransaction::{Reserve,Stage,Commit,Rollback}`, and the normal
Moon/controller/core first- and second-order movers. The thin test reuses the
existing distributed PIC controller fixture's mesh and callbacks. It supplies
one explicitly analytical, nominal-heavy-mass zero-Q charge-transfer topology.
No measured channel, cross section, recoil PDF or physical state is inferred.

Two immutable compile recipes are distinguished. Canonical Moon and the
existing eight-species fixture keep legacy constants. A separate analytical
eight-species recipe uses the exact SI elementary charge 1.602176634e-19 C,
converted to the normal species parser's legacy charge-token convention and
verified through the actual emitted C++ expression. Repository masses are
unchanged and explicitly nominal. This is not a physical Na state profile.
Actual compiled charges must agree with the explicitly supplied controller,
transaction and optional collector charge unit; a legacy-unit mismatch must
refuse before publication. Shared core constants are never changed.

Before allocation, the normal callback must run exactly once with no particle
buffer and no physical records, register semantic storage, and return. A
second invocation and late registration must refuse. The actual initializer,
SPICE furnishing, material table, mesh, weights and buffer remain production
owned. Runtime mover order is explicit; MPI layouts are one and eight ranks.

The independent random oracle is the existing HMAC-SHA256 counter definition:
seed, cohort, original global parent ID, complete genealogy, purpose and draw
index. Rank, pointer and allocation slot are excluded. Ownership clones and
MPI preserve every opaque byte. A daughter appends its canonical product
ordinal and records the already supplied product-clock draw as counter one.
Source, channel and kinematics counters remain distinct. Full genealogy
capacity is explicit resource admission, never a physical generation cutoff.

The transaction controls reserve actual private PIC slots, initialize complete
daughter records, deliberately fail after partial staging, and require all
slots, the parent bytes/list and the accepted ledger to remain unchanged after
rollback. Retrying the same immutable plan uses the same clocks and keys,
without drawing again. Capacity exhaustion must refuse before publication.
The accepted packet preserves represented number and the parent's original
physical interval despite unequal local species time steps and weights.
Normal mover continuation must finish all daughters and leave no pending
packet or double-deleted parent.

Artifacts retain command/stdout/stderr, compile input, generated headers,
library/executable/source identities, per-layout named check CSVs, semantic
key/clock evidence, result and run manifests. Comparison plots illustrate
actual versus independently expected counters/clock values; they are labelled
analytical software controls. All names and counts are checked explicitly;
exit zero or an empty summary alone cannot pass. Infrastructure failures are
ERROR; evaluated invariant failures are FAIL; missing compiler/MPI or an
unavailable authoritative physical profile is SKIPPED. No observational
scoring is performed. Every isolated package has one original wall deadline
of at most 3600 seconds including setup, builds, retries, cleanup and reports.
