# M3 conservative PIC particle backend: native verification contract

Written before fixture implementation, 2026-11-17. Supporting extension of
existing chemistry bookkeeping/production wiring coverage; **no new permanent
U/I ID, no physical charge-exchange or recoil qualification**. Formal I21–I24
and predecessor gates keep their existing status.

| Required field | Contract |
|---|---|
| ID/name | Supporting M3 native conservative particle transaction |
| Purpose | Verify that an immutable, already validated reaction plan is published through the real PIC particle allocator and actual mesh-cell lists atomically, or leaves no physical changes. |
| Production capability/entries | `Moon::Chemistry::RegisterReactionStorage`, `ReadReactionState`, `ParticleTransaction::{Reserve,Stage,Commit,Rollback}`; `ReactionExecution::{MakeEventPlan,PendingEvent::TryCommit}`; actual `PIC::ParticleBuffer` and mesh-block weight/time APIs. |
| Configuration | Clean isolated normal Moon Config.pl/make build; eight compiled atomic species Na,NA_PLUS,He,HE_PLUS,Ne,NE_PLUS,Ar,AR_PLUS; individual weights ON; species-local weight/time-step modes; MPI-only, OMP=1, tracker OFF and stationary-flux resampling OFF. The unrelated canonical velocity-distribution sampler is OFF because its physical Moon sample coordinates are outside this analytical unit-box mesh; ordinary native Init_AfterParser remains used. Particle reactions remain globally disabled. |
| Inputs | Explicitly analytical stationary zero-Q nominal-heavy-mass topology: tracked Na + bath He+ → tracked Na+ + tracked He. No cross section, probability, solar epoch, photon/electron rest-mass assumption or measured observation is supplied. Catalog uses actual compiled species masses/charges. The fixture refuses unexpected neutral/ion nominal-mass inequality instead of fitting a bath residual. |
| Independent oracle | Exact integer nuclei/charge; explicit incoming He+ reservoir contribution; represented parent number1.5; each product number1.5; actual cell list connectivity/allocation counts and immutable byte snapshots; exact declared local interval and original fresh optical clocks. Production computes the plan and executes the transaction. |
| Units | SI kg, m, m/s, J; weights and statistical corrections are particle-number quantities; charge compared with actual generated Coulomb values and integer elementary-charge registry. |
| Frame/time | Declared analytical Cartesian fixture. Local physical interval start1000s/duration10s/elapsed3s; daughters retain remaining7s. These numbers are numerical controls and are not UTC/TDB lunar epochs. |
| Deterministic/stochastic controls | Fixed explicit open-interval clock draws; no RNG in reserve/stage/commit/retry. Compare unchanged immutable plan/daughter clock values over capacity deferral and retry. No stochastic branch inference is made. |
| Procedure | Initialize normal PIC MPI/mesh topology, split at least once, distribute/allocate real blocks/cells, register state before small real particle buffer allocation. On every rank choose an actual owned block/cell. Run the cases below at1 and8 MPI ranks with separate retained artifacts. |
| Metrics/acceptance | Every declared boolean/integer identity exact. Byte invariants exact. Weight round trips compared with an explicit rounding-level bound1e-14 relative, not a physics tolerance. No product-dt factor: parent dt.5s, Na+dt.125s, Hedt.0625s while each represented product remains1.5. |
| Expected PASS | All named controls evaluated and passed on every rank; no allocated-slot leak; cell linked-list integrity and physical ledger/counts agree; actual compiled eight-species mass/charge inventory retained. |
| Expected failure modes | Capacity gives Deferred without changing parent/list/physical ledger; stale species/reference/position/opaque parent bytes or catalog identities refused; invalid retained-parent placement refused; rate-only plan cannot commit. Compiler/link/launcher/parser/artifact corruption is ERROR, completed invariant miss FAIL. |
| Required data/hashes | No external validation dataset/observations. SPICE toolkit may be a build dependency only. Retain source revision/source hashes, exact derived compile input and generated macro snapshots, library/executable hashes, compiler/MPI details, commands and run layouts. |
| Predecessors | Local reaction execution analytical controls available; successful actual native build. Formal M1/M2 and I21–I24 closure are not inferred from this supporting run. |
| Artifacts | stdout/stderr, result.json/run_manifest.json, per-rank named-check CSVs and actual compiled species inventory, build/config/probe logs, input/generated snapshots, process-budget records, binary/library/source identities. No prior meaningful output overwritten. |
| Status semantics | PASS local/native adapter controls only; missing prerequisites SKIPPED; completed invariant miss FAIL; infrastructure/runtime/partial-output failure ERROR. Never label this analytical fixture as observations or physical CEX/ionization validation. |

## Required native cases

1. Real node/block ownership and actual initialized cell list, generated
   species mass/charge values, individual-weight configuration and unequal
   source/product `GetLocalTimeStep` values.
2. Complete-packet capacity deferral with only one free slot for two products;
   parent bytes, cell head, allocated count and accepted-event ledger remain
   unchanged. Free real filler slots and retry the **same** pending plan;
   original clock decisions and physical interval remain unchanged.
3. Successful rebuilt-list transaction deletes parent once, publishes both
   products exactly once, preserves sentinel-list connectivity and increments
   accepted-event ledger once. Destructor/rollback after successful commit
   cannot remove committed products.
4. Explicit rollback after actual reserve/stage releases all private product
   slots and preserves parent bytes/list/ledger. Staged daughters are never
   in the public list before commit.
5. Parent opaque-byte mutation between stage/commit is refused. Restore the
   fixture's deliberately changed byte and retry without changing event
   decisions. Also refuse changed parent coordinates/reference weight or
   catalog state identity before reservation.
6. Product species/state index, actual compiled mass/charge, reference-weight
   correction, fresh clock, internal energy and original physical interval
   are read from actual committed native records. No daughter local-time-step
   birth factor is applied.
7. Detached continuation with `retainParent=true` leaves the parent allocated
   and byte-identical for mover-owned deletion; both products are allocated
   and unlinked. Delete the old parent exactly once through the normal PIC
   deletion API and subsequently free the two detached products. This case
   verifies ownership boundaries, not transport across MPI faces.
8. Rate-only/invalid retained-parent placement and a repeated commit are
   refused without duplicating physical publication or ledger entries.
9. Two actual transactions reserve/stage against a shared ledger. The second
   publishes first; then a one-shot C++ allocation fault forces the first's
   ledger-capacity growth to fail. This must happen **before** parent/list/
   physical-ledger publication. Rollback and identical-plan retry must work.
   The fault affects only the targeted C++ allocation, not PIC allocation
   behavior or any physical process. Missing actual fault coverage is ERROR.
10. A separate explicitly analytical one-dimensional elastic, zero-Q
    nominal-heavy-mass fixture uses parentNa velocity+1m/s and incoming bathHe+
    velocity−1m/s. Standard two-body energy/momentum conservation determines
    distinct daughter velocities. The incoming bath moments are specified
    from its declared velocity, never fitted to an observed residual. Native
    stored daughter position/velocity must match the supplied plan. This
    tests PIC field publication; it does not assert a physical Na/He CEX law,
    angular distribution, channel threshold or collision rate.

The fixture will not mock allocator success, invent daughter recoil, mutate
reference observations, or use source/preflight results as linked native PASS.
Actual mover continuation/MPI migration and physical reaction distributions
remain separate required gates.
