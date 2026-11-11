# Periodic associated-data transport

`ordered_periodic_exchange.h` is the transport helper used by
`Periodic::UpdateGhostBlockData` in `../pic_bc_periodic.cpp`. It depends only
on MPI and the C++ standard library. It handles opaque bytes, not particles,
physical units, field normalization or a new boundary condition.

`Workspace::Exchange` accepts a replicated ordered list of source/destination
ranks, record byte length, rank/size, the **supplied communicator**, pack/unpack
callbacks, and a snapshot-independence flag. Every participant must use the
same ordered table, tag, length and mode. No global collective is introduced;
a rank without participating pairs need not enter a communicator-wide barrier.

For independent pure copies it posts per-peer receives, packs sources in
pair order, sends aggregated buffers, waits for completion and applies ghost
writes in the original pair order. Messages larger than the MPI C count limit
are fragmented in same-peer/tag order. Retained memory is byte capacity only:
pair offsets, ranks and lengths are reconstructed every call, avoiding stale
plans after repartition, communicator or associated-data layout changes.

`snapshot_safe=false` preserves the streaming reference algorithm. A write
may then become a later source, or a callback may have side effects. No
pack-ahead or reordering is allowed in that mode. Invalid ranks, chunk sizes,
overflow or MPI errors reject the operation; they never become an empty PASS.
Zero-length/empty plans transfer no bytes. Workspaces require serialized use;
the PIC adapter calls this outside its particle/field worker regions.

The PIC adapter snapshots only default copy callbacks and distinct periodic
ghost/real leaves. Low-index corner serialization excludes upper block
planes, giving disjoint source/target half-open support; center nodes are
interior. Ghost-as-source, uncertain pairing and custom callbacks use the
reference path. The existing full associated-data payload and corner/center
order are unchanged, including the historical excluded upper corner planes.
Pack/unpack callback-selected smaller payloads are not implemented by this
optimization and must not be inferred from the function signature.

Optional diagnostic controls, supplied consistently to all rank executables:

- `AMPS_PERIODIC_EXCHANGE_REFERENCE=1`: force the streaming reference path.
- `AMPS_PERIODIC_EXCHANGE_DIAGNOSTICS=1`: print the first call's actual record
  bytes, traffic, message counts and local-pair count on each rank.

These are A/B instrumentation controls; they do not alter physics or gates.
The first four-rank wave call transfers 1,335,296 bytes per record and sends
about 16--17 MB per rank, aggregated into three peer messages per rank.

Independent tests: `CORE-PIC-PERFORMANCE-001`, seven transport identities at
1/2/4 ranks covering local/remote bytes, changing layouts/order, fragmentation,
alias-sensitive fallback, empty/invalid plans and split/reordered communicators.
The native four-rank fast-wave regression retains both numeric checks and all
six exact sampled-file identities. See
[`src/test/tests/pic/performance/README.md`](../../test/tests/pic/performance/README.md)
and [`test/srcFastWave/PERFORMANCE_OPTIMIZATION_PLAN.md`](../../../test/srcFastWave/PERFORMANCE_OPTIMIZATION_PLAN.md).
