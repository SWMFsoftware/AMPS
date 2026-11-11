# Complete pointer-stencil fast path

`complete_stencil.h` supplies the dependency-light host/device template
`amps_complete_stencil::Assign`. It takes eight ordered node pointers, local
IDs and weights, the periodic/domain bounds policy and an output stencil.
Positions use the caller's coordinate units; weights are dimensionless.

The helper first validates all support nodes, nonnegative weights and the
existing nonperiodic global-box test. Only then does it write eight entries.
It preserves their exact supplied order, IDs and floating-point weights,
including zero weights. It does not normalize or regroup arithmetic.
Rejection returns false without a partially modified destination. Storage
capacity of at least eight entries is a compile-time requirement.

The production users in `../pic_interpolation_routines.cpp` restrict this
path to a fully allocated **uniform** pointer stencil. Cell-center queries
with row output, adaptive geometry, missing nodes or domain clipping retain
the general implementation. Corner queries retain the historical coefficient
permutation and normalization/insertion order. Position weights are computed
for every query; no old particle coordinate or field value is reused.

Deposition additionally builds per-call block eligibility, cell geometry and
uniform-grid magnetic-node address views in `../pic_field_solver_ecsim.cpp`.
The division/product order and length-conversion convention are unchanged.
These tables are read-only during the existing worker schedule and rebuilt
before every call, so ownership/allocation/layout changes cannot reuse stale
addresses. Field values are dereferenced live. Adaptive, standalone and device
callers keep the original field-neighborhood construction through null optional
`cProcessCellData` views. GPU performance qualification is not claimed.

`PERF-STENCIL-EXACT-001`, `PERF-STENCIL-FALLBACK-001` and
`PERF-STENCIL-REFRESH-001` independently check entry order/zero weights,
no-write rejection and changed-query refresh. Existing native mesh tests
remain the physical constant/linear/AMR-interface interpolation oracle; the
four-rank fast-wave checks preserve all six sampled files byte-for-byte.
See [the performance owner](../../test/tests/pic/performance/README.md) and
[the implementation ledger](../../../test/srcFastWave/PERFORMANCE_OPTIMIZATION_PLAN.md).
