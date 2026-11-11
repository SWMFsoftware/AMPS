# General utilities: performance-preserving kernels

This directory contains shared containers, numerical tools, random streams,
communication and input support. The generic input reader has its own
[README](input/README.md). This section documents the numerical compatibility
contracts of the new CRC and sparse-multiply optimizations.

## Historical CRC compatibility

`CRC32` in `specfunc.h` retains its public API, object size/base lookup-table
storage, zero initialization, append/split-update semantics and returned
`unsigned long`. Its historical recurrence is

`R(s,b) = (s << 8) XOR T[((s >> 24) XOR b) & 255]`.

The table uses polynomial `0x04c11db7`, but updates are **not explicitly
masked to 32 bits**. LP64 upper bits are part of existing signatures. Do not
replace this implementation with a standard CRC-32, CRC-32C or hardware
instruction without proving full returned-width equivalence.

Host code shares immutable tables and expands four byte updates into four
table lookups. The surviving `s << 32` term is retained on wider unsigned
long platforms and is zero on 32-bit platforms. Byte loads avoid alignment
and endian assumptions. The constructor still initializes every object's
original base table, preserving copy semantics and host-to-device object
validity. The device branch retains its original table generation and bytewise
update; this work does not claim new CUDA qualification.

Tests under `src/test/tests/general/numerics` compare an independent bitwise
recurrence on all-byte/unaligned/long streams, split updates, scalar templates,
clear/empty updates, object copies and concurrent host construction. The native
`PIC-PBUFFER-SIGNATURE-001` checks canonical/list-inclusive particle records and
incremental accumulation without changing input bytes. Repeatability checks
and signatures remain enabled.

## Ordered sparse multiplication view

`ordered_sparse_matvec.h` provides `amps_ordered_sparse::View`. `Build(rows,
count,element)` snapshots coefficients and RHS pointers in insertion order.
`Row(i)` performs the same scalar multiply-and-add sequence; pointed-to RHS
values remain live. Empty rows return zero, and invalid row lengths, missing
RHS pointers or out-of-range row requests reject the operation.

The caller owns RHS storage and must rebuild after coefficients, topology or
storage change. The production adapter restricts the view to one fixed-operator
linear solve; direct public matrix-vector calls retain their live-data path.
No MPI communicator, preconditioner, iteration ceiling, convergence tolerance
or floating-point reduction policy is changed. See
[`srcInterface/README.md`](../../srcInterface/README.md) for lifetime details.

`PERF-SPMV-ORDER-001`, `PERF-SPMV-REFRESH-001` and
`PERF-SPMV-EMPTY-INVALID-001` use independent index-based references, including
cancellation-sensitive insertion order, live RHS updates, coefficient/storage
replacement and invalid inputs. They belong to `CORE-PIC-PERFORMANCE-001`.
