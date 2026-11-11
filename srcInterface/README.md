# AMPS component and linear-system interfaces

This tree contains the SWMF component wrappers/bridge and reusable linear-system
interfaces. `LinearSystemCornerNode.h` assembles corner-node systems and uses
the existing Fortran GMRES wrapper. The full framework coupling scope remains
documented by the manual and core interface tests; this change qualifies the
bounded CPU sparse-multiply optimization described here.

## Solve-scoped compact multiplication

For the standard host scalar kernel, `Solve` activates an instance-owned
`amps_ordered_sparse::View` for the fixed operator passed to GMRES. The first
matrix-vector callback binds fresh receiver storage and snapshots coefficient
values/RHS pointers in original row order. Subsequent callbacks read current
unknown-vector values through those pointers. The existing local input-vector
copy is retained, preserving input/output alias protection.

A new solve always refreshes the snapshot. Matrix-generation/vector-length
changes and explicit `UpdateMatrixNonZeroCoefficients` calls invalidate it.
The scope guard deactivates it immediately after the solver wrapper returns,
including C++ exception unwinding. Thus direct `MultiplyVector` calls outside
`Solve` continue to read public coefficient entries live. User mutations made
between solves are observed; a linear solve otherwise requires its declared
fixed operator. Do not reuse a compact snapshot after storage/topology changes.

The initial callback rebinds the pre-existing static communication buffers,
preventing a solve from using another matrix instance's receiver addresses.
This does not claim general reentrant/concurrent solver support. Explicit
AVX, CUDA and the separate std::thread multiplication variants keep their
existing kernels. The ordinary OpenMP row schedule reads the same immutable
compact data. No new reduction ordering, fast math, preconditioner, physical
operator, solver tolerance or iteration ceiling is introduced.

Independent sparse-kernel tests are registered under
[`CORE-PIC-PERFORMANCE-001`](../src/test/tests/pic/performance/README.md).
The full native wave comparison checks the actual Fortran callback path and
all six exact sampled-file identities. Existing ECSIM manufactured operators,
boundary cases and native MPI/core tests remain required regression gates.
New backend-specific speedups are not inferred from these scalar results.
