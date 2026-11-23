# AMPS component and linear-system interfaces

This tree contains the SWMF component wrappers/bridge and reusable linear-system
interfaces. `LinearSystemCornerNode.h` assembles corner-node systems and uses
the existing Fortran GMRES wrapper. The full framework coupling scope remains
documented by the manual and core interface tests; this change qualifies the
bounded CPU sparse-multiply optimization described here.

## Coupled SEP fixed-step clock

The SWMF/single-global/counter-OFF SEP application opts into
`AMPS2SWMF::CoupledStep` after timestep initialization. Its nominal duration
[s] is immutable: a differing `amps_impose_global_time_step_` request rejects
before writing timestep storage. The existing complete-step scheduler is
retained. If the next full endpoint exceeds the host limit, AMPS acknowledges
that limit without particle work or sampling and preserves its completed
physical time in `swmfTimeSimulation`, including across ordinary sessions.

Each accepted step binds an immutable start/duration/end window. SEP freezes
its background and native birth window at departure, then calls the coupler's
explicit particle-phase boundary. That boundary presents the end counter,
matching standalone counter-OFF ordering; SWMF PIC does not increment it again.
Only a returned particle step with successful source-window closure can mark
the phase complete. The wrapper commits the physical ledger after the whole
application step returns. An error is fatal, not a claim of particle rollback.
Legacy applications and FLUID coupling retain their existing clock paths.

This contract does not reconstruct time-dependent SWMF backgrounds or change
sampling placement/normalization. Missing complete-window background coverage
still rejects; a newer import cannot be relabelled as an older realization.
Full PT restart continuation requires separately qualified physical/sampling/
background restoration; an SC/IH restart that initializes PT is not that gate.
`bash srcSEP/test/run_birth_window_tests.sh` checks the production value contract
with strict warnings and ASan/UBSan. The original eight-rank
`SEP--focused_transport-IH-mflampa-fl` run and energy/pitch-angle comparisons
remain required numerical regression gates.

PT also forwards an explicit source-line mapping before importing geometry.
The current CON_bline producer partitions a fixed `nLon*nLat` origin grid:
`BL_init` assigns the contiguous global range to each PT process, and
`BL_iblock_to_lon_lat` identifies a line by its one-based global grid label.
The wrapper checks that decomposition and forwards each label, origin radius
[solar radii], angular bounds [rad], dimensions and the same HGR frame used by
`BL_set_grid`. No moving vertex coordinate or disconnected PIC rank is used
as a physical source identity. Registration validates the complete candidate
and rejects changed labels or namespace metadata before replacing anything.
The first explicitly fresh PT step initializes the existing source scheduler
after real labelled lines arrive and before prepopulation. PT restart/source
history is never silently reconstructed. This supports this producer's fixed
origin-grid line namespace, not arbitrary line replacement or redistribution
of an already populated native state.

The same declaration forwards the producer's `cProtonMass` [kg] and `rSun` [m].
AMPS uses these actual export units when recovering mass density and importing
coordinates; it does not assume separately compiled constants are identical.
The received HGR/HGI frame explicitly establishes a heliocentric origin for
radius-dependent scattering coefficients.

For prescribed SEP transport, reception also prepares Cartesian velocity
derivatives [s^-1] from simultaneous spatial samples on all PT processes.
This uses the original PT communicator before PIC's disconnected partition;
particle workers consume the prepared immutable tensors without collectives.
Canonical producer-line/Lagrangian-node ordering makes stencil ties independent
of decomposition. Invalid packets, inconsistent epochs, deficient stencils and
replacement of an already consumed packet at the same epoch reject. Seed-only
bootstrap packets do not constitute a three-dimensional medium. The gradient
method and numerical controls enter the background configuration fingerprint.
This spatial reconstruction supplies no temporal interpolation contract.

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

PT also forwards the particle coordinate convention, separately from the HGR/HGI
vector basis. Its CON_bline material-node contract maps to
`plasma-comoving-v1`; changing it on a live population rejects. Material prefix
changes remap particle indices by the change in first producer label and refresh
Cartesian views from the new geometry. Removed segments retain absorbing
boundaries with count/weight diagnostics. Noncontiguous labels require a separate
explicit remeshing contract. Generic Eulerian PIC coupling does not enter this
field-line import path.

Imported material lines register their stable producer identity and geometry
epoch with the native container before sampling. SWMF owns their positions:
the local callback motion API must not also apply solar-wind advection. Import
updates the native geometry generation and first material label, then rebuilds
the MPI segment-ownership table after all line changes so exchange cannot use
retired segment pointers. See [the material-container contract](../src/pic/README.md)
for per-segment exposure normalization and `.sampling` output metadata. This
registration does not qualify moving-topology checkpoint continuation; the old
fixed-line codec rejects that unsupported state explicitly.
