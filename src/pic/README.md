# AMPS interpolation routines

This directory contains the AMPS particle-in-cell interpolation infrastructure. The principal implementation is in `pic_interpolation_routines.cpp`; the public types and entry points are declared in `pic.h` under `PIC::InterpolationRoutines`.

## Order-preserving performance work (2026-11-10)

The fast-wave function profile led to four bounded optimizations. Periodic
associated-data exchange now aggregates independent copies per peer and
retains ordered streaming for custom/uncertain cases; the exact payload and
ghost-write order remain unchanged. Fully allocated uniform pointer stencils
use a guarded complete-entry path. Deposition reconstructs per-call block
geometry/eligibility and uniform-grid field-node address views, preserving
coordinate units, arithmetic order and the existing worker schedule.

Particle signatures still hash exactly the canonical compiled record bytes;
single-copy canonicalization, descriptor capacity reuse and a full-width-
compatible host CRC implementation reduce work without disabling repeatability
or list checks. The host scalar sparse-multiply adapter uses a compact view
bounded to one fixed-operator solve. Iteration limits, solver tolerances,
physical operators, weights, seeds, ranks and sampled products are unchanged.

Detailed contracts, fallback domains and validation commands:

- [Periodic transport](parallel/README.md).
- [Complete stencils and deposition views](interpolation/README.md).
- [CRC and compact arithmetic](../general/README.md).
- [Linear-system snapshot lifecycle](../../srcInterface/README.md).
- [Independent performance-kernel cases](../test/tests/pic/performance/README.md).
- [Implementation ledger and measured wave results](../../test/srcFastWave/STATUS.md).

Adaptive interpolation, missing/unallocated support and row-only stencils
keep their general paths; evolving field values are never cached across a
step. Explicit CUDA/AVX/std::thread sparse variants retain their previous
kernels. Backend qualification is not inferred from a scalar timing gain.

For the application-facing contract for associated center/corner storage,
Tecplot variable/data callbacks, MPI ownership, sampled-data buffers, and
temporary output-node interpolation, see the core developer guide in
[`../README.md`](../README.md).

## Sampling-window lifetime and visualization flags

`PIC::Sampling::Sampling()` exchanges completed cell samples and switches the
completed/collecting offsets at a sample boundary. It then clears the new
collecting storage exactly once on every locally owned cell, across all
compiled species. The just-completed storage remains readable by output and
external completed-sample callbacks. This reset is independent of
`SaveOutputDataFile`, the minimum output sample length, and visualization
cadence. The per-species flag controls file creation only. Tying reset to that
flag caused alternating previous windows to accumulate when every flag was
false, including the srcMoon I11 duration pilot.

The production-path regression and constant-window oracle are documented in
[`../../srcMoon/test/integration-tests/I11/SAMPLING_CONTRACT.md`](../../srcMoon/test/integration-tests/I11/SAMPLING_CONTRACT.md).
The fixture tests three consecutive windows with all species files disabled,
only the first enabled, and all enabled. Its callback checks completed storage
before the next collection, and it checks zero collecting storage afterward.
Scientific sampling values and the normal enabled-output regression must both
remain qualified; a storage-fixture PASS alone is not physical validation.

## Particle-mover selection

Ordinary host-side single-particle advances now pass through
`PIC::Mover::DispatchParticleMover(ptr, dt, node)`. The compile-time switch
`_PIC_PARTICLE_MOVER_LEGACY_SETTINGS_` is defined next to the historical mover
macro in `picGlobal.dfn` and defaults to `_PIC_MODE_ON_`. In that mode the
dispatcher expands `_PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_`, preserving
existing generated configurations exactly.

An application that needs input-file selection must set the switch to
`_PIC_MODE_OFF_` in its generated configuration and install a callback before
the first particle advance:

```cpp
PIC::Mover::SetUserDefinedParticleMover(MyMover);
```

`MyMover` has type `PIC::Mover::fSpeciesDependentParticleMover`:

```cpp
int MyMover(long int particle, double dt,
            cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
```

The public pointer `PIC::Mover::UserDefinedParticleMover` defaults to
`PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder`.
The setter rejects null, so pointer mode cannot defer an empty-callback error
until a populated timestep. The callback inherits the legacy ownership
contract: it must perform any deletion, destination-list insertion, or
migration staging required by the configured particle-list mode. PIC discards
the integer result at these dispatch sites, as the prior macro call sites did.

The dispatcher covers the normal block/cell/particle loops, field-line loops,
open-flow remainder advance, and initialization mode `_MOVE_`. Specialized
boundary-injection mover macros and accelerator-specific bulk kernels retain
their own established selection paths. Changing the switch does not select a
physics model by itself; model parsing and callback registration are
application responsibilities.

## T4/T5 validated physics seams

The configured core-test application supplies uniform SI electric and
magnetic fields through generated hooks. The default relativistic-Boris hook
still performs the historical CPLR reads, so applications not selecting the
test provider are unchanged. `PIC::GYROKINETIC::EvaluateRhsFromFields` exposes
the local first-order guiding-center equations for provider validation. It
returns typed neutral, zero-field, non-finite, and invalid-magnetic-moment
states; it does not claim that gyro-radius/gradient scale separation is valid
for a particle trajectory.

Boundary/source helpers now centralize two distributions used by production:
the exponential waiting time `-log(U)/lambda` and the uniform-ball radius
`R*cbrt(U)`. Non-positive/non-finite rates return an infinite no-event time;
invalid radii return NaN. The corresponding one-rank tests cover their
analytical moments/CDF but do not establish distributed source ownership.

The collision generator records whether HS, VHS, or VSS was compiled, and the
elastic scattering dispatcher uses the matching angular law. VHS/HS are
isotropic; VSS uses `cos(chi)=2*U^(1/alpha)-1` with a declared arithmetic
binary alpha. All five cell selectors call this dispatcher. NTC and MF report
individual-record-weight support; their simplified forms require fixed record
weight, and Yinsi reports individual weights unimplemented. Collision trial
limiting uses a dimensionless per-model-particle cap and an unbiased sample
multiplier; a cell-local selector never overwrites the global threshold.

The LB rotational collision number is read from the generated species table.
One rotating partner or equal values for two rotating partners are
unambiguous; unequal rotating-pair values are rejected until a physical
mixing law is selected. RT is the validated channel. VT remains explicitly
unverified/terminating and VV is not implemented. Exact covered tests and
remaining T4/T5 boundaries are listed in
[`../../test/srcCoreUnitTests/README.md`](../../test/srcCoreUnitTests/README.md).

## Optional event inside a first-order trajectory

The first-order force-free trace mover has a default-inactive, three-callback
interface for an application event that must occur at an exact point inside a
particle trajectory. The generic core does not define event physics. All
three pointers begin as `NULL`, and an application must install either the
complete set or none before `PIC::Mover::Init()`:

```cpp
PIC::Mover::QueryTrajectoryEvent;
PIC::Mover::AdvanceTrajectoryEvent;
PIC::Mover::ProcessTrajectoryEvent;
```

`QueryTrajectoryEvent` is a side-effect-free preview. Given the current
species, owner block, local timestep, and maximum geometrically allowed
segment, it returns either `TrajectoryEventNone` or `TrajectoryEventPending`
with a finite event time in that segment. The mover compares this time with
block faces and internal surfaces. A strictly earlier event wins. Exact
equality gives the geometric boundary precedence; the event state is copied
through the normal face transfer and is processed as a zero-time event by the
next local or MPI continuation. Equality at an ordinary timestep endpoint has
no competing boundary and remains an event.

After moving a positive interval, the core calls `AdvanceTrajectoryEvent`
before any face-resampling clone or MPI exchange can copy the particle. The
application must commit only the elapsed event-clock state here. A useful
local-timestep implementation stores a dimensionless remaining event fraction
and subtracts `elapsedTime/localTimeStep`; the existing continuation machinery
then preserves that fraction while changing timestep between blocks.

At the event point, `ProcessTrajectoryEvent` receives mutable position,
velocity, species, and remaining mover time. It may update application bytes,
change species/velocity, and rescale the remaining interval for a new species.
It returns `TrajectoryEventContinue` or `TrajectoryEventDelete`. The callback
must not allocate, delete, or link particles: the caller must preallocate any
branch particles before movement, and the mover remains the sole owner of
particle-buffer deletion and destination-list insertion. A continuing event
must leave a valid allocated particle, matching particle-buffer species,
finite phase space and remaining time, and—when compiled—positive finite
individual correction. The event marker must be cleared before returning;
the mover re-queries a zero interval to prevent a forgotten marker from
causing an infinite loop.

This interface is currently implemented only by
`UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection`, the
first-order spatial mover used with local-timestep face continuation. An
application must reject the mode with second-order, Boris, guiding-center,
field-line, CUDA, or hybrid movers until those paths implement and validate
equivalent semantics. With all pointers `NULL`, no callback or random draw is
made and the historical trajectory path is retained.

## Local-timestep block-face flux resampling

The opt-in `_PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_` makes the trace mover
stop exactly at each AMR block face and complete the unused local-step fraction
through exchange/continuation sweeps. It supports a species-dependent local
timestep with either a supported spatially global particle weight or a
species-dependent particle weight stored in each block. A positive finite
individual correction is retained by the parent and every clone: it cancels
from the face multiplicity but remains in each represented physical moment.
CUDA, field-line, and periodic configurations are rejected at compile time
because those combinations do not implement this continuation contract.

For a particle moving from block `old` to block `new`, one model particle
carries physical flux `W/dt`. The destination multiplicity is drawn by
residual resampling with mean

```text
mu = (W_old/dt_old)/(W_new/dt_new)
   = W_old*dt_new/(W_new*dt_old).
```

The integer result is `floor(mu)` plus one with probability
`mu-floor(mu)`. For `mu<1`, deletion therefore has probability `1-mu`; for
`mu>1`, the original is retained and the remaining destination multiplicity is
created by cloning. Position and velocity are unchanged, so number, mass,
momentum, and kinetic-energy flux are conserved in expectation. The unused
trajectory interval is separately scaled by `dt_new/dt_old`, preserving its
dimensionless fraction of a local step; particle weight does not enter that
kinematic scaling.

Block-local weights are not stored in each particle record. They travel with
the AMR block and are read from the source and destination blocks at the face;
normal block migration already serializes the block-associated data. The
optional per-particle continuation fraction is part of the dynamically sized
particle record and is therefore copied by cloning, MPI exchange, and restart
packing. Applications must initialize valid positive finite weights and
timesteps on every allocated block copy before moving particles.

## Species-local particle restart

`PIC::Restart::SaveParticleData()` supports
`_SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_`. The local format stores:

- the global minimum-weight cache used by core bookkeeping;
- the application-visible block-data byte length; and
- the complete block-associated region beginning at
  `UserAssociatedDataOffset` for every owner block, including blocks with no
particles.

The last rule is essential: an empty block still owns local timesteps, local
particle weights, and application state that controls the first particle to
enter it. Populated records then append the existing per-cell particle counts
and particle payloads. The enveloped global-weight identity remains revision
2. Local-weight files use a distinct revision-3 layout identity including
weight mode and block-slice length, inside the checksum-protected AMPSPF1
envelope and atomic publication workflow. Unversioned local-weight streams
are rejected because their additional state has no authenticated layout.

On restore, each current owner consumes its block record. AMPS immediately
calls `ParallelBlockDataExchange()` before returning from the local-weight
restart, publishing restored W/dt and application bytes to boundary ghosts
before any face mover can use them. Restart/checksum tree recursion uses
`1<<DIM`; a hard-coded eight-child traversal is invalid in DIM1 and DIM2.

The core validates the application block-slice length but cannot identify its
semantic record order. Applications that persist structured state in this
slice must include their own layout version or fingerprint. The generic comet
systematic-residual state stores a 53-bit reaction/product layout fingerprint
in each initialized block and rejects a mismatched continuation.

Required validation for changes to this format includes empty and populated
blocks, owner-to-ghost publication, particle checksum continuation, state
migration across MPI owners, and DIM1 plus an unrelated DIM3 global-weight
restart. A local-weight file made by an executable with a different block-data
length is rejected; restart files are not a portable interchange format.

## 1. Interpolation representations

AMPS now provides two complementary representations of a cell-centered interpolation stencil.

### 1.1 Legacy pointer stencil

```cpp
PIC::InterpolationRoutines::CellCentered::cStencil
```

The legacy stencil stores:

- a pointer to an allocated `PIC::Mesh::cDataCenterNode`;
- the local center-node number; and
- the interpolation weight.

This representation is used throughout the existing particle movers, couplers, and field solvers. It requires the center node to be allocated on the current MPI rank. Cells for which no center-node pointer is available are omitted, and the existing fallback or normalization behavior is retained.

### 1.2 Decomposition-independent row stencil

```cpp
PIC::InterpolationRoutines::cRowStencil
```

The row stencil stores, for every physical cell participating in interpolation:

- `Element[s].node`: the AMR tree node that owns the physical cell;
- `Element[s].i`, `Element[s].j`, and `Element[s].k`: interior cell indices local to that tree node; and
- `Element[s].Weight`: the interpolation weight.

The row stencil does not contain a `cDataBlockAMR` or `cDataCenterNode` pointer. It is therefore valid even when the corresponding block is not allocated on the current MPI rank. The AMR tree geometry is used to canonicalize ghost indices to the actual owning leaf node and interior `(i,j,k)` indices.

This representation is intended for data sets that are globally replicated in a compact form, such as magnetic and electric fields used for cutoff-rigidity trajectory tracing. It avoids globally replicating every AMPS cell state vector.

## 2. Cell-centered constant interpolation

`PIC::InterpolationRoutines::CellCentered::Constant::InitStencil()` locates the cell containing the interpolation point and returns a one-element pointer stencil with weight one.

The linear interpolation implementation uses this routine as its legacy fallback when a complete locally allocated linear stencil cannot be constructed. The fallback behavior has not been changed.

## 3. Cell-centered linear interpolation

The public entry point is:

```cpp
PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(
    double *x,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,
    PIC::InterpolationRoutines::CellCentered::cStencil& Stencil,
    PIC::InterpolationRoutines::cRowStencil *RowStencil = NULL);
```

When `RowStencil == NULL`, the call is backward compatible with the previous implementation. Existing call sites do not need to change.

When a row stencil is requested, AMPS builds both representations in parallel:

- the legacy pointer stencil preserves the previous local-allocation behavior;
- the row stencil contains the complete geometric interpolation row, including cells on remote, unallocated blocks.

A row-only convenience overload is also available:

```cpp
PIC::InterpolationRoutines::cRowStencil row;
PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,row);
```

The convenience overload internally creates a temporary legacy stencil so that all interpolation decisions remain centralized in the same implementation.

### 3.1 Trilinear interpolation inside a uniform-resolution neighborhood

For a point with block-local continuous coordinates `(iLoc,jLoc,kLoc)`, the lower cell-center indices are

```text
i0 = floor(iLoc - 0.5)
j0 = floor(jLoc - 0.5)
k0 = floor(kLoc - 0.5)
```

with special handling below the first center. The eight support cells are `(i0+di,j0+dj,k0+dk)`, where `di`, `dj`, and `dk` are zero or one.

The weights are the standard tensor-product trilinear coefficients:

```text
wx(0) = 1-xi, wx(1) = xi
wy(0) = 1-eta, wy(1) = eta
wz(0) = 1-zeta, wz(1) = zeta

w(di,dj,dk) = wx(di)*wy(dj)*wz(dk)
```

The legacy pointer stencil adds an entry only when the corresponding center-node pointer exists. The row stencil always attempts to add the geometric cell. A support index that is a block ghost index is converted to a physical center coordinate, resolved through the replicated AMR tree, and stored as an interior index of the owning tree node.

At a non-periodic global boundary, support centers outside the global box are omitted and the retained row weights are normalized. For periodic meshes, support coordinates are wrapped into the primary global box before the owning node is resolved.

### 3.2 AMR coarse/fine interfaces

Near a refinement interface, AMPS constructs interpolation on a logical coarse-cell-center lattice.

For each of the eight logical coarse support centers:

- if the center is covered by a block at the reference refinement level, one physical cell is used;
- if the center is covered by a block one level finer, the logical value is represented by the average of the corresponding `2 x 2 x 2` fine cells.

The outer trilinear coefficient is then applied to the logical value. A fine cell therefore receives

```text
logical trilinear weight * 1/8
```

The row stencil records the actual owning node and fine-cell indices without requiring `StencilNode->block` to be allocated.

The legacy pointer stencil still uses the original center-node and ghost-cell access strategy. Duplicate physical center-node pointers are merged exactly as in the previous `cStencilGeneric::Add()` path.

### 3.3 Continuous transition on the fine side

When a point lies between one-half and one fine-cell width from a coarse/fine interface, AMPS blends the coarse-lattice and fine-grid stencils:

```text
alpha = (dmin - 0.5)/0.5
result = (1-alpha)*coarse_stencil + alpha*fine_stencil
```

The same blending coefficients are applied to both the legacy stencil and the row stencil. Duplicate row entries are merged by `(node,i,j,k)`.

### 3.4 Native AMPS and SWMF interpolation backends

The native AMPS backend builds trilinear and multiblock stencils directly from the AMR tree.

The SWMF backend receives weights and cell identifiers from the external interpolation library. The returned AMR node and cell indices are copied into the row stencil before local block allocation is checked. The existing OpenMP restriction of the SWMF backend remains unchanged.

## 4. `cRowStencil` interface

### 4.1 Stored data

```cpp
int row.Length;
row.Element[s].node;
row.Element[s].i;
row.Element[s].j;
row.Element[s].k;
row.Element[s].Weight;
```

Entries identifying the same `(node,i,j,k)` cell are merged automatically by `AddCell()`.

### 4.2 Coordinate lookup

```cpp
double xCell[3];
row.GetElementCoordinate(s,xCell);
```

`GetElementCoordinate()` computes the physical cell-center coordinate from the tree-node bounds and local indices. It never accesses `node->block`, so it works for remote cells.

### 4.3 Removing one point and renormalizing

```cpp
row.RemoveElement(s);
```

By default, the remaining weights are immediately renormalized to unit sum. The method returns `false` for an invalid index and `true` after a successful removal.

For removing several entries, it is more efficient and numerically cleaner to defer normalization:

```cpp
for (int s=row.Length-1;s>=0;s--) {
  double xCell[3];
  row.GetElementCoordinate(s,xCell);

  if (ExcludeCell(xCell)) row.RemoveElement(s,false);
}

row.Normalize();
```

Iterating backward avoids skipping an element when the array is compacted after removal.

### 4.4 Weight operations

```cpp
row.GetWeightSum();
row.Normalize();
row.MultiplyScalar(a);
row.Add(&otherRow);
row.flush();
```

`Add()` merges duplicate physical cells. `Normalize()` leaves an empty or zero-sum row unchanged.

## 5. Example: interpolation from a compact global magnetic field

Assume the application stores only the global magnetic field, indexed by AMR node and local interior cell indices:

```cpp
PIC::InterpolationRoutines::CellCentered::cStencil localStencil;
PIC::InterpolationRoutines::cRowStencil row;

double x[3]={x0,x1,x2};
auto *node=PIC::Mesh::mesh->findTreeNode(x);

PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(
    x,node,localStencil,&row);

double B[3]={0.0,0.0,0.0};

for (int s=0;s<row.Length;s++) {
  const auto &e=row.Element[s];
  const double *Bcell=GlobalField.GetMagneticField(e.node,e.i,e.j,e.k);

  for (int d=0;d<3;d++) B[d]+=e.Weight*Bcell[d];
}
```

The global field container is application-specific. Within one MPI process, the AMR node pointer can be used directly as a lookup key. If field identifiers are serialized or exchanged between processes, use a stable AMR node identifier rather than the numerical pointer value because replicated tree nodes have different addresses in different processes.

## 6. Example: reject invalid global-field points

A global field may omit cells inside a body, outside a validity mask, or beyond a data-product boundary. The row can be filtered before interpolation:

```cpp
PIC::InterpolationRoutines::cRowStencil row;
PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,row);

for (int s=row.Length-1;s>=0;s--) {
  double xCell[3];
  row.GetElementCoordinate(s,xCell);

  const auto &e=row.Element[s];
  if (!GlobalField.IsValid(e.node,e.i,e.j,e.k,xCell)) {
    row.RemoveElement(s,false);
  }
}

row.Normalize();
```

If all points are removed, `row.Length` is zero and the caller must select an application-specific fallback.

## 7. Compatibility behavior

The modification is additive.

- Existing three-argument calls to `Linear::InitStencil()` compile unchanged because the row pointer defaults to `NULL`.
- Existing calls through the thread-local `StencilTable` wrappers are unchanged.
- When no row is requested, the legacy allocation checks, missing-cell behavior, constant fallback, and weight normalization path remain active.
- When a row is requested for a point whose containing block is remote, the row is produced from tree geometry and the legacy stencil is returned empty because no valid local center-node pointers exist.
- When the containing block is local but one or more support center nodes are unavailable, the legacy stencil retains its previous constant-fallback behavior while the row retains the complete linear interpolation geometry.

## 8. Current assumptions

The AMR multiblock interpolation algorithm retains its existing assumptions:

- adjacent refinement levels use the AMPS `2:1` refinement relation;
- block dimensions are compatible with the `2 x 2 x 2` fine-cell averaging rule;
- the AMR tree is available on every rank;
- `IsUsedInCalculationFlag` identifies active leaf blocks;
- the SWMF external interpolation library obeys its existing block-index contract.

The row stencil capacity is `2*nMaxStencilLength`, allowing the current maximum 64-entry multiblock row to be blended with an additional fine-grid stencil without changing the legacy stencil capacity.

## 9. Application use of particle splitting and merging

The generic linked-list helpers are declared in `pic.h` and implemented in
`pic_particle_spliting.cpp`:

```cpp
PIC::ParticleSplitting::MergeParticleList(spec, firstParticle, target);
PIC::ParticleSplitting::SplitParticleList(spec, firstParticle, target);
```

They operate on one AMPS particle list and one species. The implementation now
returns without mutation for a non-positive target, an empty merge list, or a
split list containing fewer than two particles. These guards are required:
the bin algorithms otherwise have no valid `particles[0]` or iterator
successor. A caller must still invoke them only when particle movement/list
exchange is quiescent and must pass the actual cell-list head by reference.

The complete no-op contract at the public entry points is:

| Routine | Condition | Reason no mutation is permitted |
|---|---|---|
| both | `target <= 0` | The requested population is invalid/no-op, and a negative `int` must not be converted to a large vector capacity or unsigned count. |
| merge | selected population is empty or already `<= target` | There is no reduction to perform and no valid first velocity record. |
| merge | selected population contains fewer than three records | The conservative implementation is specifically a 3-to-2 operation. |
| split | selected population contains fewer than two records | The conservative implementation is specifically a 2-to-3 operation. |
| both | no velocity bin contains a viable tuple at the current resolution | The routine coarsens the binning before retrying; it never advances an end iterator. |

These are core precondition guards, not an alternative population-control
policy. A positive target is still only a stopping bound: the tuple size,
weight threshold, and conservation constraints can prevent an exact requested
count from being reached.

These routines conserve the generic AMPS record they know about; they cannot
infer every application's phase-space representation. In particular:

- `CloneParticle` preserves AMPS `next`/`prev` links but copies every
  application-requested extension byte. An application storing a unique ID,
  random-stream state, collision history, or optical depth must assign valid
  child metadata after cloning.
- The current generic merge reconstructs Cartesian velocity from weighted
  momentum and a weighted nonrelativistic `v^2` energy moment. Relativistic
  applications must not interpret that as conservation of
  `sqrt(p^2 c^2 + m^2 c^4)-m c^2`.
- Generic random directions use the AMPS process RNG. Applications requiring
  rank/order/restart invariant histories need a semantic keyed-random callback
  or an application-owned reconstruction.
- If an application stores momentum/pitch separately from Cartesian velocity,
  all representations must be updated together after a merge.

The recommended core extension is a resampling-policy interface supplied by
the application. The policy should provide (1) a split-extension callback,
(2) a merge reconstruction callback with the application's energy definition,
and (3) an optional deterministic/keyed direction source. The existing list
allocation/deletion and bin-selection code can then remain shared without
silently imposing a nonrelativistic or metadata-free model. Until that
interface exists, applications with such state should disable the automatic
`ParticleSplitting::Mode` and call an application-aware boundary controller
that uses `GetNewParticle`, `CloneParticle`, and `DeleteParticle` only for list
mechanics.

## Read-only accepted-particle event diagnostics

`PIC::ParticleEvents` declares a null-by-default observer in `pic.h`, with its
storage and particle-weight resolver in `pic.cpp`. Exosphere records accepted
emissions only after allocation/rejection checks, adsorption with the whole
removed gas weight, reflection as an informational event, actual chemical
parent losses/daughter creations, and signed prescribed reservoir supply.
The selected Moon first/second-order mover branches record outer-domain escape
and known numerical removals before releasing particle storage. Migration
between MPI ranks is not a physical event. The default observer changes no
RNG draws, particle data, output flags, force/source law or sampling buffers.

Applications own the observer's scope and counters. The supporting srcMoon
ledger currently requires neutral Na, a sphere and its declared global-weight
configuration; it independently sums owned real particle lists and records
counts/atoms, never just rates. Uninstrumented paths show up as unexplained
residuals. This is not a claim that all PIC movers/species/GPU configurations
have complete event instrumentation. See
`srcMoon/test/integration-tests/I11/BUDGET_CONTRACT.md` for the exact oracle,
MPI reductions, weighted roundoff bound, units, timing and tests.

The November 12 extension appends `SurfaceDeposit` and `SurfaceRelease`
without changing existing event numbers. These describe represented atoms in
the surface ledger and are informational for the older Na-only gas-count
identity. `DomainExitClassifier` is optional: the application classifies an
actual boundary exit from its current position/velocity before deletion; a
null classifier keeps the previous Escape event. Moon uses the existing I09
instantaneous lunar-energy definition, distinguishing bound finite-box loss
from positive-energy lunar escape. This is not ultimate three-body escape.
`MoonM2Budget.cpp` independently sums owned Na/Ar gas lists and logical surface
stocks counted once; its supporting scope does not cover every PIC mover.

`cMovingInjectionEpochScope` temporarily supplies a source departure offset
[s] during source placement and shortened boundary injection. Outside that
scope its NaN sentinel selects the ordinary full-step trajectory. An injection
scope is a complete no-op for that sentinel: it neither reads nor writes the
shared departure context, preserving parallel legacy injectors without moving
geometry. Finite-offset moving injection retains its documented threading
restriction. The selected
force SecondOrder entry now routes to its boundary tracer when the moving
surface callback exists. Collision search uses that mover's constant initial
acceleration position extension, and hit-velocity evaluation uses a shortened
midpoint force evaluation. This is still a second-order numerical trajectory;
time/space convergence is required. Static-boundary and null-callback modes
retain their existing dispatch.

## 10. CCMC nightly trajectory receipt

`PIC::CCMC::TraceParticles()` asks the particle tracker to write one final
Tecplot trajectory file per compiled species. In nightly-test mode, rank zero
then concatenates those already-complete files into
`test_CCMC-Individual_Trajectories.dat` for the legacy reference comparison.
This last operation is file packaging only: MPI collection, trajectory
selection/thinning, coordinates, velocities, timestamps and particle state
have already been finalized by `PIC::ParticleTracker::OutputTrajectory()`.

The merge is implemented inside `pic_ccmc.cpp`; it does not invoke a shell.
Input names are reconstructed from the same trajectory base, numerical species
index and chemical symbol used by `CreateTrajectoryOutputFiles()`, then sorted
lexically before copying. Consequently, a stale file belonging to a species
that is not compiled into the current executable cannot enter the receipt.
Paths may contain spaces or shell metacharacters. The configured output
directory retains its existing `_MAX_STRING_LENGTH_PIC_` contract, but adding
the input/output suffixes now uses dynamically sized strings instead of
overflowing a second fixed-size command buffer.

Rank zero writes a bounded-buffer, byte-for-byte concatenation to
`test_CCMC-Individual_Trajectories.dat.tmp`. Every input-open, read, write and
close operation is checked. Only after the complete temporary file closes
successfully is it renamed over the public receipt in the same directory. On
failure, the temporary file is removed and AMPS terminates; an older complete
receipt is not truncated or presented as new evidence. The same-directory
rename is atomic on the POSIX filesystems supported by native AMPS runs.

The formatting regression can be checked without running the physics case by
compiling `pic_ccmc.cpp` with the configured generated headers and
`-Wformat=2 -Wformat-overflow=2`. End-to-end qualification remains the
`CCMC-Individual_Trajectories` nightly case, whose reference comparison checks
the merged bytes. A native rebuild must follow the repository preflight:
confirm the Mars2 root and that no build/test process is active, remove only
the root `build`, regenerate the selected application and production hooks,
then compile through the top-level workflow with `-j16`.

## 11. Optional pre-exchange and transactional photochemistry hooks

Two nullable interfaces allow an application to control particle-buffer
capacity without embedding its population policy in the AMPS core.

`PIC::Parallel::BeforeParticleDataExchange` is called collectively immediately
before either ordinary exchange implementation allocates receive-side particle
slots. The default is `NULL`, so existing applications retain their previous
behavior. An installed callback may inspect and resample application particle
populations and may use MPI collectives, but it must not call
`ExchangeParticleData()` recursively. The location of this callback is part of
its contract: moving it after the buffered/unbuffered dispatch would be too
late to protect receive allocation.

`PIC::ChemicalReactions::PhotolyticReactions` also provides a complete
five-callback transactional batch interface:

- `BeginBatchedModel` prepares per-stage application state;
- `IsBatchedParticleEligible` distinguishes deferred parents from completed
  parents and newly created daughters;
- `GetBatchedParticleLimit` supplies the maximum parents attempted in the
  current pass;
- `BatchedReactionProcessor` commits one complete event or reinserts an
  unchanged deferred parent and returns `false`; and
- `FinishBatchedPass` is the collective-safe checkpoint at which the
  application can create headroom before another pass.

All five pointers must be either null or installed. A partial set terminates
with an error. With a complete set, the core rebuilds cell lists on each pass,
keeps ineligible particles without processing them again, globally counts
remaining eligible parents, and repeats until the list is drained. The
application owns event atomicity and eligibility metadata; the generic core
does not assume a reaction network, statistical-weight law, or buffer
watermark. The historical reaction processor and single-pass behavior are
unchanged when all pointers are null.

The current batch counters and application checkpoint protocol are validated
for MPI compilation only. A hybrid/OpenMP build explicitly rejects an
installed batch interface. Before a future master integration enables hybrid
use, counters and list rebuilding must be made thread-safe and tested with
multiple OpenMP threads. Merge `pic.h` and `pic_photolytic_reactions.cpp`
together; merging declarations without the null definitions/driver, or the
driver without the callback contract, creates an incomplete core API. The
pre-exchange hook in `pic_parallel.cpp`/`pic.h` is a separate nullable core
primitive and should be reviewed as a second semantic unit.

## 12. Reusable NTC per-species particle buckets

The general `ModelCellCollisions_ntc` path used by species-local weights and
timesteps no longer allocates two scratch arrays for every occupied cell or
rescans the cell linked list once for every species. It owns one
`thread_local std::vector<cParticleDataList>` per compiled species. For each
cell, one linked-list traversal fills those buckets and accumulates the
species particle count, minimum individual weight correction, and weight sum.
`clear()` empties a bucket without releasing its capacity, so later cells and
timesteps reuse the allocation. `thread_local` is required: making this state
process-global would introduce a race in OpenMP builds.

The optimization does not change NTC physics. Particle order within a species
is the original linked-list order; the fixed 100-probe majorant construction,
pair-loop order, acceptance logic, and random-number calls are unchanged.
Velocity changes are retained in the common bucket, so a species used by a
later pair observes accepted collisions from earlier pairs exactly as it did
when the old implementation repopulated its scratch array from the particle
buffer. Flags are reset before a species becomes the first partner, and before
a distinct second partner, then accepted updates are written back at the same
pair/outer-loop boundaries as before.

This is intentionally limited to the general NTC routine. Specialized
uniform-weight collision paths were not changed and should not be folded into
this implementation during master integration without separate deterministic
and statistical validation.

## 13. Individual corrections with local face transfer and general NTC

Stage-4 comet work adds a narrowly scoped correction-aware path to the generic
core. The feature remains selected by existing compile-time configuration and
does not change correction-off builds.

### Local-timestep face transfer

`pic_mover.cpp` permits individual corrections when local-timestep flux
resampling is enabled. The supported policy holds a particle's individual
correction `c` fixed across a block face. Because the physical rate is
`W*c/dt`, `c` cancels from the residual split/roulette multiplicity:

```text
mu = W_old*dt_new/(W_new*dt_old).
```

The original particle and all clones retain the complete particle record,
including `c`, rotational energy, and vibrational state. The mover rejects a
non-positive or non-finite correction before changing a particle. Its opt-in
validation diagnostic uses `W*c/dt` for number, mass, momentum, kinetic
energy, and internal energy; the variance terms contain the corresponding
squared outgoing rate. The diagnostic is not a production sampler.

### General NTC collisions

For species bucket `s`, particle `i` has physical carrier
`q_i=W_s*c_i`. `ModelCellCollisions_ntc()` builds strictly increasing
correction prefix sums and selects each trial particle with probability
`q_i/sum(q)`. Cross-species trials use

```text
N_trial = sum_i(q_i) sum_j(q_j) sigma_g,max dt_max / (q_min V),
```

with the usual one-half factor for one-species unordered pairs. A selected
particle is updated with probability `q_min/q_i` times its local-timestep
ratio. The implementation validates that every probability is finite and in
`[0,1]`, allowing only a roundoff-sized clamp. Same-species draws are
independent; a self pair is discarded rather than redrawn, preserving the
ordered-pair normalization and preventing the historical one-particle loop.
The collision limiter uses `ancoll/(N0*N1)` with explicit parentheses.

`cNTCValidationDiagnostic` in `pic.h` is disabled by default. A focused MPI
validator may enable it to count trials, accepted collisions, self-pair
rejections, particle updates, invalid probabilities, and probability extrema.
Enabling that process-local diagnostic in a hybrid build fails explicitly
until a thread-safe reduction exists.

The simplified NTC kernel is not correction-aware. In a correction-on build
it terminates at entry with a nonzero exit code and directs the caller to the
general NTC implementation. Application setup must similarly reject MF or
any unvalidated collision implementation; the core does not reinterpret
those models.

### Qualification and merge constraints

The comet validator covers corrections `1/16, 1/4, 1, 4, 16`, all-correction-
one identity cases, local/remote/multiple faces, internal energy, adaptation,
sampling, collision relaxation, boundary handling, restart, migration, and
chemistry at 1, 4, and 12 MPI ranks. Evidence is listed in
`srcHartley2/CHEMISTRY.md` and `srcHartley2/STATUS.md`.

For a later master integration, treat the core work as two reviewable semantic
units:

1. merge the `pic.h` flux-diagnostic fields with `pic_mover.cpp`; and
2. merge the `pic.h` NTC diagnostic interface with
   `pic_particle_collision_models.cpp`.

Do not merge only the declarations, weaken the simplified-NTC failure, or
broaden the supported mover/collision matrix during conflict resolution.
Retain default-off behavior and repeat correction-off plus correction-on
1/4/12-rank tests after a clean regenerated build. CUDA, hybrid diagnostic
reduction, field-line/periodic movers, and non-radial correction-on chemistry
remain outside this qualification.
