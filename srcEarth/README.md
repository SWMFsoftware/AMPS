# SEP-in-geospace model: compact global fields in Mode3D

## Purpose

Mode3D backward calculations trace an energetic particle through the complete AMR domain. The AMPS mesh uses normal MPI domain decomposition: every rank has the global AMR tree, but only the owner domain and a limited neighbor layer have allocated `cDataBlockAMR` objects. A trajectory assigned to one MPI rank can enter any leaf of the global tree, so field interpolation must not depend on whether that leaf's block is allocated locally.

The previous implementation solved this by allocating every used AMR block on every MPI rank, gathering the magnetic field, and populating all replicated interior and ghost center nodes. That approach was correct for field access but replicated complete AMPS cell-associated state vectors and ghost storage. Its memory cost was much larger than the physical B field required by cutoff-rigidity and density/flux backtracing.

Mode3D now keeps the standard distributed AMPS mesh and replicates only compact cell-centered magnetic and electric field arrays. AMPS linear interpolation is evaluated with the decomposition-independent row stencil introduced in `src/pic`.

## Modified files

- `3d/GlobalMagneticField.h`
- `3d/GlobalMagneticField.cpp`
- `3d/CutoffRigidityMode3D.cpp`
- `3d/CutoffRigidityMode3D.h`
- `3d/Mode3D.cpp`
- `3d/DensityMode3D.cpp`
- `3d/DensityMode3D.h`
- `3d/Mode3DParallel.cpp`
- `3d/Mode3DParallel.h`
- `3d_forward_swmf/Mode3DForwardSWMF.cpp`
- `3d_forward_swmf/Mode3DForwardSWMF.h`

## Global cell indexing

AMPS replicates the AMR tree topology on every MPI rank. The temporary tree-node member `Temp_ID` is therefore used as a dense global leaf index.

Before each field snapshot is assembled, `GlobalMagneticField.cpp` recursively resets `Temp_ID` to `-1` on every tree node. It then traverses the tree in deterministic child order and assigns consecutive IDs to leaves satisfying `IsUsedInCalculationFlag`.

For a block containing `Nx`, `Ny`, and `Nz` interior cells, the scalar global cell index is

```text
cellsPerBlock = Nx * Ny * Nz
localCell     = i + Nx * (j + Ny*k)
globalCell    = node->Temp_ID * cellsPerBlock + localCell
```

The vector components are stored at `3*globalCell + component`.

`Temp_ID` is shared AMPS scratch storage. No other algorithm may overwrite it between compact-field assembly and completion of the corresponding backward calculation. Field evaluation checks every row-stencil node ID and terminates with a diagnostic if an invalid or stale `Temp_ID` is detected.

## Compact arrays

The implementation owns three process-local arrays that are identical on all MPI ranks after assembly:

```cpp
std::vector<double> GlobalMagneticField_; // 3 * NglobalCells
std::vector<double> GlobalElectricField_; // 3 * NglobalCells
std::vector<int>    GlobalCellPresence_;  // 1 * NglobalCells
```

Only physical interior cells are stored. There are no global ghost cells and no replicated AMPS data blocks.

The persistent field memory per MPI rank is approximately

```text
B: 24 * NglobalCells bytes
E: 24 * NglobalCells bytes
presence: 4 * NglobalCells bytes
```

The presence array is retained for defensive validation. Every used physical cell must be contributed by exactly one owner rank.

## Assembly procedure

The main entry point is

```cpp
Earth::Mode3D::GlobalMagneticField::AssembleCellCenteredFieldsForCutoff(
    diagnosticTag,
    magneticFieldDataOffset,
    electricFieldDataOffset,
    plasmaVelocityDataOffset,
    verbose);
```

The procedure is:

1. Mark the previous snapshot unavailable.
2. Reset `Temp_ID` throughout the global AMR tree.
3. Assign deterministic dense IDs to all used leaves.
4. Allocate compact B, E, and presence arrays.
5. Visit only leaves for which `node->Thread == PIC::ThisThread`.
6. Pack each allocated owner interior cell into its `Temp_ID`-based global location.
7. Use chunked in-place `MPI_Allreduce` operations to replicate B, E, and the presence counts.
8. Verify that every cell has exactly one owner contribution.
9. Mark the snapshot ready and keep the arrays read-only during trajectory tracing.

No call to `AllocateBlock()` is made by the compact-field module.

### Standalone Mode3D

`Mode3D::InitMeshFields()` evaluates the configured standalone B and E models on owner-rank AMPS blocks. `Mode3DPrepareMagneticFieldSnapshot()` then calls

```cpp
const long int eOffset =
    PIC::CPLR::DATAFILE::Offset::ElectricField.active ?
    Earth::Mode3D::GlobalMagneticField::DataFileElectricFieldDataOffset() : -1;

Earth::Mode3D::GlobalMagneticField::AssembleCellCenteredFieldsForCutoff(
    "Mode3D",
    Earth::Mode3D::GlobalMagneticField::DataFileMagneticFieldDataOffset(),
    eOffset,
    -1,
    true);
```

When the DATAFILE electric-field slot is inactive, the compact E array is explicitly zero.

### SWMF-coupled Mode3D

SWMF stores magnetic field and plasma velocity in separate coupler-buffer locations. It does not store an independent electric-field vector. `PrepareGlobalSWMFCoupledMagneticFieldForCutoff()` therefore calls

```cpp
Earth::Mode3D::GlobalMagneticField::AssembleCellCenteredFieldsForCutoff(
    "Mode3DForwardSWMF",
    PIC::CPLR::SWMF::MagneticFieldOffset,
    -1,
    PIC::CPLR::SWMF::BulkVelocityOffset,
    true);
```

For each owner cell, the compact electric field is calculated as

```text
E = -v x B
```

Both offsets and completion of the first SWMF coupling receive are checked before assembly.

## Row-stencil interpolation

The public field evaluators are

```cpp
bool InterpolateMagneticField(const double* x,cAMRNode* node,double* B);
bool InterpolateElectricField(const double* x,cAMRNode* node,double* E);
```

They build a stack-local row stencil:

```cpp
PIC::InterpolationRoutines::cRowStencil row;
PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,row);
```

Each row element contains

```cpp
row.Element[s].node;
row.Element[s].i;
row.Element[s].j;
row.Element[s].k;
row.Element[s].Weight;
```

The owning block can be unallocated on the current MPI rank. The field value is obtained directly from the compact array:

```cpp
for (int s=0;s<row.Length;s++) {
  const auto& e=row.Element[s];

  const long int localCell =
      e.i + _BLOCK_CELLS_X_*(e.j + _BLOCK_CELLS_Y_*e.k);

  const long int globalCell =
      e.node->Temp_ID*cellsPerBlock + localCell;

  B[0] += e.Weight*GlobalB[3*globalCell+0];
  B[1] += e.Weight*GlobalB[3*globalCell+1];
  B[2] += e.Weight*GlobalB[3*globalCell+2];
}
```

The row stencil preserves the AMPS cell-centered linear interpolation geometry, including:

- ordinary eight-cell trilinear interpolation;
- same-level block-boundary interpolation;
- coarse/fine multiblock reconstruction;
- averaging of the corresponding `2 x 2 x 2` fine cells for a logical coarse support point;
- continuous coarse/fine transition blending.

Unlike the legacy pointer stencil, a remote cell is not omitted merely because its center-node object is unavailable locally.

## Updated Mode3D field evaluator

`cMode3DMeshFieldEval::GetB_T()` now:

1. finds the containing leaf in the global AMR tree using its private `lastNode_` hint;
2. accepts a valid leaf even when `node->block == NULL`;
3. calls `GlobalMagneticField::InterpolateMagneticField()`;
4. returns the interpolated global-array value.

The evaluator no longer branches between DATAFILE and SWMF pointer-stencil access. The source-specific work is performed once during snapshot assembly. During particle tracing the field path is identical for standalone and coupled calculations.

The optional direct analytic-field debug path remains available and bypasses the compact arrays.

## Error handling

The compact-field implementation treats the following conditions as fatal consistency errors:

- an owner-cell gather misses a used interior cell;
- more than one rank contributes the same used interior cell;
- a row-stencil entry contains a non-interior index;
- a row-stencil node has an invalid `Temp_ID`;
- a requested global cell is absent from the assembled snapshot.

A missing row point is not removed and the remaining interpolation weights are not renormalized. Doing so would change the field and could make cutoff results depend on MPI decomposition. Row-point removal remains available to external applications through `cRowStencil`, but it must be an explicit physical masking decision rather than an automatic response to missing distributed data.

A field interpolation call returns `false` only when:

- no compact snapshot has been assembled;
- the requested point is outside the used AMR tree; or
- AMPS cannot construct an interpolation row.

## Compatibility

The old function

```cpp
MaterializeCellCenteredMagneticFieldForCutoff(tag,bOffset,verbose)
```

is retained as a B-only wrapper. It now assembles compact B plus a zero E array.

The old debug signature

```cpp
RedefineAllAllocatedMagneticField(...)
```

is also retained so existing source code compiles. Its implementation delegates to `RedefineGlobalMagneticField()` and intentionally does not allocate remote blocks.

Existing forward Monte Carlo transport continues to use the normal distributed AMPS cell buffers and is not redirected through these compact arrays. The change is specific to Mode3D backward products that require arbitrary global field access.

## Recommended validation

The following comparisons should be run before removing any legacy reference branch from a development tree:

1. **Uniform mesh:** compare old pointer-stencil and new row-stencil B at random positions.
2. **Same-level block faces:** sample both sides and exactly on the transition region.
3. **AMR coarse/fine faces, edges, and corners:** compare every physical row entry and final B.
4. **MPI decomposition invariance:** run the same snapshot with 1, 2, 4, and 8 ranks and compare fields and cutoff products bitwise where reduction order permits, otherwise to roundoff.
5. **Standalone DIPOLE:** compare interpolated B and cutoff maps with the former replicated-block implementation.
6. **SWMF snapshot:** compare compact B and derived E against direct owner-cell values.
7. **Time series:** verify that `Temp_ID` is reset and arrays are rebuilt for every field snapshot.
8. **Memory scaling:** confirm that increasing MPI rank count does not allocate additional global AMPS blocks and that per-rank growth is limited to compact B/E/presence arrays.

## Mode3D DIPOLE magnetic-field interpolation error statistics

For a mesh-backed `FIELD_MODEL DIPOLE` calculation, the exact analytic field is known at
all trajectory sample coordinates. Mode3D now uses this reference to measure the error of
the actual compact-array and row-stencil field value returned to the particle mover.

For each successful mesh-field evaluation, the diagnostic records

```text
relative_error = |B_mesh - B_dipole| / |B_dipole|
```

The statistics are sample-weighted. A coordinate visited repeatedly by the trajectory
integrator contributes repeatedly because those are the field determinations that affect
the calculated trajectory. Every field evaluator keeps a private count, error sum,
maximum error, and maximum-error coordinate. It merges those values into a rank-local
accumulator only when the evaluator is destroyed, avoiding a mutex in the high-frequency
field-evaluation path. At the end of cutoff or density/flux processing, MPI reductions
produce the global sample count, mean relative error, and maximum relative error.
`MPI_MAXLOC` identifies the rank containing the maximum and that rank broadcasts the
associated location.

Rank zero prints a summary similar to

```text
========== Mode3D DIPOLE magnetic-field interpolation error ==========
Calculation                 : Mode3D cutoff rigidity
Definition                  : |B_mesh-B_dipole|/|B_dipole|
Number of field samples     : ...
Sample mean relative error  : ...
Maximum relative error      : ...
Max-error location [m]      : x y z
Max-error location [km]     : x y z
Max-error location [Re]     : x y z
=======================================================================
```

The diagnostic is automatically enabled only when

- the selected field model is `DIPOLE`; and
- Mode3D uses mesh-backed field evaluation.

It is disabled for `-mode3d-field-eval ANALYTIC`, because that path directly returns the
reference field and would only report a trivial zero interpolation error. A new statistics
window is started independently for each cutoff and density/flux calculation, including
each separately processed time snapshot.
