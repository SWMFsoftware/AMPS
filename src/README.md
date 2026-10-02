# AMPS core guide: application-defined mesh output

This document describes the supported AMPS core interfaces for adding
application-owned quantities to mesh data files.  It is intentionally
application-neutral: a new AMPS application should be able to follow this
contract without copying initialization or output code from another
application.

The guide applies to data written by
`PIC::Mesh::mesh->outputMeshDataTECPLOT(...)` and the cut-cell-aware
`PIC::Mesh::mesh->OutputDistributedDataTECPLOT(...)`. These writers emit AMR
geometry, built-in particle moments, registered node quantities, and the
block-local particle time step and statistical weight. Their boundary geometry
and MPI print protocols differ; section 8 explains how to select them. The
related `outputMeshTECPLOT(...)` path is geometry-only and does **not** exercise
the data callbacks described here. Registering a spherical body or correcting
its cell measures alone does not make the whole-brick writer clip its surface.

The public declarations are in `pic/pic.h`.  The associated-data layout and
center-node callback dispatch are implemented in `pic/pic_mesh.cpp`, and the
Tecplot traversal is implemented in `meshAMR/meshAMRgeneric.h`.

## 1. Choose the physical meaning before choosing an API

AMPS supports three different classes of cell data.  They have different
lifetime, normalization, and restart requirements and must not be mixed.

| Data class | Examples | Storage request | Output rule |
|---|---|---|---|
| Persistent/static center state | magnetic field, fluid state, a model-valid flag | `PIC::IndividualModelSampling::RequestStaticCellData` | Initialize once per accepted state generation and interpolate it for output vertices. |
| Accumulated sampled state | particle count, a velocity moment, a reaction-rate sum | `PIC::IndividualModelSampling::RequestSamplingData` or `PIC::Datum` | Accumulate in the collecting sample, read the completed sample, and normalize exactly once. |
| Derived display state | field magnitude, pressure computed from stored primitives | no additional storage when all inputs already exist | Compute from accepted/interpolated primitive values in the print callback. |

Corner-node persistent state uses
`PIC::IndividualModelSampling::RequestStaticCellCornerData`.  Block state has a
different ABI: `PIC::Mesh::cDataBlockAMR` currently prints the compiled local
time step and particle weight.  There is no general application callback
registry for arbitrary block-output columns.  Prefer center- or corner-node
associated data for an application diagnostic rather than patching an
application-specific field into `cDataBlockAMR`.

Do not use a sampling buffer for an immutable background field, and do not use
a static field for a quantity that must be cleared at each sampling-window
transition.  A derived output callback must not rerun a physical model or
invent a replacement when source data are unavailable; it should report the
already accepted state and an explicit validity indicator.

## 2. Required initialization order

The associated-data layout becomes an ABI as soon as AMPS calculates its final
offsets.  A typical application must use this order:

```cpp
void amps_init_mesh() {
  // Initializes MPI, the AMPS mesh object, core output callbacks, and the
  // corner-data request registry.  Core models enabled at compile time can
  // also register their storage and output callbacks here.
  PIC::Init_BeforeParser();

  // Registers this application's storage requests and output callbacks.
  // This function must be idempotent or guarded against duplicate calls.
  MyApplication::Init_BeforeParser();

  // Invokes every registered request, assigns final byte offsets, and fixes
  // the completed/collecting sample-buffer locations.  No request may be
  // added after this call.
  PIC::Mesh::initCellSamplingDataBuffer();

  // Build/decompose the mesh and allocate blocks only after the final center,
  // corner, and block associated-data sizes are known.
  BuildDecomposeAndAllocateMesh();

  // Populate every owner-local physical node.  If interpolation can touch
  // ghost values, exchange the completed state before any data output.
  MyApplication::InitializeOwnerLocalCellState();
  PIC::Mesh::mesh->ParallelBlockDataExchange();
}
```

The important boundaries are:

1. Call `PIC::Init_BeforeParser()` before registering an application corner
   request because that routine constructs
   `RequestStaticCellCornerData`.
2. Register all center/corner static and sampled requests before
   `PIC::Mesh::initCellSamplingDataBuffer()`.
3. Freeze the layout before block allocation.  An offset returned afterward
   would not exist in already allocated node buffers.
4. Initialize owner-local nodes after allocation, and exchange complete state
   before interpolation or output.
5. Do not emit a data-bearing initialization file until background state,
   particle time steps, particle weights, application storage, and ghost cells
   all represent the same completed initialization state.

An application may parse and validate its run configuration before this
sequence, but it must not store node pointers or write associated buffers
before the mesh storage has been allocated.

## 3. Reserving persistent center-node storage

A request callback receives the next available byte offset and returns the
number of bytes it consumes.  Keep every application field relative to the
single base offset returned to the application.  This makes the layout easy to
audit and avoids unrelated absolute offsets scattered through the code.

```cpp
#include "pic.h"

#include <cstddef>
#include <cstring>

namespace MyApplication::CellState {

struct Layout {
  int baseOffset = -1;
  int magneticFieldOffset = -1;  // Three doubles: Bx, By, Bz [T].
  int validOffset = -1;          // One double: exactly 0.0 or 1.0.
  int byteCount = 0;
};

inline Layout layout;

int RequestStaticCellData(int nextFreeOffset) {
  // AMPS supplies byte offsets, not an alignment promise.  Store byte
  // positions here and use memcpy at the access boundary below.
  layout = Layout{};
  layout.baseOffset = nextFreeOffset;
  layout.magneticFieldOffset = layout.byteCount;
  layout.byteCount += 3 * static_cast<int>(sizeof(double));
  layout.validOffset = layout.byteCount;
  layout.byteCount += static_cast<int>(sizeof(double));
  return layout.byteCount;
}

template <class T>
void Store(PIC::Mesh::cDataCenterNode* node, int relativeOffset,
           const T& value) {
  std::memcpy(node->GetAssociatedDataBufferPointer() + layout.baseOffset +
                  relativeOffset,
              &value, sizeof(T));
}

template <class T>
T Load(PIC::Mesh::cDataCenterNode* node, int relativeOffset) {
  T value{};
  std::memcpy(&value,
              node->GetAssociatedDataBufferPointer() + layout.baseOffset +
                  relativeOffset,
              sizeof(T));
  return value;
}

}  // namespace MyApplication::CellState
```

Register the allocator once, after `PIC::Init_BeforeParser()` and before the
layout-freeze call:

```cpp
PIC::IndividualModelSampling::RequestStaticCellData.push_back(
    MyApplication::CellState::RequestStaticCellData);
```

The `memcpy` accessors are deliberate.  An associated-data byte offset is not
an application-level guarantee that a `reinterpret_cast<double*>` is aligned.
They also make the exact serialized type and component count visible at the
only access boundary.

Use a layout version or fingerprint when a restart file or external coupling
product contains the application slice.  A changed component order, type,
unit, or size is a changed data ABI even if the C++ code still compiles.

## 4. The center-node output callback contract

A complete application-defined center-node output procedure normally
registers three matching callbacks:

```cpp
PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
```

Registration order is observable output order.  If several modules register
callbacks, their variable-list and data vectors must be populated in the same
module order.  If a variable callback appends `N` Tecplot names, the
corresponding data callback must append exactly `N` numeric values for every
`DataSetNumber`.

### 4.1 Variable names and units

The variable-list callback runs while the selected writer creates its header:
on rank zero for the single-file writer, and on each rank for distributed
fragments. It must produce the same column list on every rank and contain no
MPI collective or provider update.
Append a leading comma before every new variable and put an unambiguous unit
in each physical quantity's name:

```cpp
void PrintVariableList(FILE* output, int dataSetNumber) {
  (void)dataSetNumber;
  std::fprintf(output,
      ", \"B_x_T\", \"B_y_T\", \"B_z_T\", \"background_valid\"");
}
```

The core center-node column order is:

1. active built-in sampled data;
2. active built-in derived data;
3. `PIC::Mesh::PrintVariableListCenterNode` callbacks, in registration order;
4. optional drift velocity;
5. `PIC::IndividualModelSampling::PrintVariableList` callbacks;
6. optional core SI-converted columns.

`PIC::Mesh::cDataCenterNode::PrintData()` emits numeric values in the same
order.  Never add a header in one registry and its values in another unless
that placement is intentional and tested.

### 4.2 MPI-safe data printing

For `outputMeshDataTECPLOT`, AMPS invokes the callback on the rank that owns
the node and on rank zero, which owns the output stream. For a remote node,
rank zero receives the prepared values; it must not dereference the remote
node's associated buffer. The owner sends values through the provided
`CMPI_channel` and does not write to the shared `FILE*`.

For `OutputDistributedDataTECPLOT`, the callback receives `pipe == nullptr`.
Each rank prints its owner-local nodes to its own fragment stream, including
on nonzero ranks. The surrounding writer assembles those fragments afterward.
A null channel therefore means direct/local output, not necessarily a
single-rank run. The following callback supports both protocols.

```cpp
#include <array>
#include <cmath>
#include <cstdio>

namespace MyApplication::Output {

constexpr std::size_t kValueCount = 4;

void PrintData(FILE* output, int dataSetNumber, CMPI_channel* pipe,
               int centerNodeThread,
               PIC::Mesh::cDataCenterNode* centerNode) {
  (void)dataSetNumber;
  std::array<double, kValueCount> values{};

  // With a channel, only centerNodeThread owns authoritative bytes.  A null
  // channel is direct/local output, including distributed MPI fragments.
  const bool ownsNode =
      pipe == nullptr || pipe->ThisThread == centerNodeThread;

  if (ownsNode) {
    LoadAlreadyInitializedOutputValues(centerNode, values.data());

    // Fail before serialization when a valid cell contains a non-finite
    // physical value.  Missing state should use an explicit valid flag and a
    // documented finite placeholder, never NaN as an availability protocol.
    for (double value : values) {
      if (!std::isfinite(value)) ReportInvalidOutputValueAndAbort();
    }
  }

  if (PIC::ThisThread == 0 || pipe == nullptr) {
    if (pipe != nullptr && centerNodeThread != 0) {
      pipe->recv(values.data(), static_cast<int>(values.size()),
                 centerNodeThread);
    }

    for (double value : values) std::fprintf(output, "%e ", value);
  } else {
    // The mesh writer calls this branch only on the remote owner.  The channel
    // was opened by the surrounding AMPS writer and is directed to rank zero.
    pipe->send(values.data(), static_cast<int>(values.size()));
  }
}

}  // namespace MyApplication::Output
```

Do not call `MPI_Barrier`, `MPI_Allreduce`, or any other collective from a
per-node print callback.  Ranks can be at different positions in the mesh
traversal, so a collective can deadlock.  Perform global completeness checks,
range checks, and halo exchange before entering the writer.

### 4.3 Interpolation into temporary output nodes

Tecplot FEBRICK and boundary-tetrahedron vertices do not generally coincide
with physical center nodes. AMPS constructs a temporary `cDataCenterNode`, interpolates the
contributing physical center nodes into it, and then prints the temporary
node.  The core automatically handles built-in sampled data, but it cannot
infer the meaning of application-owned bytes.  A static center field therefore
needs an `InterpolateCenterNode` callback.

Internal surfaces and inactive corridor neighbors can leave an output vertex
with no physical interpolation donors. `sourceCount == 0` is a valid empty
stencil: initialize the entire application slice to finite placeholders and
mark it unavailable, without reading either source array. Do not skip the
subsequent print row, especially in the channel-based path where every owner
must still send the expected number of columns. A negative count, null
destination, or missing source arrays for a positive count is a malformed
request. Validate those cases separately; finite placeholders do not authorize
access to an unallocated or unfrozen destination buffer.

```cpp
void Interpolate(PIC::Mesh::cDataCenterNode** sourceNodes,
                 double* coefficients, int sourceCount,
                 PIC::Mesh::cDataCenterNode* destination) {
  std::array<double, 3> interpolatedB{0.0, 0.0, 0.0};
  bool allContributingSourcesAreValid = sourceCount > 0;

  for (int source = 0; source < sourceCount; ++source) {
    const double coefficient = coefficients[source];
    const auto magneticField = LoadMagneticField(sourceNodes[source]);
    const double valid = LoadValidFlag(sourceNodes[source]);

    for (int component = 0; component < 3; ++component) {
      interpolatedB[component] += coefficient * magneticField[component];
    }

    // The stored flag is required to be exactly 0.0 or 1.0.  A zero-weight
    // stencil entry does not contribute to the destination.  Every source
    // with nonzero weight must be valid for this conservative rule to pass.
    if (coefficient != 0.0 && valid != 1.0)
      allContributingSourcesAreValid = false;
  }

  StoreMagneticField(destination, interpolatedB);

  // A validity flag is categorical state, not a continuous physical field.
  // This example marks the destination valid only when all effective stencil
  // weight came from valid sources.  Applications must document their own
  // physically appropriate discrete rule.
  StoreValidFlag(destination, allContributingSourcesAreValid ? 1.0 : 0.0);
}
```

Interpolate authoritative primitives and derive display quantities afterward.
For example, interpolate vector components or directional wave variances, then
derive magnitudes or total energy density in `PrintData`.  Directly
interpolating a nonlinear derived value can describe a different physical
quantity.

Do not linearly interpolate identifiers, ownership values, enum encodings,
boolean flags, or counters.  Define a categorical rule or recompute them from
the interpolated primitive state.

## 5. Populating and exchanging static state

Static data should be initialized in one owner-local pass after mesh block
allocation.  For each physical center node:

1. evaluate or import the authoritative model at the cell-center coordinate;
2. validate units, finite values, and physical bounds;
3. store all components of the accepted state;
4. store an explicit validity/status value;
5. optionally read back critical fields to catch offset or layout errors.

Zero-initialize padding and unused nodes, but do not let those zeros count as
successful physical initialization.  Track the number of expected and valid
physical cells and reduce those counts outside the writer.  For quantities
that must be positive, also reduce a positive-value count and finite minimum
and maximum; a zero-only test cannot detect a missing write or missing
interpolation callback.

After every owner-local node has a complete state, call:

```cpp
PIC::Mesh::mesh->ParallelBlockDataExchange();
```

Exchange a logically complete generation, not individual fields one at a
time.  A temporary output stencil can cross a block or MPI boundary, so stale
ghost bytes otherwise appear as artificial zero seams in an otherwise valid
Tecplot file.

If the background changes during a run, publish it transactionally: prepare
and validate the next generation, update every owner-local node, exchange the
whole generation, and only then allow movers or output to read it.  Output
callbacks must remain read-only.

## 6. Sampled application data

Sampled quantities are accumulated over a sampling window and are distinct
from persistent static state.  A legacy model can reserve raw sampling bytes
with:

```cpp
PIC::IndividualModelSampling::RequestSamplingData.push_back(
    RequestSamplingData);
```

`RequestSamplingData` is invoked only in builds with
`_PIC_SAMPLING_MODE_ == _PIC_MODE_ON_`.  Keep every sampling offset at an
invalid sentinel such as `-1` until the request actually runs, and make the
sampler and output callbacks handle a disabled sampling build explicitly.

The request receives an offset within one sample set.  AMPS then places that
sample set at:

- `PIC::Mesh::collectingCellSampleDataPointerOffset` while the current window
  is being accumulated; and
- `PIC::Mesh::completedCellSampleDataPointerOffset` while the previous
  completed window is being reported.

The two offsets can be identical when previous-cycle storage is disabled.  Do
not assume that they differ, and do not print from the collecting set merely
because it currently contains nonzero values.

An application using the `PIC::IndividualModelSampling` callback family should
register its corresponding variable, interpolation, and data callbacks in
matching order:

```cpp
PIC::IndividualModelSampling::SamplingProcedure.push_back(
    AccumulateOneSamplingStep);
PIC::IndividualModelSampling::PrintVariableList.push_back(PrintSampleNames);
PIC::IndividualModelSampling::InterpolateCenterNodeData.push_back(
    InterpolateCompletedSample);
PIC::IndividualModelSampling::PrintSampledData.push_back(PrintCompletedSample);
```

Keep these three lists one-to-one.  The current center-node interpolation path
indexes the model-sampling interpolation callbacks alongside the registered
model-sampling variable callbacks.

Normalization is part of the quantity's definition and must occur exactly
once.  Depending on what was accumulated, output can require division by
`PIC::LastSampleLength`, physical cell measure, total macroparticle weight, or
a combination.  Document the accumulator equation and the final units next to
the code.  An empty sampling window and a window containing zero physical
particles are different states when the distinction matters; expose a finite
count/presence field rather than using NaN.

Sampling updates executed by OpenMP workers must use the synchronization or
thread-private reduction strategy appropriate to the accumulator.  The print
and interpolation callbacks read the completed set and must not modify the
collecting set.

For new built-in-style moments, prefer the typed `PIC::Datum` infrastructure
when it represents the required weighting and normalization.  Use raw request
callbacks only when the quantity's semantics do not fit that infrastructure.

## 7. Corner-node data

Persistent corner data uses the same byte-allocation principle but a separate
registry and print callback pair:

```cpp
PIC::IndividualModelSampling::RequestStaticCellCornerData->push_back(
    RequestStaticCornerData);
PIC::Mesh::PrintVariableListCornerNode.push_back(PrintCornerVariableList);
PIC::Mesh::PrintDataCornerNode.push_back(PrintCornerData);
```

Register only after `PIC::Init_BeforeParser()` has created
`RequestStaticCellCornerData`, and before
`PIC::Mesh::initCellSamplingDataBuffer()` invokes it.

The public core API has no general corner-node interpolation callback vector.
Whole-cell output reports actual corner-node state, using the selected MPI
protocol described above. Therefore initialize and exchange every corner
value that the selected writer can visit.  Do not register a center-node
interpolator and expect it to populate corner storage.
The current boundary-tetrahedron printer does not append the corner-data
callbacks; see the explicit row-layout limitation in section 8.5 before adding
corner columns to a cut-cell product.

## 8. Dataset/species selection and output calls

### 8.1 Public writers and their geometry

| Function on `PIC::Mesh::mesh` | Product | Internal-body treatment | When to use it |
|---|---|---|---|
| `outputMeshTECPLOT(fileName)` | Geometry-only Cartesian mesh; no particle, background, or application data columns | Outputs whole cells of retained leaves; does not tessellate the sphere surface | Early mesh/refinement diagnostics, including before block allocation |
| `outputMeshDataTECPLOT(fileName, dataSetNumber)` | Data-bearing whole-cell mesh; FEBRICK in 3-D | Does not call the cut-cell tessellator; a partially solid cell remains a whole brick | Whole-cell diagnostics when curved-boundary clipping is not required |
| `OutputDistributedDataTECPLOT(fileName, true, dataSetNumber)` | Rank-local whole-cell and boundary-tetrahedron zones, optionally assembled into one file | Omits removed cells and invokes `GetCutcellTetrahedronMesh` for boundary cells | Initialized 3-D data with a resolved internal sphere or supported cut surface |

Both data-bearing writers are collective operations: every MPI rank in the
mesh communicator must call them in the same dataset/file order. The
distributed writer also requires allocated blocks and initialized associated
data. It is not a replacement for a preallocation octree diagnostic.

### 8.2 Dataset/species identity

The integer is passed unchanged to the node callbacks and block printers.
AMPS built-in particle moments and block time-step/weight accessors interpret
it as the compiled species index. When those columns are enabled, pass an
index in `[0, PIC::nTotalSpecies)` and normally emit a separate file per
species. Application callbacks must not reinterpret that same value
incompatibly. Although the whole-cell API declares a default dataset of `-1`,
that default is unsuitable when the block printers access species numerics.

For example, write a whole-cell diagnostic with an explicit species:

```cpp
PIC::Mesh::mesh->outputMeshDataTECPLOT(fileName, dataSetNumber);
```

The compiled species table is fixed when the application is configured and
built.  Runtime input must not add, remove, reorder, or relabel species.  Use
`PIC::nTotalSpecies` for the runtime table length and
`PIC::MolecularData::GetChemSymbol(...)`, `GetMass(...)`, and
`GetElectricCharge(...)` to describe the selected compiled entry.

Use `outputMeshTECPLOT(...)` only when geometry is intentionally sufficient.
It is not evidence that center state, background fields, sampling buffers,
time steps, or particle weights were initialized correctly.

### 8.3 Cut-cell initialization output

Register the internal body after AMPS base initialization but **before**
`mesh->init()` creates the tree. `Sphere::RegisterInternalSphere()` already
calls `mesh->RegisterInternalBoundary()`; do not register its returned
descriptor a second time. Analytic sphere support requires
`_INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_`; applications using
the particle-interaction callback also require the user-defined sphere mode.
These are compile-time AMPS settings, not runtime plotting options.

After layout freeze, mesh construction/decomposition, block allocation,
`InitCellMeasure()`, owner-local field population, valid species numerics and
halo exchange, write the data collectively:

```cpp
// Every MPI rank executes this block. The caller has already validated the
// complete state and exchanged ghosts; no initialization occurs in callbacks.
MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
PIC::Mesh::mesh->SetAssembleDistributedOutputFileFlag(true);
PIC::Mesh::mesh->OutputDistributedDataTECPLOT(fileName, true, speciesIndex);
MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
```

The `true` argument requests data columns; it is not the assembly flag. The
separate setter enables assembly into `fileName`. It changes a persistent mesh
setting, so a caller that needs a different assembly policy later must restore
that policy before its next output. All ranks must use the same policy.

The writer visits active owner-local blocks, classifies cells against the
registered body, prints complete cells as bricks, and appends tetrahedra for
resolved boundary cuts. It gathers zone counts and rank zero assembles the
fragments. The path must be writable on every rank and the fragments readable
by rank zero. With assembly enabled, intermediate fragments are removed after
assembly; with it disabled, the current core retains separate `.variables`,
`.header`, `.data`, `.connectivity`, and optional tetrahedron fragments. Those
pieces require assembly before use as a complete multi-zone dataset.

srcSEP3D uses this path in `WriteInitializationDataTecplotAfterBackground()`
for `sep3d-initialization-data.dat`, independently for each compiled species.
Its earlier `sep3d-initialization-mesh.dat` remains an octree diagnostic. A
visible brick across the solar interior in that early file does not prove
that the absorbing sphere is missing; inspect the later cut-cell data product.

### 8.4 Automatic sampling output versus explicit calls

The macroscopic Tecplot sampling path in `pic/pic.cpp` selects its writer
using `_PIC_OUTPUT_MODE_`, defined in the configured `picGlobal.dfn`:

| Compile-time mode | Automatic sampling writer |
|---|---|
| `_PIC_OUTPUT_MODE_SINGLE_FILE_` | `outputMeshDataTECPLOT(fileName, species)` |
| `_PIC_OUTPUT_MODE_DISTRIBUTED_FILES_` | `OutputDistributedDataTECPLOT(fileName, true, species)` |
| `_PIC_OUTPUT_MODE_OFF_` | No automatic macroscopic mesh data file |

To select cut-cell output for ordinary sampling, configure/build with:

```cpp
#define _PIC_OUTPUT_MODE_ _PIC_OUTPUT_MODE_DISTRIBUTED_FILES_
```

Ensure the generated `build/pic/picGlobal.dfn` contains the setting and rebuild
the relevant AMPS objects. This setting is read at compilation. An explicit
application call to either writer bypasses the selector: changing the macro
does not redirect a hard-coded call, and the OFF mode does not suppress one.
The word "distributed" identifies the rank-local output path; assembled
output can still be a single final file.

### 8.5 Helper functions and current core limitations

These functions implement the public paths; applications should normally
call the complete writers above rather than enter a traversal halfway through.

| Core helper | Role |
|---|---|
| `outputMeshTECPLOTinternal(...)` | Sets the header/counts and coordinates the whole-cell node/connectivity traversal; public wrappers preserve mesh modification flags around it |
| `outputMeshTECPLOT_BlockCornerNode_BlockConnectivityList(...)` | Whole-cell node and connectivity passes, including channel-based owner/rank-zero data transfer |
| `GetCutcellTetrahedronMesh(...)` | Classifies one cell and appends its boundary tetrahedra; it does not print data or modify the physical cell measure |
| `PrintTetrahedronMeshData(...)` | Writes tetrahedron vertices, interpolated center state, block numerics, and connectivity to the rank's tetrahedron fragment |
| `PrintTetrahedronMesh(...)` | Geometry-only output for an already constructed tetrahedron list; it does not construct the AMR cut-cell list |
| `cDataCenterNode::Interpolate(...)` | Fills a temporary output node through built-in and registered application interpolators |
| Node `PrintVariableList(...)` / `PrintData(...)` | Declare columns and serialize node values in matching registration order |
| `cDataBlockAMR::PrintVariableList(...)` / `PrintData(...)` | Append species-local time step and particle-weight columns |

The current core implementation has three relevant restrictions:

- Use `PrintMeshData=true` for `OutputDistributedDataTECPLOT`. Its present
  whole-cell/tetrahedron row writers still emit data columns when the flag is
  false, while the header omits them; that is not a valid geometry-only path.
- `PrintTetrahedronMeshData` currently prints center and block fields but does
  not serialize the corner-node callback columns included in the common
  header. If an application registers additional corner columns, correct that
  tetrahedron row layout or omit those columns from the product, and verify
  header/row agreement for **both** zone types before using the result.
- The analytic-sphere tests in `meshAMR/cut_cell.hpp` currently assume the
  sphere is centered at the coordinate origin. A nonzero body origin needs a
  corresponding core correction before this path can be trusted. Placing the
  Sun near a box corner by changing box bounds, while retaining the physical
  Sun at `(0,0,0)`, satisfies the existing assumption.

Boundary tetrahedra approximate curved geometry with straight faces. Resolve
the body with the surface mesh resolution and assess refinement convergence;
zeroing interior measures is a physical-volume correction, not a substitute
for the selected writer's geometry classification.

## 9. Registration and output invariants

Every application should enforce these invariants:

- Register each allocator and callback exactly once.  Duplicate registration
  duplicates storage and columns and can desynchronize headers from rows.
- Keep the exact number and order of names and values identical.
- Include units in names and define the frame/basis for vectors.
- Store finite placeholders plus an explicit validity/presence/status column;
  do not encode missing data as NaN or infinity.
- Do not mutate providers, the runtime clock, sampling windows, particles, or
  mesh state from an output callback.
- Do not perform per-node collectives.
- Do not read a remote node buffer on rank zero.
- Do not cast unverified byte offsets to aligned C++ types.
- Do not append storage after layout freeze or block allocation.
- Interpolate every application-owned center field needed at output vertices.
- Exchange ghost state after a complete update and before output.
- Treat a callback ABI or associated-data layout change as a restart/output
  compatibility change.

## 10. Validation checklist for a new application

Before accepting a new output procedure:

1. **Layout test:** verify offsets are assigned exactly once, remain within the
   reserved slice, and are frozen before block allocation.
2. **Header/row test:** parse a produced Tecplot file and prove that every data
   row has exactly the number of values declared in `VARIABLES`.
3. **Finite-value test:** reject `nan`, `inf`, and `-inf`; separately verify the
   meaning of every validity/presence column.
4. **Nonzero fixture:** initialize a field with varying, physically meaningful
   nonzero values.  A zero-only fixture cannot detect missing initialization or
   interpolation.
5. **Interpolation test:** test constant, linear-gradient, and invalid-stencil
   cases independently of the full application.
6. **Single-rank run:** inspect units, signs, ranges, dataset selection, local
   time steps, and particle weights.
7. **Multi-rank run:** compare against the single-rank product and inspect MPI
   and AMR interfaces for zero or discontinuous seams.
8. **Empty-particle run:** confirm zero-particle cells produce finite sampled
   values and an explicit absence/count signal.
9. **Sampling-window test:** confirm the writer reads the completed buffer and
   applies time/volume/weight normalization once.
10. **Restart test:** when the field is restart-relevant, verify byte-layout
    compatibility checks and exact restoration before output.
11. **Internal-boundary output test:** verify the body is registered before
    mesh construction and the selected writer produces boundary-tetrahedron
    zones for a fixture with resolved cuts. Verify fully interior cells are
    absent and compare the represented excluded volume under refinement.
12. **Both MPI print protocols:** exercise channel-based owner/rank-zero rows
    and null-channel owner-local fragments, including zero-donor vertices.
    Parse every brick and tetrahedron zone separately for matching column
    counts, species identity, and finite data.

For an initialization-only mode, success should mean more than “the mesh was
built.”  Require complete owner-cell initialization, halo exchange, finite
background state, valid species time steps and weights, and successful
generation of at least one data-bearing file.

## 11. Core implementation examples

The following files illustrate parts of the contract.  They are examples of
the core ABI, not templates whose physical assumptions should be copied:

| File | Relevant behavior |
|---|---|
| `pic/pic_mesh.cpp` | layout freeze, completed/collecting sample offsets, callback dispatch, and built-in interpolation |
| `pic/pic_datafile.cpp` | registration and interpolation of persistent background variables |
| `pic/pic_background_atmosphere.cpp` | a raw sampled-data request plus MPI-safe printing and interpolation |
| `pic/pic_stopping_power.cpp` | another model-specific sampled output path |
| `models/dust/Dust.cpp` | application/model registration of the center-node callback trio |
| `models/exosphere/Exosphere.cpp` | sampled storage and output registration |
| `meshAMR/meshAMRgeneric.h` | construction of temporary output nodes and distributed Tecplot traversal |
| `meshAMR/cut_cell.hpp` | cut-cell classification/tessellation and tetrahedron vertex/data printing |
| `pic/pic.cpp` | automatic sampling writer selection through `_PIC_OUTPUT_MODE_` |
| `general/mpichannel.h` | `CMPI_channel` send/receive semantics used by output callbacks |

Read the current declarations in `pic/pic.h` when adding a new application;
the signatures there are authoritative.  An application README should then
document its own storage layout, physical equations, units, validity policy,
initialization boundary, normalization, and acceptance tests in addition to
referencing this core lifecycle contract.

## 12. Runtime model fields and AMPS publication

`srcSEP3D/background/BackgroundProvider` (declared in `bg_provider.h`) separates
model evaluation from AMPS storage. A source prepares one frozen epoch and
returns SI plasma/IMF samples, derivatives, coverage and a configuration
fingerprint. `srcSEP3D/runtime/background_factory.h` registers named future
models; `srcSEP3D/main_lib.cpp` owns the common publication and MPI boundary.

Evaluate into candidate snapshots, validate ALL owner fields and associated
turbulence, and collectively agree on identity before changing live buffers.
Write application storage and allocated native DATAFILE current/next slots
from the same sample. Exchange associated data and complete enabled derived
guiding-center fields before committing the generation and resuming movers.
Both native slots deliberately contain one epoch: cadence-frozen fields must
not interpolate an updated B against stale U. Empty owner ranks join the same
collective while global physical coverage remains mandatory.

SWCME now uses this path for evolving U, vector B, plasma state and gradients.
`E=-U×B` and allocated current/electron-pressure fields are updated with it.
The internal solar sphere and cut-cell Tecplot output contract described above
remain independent of the background source. For the adapter implementation,
model limitations and adding a provider, see
[`SWCME_MESH_BACKGROUND_CHANGES.md`](../SWCME_MESH_BACKGROUND_CHANGES.md).

The [background module](../srcSEP3D/background/README.md) specifies sample
units, vector-Jacobian conventions and transactional evaluation. The
[runtime module](../srcSEP3D/runtime/README.md) specifies factory registration
and the collective publication sequence. The adapter validates byte offsets
and component counts before writing, including optional native fields; adding
a stored field requires updating both the mapping and that layout check.
Model code must stay independent of PIC buffers and MPI ownership.

Publication evidence reads actual allocated bytes without a second field fill
or halo exchange. Owner buffers are checked for every physical cell; remote
coverage is one physical representative per received active block. This is a
direct-copy/communication check, separate from physical convergence or a
spacecraft fit. The absorbing solar sphere and distributed cut-cell writer do
not change when a different background provider is selected.
