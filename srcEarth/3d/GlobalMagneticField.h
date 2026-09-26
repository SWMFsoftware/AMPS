#ifndef _SRC_EARTH_3D_GLOBALMAGNETICFIELD_H_
#define _SRC_EARTH_3D_GLOBALMAGNETICFIELD_H_

//======================================================================================
// GlobalMagneticField.h
//======================================================================================
//
// Compact, decomposition-independent cell-centered field storage for Mode3D backward
// trajectory calculations.
//
// AMPS keeps the AMR tree topology on every MPI rank, but allocates cDataBlockAMR and
// cDataCenterNode objects only for the rank-local domain and a limited neighbor layer.
// A cutoff-rigidity trajectory can cross the complete magnetosphere, so field lookup
// cannot depend on the local block allocation.  Replicating every AMPS block solves
// that problem but also replicates every cell-associated state vector and all ghost
// storage, which is prohibitively expensive for realistic SWMF meshes.
//
// This module instead replicates compact B/E arrays and, for a Step-9 SWMF snapshot,
// the bulk velocity used to define the optional ideal-MHD electric field. Every leaf is
// assigned a deterministic dense index in node->Temp_ID.  A physical interior cell is
// addressed by
//
//   globalCell = node->Temp_ID * (Nx*Ny*Nz)
//              + i + Nx*(j + Ny*k).
//
// Field interpolation uses PIC::InterpolationRoutines::cRowStencil.  A row-stencil
// element contains the owning tree node, interior (i,j,k), and interpolation weight,
// and therefore remains valid even when the corresponding node->block is NULL on the
// current rank.
//======================================================================================

#include "pic.h"
#include "../util/FieldProvider.h"

#include <memory>
#include <string>

namespace Earth {
namespace Mode3D {
namespace GlobalMagneticField {

// Public node alias used by the row-stencil field-evaluation interface.
typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cAMRNode;

// Diagnostic summary returned after a compact global field snapshot is assembled.
struct MaterializationStats {
  long int usedLeafBlocks;
  long int ownerInteriorCells;
  long int expectedInteriorCells;
  long int missingInteriorCells;
  long int duplicateInteriorCells;
  long int magneticFieldBytes;
  long int electricFieldBytes;
  long int plasmaVelocityBytes;
  bool electricFieldReadFromBuffer;
  bool electricFieldDerivedFromVelocity;
  bool plasmaVelocityAvailable;
  bool ownerCellParityValidated;
  std::string snapshotId;
  std::string contentFingerprint;
  std::string meshRevision;

  MaterializationStats() :
    usedLeafBlocks(0), ownerInteriorCells(0), expectedInteriorCells(0),
    missingInteriorCells(0), duplicateInteriorCells(0),
    magneticFieldBytes(0), electricFieldBytes(0), plasmaVelocityBytes(0),
    electricFieldReadFromBuffer(false),
    electricFieldDerivedFromVelocity(false),plasmaVelocityAvailable(false),
    ownerCellParityValidated(false) {}
};

// Return absolute DATAFILE offsets relative to
// cDataCenterNode::GetAssociatedDataBufferPointer().
long int DataFileMagneticFieldDataOffset();
long int DataFileElectricFieldDataOffset();

// Assemble compact global B and E arrays on every MPI rank.
//
// magneticFieldDataOffset:
//   Required absolute offset of three consecutive double-precision B components.
//
// electricFieldDataOffset:
//   Optional absolute offset of three consecutive E components.  Pass -1 when E is
//   not stored directly.
//
// plasmaVelocityDataOffset:
//   Optional absolute offset of three consecutive plasma-velocity components. The
//   values are always retained for SWMF export/replay. When metadata explicitly marks
//   E available and electricFieldDataOffset is negative, E is derived cell-by-cell as
//   E=-v x B. Released Phase-1 SWMF products leave E unavailable (magnetic-only);
//   derived E is an explicitly experimental option until its later validation gates.
//
// When neither E source is available, a valid zero electric field is stored.  The
// routine resets node->Temp_ID over the complete tree before assigning new dense IDs;
// callers must not overwrite Temp_ID while the snapshot is in use.
MaterializationStats AssembleCellCenteredFieldsForCutoff(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    long int electricFieldDataOffset=-1,
    long int plasmaVelocityDataOffset=-1,
    bool verbose=true);

// Step-3 overload: publish validated provenance with the same generation as the
// compact arrays.  The legacy overload above remains source-compatible, while all
// production standalone/SWMF callers provide physical metadata through this form.
MaterializationStats AssembleCellCenteredFieldsForCutoff(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    long int electricFieldDataOffset,
    long int plasmaVelocityDataOffset,
    const Earth::Field::SnapshotMetadata& metadata,
    bool verbose=true,
    double sourceSimulationTime_s=0.0);

// Backward-compatible B-only entry point.  Existing call sites can continue using this
// function; it now creates compact arrays and a zero E array instead of allocating and
// populating nonlocal AMR blocks.
MaterializationStats MaterializeCellCenteredMagneticFieldForCutoff(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    bool verbose=true);

// Interpolate a compact global field with an AMPS decomposition-independent row
// stencil.  Return false only when no snapshot is ready, the point is outside the used
// AMR tree, or no interpolation row can be constructed.  A malformed Temp_ID or a
// missing row cell is a fatal consistency error because silently dropping that cell
// would change the physical interpolation result.
bool InterpolateMagneticField(const double* x,cAMRNode* node,double* B);
bool InterpolateElectricField(const double* x,cAMRNode* node,double* E);

// Direct cell access is useful for diagnostics and unit tests.  The indices must be
// interior indices of the supplied owning leaf node.
bool GetCellCenteredMagneticField(cAMRNode* node,int i,int j,int k,double* B);
bool GetCellCenteredElectricField(cAMRNode* node,int i,int j,int k,double* E);
bool GetCellCenteredPlasmaVelocity(cAMRNode* node,int i,int j,int k,double* velocity);

bool GlobalFieldsReady();
long int GlobalCellCount();
void ClearGlobalFields();

// Read-only metadata and provider views of the published compact generation.  A
// snapshot view captures its generation; after reassembly or ClearGlobalFields(), its
// Sample() method returns STALE_EPOCH instead of reading replacement arrays.
const Earth::Field::SnapshotMetadata& CurrentSnapshotMetadata();
std::shared_ptr<const Earth::Field::IFieldSnapshot> CurrentSnapshot();
std::shared_ptr<Earth::Field::IFieldProvider> CurrentFieldProvider();

// Freeze/release one published generation around a product batch. Mutating operations
// fail while frozen, so a late field update cannot replace arrays underneath active
// trajectories. In SWMF execution the coupler scheduler naturally queues the next
// receive until the callback returns; these guards make that serialization explicit.
void BeginFrozenFieldBatch(const std::string& expectedSnapshotId);
void EndFrozenFieldBatch(const std::string& expectedSnapshotId);
bool FrozenFieldBatchActive();

// Roadmap Step 9 live/export/replay bridge. Export reads only the frozen compact
// arrays, never mutable coupler storage. Import validates schema, SI/GSM declarations,
// epoch, topology, mesh revision, cell centres, content identity, and the requested
// magnetic-only/experimental-E mode before publishing a normal field generation.
std::string ExportCurrentSWMFSnapshot(const std::string& fileName);
MaterializationStats ImportSWMFSnapshot(
    const std::string& fileName,
    const std::string& expectedEpochUTC,
    bool enableExperimentalDerivedElectric,
    bool verbose=true);
bool PlasmaVelocityAvailable();
std::string CurrentContentFingerprint();
std::string CurrentMeshRevision();

// Replace the compact global magnetic field by values generated from a coordinate
// callback.  No AMPS block allocation is performed.  The electric field is reset to
// zero because it would generally be inconsistent with the replacement B field.
long int RedefineGlobalMagneticField(
    const char* diagnosticTag,
    void (*fieldCallback)(double*,double*),
    bool verbose=true);

// Compatibility wrapper for the previous replicated-block debug API.  The data offset
// and allocateMissingBlocks arguments are retained so old callers compile, but the
// implementation intentionally updates only the compact global array and never
// allocates nonlocal blocks.
long int RedefineAllAllocatedMagneticField(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    void (*fieldCallback)(double*,double*),
    bool allocateMissingBlocks=true,
    bool verbose=true);

} // namespace GlobalMagneticField
} // namespace Mode3D
} // namespace Earth

#endif // _SRC_EARTH_3D_GLOBALMAGNETICFIELD_H_
