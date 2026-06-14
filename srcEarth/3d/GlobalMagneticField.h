#ifndef _SRC_EARTH_3D_GLOBALMAGNETICFIELD_H_
#define _SRC_EARTH_3D_GLOBALMAGNETICFIELD_H_

//======================================================================================
// GlobalMagneticField.h
//======================================================================================
//
// Shared support for backward cutoff-rigidity calculations that are executed with a
// distributed AMPS mesh.
//
// In the old standalone 3-D cutoff path each MPI rank called
//
//   PIC::InitMPI(/*independentDomainMode=*/true)
//
// and therefore built a complete private copy of the AMR domain.  That made field
// interpolation simple, but it bypassed the normal AMPS parallel-domain layout and
// duplicated all mesh/block memory on every rank from the start of the run.
//
// The SWMF-coupled cutoff path later implemented a better approach: keep the normal
// AMPS distributed mesh during initialization, then immediately before cutoff tracing
// materialize a read-only global magnetic-field snapshot on every MPI rank.  The
// generic helper declared here is that data-management operation factored out so it
// can be reused by both
//
//   * standalone 3-D cutoff: source/destination is the DATAFILE magnetic-field slot;
//   * SWMF-coupled cutoff:   source/destination is the SWMF coupler magnetic-field slot.
//
// The helper does not know how the magnetic field was produced.  It only needs the
// byte offset, relative to cDataCenterNode::GetAssociatedDataBufferPointer(), where
// the three double precision B components are stored.
//======================================================================================

namespace Earth {
namespace Mode3D {
namespace GlobalMagneticField {

// Compact diagnostic summary returned by MaterializeCellCenteredMagneticFieldForCutoff().
// The values are useful for rank-0 log messages and for debugging incomplete gathers.
struct MaterializationStats {
  long int usedLeafBlocks;
  long int ownerInteriorCells;
  long int expectedInteriorCells;
  long int newlyAllocatedBlocksMin;
  long int newlyAllocatedBlocksMax;
  long int missingGhostCellsMin;
  long int missingGhostCellsMax;

  MaterializationStats() :
    usedLeafBlocks(0), ownerInteriorCells(0), expectedInteriorCells(0),
    newlyAllocatedBlocksMin(0), newlyAllocatedBlocksMax(0),
    missingGhostCellsMin(0), missingGhostCellsMax(0) {}
};

// Return the absolute DATAFILE magnetic-field offset in a cell-associated data buffer.
// This is the offset corresponding to the layout used by Mode3D::InitMeshFields().
long int DataFileMagneticFieldDataOffset();

// Build a replicated, read-only cell-centered magnetic-field snapshot on every MPI rank.
//
// Parameters:
//   diagnosticTag
//     Short text used in log/error messages, e.g. "Mode3D" or "Mode3DForwardSWMF".
//
//   magneticFieldDataOffset
//     Absolute byte offset from cell->GetAssociatedDataBufferPointer() to Bx.  By and
//     Bz must follow immediately as the next two doubles.  This intentionally supports
//     both DATAFILE and SWMF-coupler storage without branching inside the helper.
//
//   verbose
//     If true, rank 0 prints a one-line summary.
//
// Operation:
//   1. assign deterministic global Temp_ID values to all used AMR leaves;
//   2. allocate missing leaf blocks locally so every rank has every used leaf;
//   3. pack only owner-rank interior-cell B values into a dense array;
//   4. MPI_Allreduce the dense B array and an integer presence mask;
//   5. fill every allocated block, including ghost cells, from the global cache.
//
// After this call, the usual AMPS cell-centered interpolation routines can be used
// unchanged during cutoff trajectory tracing even though the run was initialized with
// the normal distributed MPI domain decomposition.
MaterializationStats MaterializeCellCenteredMagneticFieldForCutoff(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    bool verbose=true);

// Debug/utility path: replace B in all currently used leaves with a callback-generated
// field.  This is used by the SWMF bridge's analytic-dipole debug override.  The same
// block-allocation logic is available here so the override can operate on a fully
// replicated cutoff mesh when requested.
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
