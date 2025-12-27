

#include "../pic.h"

/**
 * @file GetStencil_Legacy_Wrapper.cpp
 *
 * @brief Legacy-to-new RHS support wrapper for ECSIM stencil construction.
 *
 * This wrapper allows code that expects the *new* stencil interface (semantic RHS support vectors)
 * to still reuse the *legacy* GetStencil() implementation that fills fixed-size RHS support tables:
 *   - RhsSupportTable_CornerNodes + RhsSupportLength_CornerNodes
 *   - RhsSupportTable_CenterNodes + RhsSupportLength_CenterNodes
 *
 * The wrapper:
 *   1) Allocates local (thread-local) legacy RHS support tables and local lengths.
 *   2) Calls the legacy GetStencil() to populate the matrix row, rhs constant term, and legacy support tables.
 *   3) Copies the local legacy tables into the semantic vectors support_corner_vector and support_center_vector.
 *   4) Optionally mirrors the local legacy results back into the caller-provided legacy arrays/lengths
 *      (useful while the codebase is transitioning).
 *
 * Notes:
 *   - The thread_local scratch tables avoid large stack allocations and are safe for OpenMP/MPI ranks.
 *   - LEGACY_MAX_* must be >= the maximum number of legacy RHS support entries emitted by GetStencil().
 *   - Error handling uses AMPS-style exit(__LINE__,__FILE__,msg).
 */

#include <vector>
#include <array>
#include <cstdio>   // std::snprintf

using std::vector;

//------------------------------------------------------------------------------
// Forward declaration of the LEGACY GetStencil() implementation.
// This function must exist in the build (e.g., in pic_field_solver_ecsim.cpp).
//
// Wrapper requested: uses local legacy tables, calls legacy GetStencil(), then
// populates semantic vectors from the local tables.
//------------------------------------------------------------------------------

namespace PIC {
namespace FieldSolver {
namespace Electromagnetic {
namespace ECSIM {
void GetStencil_Legacy_Wrapper(
    int i,int j,int k,int iVar,
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,
    int& NonZeroElementsFound,
    double& rhs,
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable* RhsSupportTable_CornerNodes,
    int& RhsSupportLength_CornerNodes,
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable* RhsSupportTable_CenterNodes,
    int& RhsSupportLength_CenterNodes,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,
    vector<cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable>& support_corner_vector,
    vector<cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable>& support_center_vector
) {
  // Local type aliases for readability
  using LS       = cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>;
  using RhsEntry = LS::cRhsSupportTable;

  // Maximum legacy RHS support sizes.
  // These must be at least as large as the maximum produced by the legacy GetStencil().
  // If you know tighter bounds, you can reduce them.
  constexpr int LEGACY_MAX_CORNER = 4096;
  constexpr int LEGACY_MAX_CENTER = 4096;

  // Thread-local scratch tables prevent large stack allocations and allow OpenMP safety.
  static thread_local std::array<RhsEntry, LEGACY_MAX_CORNER> localRhsSupportTable_CornerNodes;
  static thread_local std::array<RhsEntry, LEGACY_MAX_CENTER> localRhsSupportTable_CenterNodes;

  // Local lengths returned by the legacy GetStencil()
  int localCornerLen = 0;
  int localCenterLen = 0;

  // ---------------------------------------------------------------------------
  // 1) Call the legacy GetStencil() using the local scratch tables.
  //    This fills:
  //      - MatrixRowNonZeroElementTable + NonZeroElementsFound
  //      - rhs (constant part)
  //      - localRhsSupportTable_* + local*Len
  // ---------------------------------------------------------------------------
  GetStencil_Legacy(
      i,j,k,iVar,
      MatrixRowNonZeroElementTable, NonZeroElementsFound, rhs,
      localRhsSupportTable_CornerNodes.data(), localCornerLen,
      localRhsSupportTable_CenterNodes.data(), localCenterLen,
      node
  );

  // ---------------------------------------------------------------------------
  // 2) Validate the returned lengths to ensure we did not overflow our scratch.
  // ---------------------------------------------------------------------------
  if (localCornerLen < 0 || localCornerLen > LEGACY_MAX_CORNER) {
    char msg[256];
    std::snprintf(msg, sizeof(msg),
      "GetStencil_Legacy_Wrapper: localCornerLen=%d exceeds LEGACY_MAX_CORNER=%d",
      localCornerLen, LEGACY_MAX_CORNER);
    exit(__LINE__, __FILE__, msg);
  }

  if (localCenterLen < 0 || localCenterLen > LEGACY_MAX_CENTER) {
    char msg[256];
    std::snprintf(msg, sizeof(msg),
      "GetStencil_Legacy_Wrapper: localCenterLen=%d exceeds LEGACY_MAX_CENTER=%d",
      localCenterLen, LEGACY_MAX_CENTER);
    exit(__LINE__, __FILE__, msg);
  }

  // ---------------------------------------------------------------------------
  // 3) Populate the semantic RHS support vectors from the local legacy tables.
  //    The new UpdateRhs() path should consume these vectors.
  // ---------------------------------------------------------------------------
  support_corner_vector.assign(
      localRhsSupportTable_CornerNodes.data(),
      localRhsSupportTable_CornerNodes.data() + localCornerLen
  );

  support_center_vector.assign(
      localRhsSupportTable_CenterNodes.data(),
      localRhsSupportTable_CenterNodes.data() + localCenterLen
  );

  // ---------------------------------------------------------------------------
  // 4) Optional: mirror results back into the caller-provided legacy arrays/lengths.
  //    This is useful if any remaining code still inspects the legacy members.
  // ---------------------------------------------------------------------------
  RhsSupportLength_CornerNodes = localCornerLen;
  RhsSupportLength_CenterNodes = localCenterLen;

  if (RhsSupportTable_CornerNodes != nullptr) {
    for (int n = 0; n < localCornerLen; ++n) {
      RhsSupportTable_CornerNodes[n] = localRhsSupportTable_CornerNodes[n];
    }
  }

  if (RhsSupportTable_CenterNodes != nullptr) {
    for (int n = 0; n < localCenterLen; ++n) {
      RhsSupportTable_CenterNodes[n] = localRhsSupportTable_CenterNodes[n];
    }
  }
}
}}}}
