
#include "indexer_4d.h"

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

/**
 * Assemble one matrix row and RHS support for the ECSIM semi-implicit E-field solve.
 *
 * This routine is invoked by the corner-node linear system builder (cLinearSystemCornerNode)
 * for a single unknown: the electric-field component E_{iVar} at the *corner node* (i,j,k).
 *
 * ECSIM context (semi-implicit Maxwell + particle response)
 * ---------------------------------------------------------
 * AMPS stores:
 *   • E on corner nodes (unknowns of this linear system)
 *   • B on center nodes (appears in RHS through curl(B))
 *   • J and the particle "mass matrix" response on corner nodes (from UpdateJMassMatrix()).
 *
 * The ECSIM update solves for E^{n+θ} (θ∈[0,1], typically 0.5) from a discretized Ampère–Maxwell law
 * where the particle response is treated implicitly through the mass matrix. In operator form:
 *
 *   A(E^{n+θ}) = b(E^n, B^n, J^n, ...)
 *
 * where A contains:
 *   • the identity term
 *   • discrete curl-curl / grad-div pieces (coming from ∇×B and the semi-implicit Faraday coupling)
 *   • + (4π θ Δt) * (mass-matrix operator) contributions
 *
 * Output contract
 * ---------------
 * This function fills:
 *   • MatrixRowNonZeroElementTable / NonZeroElementsFound : sparse stencil structure for the row of A
 *   • rhs                                              : the *constant* part of b
 *   • RhsSupportTable_* / RhsSupportLength_*            : legacy RHS representation (kept for compatibility)
 *   • support_corner_vector / support_center_vector     : semantic RHS representation used by the new UpdateRhs()
 *
 * Semantic RHS representation (new path)
 * -------------------------------------
 * The RHS is represented as a dot-product over a small set of sampled quantities:
 *
 *   rhs_total = rhs + Σ_s ( coeff_s * sample_s )
 *
 * where each sample_s is a *typed* (semantic) reference to:
 *   • Corner(E component) or Corner(J component)
 *   • Center(B component)
 *
 * UpdateRhs() later evaluates these semantic entries using the correct time-level offsets
 * (CurrentEOffset / CurrentBOffset) without relying on hard-coded array sizes like _PIC_STENCIL_NUMBER_.
 *
 * NOTE: For compatibility with older code paths and optional self-tests, the function signature
 * still includes the legacy RHS arrays; the semantic vectors are the authoritative output for
 * the refactored RHS path.
 */

void PIC::FieldSolver::Electromagnetic::ECSIM::GetStencil_import_cStencil(int i,int j,int k,int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,
			     cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable* RhsSupportTable_CornerNodes,int& RhsSupportLength_CornerNodes,
			     cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable* RhsSupportTable_CenterNodes,int& RhsSupportLength_CenterNodes, 
			     cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,                                                                           vector<cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable>& support_corner_vector,
              vector<cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable>& support_center_vector) {
  
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  // ---------------------------------------------------------------------------
  // (1) Initialize semantic RHS support containers.
  //     In the ECSIM linear solve, the RHS consists of:
  //       • a constant part (stored in `rhs`)
  //       • plus linear combinations of sampled E, B, J values from neighboring nodes.
  //     The new path stores those samples in the semantic vectors below.
  // ---------------------------------------------------------------------------

  support_corner_vector.clear();
  support_center_vector.clear();

//  GetStencil(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,RhsSupportTable_CornerNodes,RhsSupportLength_CornerNodes,RhsSupportTable_CenterNodes,RhsSupportLength_CenterNodes,node);
//  return;
  
  // No.0-No.26  stencil Ex
  // No.27-No.53 stencil Ey
  // No.54-No.80 stencil Ez
  
  // i+indexadd[ii](ii:0,1,2), j+indexAdd[jj](jj:0,1,2), k+indexAdd[kk](kk:0,1,2)
  // index number: ii+3*jj+9*kk 
  // No.0: i,j,k No.1 i-1,j,k No.2 i+1,j,k
  // No.3: i,j-1,k No.4:i-1,j-1,k No.5 i+1,j-1,k 
 
  // double cLighgt, dt;
  double x[3];
  //power of 3 array created
  //for some pgi compilers that cannot convert result of pow() from double to int correctly
  int pow3_arr[3]={1,3,9};
  int index[3] = {i,j,k};
  int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
  double dx[3],coeff[3],coeffSqr[3],coeff4[3]; 

  // ---------------------------------------------------------------------------
  // (4) Stencil indexing (do NOT assume _PIC_STENCIL_NUMBER_).
  //     The legacy implementation relied on fixed-size arrays sized by _PIC_STENCIL_NUMBER_.
  //     Here we use cIndexer4D to map (di,dj,dk,component) -> compact integer slot.
  //     Only slots that are actually referenced are created/used in the semantic vectors.
  // ---------------------------------------------------------------------------
  //init indexers 
  //static 
	  cIndexer4D indexer[3]; //iVar=0...2 

//  static cIndexer4D indexer_rhs;

  if (indexer[0].isInitialized()==false) {
    //init indexer

    int min_limit[4]={-5,-5,-5,0},max_limit[4]={5,5,5,2};

    for (int i=0;i<3;i++) indexer[i].init(min_limit,max_limit);

 //   indexer_rhs.init(min_limit,max_limit);
  }

for (int i=0;i<3;i++) indexer[i].reset();

  // Canonical mapping from physical offsets (di,dj,dk) and unknown component
  // (iVarIndex: 0->Ex, 1->Ey, 2->Ez) into the compact element slot *for this row iVar*.
  // This keeps iElement <-> (di,dj,dk,iVarIndex) consistent across all terms.
  auto map_iElement = [&](int di, int dj, int dk, int iVarIndex)->int {
    return indexer[iVar].get_idx(di, dj, dk, iVarIndex);
  };

  //------------------------------------------------------------------------------
  // ensure_curlcurl_stencil_initialized
  //------------------------------------------------------------------------------
  // PURPOSE
  //   One-shot initialization wrapper for the 4th-order curl–curl(E) stencil.
  //   On the *first* call, it invokes the provided initializer
  //
  //     PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::FourthOrder::
  //     InitCurlCurlEStencils(cCurlCurlEStencil* CurlCurl, double dx, double dy, double dz)
  //
  //   with the given spacings (dx,dy,dz) and caches the result in a static
  //   object. On subsequent calls, it *does not* re-initialize; it simply copies
  //   the cached stencil into the caller-provided `out`.
  //
  // NOTES
  //   - This is intentionally minimal: a single static boolean gate. If your
  //     grid spacings can change at runtime, extend this to also cache (dx,dy,dz)
  //     and rebuild when they differ. As written, it initializes once and reuses.
  //   - Copying `s_stencil` into `out` relies on `cStencil`/`cCurlCurlEStencil`
  //     being trivially copyable or having proper copy semantics.
  //   - Not thread-safe. If used from multiple threads, guard the first-call
  //     section with a mutex.
  //
  // Example usage:
  // cCurlCurlEStencil *curlcurl;
  // curlcurl=get_curlcurl_stencil(dx, dy, dz);
  //------------------------------------------------------------------------------
  auto get_grad_div_stencil =
    [&]() -> PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cGradDivEStencil* {
      using PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::FourthOrder::InitGradDivEBStencils;

      static bool               s_inited   = false;
      static PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cGradDivEStencil s_stencil[3];          // cached once

      if (!s_inited) {
        InitGradDivEBStencils(s_stencil, 1.0,1.0,1.0);  // build exactly once
        s_inited = true;
      }

      // hand back the cached stencil to the caller
      return s_stencil;
    };

  auto get_laplacian_stencil =
    [&]() -> PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cLaplacianStencil* {
      using PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::FourthOrder::InitLaplacianStencil;

      static bool               s_inited   = false;
      static PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cLaplacianStencil  s_stencil;          // cached once

      if (!s_inited) {
        InitLaplacianStencil(&s_stencil, 1.0,1.0,1.0);  // build exactly once
        s_inited = true;
      }

      // hand back the cached stencil to the caller
      return &s_stencil;
    };

  //------------------------------------------------------------------------------
  // init_metrics_and_time_factors
  //------------------------------------------------------------------------------
  // PURPOSE
  //   Compute per-dimension grid metrics and ECSIM time-centering factors used
  //   throughout stencil assembly:
  //     • dx[d]      : physical cell size along axis d (after unit conversion)
  //     • x[d]       : coordinate of the current corner index along axis d
  //     • coeff[d]   : θ c Δt / dx[d]           (single-derivative scale factor)
  //     • coeffSqr[d]: (θ c Δt / dx[d])^2       (second-derivative scale factor)
  //     • coeff4[d]  : (θ c Δt / dx[d]) / 4     (used for face-avg in curl(B))
  //
  // INPUTS (must exist in the enclosing scope):
  //   node   → current AMR tree node (provides xmin/xmax bounds)
  //   nCell  → {NX, NY, NZ} number of cells per block in each dimension
  //   index  → {i, j, k} integer corner indices within the block
  //   cDt    → c * Δt  (speed of light times timestep)
  //   theta  → ECSIM time-centering parameter
  //   length_conv → optional unit conversion factor for lengths
  //
  // OUTPUTS (filled by this lambda; exist in enclosing scope):
  //   dx[3], x[3], coeff[3], coeffSqr[3], coeff4[3]
  //
  // NOTES
  //   • x[d] is evaluated at the *corner* location corresponding to (i,j,k).
  //   • coeff and coeffSqr set the correct scaling for first/second derivatives
  //     in discrete operators (e.g., grad, div, Laplacian, grad–div).
  //   • coeff4 is exactly the 1/4 factor used when averaging 4 center samples
  //     on a face to approximate curl(B) with a compact stencil.
  //------------------------------------------------------------------------------
  auto init_metrics_and_time_factors = [&]() {
    for (int iDim = 0; iDim < 3; ++iDim) {
      // Raw cell size in node coordinates
      dx[iDim] = (node->xmax[iDim] - node->xmin[iDim]) / nCell[iDim];

      // Apply unit conversion (if any)
      dx[iDim] *= length_conv;

      // Corner coordinate along this axis at index[iDim]
      x[iDim] = node->xmin[iDim]
            + index[iDim] * (node->xmax[iDim] - node->xmin[iDim]) / nCell[iDim];

      // Time-centering / metric factors:
      //   coeff    = θ c Δt / Δx
      //   coeffSqr = (θ c Δt / Δx)^2
      //   coeff4   = (θ c Δt / Δx) / 4 (face-average factor for curl(B))
      coeff[iDim]    = (cDt / dx[iDim]) * theta;
      coeffSqr[iDim] = coeff[iDim] * coeff[iDim];
      coeff4[iDim]   = 0.25 * coeff[iDim];
    }
  };

  init_metrics_and_time_factors();


  //------------------------------------------------------------------------------
  // boundary_short_circuit
  //------------------------------------------------------------------------------
  // PURPOSE
  //   Handle non-periodic physical boundaries for the ΔE solve at corner (i,j,k).
  //   If (i,j,k) is a boundary corner, we enforce the Dirichlet condition ΔE=0 by
  //   emitting a *single diagonal row* and skipping all other assembly. This
  //   makes the row:
  //        A_00 := 1  (stored via PARAMETER[0] = 1.0; numeric value filled later)
  //        rhs  := 0
  //   and sets lengths of RHS support tables to zero.
  //
  // WHAT IT WRITES (index 0 only):
  //   MatrixRowNonZeroElementTable[0].{i,j,k}            ← (i,j,k)
  //   MatrixRowNonZeroElementTable[0].iVar               ← iVar (row/equation component)
  //   MatrixRowNonZeroElementTable[0].MatrixElementValue ← 0.0      (value is filled later)
  //   MatrixRowNonZeroElementTable[0].MatrixElementParameterTable[0] ← 1.0 (unit diagonal)
  //   MatrixRowNonZeroElementTable[0].MatrixElementParameterTableLength ← 1
  //   MatrixRowNonZeroElementTable[0].MatrixElementSupportTableLength   ← 0  (no supports)
  //   MatrixRowNonZeroElementTable[0].MatrixElementSupportTable[0]      ← NULL
  //   MatrixRowNonZeroElementTable[0].Node ← node
  //
  // SIDE EFFECTS:
  //   NonZeroElementsFound      ← 1 (or 0 for the "right" boundary corner case)
  //   rhs                       ← 0
  //   RhsSupportLength_CenterNodes ← 0
  //   RhsSupportLength_CornerNodes ← 0
  //
  // RETURN VALUE:
  //   true  → boundary handled (caller should `return;` from GetStencil)
  //   false → not a boundary (caller continues with normal assembly)
  //
  // NOTES:
  //   • Active only when periodic BCs are OFF: _PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_
  //   • The special case `isRightBoundaryCorner(x,node)` mirrors the original code’s
  //     treatment where the single entry is suppressed by setting NonZeroElementsFound=0.
  //------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // (3) Boundary handling for the linear system row.
  //     For physical (non-periodic) boundaries, ECSIM imposes a simplified row
  //     (often Dirichlet-like) to avoid referencing off-domain nodes.
  //     This block mirrors the legacy GetStencil() behavior: if the corner lies on
  //     a boundary, assemble a minimal row and return early.
  // ---------------------------------------------------------------------------
  auto boundary_short_circuit = [&]() -> bool {
    if (_PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_OFF_)
      if (isBoundaryCorner(x, node)) {
        MatrixRowNonZeroElementTable[0].i = i;
        MatrixRowNonZeroElementTable[0].j = j;
        MatrixRowNonZeroElementTable[0].k = k;

        MatrixRowNonZeroElementTable[0].MatrixElementValue = 0.0; // numeric set later by updater
        MatrixRowNonZeroElementTable[0].iVar = iVar;

        // One-parameter design: put 1.0 on PARAMETER[0] to encode the unit diagonal.
        MatrixRowNonZeroElementTable[0].MatrixElementParameterTable[0] = 1.0;
        MatrixRowNonZeroElementTable[0].MatrixElementParameterTableLength = 1;

        // No supports for boundary diagonal.
        MatrixRowNonZeroElementTable[0].MatrixElementSupportTableLength = 0;
        MatrixRowNonZeroElementTable[0].MatrixElementSupportTable[0] = NULL;

        MatrixRowNonZeroElementTable[0].Node = node;

        // Set row sizes / rhs
        NonZeroElementsFound = 1;
        rhs = 0.0;
        RhsSupportLength_CenterNodes = 0;
        RhsSupportLength_CornerNodes = 0;

        // Preserve original behavior: suppress even this entry on the "right" boundary.
        if (isRightBoundaryCorner(x, node)) NonZeroElementsFound = 0;

        return true; // caller should `return;`
      }
    return false; // not handled; proceed with normal stencil assembly
  };

 if (boundary_short_circuit()) return;



  if (!initMassMatrixOffsetTable) computeMassMatrixOffsetTable(); 

  const int indexAddition[3] = {0,-1,1};
  const int reversed_indexAddition[3] = {1,0,2};  //table to determine ii,jj,kk from i,j,k 

  char * NodeDataOffset = node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;


//====================================================================
  // Populate semantic RHS support entries for the **mass-matrix dot product**
  //====================================================================
  //
  // Legacy RHS form (inside old UpdateRhs()):
  //   The RHS included a term that (schematically) looked like
  //
  //     res += (-4*pi*dtTotal*theta) * [  Σ_{ii=0..26}  (E(neib_ii) · M_row(ii))  ] * E_conv
  //
  //   where:
  //     - E(neib_ii) is the electric field sampled from the 3×3×3 neighbor
  //       corner nodes around (i,j,k),
  //     - M_row(ii) are the stored mass-matrix coefficients for the current
  //       row component p=iVar, and
  //     - the compact mass matrix lives in the **row corner node** buffer.
  //
  // New goal:
  //   The new UpdateRhs() should contain only two simple loops over support
  //   containers (corner + center). Therefore, we *encode the dot product*
  //   as explicit semantic support entries of Quantity::MassMatrix and store
  //   the constant prefactor (-4*pi*dtTotal*theta) in CoefficientNEW.
  //
  // Strategy used here:
  //   - We append 81 entries into the caller-provided 'support_corner_vector':
  //       27 neighbors × 3 E-components (Ex,Ey,Ez).
  //   - Each entry represents one scalar term:  Mcoeff(aux_index) * Ecomp(neighbor)
  //     (SampleRhsScalar implements this when quantity==MassMatrix).
  //   - UpdateRhs() then accumulates:
  //        res += s.CoefficientNEW * SampleRhsScalar(s,iVar)
  //     where CoefficientNEW already includes the requested constant factor.
  //
  // IMPORTANT:
  //   This block only ADDS population of the new semantic support vector.
  //   It does NOT change any legacy RHS support arrays (RhsSupportTable_*),
  //   and it does NOT remove or modify any existing stencil/matrix assembly.
  //====================================================================
  {
    // Type alias for readability (matches the template used in this function signature).
    using RhsEntry = cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,_PIC_STENCIL_NUMBER_,_PIC_STENCIL_NUMBER_+1,16,1,1>::cRhsSupportTable;

    // The requested design requires that "-4*pi*dtTotal*theta" is stored in the
    // support table (decouples RHS evaluation from the rest of the field solver).
    //
    // NOTE:
    //   The legacy mass-matrix RHS also multiplied by E_conv, so we fold E_conv
    //   here to preserve the original scaling exactly.
    const double mmScale = (-4.0*Pi*dtTotal*theta) * E_conv;

    // The mass-matrix buffer is stored on the *row* corner node (i,j,k).
    // We keep an explicit pointer so SampleRhsScalar() can access the correct
    // MM buffer even when E is sampled from a neighbor corner.
    PIC::Mesh::cDataCornerNode *mmOwnerCorner =
      node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));

    // We are going to push 81 entries (27 neighbors × 3 components).
    // Reserve upfront to avoid repeated reallocations.
    support_corner_vector.reserve(support_corner_vector.size() + 81);

    // Build entries for:
    //   Σ_{q∈{x,y,z}} Σ_{neighbor∈3×3×3}  [ E_q(neighbor) * M_{p=iVar,q}(neighbor) ]
    //
    // Here:
    //   - q selects the E component (0=Ex,1=Ey,2=Ez),
    //   - neighbor_linear = ii + 3*jj + 9*kk selects the neighbor corner in the 3×3×3 block,
    //   - legacyElem = q*27 + neighbor_linear matches the historic "component block" ordering,
    //   - MassMatrixOffsetTable[iVar][legacyElem] returns the scalar index into the MM buffer.
    for (int q = 0; q < 3; q++) {
      for (int kk = 0; kk < 3; kk++) {
        for (int jj = 0; jj < 3; jj++) {
          for (int ii = 0; ii < 3; ii++) {

            // Neighbor offsets in {-1,0,+1} using the same ordering as the legacy code.
            const int di = indexAddition[ii];
            const int dj = indexAddition[jj];
            const int dk = indexAddition[kk];

            // Neighbor corner node that provides E_q(neighbor).
            PIC::Mesh::cDataCornerNode *neighborCorner =
              node->block->GetCornerNode(_getCornerNodeLocalNumber(i+di,j+dj,k+dk));

            // Linear neighbor index in the compact 3×3×3 neighborhood.
            const int neighbor_linear = ii + 3*jj + 9*kk;

            // Historic component-block layout:
            //   [0..26]  -> Ex neighbors
            //   [27..53] -> Ey neighbors
            //   [54..80] -> Ez neighbors
            const int legacyElem = q*27 + neighbor_linear;

            // Scalar index into the MM buffer for row component p=iVar.
            // The MM itself lives on the row node mmOwnerCorner.
            const int mmScalarIndex = MassMatrixOffsetTable[iVar][legacyElem];

            // Create and populate the semantic entry.
            //
            // For quantity==MassMatrix, SampleRhsScalar() is expected to:
            //   - sample E_q(neighborCorner) at the current E time level, and
            //   - multiply by the MM coefficient M[mmScalarIndex] stored at mmOwnerCorner,
            // returning one scalar dot-product term.
            RhsEntry s;
            s.SetCornerMassMatrix(mmScale, neighborCorner, mmOwnerCorner, mmScalarIndex, (unsigned char)q);

            // Keep legacy fields benign. The legacy RHS path uses only
            // (Coefficient, AssociatedDataPointer) and does not look at CoefficientNEW.
            // We do not insert these semantic entries into the legacy arrays, so they
            // should not affect the legacy UpdateRhs() behavior.
            s.Coefficient           = 0.0;
            s.AssociatedDataPointer = NULL;

            support_corner_vector.push_back(s);
          }
        }
      }
    }

    // NOTE:
    //   We intentionally do NOT modify RhsSupportLength_CornerNodes here in order
    //   to avoid breaking any existing legacy RHS layout assumptions. The mass-matrix
    //   semantic entries are added ONLY to support_corner_vector for consumption by
    //   the new UpdateRhs() path.
  }

  /*
  emit_compact_mass_structure
  -----------------------------------------
  PURPOSE
    Populate the 3×3×3 *compact* (radius-1) coupling structure of the ECSIM
    ΔE system’s **mass matrix term** for a single row (equation component p=iVar)
    at corner (i,j,k). For each target/column component q∈{Ex,Ey,Ez} and each
    neighbor corner (i+αx, j+αy, k+αz), α•∈{−1,0,+1}, the lambda:
      1) records WHICH column unknown it couples to (indices + component),
      2) allocates one PARAMETER slot (kept at length=1, value set later),
      3) stores one SUPPORT pointer to the mass/metric entry M_{pq} for this tap.

  MATHEMATICAL ROLE
    This is the structural encoding of the compact product (per row component p):
      \[
        (M\,\Delta\mathbf{E})_p(i,j,k)
        = \sum_{q\in\{x,y,z\}}
          \sum_{\alpha_x,\alpha_y,\alpha_z\in\{-1,0,1\}}
          M_{pq}(i,j,k;\boldsymbol{\alpha})\;
          \Delta E_q(i+\alpha_x,\;j+\alpha_y,\;k+\alpha_z).
      \]
    In the full operator
      \[
        \big[I + (\theta c\Delta t)^2(\nabla\nabla\!\cdot - \nabla^2 I)
          + 4\pi\,\theta\,\Delta t\,M\big]\;\Delta\mathbf{E}=\mathbf{r},
      \]
    this mass product is later scaled by \(4\pi\,\theta\,\Delta t\) and added to
    the numeric matrix entry during the **matrix-update** pass.

  WHAT GETS WRITTEN PER NONZERO (iElement)
    • Column addressing:
        .i,.j,.k   ← (i+αx, j+αy, k+αz),
        .iVar      ← q (the column component).
    • Value placeholder:
        .MatrixElementValue ← 0.0 (the actual number is computed later).
    • PARAMETER (for I + curl–curl part; original design = one scalar):
        .MatrixElementParameterTable[0] ← 0.0;
        .MatrixElementParameterTableLength ← 1;
      The matrix-update routine will overwrite/accumulate this to form the
      compact \(I + (\nabla\nabla\!\cdot - \nabla^2 I)\) contribution.
    • SUPPORT (for mass matrix M; original design = one pointer):
        .MatrixElementSupportTableLength ← 1;
        .MatrixElementSupportTable[0] ←
            (double*)NodeDataOffset + MassMatrixOffsetIndex
          + MassMatrixOffsetTable[iVar][iElement];
      This is an opaque address to the stored coefficient
      \(M_{pq}(i,j,k;\boldsymbol{\alpha})\). The update pass dereferences it and
      adds \( (4\pi\,\theta\,\Delta t)\times M_{pq} \) to the numeric entry.

  LAYOUT / ORDERING
    • 27 neighbors per component, 3 components → **81** entries per interior row.
      Linear index: iElement = q*27 + (ii + 3*jj + 9*kk), with ii,jj,kk∈{0,1,2}.
    • (p=iVar) identifies the row component; (q=iVarIndex) identifies the column
      component block. The per-tap mass offset uses both: MassMatrixOffsetTable[p][iElement].

  NOTES / GOTCHAS
    • Do not assign numeric coefficients here—leave MatrixElementValue at 0.0.
      Numbers are produced by the matrix-update routine using PARAMETER[0] & SUPPORT[0].
    • Pointer arithmetic assumes NodeDataOffset is aligned to “double*” slots for this
      region. If your build uses byte offsets, switch to `char*` arithmetic consistently.
    • Boundary/periodic remapping and row compaction are performed elsewhere.
*/
  auto emit_compact_mass_structure = [&]() {
    for (int iVarIndex = 0; iVarIndex < 3; iVarIndex++) {
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
            int iNode = i + indexAddition[ii];
            int jNode = j + indexAddition[jj];
            int kNode = k + indexAddition[kk];

            // Linear slot: 27 per component, 3 components total → 81 entries
            int iElement = iVarIndex * 27 + ii + jj * 3 + kk * 9;
	    int legacyElem = iVarIndex * 27 + ii + jj * 3 + kk * 9;

            const int di = indexAddition[ii];
            const int dj = indexAddition[jj];
            const int dk = indexAddition[kk];
            //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	    iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex); 


            MatrixRowNonZeroElementTable[iElement].i = iNode;
            MatrixRowNonZeroElementTable[iElement].j = jNode;
            MatrixRowNonZeroElementTable[iElement].k = kNode;

            // already initialized in LinearSystemSolver->h
            MatrixRowNonZeroElementTable[iElement].MatrixElementValue = 0.0;
            MatrixRowNonZeroElementTable[iElement].iVar = iVarIndex;

            // Original design: exactly one parameter per element
            MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0] = 0.0;
            MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength = 1;

            // Exactly one support pointer per element (mass/metric slot)
            MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength = 1;

            // Pointer math as in the original codebase:
            MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTable[0] =
              (double*)NodeDataOffset
              + MassMatrixOffsetIndex
              + MassMatrixOffsetTable[iVar][legacyElem];
          }
        }
      }
    }
  };

  // invoke it
  emit_compact_mass_structure();
    

  int indexOffset[5] = {0,-1,1,-2,2};

  int reversed_indexOffset[5]={3,1,0,2,4}; //table to convert i,j,k ->ii,jj,kk

if (_PIC_STENCIL_NUMBER_==375) {
  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){
    int cntTemp =0;
    
    for (int kk=0;kk<5;kk++){
      for (int jj=0;jj<5;jj++){	
	for (int ii=0;ii<5;ii++){       
	  
	  if (ii<3 && jj<3 && kk<3) continue;
          int iNode = i+indexOffset[ii];
          int jNode = j+indexOffset[jj];
          int kNode = k+indexOffset[kk];
          int iElement = 81+iVarIndex*98+cntTemp;

      {
        const int di = indexOffset[ii];
        const int dj = indexOffset[jj];
        const int dk = indexOffset[kk];
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

	  cntTemp++;
          MatrixRowNonZeroElementTable[iElement].i=iNode;
          MatrixRowNonZeroElementTable[iElement].j=jNode;
          MatrixRowNonZeroElementTable[iElement].k=kNode;

          MatrixRowNonZeroElementTable[iElement].MatrixElementValue=0.0;
          MatrixRowNonZeroElementTable[iElement].iVar=iVarIndex;
          MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]=0.0;
          MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTableLength=1;
          MatrixRowNonZeroElementTable[iElement].MatrixElementSupportTableLength = 0;
        }
      }
    }
  }
}


  //laplacian
  auto *laplacian=get_laplacian_stencil();

  for (int idim=0;idim<3;idim++) {
    cStencil::cStencilData *st=LaplacianStencil+idim;

    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1]; 
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int index=ii+jj*3+kk*9;
      int iElement = index + iVar*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVar); 
	iElement=indexer[iVar].get_idx(di, dj, dk, iVar);
      }

      //minus laplacian
      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]-=st->Data[it].a*coeffSqr[idim]; 
    }
  }
  

  //plus self
  //
  
int idx1;
      {
        const int di = 0;
        const int dj = 0;
        const int dk = 0;
        //indexer[iVar].check_and_set(27*iVar,di, dj, dk, iVar);
	idx1=indexer[iVar].get_idx(di, dj, dk, iVar);
      }
  
  MatrixRowNonZeroElementTable[idx1].MatrixElementParameterTable[0]+=1;


  //plus graddiv E
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++) {
    cStencil::cStencilData *st=&GradDivStencil[iVar][iVarIndex];
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1];
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]+=(1-corrCoeff)*st->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }


    //add contribution from a larger stencil (divE correction)
    for (int it=0;it<st375->Length;it++) {
      if ((st375->Data[it].i<-1)||(st375->Data[it].i>1) || (st375->Data[it].j<-1)||(st375->Data[it].j>1) ||(st375->Data[it].k<-1)||(st375->Data[it].k>1) ) continue;

      int ii=reversed_indexAddition[st375->Data[it].i+1];
      int jj=reversed_indexAddition[st375->Data[it].j+1];
      int kk=reversed_indexAddition[st375->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]+=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }
  }


int OrderingOffsetTable[5][5][5];

if (_PIC_STENCIL_NUMBER_==375) {   
  int cntTemp = 0;

  for (int kk=0;kk<5;kk++){
    for (int jj=0;jj<5;jj++){
      for (int ii=0;ii<5;ii++){
        if (ii<3 && jj<3 && kk<3) continue;

        OrderingOffsetTable[ii][jj][kk]=cntTemp;
        cntTemp++;
      }
    }
  }

  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st375->Length;it++) {

      int ii=reversed_indexOffset[st375->Data[it].i+2];
      int jj=reversed_indexOffset[st375->Data[it].j+2];
      int kk=reversed_indexOffset[st375->Data[it].k+2];

      if (ii<3 && jj<3 && kk<3) continue;
      int iElement = 81+iVarIndex*98+OrderingOffsetTable[ii][jj][kk];
      int nodeIndex= 27+OrderingOffsetTable[ii][jj][kk];

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
      }

      MatrixRowNonZeroElementTable[iElement].MatrixElementParameterTable[0]+=corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex];
    }
  }
}


  //find corners outside the boundary
  vector<int> pointLeft;
  int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

  for (int ii=0;ii<_PIC_STENCIL_NUMBER_;ii++) {
    MatrixRowNonZeroElementTable[ii].Node=node;

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * nodeTemp = node;
    int i0,j0,k0; //for test
    i0 = MatrixRowNonZeroElementTable[ii].i;
    j0 = MatrixRowNonZeroElementTable[ii].j;
    k0 = MatrixRowNonZeroElementTable[ii].k;

    if ((i0>=iMax) && (nodeTemp!=NULL)) {
      i0-=_BLOCK_CELLS_X_;
      nodeTemp=nodeTemp->GetNeibFace(1,0,0,PIC::Mesh::mesh);
    }
    else if (i0<0 && nodeTemp!=NULL) {
      i0+=_BLOCK_CELLS_X_;
      nodeTemp=nodeTemp->GetNeibFace(0,0,0,PIC::Mesh::mesh);
    }

    if ((j0>=jMax) && (nodeTemp!=NULL)) {
      j0-=_BLOCK_CELLS_Y_;
      nodeTemp=nodeTemp->GetNeibFace(3,0,0,PIC::Mesh::mesh);
    }
    else if (j0<0 && nodeTemp!=NULL) {
      j0+=_BLOCK_CELLS_Y_;
      nodeTemp=nodeTemp->GetNeibFace(2,0,0,PIC::Mesh::mesh);
    }

    if ((k0>=kMax) && (nodeTemp!=NULL)) {
      k0-=_BLOCK_CELLS_Z_;
      nodeTemp=nodeTemp->GetNeibFace(5,0,0,PIC::Mesh::mesh);
    }
    else if (k0<0 && nodeTemp!=NULL) {
      k0+=_BLOCK_CELLS_Z_;
      nodeTemp=nodeTemp->GetNeibFace(4,0,0,PIC::Mesh::mesh);
    }

    double xlocal[3];
    int indlocal[3]={MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k};
    int indexG_local[3];
    bool isFixed=false;

    for (int idim=0; idim<3; idim++) {
      xlocal[idim]=MatrixRowNonZeroElementTable[ii].Node->xmin[idim]+indlocal[idim]*dx[idim];
    }
     
    pointLeft.push_back(ii);

    if (MatrixRowNonZeroElementTable[ii].Node==NULL){
      pointLeft.pop_back();
      continue;
    }
    else if (MatrixRowNonZeroElementTable[ii].Node->IsUsedInCalculationFlag==false) {
      pointLeft.pop_back();
      continue;
    }
    
    if (nodeTemp==NULL){
      pointLeft.pop_back();
      continue;
    }
    else if (nodeTemp->IsUsedInCalculationFlag==false){
      pointLeft.pop_back();
      continue;
    }

    if (isBoundaryCorner(xlocal,node)) pointLeft.pop_back();
  }

  for (int ii=0; ii<pointLeft.size();ii++){
    int copyFrom = pointLeft[ii];

    if (ii!=copyFrom){
      MatrixRowNonZeroElementTable[ii].i=MatrixRowNonZeroElementTable[copyFrom].i;
      MatrixRowNonZeroElementTable[ii].j=MatrixRowNonZeroElementTable[copyFrom].j;
      MatrixRowNonZeroElementTable[ii].k=MatrixRowNonZeroElementTable[copyFrom].k;

      MatrixRowNonZeroElementTable[ii].MatrixElementValue= MatrixRowNonZeroElementTable[copyFrom].MatrixElementValue;
      MatrixRowNonZeroElementTable[ii].iVar=MatrixRowNonZeroElementTable[copyFrom].iVar;

      MatrixRowNonZeroElementTable[ii].MatrixElementParameterTable[0]=MatrixRowNonZeroElementTable[copyFrom].MatrixElementParameterTable[0];
      MatrixRowNonZeroElementTable[ii].MatrixElementParameterTableLength=MatrixRowNonZeroElementTable[copyFrom].MatrixElementParameterTableLength;
      MatrixRowNonZeroElementTable[ii].MatrixElementSupportTableLength = MatrixRowNonZeroElementTable[copyFrom].MatrixElementSupportTableLength;
         
      MatrixRowNonZeroElementTable[ii].MatrixElementSupportTable[0]=MatrixRowNonZeroElementTable[copyFrom].MatrixElementSupportTable[0];
    }

  }
  
  NonZeroElementsFound=pointLeft.size();
  
  //NonZeroElementsFound=81;

  //  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){

  // ---------------------------------------------------------------------------
  // PART 2a: BUILD CORNER-NODE RHS SUPPORT FOR E^n (semantic-vector path)
  // ---------------------------------------------------------------------------
  // Goal:
  //   Populate `support_corner_vector` directly with semantic Corner/E entries,
  //   without relying on (or copying from) `RhsSupportTable_CornerNodes`.
  //
  // Legacy compatibility:
  //   We still *initialise/update* `RhsSupportTable_CornerNodes` and
  //   `RhsSupportLength_CornerNodes` because older code paths may still read them.
  //   However, the semantic RHS evaluation path (UpdateRhs v2) uses ONLY the
  //   support vectors.
  //
  // IMPORTANT:
  //   Do NOT assume the RHS corner support length equals `_PIC_STENCIL_NUMBER_`.
  //   We use the runtime `indexer[iVar]` mapping and only create entries that are
  //   actually referenced.
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // (5) Semantic RHS support helpers: corner-node samples (E and J).
  //
  //     In Ampère–Maxwell, the RHS contains (1-θ) contributions built from known
  //     fields at time level n, e.g. E^n and J^n, while the θ-part is moved to the
  //     left-hand side to form the matrix operator.
  //
  //     We represent each required corner sample with a unique semantic entry and
  //     accumulate coefficients into it. The uniqueness key is the runtime indexer slot.
  // ---------------------------------------------------------------------------

  // Map: iElement (indexer idx) -> position in support_corner_vector
  std::map<int,int> cornerE_pos;

  // Ensure a semantic Corner/E entry exists for (di,dj,dk,comp) and return it.
  auto ensure_corner_E = [&](int di,int dj,int dk,unsigned char comp) -> RhsEntry& {
    const int iElement = indexer[iVar].get_idx(di,dj,dk,(int)comp);

    auto it = cornerE_pos.find(iElement);
    if (it!=cornerE_pos.end()) return support_corner_vector[it->second];

    PIC::Mesh::cDataCornerNode *cn =
      node->block->GetCornerNode(_getCornerNodeLocalNumber(i+di,j+dj,k+dk));

    char *pnt   = (cn!=NULL) ? cn->GetAssociatedDataBufferPointer() : NULL;
    char *assoc = (pnt!=NULL) ? pnt + PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset : NULL;

    // Semantic entry for new RHS evaluation path
    RhsEntry e;
    e.Coefficient = 0.0;                  // legacy field kept consistent
    e.AssociatedDataPointer = assoc;       // legacy pointer path (optional)
    e.SetCornerE(0.0, (assoc!=NULL) ? cn : NULL, comp); // sets CoefficientNEW
    e.idx = iElement;

    support_corner_vector.push_back(e);
    const int pos = (int)support_corner_vector.size()-1;
    cornerE_pos[iElement]=pos;

    // Legacy array init (kept for compatibility / self-checks)
    RhsSupportTable_CornerNodes[iElement].Coefficient = 0.0;
    RhsSupportTable_CornerNodes[iElement].CoefficientNEW = 0.0;
    RhsSupportTable_CornerNodes[iElement].AssociatedDataPointer = assoc;
    RhsSupportTable_CornerNodes[iElement].SetCornerE(0.0, (assoc!=NULL) ? cn : NULL, comp);

    return support_corner_vector[pos];
  };

  // Accumulate a delta coefficient into BOTH the semantic vector entry and the
  // legacy array entry.
  auto add_corner_E_coeff = [&](int di,int dj,int dk,unsigned char comp,double delta) -> void {
    RhsEntry &e = ensure_corner_E(di,dj,dk,comp);

    // Keep both coefficient fields consistent (UpdateRhs v2 uses CoefficientNEW).
    e.Coefficient    += delta;
    e.CoefficientNEW += delta;

    const int idxLocal = e.idx;
    RhsSupportTable_CornerNodes[idxLocal].Coefficient    += delta;
    RhsSupportTable_CornerNodes[idxLocal].CoefficientNEW += delta;
  };

  // Pre-create ALL Corner/E entries referenced by the configured stencil for this row.
  // 3×3×3 block (always present)
  for (unsigned char comp=0; comp<3; comp++) {
    for (int kk=0; kk<3; kk++) for (int jj=0; jj<3; jj++) for (int ii=0; ii<3; ii++) {
      ensure_corner_E(indexAddition[ii], indexAddition[jj], indexAddition[kk], comp);
    }
  }

  // 5×5×5 extension (only for the larger ECSIM stencil)
  if (_PIC_STENCIL_NUMBER_==375) {
    for (unsigned char comp=0; comp<3; comp++) {
      for (int kk=0; kk<5; kk++) for (int jj=0; jj<5; jj++) for (int ii=0; ii<5; ii++) {
        if (ii<3 && jj<3 && kk<3) continue;
        ensure_corner_E(indexOffset[ii], indexOffset[jj], indexOffset[kk], comp);
      }
    }
  }

  //laplacian
  for (int idim=0;idim<3;idim++) {
    cStencil::cStencilData *st=LaplacianStencil+idim;

    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1];
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int index=ii+jj*3+kk*9;
      int iElement = index + iVar*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVar);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVar);
              // plus laplacian
        add_corner_E_coeff(di,dj,dk,(unsigned char)iVar, st->Data[it].a*coeffSqr[idim]);
      }
    }
  }


  //minus graddiv E
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    cStencil::cStencilData *st=&GradDivStencil[iVar][iVarIndex];


    for (int it=0;it<st->Length;it++) {
      int ii=reversed_indexAddition[st->Data[it].i+1];
      int jj=reversed_indexAddition[st->Data[it].j+1];
      int kk=reversed_indexAddition[st->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st->Data[it].i;
        const int dj = st->Data[it].j;
        const int dk = st->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
              add_corner_E_coeff(di,dj,dk,(unsigned char)iVarIndex, -(1-corrCoeff)*st->Data[it].a*coeff[iVar]*coeff[iVarIndex]);
      }
    }
  }

if (_PIC_STENCIL_NUMBER_==375) {
  for (int iVarIndex=0;iVarIndex<3;iVarIndex++){
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st375->Length;it++) {

      if ((st375->Data[it].i<-1)||(st375->Data[it].i>1) || (st375->Data[it].j<-1)||(st375->Data[it].j>1) ||(st375->Data[it].k<-1)||(st375->Data[it].k>1) ) continue;

      int ii=reversed_indexAddition[st375->Data[it].i+1];
      int jj=reversed_indexAddition[st375->Data[it].j+1];
      int kk=reversed_indexAddition[st375->Data[it].k+1];

      int nodeIndex=ii+jj*3+kk*9;
      int iElement = nodeIndex + iVarIndex*27;

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
              add_corner_E_coeff(di,dj,dk,(unsigned char)iVarIndex, -corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex]);
      }
    }
  }


  for (int iVarIndex=0; iVarIndex<3; iVarIndex++){
    cStencil::cStencilData *st375=&GradDivStencil375[iVar][iVarIndex];

    for (int it=0;it<st375->Length;it++) {
      int ii=reversed_indexOffset[st375->Data[it].i+2];
      int jj=reversed_indexOffset[st375->Data[it].j+2];
      int kk=reversed_indexOffset[st375->Data[it].k+2];

      if (ii<3 && jj<3 && kk<3) continue;

      int iElement = 81+iVarIndex*98+OrderingOffsetTable[ii][jj][kk];
      int nodeIndex = 27+OrderingOffsetTable[ii][jj][kk];

      {
        const int di = st375->Data[it].i;
        const int dj = st375->Data[it].j;
        const int dk = st375->Data[it].k;
        //indexer[iVar].check_and_set(iElement,di, dj, dk, iVarIndex);
	iElement=indexer[iVar].get_idx(di, dj, dk, iVarIndex);
              add_corner_E_coeff(di,dj,dk,(unsigned char)iVarIndex, -corrCoeff*st375->Data[it].a*coeff[iVar]*coeff[iVarIndex]);
      }
    }
  }

}

int idx=indexer[iVar].get_index_count(); //_PIC_STENCIL_NUMBER_  
int idx2;
      {
        const int di = 0;
        const int dj = 0;
        const int dk = 0;
        //indexer[iVar].check_and_set(27*iVar,di, dj, dk, iVar);
	idx2=indexer[iVar].get_idx(di, dj, dk, iVar);
      }

  RhsSupportTable_CornerNodes[idx].AssociatedDataPointer=RhsSupportTable_CornerNodes[idx2].AssociatedDataPointer;
  RhsSupportTable_CornerNodes[idx].Coefficient=-4*Pi*dtTotal*theta;

  RhsSupportTable_CornerNodes[idx].SetCornerJ(-4*Pi*dtTotal*theta, RhsSupportTable_CornerNodes[idx2].corner, (unsigned char)iVar);
 
  RhsSupportLength_CornerNodes=idx+1;

  // ---------------------------------------------------------------------------
  // (8) Add the current density J term to the RHS support.
  //
  //     In Gaussian units, Ampère’s law contributes a -4π J term. Under ECSIM time-centering,
  //     the explicit part of J (at time level n) is included on the RHS with coefficient
  //     (-4π θ Δt) in this formulation (matching the legacy code’s convention).
  // ---------------------------------------------------------------------------
  // Append the J term to the *semantic* corner support vector (do not copy from the legacy array).
  {
    PIC::Mesh::cDataCornerNode *rowCorner = node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));
    char *pnt   = (rowCorner!=NULL) ? rowCorner->GetAssociatedDataBufferPointer() : NULL;
    char *assoc = (pnt!=NULL) ? pnt + PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset : NULL;

    RhsEntry jEntry;
    jEntry.Coefficient = -4*Pi*dtTotal*theta;
    jEntry.AssociatedDataPointer = assoc;
    jEntry.SetCornerJ(-4*Pi*dtTotal*theta, (assoc!=NULL) ? rowCorner : NULL, (unsigned char)iVar);
    jEntry.idx = idx;

    support_corner_vector.push_back(jEntry);
  }


  //Ex^n,Ey^n,Ez^n
  rhs=0.0;


  int indexAdditionB[2] = {-1,0};
  // ---------------------------------------------------------------------------
  // PART 3: BUILD CENTER-NODE RHS SUPPORT (curl B)  [semantic-vector path]
  // ---------------------------------------------------------------------------
  // Goal:
  //   Populate `support_center_vector` directly with semantic Center/B entries,
  //   WITHOUT relying on (or copying from) `RhsSupportTable_CenterNodes`.
  //
  // Notes:
  //   - We aggregate duplicates by (center-node pointer, B-component) so the
  //     semantic RHS is a compact dot-product.
  //   - The legacy center RHS table is not populated here anymore; the signature
  //     is preserved for compatibility during the transition.
  // ---------------------------------------------------------------------------

  // (Optional) silence unused-parameter warnings in builds that still pass these.
  (void)RhsSupportTable_CenterNodes;

  // ---------------------------------------------------------------------------
  // (6) Semantic RHS support helpers: center-node samples (B for curl(B)).
  //
  //     The discrete Ampère–Maxwell RHS includes a curl(B^n) term.
  //     B is stored at *cell centers*, so we build support_center_vector entries
  //     that point to specific center nodes and B components (Bx/By/Bz).
  //
  //     As with corner E, we ensure uniqueness and accumulate coefficients rather
  //     than relying on a fixed-size legacy table.
  // ---------------------------------------------------------------------------
  struct CenterBKey {
    PIC::Mesh::cDataCenterNode* cn;
    unsigned char comp; // 0=Bx,1=By,2=Bz
    bool operator<(const CenterBKey& o) const {
      if (cn != o.cn) return cn < o.cn;
      return comp < o.comp;
    }
  };

  std::map<CenterBKey,int> centerB_pos;

  auto ensure_center_B = [&](PIC::Mesh::cDataCenterNode* cn, unsigned char bcomp) -> RhsEntry& {
    CenterBKey key{cn, bcomp};
    auto it = centerB_pos.find(key);
    if (it != centerB_pos.end()) return support_center_vector[it->second];

    // Keep legacy pointer field consistent (even though semantic sampling should not require it).
    char *pnt   = (cn!=NULL) ? cn->GetAssociatedDataBufferPointer() : NULL;
    char *assoc = (pnt!=NULL) ? pnt + PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset : NULL;

    RhsEntry e;
    e.Coefficient = 0.0;
    e.AssociatedDataPointer = assoc;
    e.SetCenterB(0.0, cn, bcomp);  // also initializes CoefficientNEW
    e.idx = -1; // no indexer meaning for center taps

    support_center_vector.push_back(e);
    const int pos = (int)support_center_vector.size() - 1;
    centerB_pos[key] = pos;
    return support_center_vector[pos];
  };

  auto add_center_B_coeff = [&](PIC::Mesh::cDataCenterNode* cn, unsigned char bcomp, double delta) -> void {
    RhsEntry& e = ensure_center_B(cn, bcomp);
    e.Coefficient    += delta;  // legacy field (kept consistent)
    e.CoefficientNEW += delta;  // semantic UpdateRhs uses CoefficientNEW
  };

  // curlB contributions:
  //   iVar==0: rhs += dBz/dy - dBy/dz
  //   iVar==1: rhs += dBx/dz - dBz/dx
  //   iVar==2: rhs += dBy/dx - dBx/dy
  rhs = 0.0;

  if (iVar==0) {
    // + dBz/dy  (Bz at j and j-1)
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[ii], j,   k+indexAdditionB[jj]));
      add_center_B_coeff(center, /*Bz*/2, +coeff4[1]);
    }
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[ii], j-1, k+indexAdditionB[jj]));
      add_center_B_coeff(center, /*Bz*/2, -coeff4[1]);
    }

    // - dBy/dz  (By at k and k-1)
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k));
      add_center_B_coeff(center, /*By*/1, -coeff4[2]);
    }
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k-1));
      add_center_B_coeff(center, /*By*/1, +coeff4[2]);
    }
  }

  if (iVar==1) {
    // + dBx/dz (Bx at k and k-1)
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k));
      add_center_B_coeff(center, /*Bx*/0, +coeff4[2]);
    }
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[ii], j+indexAdditionB[jj], k-1));
      add_center_B_coeff(center, /*Bx*/0, -coeff4[2]);
    }

    // - dBz/dx (Bz at i and i-1)
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i,   j+indexAdditionB[jj], k+indexAdditionB[ii]));
      add_center_B_coeff(center, /*Bz*/2, -coeff4[0]);
    }
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i-1, j+indexAdditionB[jj], k+indexAdditionB[ii]));
      add_center_B_coeff(center, /*Bz*/2, +coeff4[0]);
    }
  }

  if (iVar==2) {
    // + dBy/dx (By at i and i-1)
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i,   j+indexAdditionB[jj], k+indexAdditionB[ii]));
      add_center_B_coeff(center, /*By*/1, +coeff4[0]);
    }
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i-1, j+indexAdditionB[jj], k+indexAdditionB[ii]));
      add_center_B_coeff(center, /*By*/1, -coeff4[0]);
    }

    // - dBx/dy (Bx at j and j-1)
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[jj], j,   k+indexAdditionB[ii]));
      add_center_B_coeff(center, /*Bx*/0, -coeff4[1]);
    }
    for (int ii=0; ii<2; ii++) for (int jj=0; jj<2; jj++) {
      auto *center = node->block->GetCenterNode(
        _getCenterNodeLocalNumber(i+indexAdditionB[jj], j-1, k+indexAdditionB[ii]));
      add_center_B_coeff(center, /*Bx*/0, +coeff4[1]);
    }
  }
}


