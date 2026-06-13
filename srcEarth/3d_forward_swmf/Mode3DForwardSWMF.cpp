//======================================================================================
// Mode3DForwardSWMF.cpp
//======================================================================================
//
// SWMF-coupled initialization bridge for the Earth 3-D forward energetic-particle
// branch.
//
// DESIGN — TWO-PHASE INITIALIZATION
// -----------------------------------
// The original single amps_init() was called at the END of main_lib.cpp::amps_init(),
// after both amps_init_mesh() and amps_init() had already executed.  That ordering
// caused six initialization gaps relative to the standalone 3d_forward path:
//
//   Gap 1  Earth::ModelMode not set before amps_init_mesh()
//          → main_lib.cpp took the cutoff-rigidity branch instead of the
//            BoundaryInjection/MAIN_LIB_GEO branch.
//   Gap 2  cDensity3D::ConfigureEnergyGrid() called too late
//          → PIC::Mesh::initCellSamplingDataBuffer() allocated per-cell buffers
//            using the wrong (compile-time default) nEnergyBins.
//   Gap 3  Earth::Init_AfterParser() never called
//          → only reachable via MAIN_LIB_GEO::amps_init_mesh(), which was
//            never entered because of Gap 1.
//   Gap 4  BoundingBoxInjection::InitDirectionIMF() not called
//          → main_lib.cpp::amps_init() calls it in the BoundaryInjectionMode
//            branch, but that branch was skipped because of Gap 1.
//   Gap 5  PIC::Mover::ProcessOutsideDomainParticles left pointing at the
//          cutoff-rigidity handler set by the cutoff-rigidity amps_init_mesh()
//          path; should be null for the forward mode.
//   Gap 6  BoundingBoxInjection::SetPrm() not called before InitDirectionIMF().
//
// The fix splits initialization into two hooks:
//
//   amps_pre_init()  — called BEFORE amps_init_mesh() in main_lib.cpp
//      Fixes Gaps 1, 2, 3, 4, 6 by setting the mode and energy grid before
//      the mesh is built, and by registering prm before InitDirectionIMF().
//
//   amps_init()      — called at the END of main_lib.cpp::amps_init()
//      Fixes Gap 5.  Continues to do all post-mesh forward-mode setup,
//      reusing the prm parsed by amps_pre_init() (no second file read).
//
// The only intended physical difference from standalone 3d_forward remains the
// source of background B/E fields: in _PIC_COUPLER_MODE__SWMF_ builds they come
// from PIC::CPLR rather than Tsyganenko/DATAFILE.
//
//======================================================================================

#include "Mode3DForwardSWMF.h"

#include "../3d_forward/Mode3DForward.h"
#include "../3d_forward/ForwardParticleMovers.h"
#include "../3d_forward/Density3D.h"
#include "../3d_forward/SphereFlux3D.h"
#include "../3d/Mode3D.h"
#include "../3d/CutoffRigidityMode3D.h"
#include "../gridless/DipoleInterface.h"
#include "../boundary/spectrum.h"
#include "../util/amps_param_parser.h"
#include "../Earth.h"

#include "pic.h"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <cstring>

namespace Earth {
namespace Mode3DForwardSWMF {

//======================================================================================
// Shared parsed parameters
//======================================================================================
// Populated once by amps_pre_init() and consumed read-only by amps_init().
// Using a file-scope static avoids re-parsing AMPS_PARAM.in and guarantees both
// hooks operate on an identical parameter set.
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
static EarthUtil::AmpsParam s_prm;
static bool                 s_pre_initialized = false;
static long int             s_cutoff_call_counter = 0;
static double               s_cutoff_elapsed_time_s = 0.0;
#endif


namespace {

void AnalyticDipoleMagneticField_(double *xIn,double *bOut) {
  Earth::GridlessMode::Dipole::GetB_Tesla(xIn,bOut);
}

bool IsCutoffTarget_(const EarthUtil::AmpsParam& prm) {
  // The parser default is CUTOFF_RIGIDITY for standalone backward
  // compatibility.  In coupled SWMF runs, require an explicit CALC_TARGET so
  // old 3d_forward input files that omit #CALCULATION_MODE keep their
  // historical forward-injection behavior.
  return prm.calc.targetExplicit &&
         EarthUtil::ToUpper(prm.calc.target) == "CUTOFF_RIGIDITY";
}

bool IsForwardTarget_(const EarthUtil::AmpsParam& prm) {
  return !IsCutoffTarget_(prm);
}

void ApplyParsedDomainForCutoff_(const EarthUtil::AmpsParam& prm) {
  Earth::Mode3D::ParsedDomainActive = true;
  Earth::Mode3D::ParsedDomainMin[0] = prm.domain.xMin * 1000.0;
  Earth::Mode3D::ParsedDomainMin[1] = prm.domain.yMin * 1000.0;
  Earth::Mode3D::ParsedDomainMin[2] = prm.domain.zMin * 1000.0;
  Earth::Mode3D::ParsedDomainMax[0] = prm.domain.xMax * 1000.0;
  Earth::Mode3D::ParsedDomainMax[1] = prm.domain.yMax * 1000.0;
  Earth::Mode3D::ParsedDomainMax[2] = prm.domain.zMax * 1000.0;
}

std::string FormatCutoffOutputSuffix_(long int callIndex,double tSim_s) {
  std::ostringstream ss;
  ss << ".swmf_n" << std::setw(6) << std::setfill('0') << callIndex
     << "_t" << std::setw(12) << std::setfill('0')
     << std::fixed << std::setprecision(3) << std::max(0.0,tSim_s) << "s";
  return ss.str();
}

double CoupledCutoffOutputDtSeconds_(const EarthUtil::AmpsParam& prm) {
  // amps_time_step() does not receive an explicit SWMF time argument.  Use the
  // configured temporal field-update cadence when present; otherwise fall back
  // to the AMPS particle time step if it is available.  The call index is always
  // included in the file name, so outputs remain unique even if this fallback is
  // zero.  A future SWMF-side hook can update s_cutoff_elapsed_time_s directly
  // before calling amps_time_step() if exact MHD time stamps are required.
  if (prm.temporal.fieldUpdateDt_min > 0.0) return prm.temporal.fieldUpdateDt_min * 60.0;

  if (PIC::nTotalSpecies > 0) {
    const double dt = PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
    if (std::isfinite(dt) && dt > 0.0) return dt;
  }

  return 0.0;
}

//======================================================================================
// SWMF-COUPLED CUTOFF: GLOBAL CELL-CENTERED MAGNETIC FIELD MATERIALIZATION
//======================================================================================
//
// Normal SWMF coupling in AMPS is memory-local:
//   * on an MPI rank, AMPS allocates blocks owned by that rank plus the usual
//     domain-boundary/ghost-layer blocks;
//   * RecieveCenterPointData() stores the coupled magnetic field only in those
//     allocated blocks, at
//       cell->GetAssociatedDataBufferPointer() + PIC::CPLR::SWMF::MagneticFieldOffset;
//   * blocks owned by other MPI ranks normally have node->block == NULL on this rank.
//
// That is efficient for local particle transport, but it is not enough for backward
// cutoff tracing.  A single trajectory can move through any region of the global AMR
// tree, and the Mode3D cutoff evaluator expects the magnetic field to be available
// wherever findTreeNode(x,...) lands on the local process.
//
// The helpers below convert the distributed SWMF field into a replicated, read-only
// mesh field immediately before each cutoff snapshot:
//
//   1. Walk the replicated AMR tree and assign each used leaf block a deterministic
//      global dense index in node->Temp_ID.  The AMR tree topology is the same on all
//      ranks, so this traversal produces identical IDs everywhere.
//
//   2. Allocate any missing leaf blocks on every rank.  This makes node->block valid
//      for all used leaves on every process, which is what the cutoff tracer needs.
//
//   3. On each rank, pack only the interior cells of leaf blocks whose owner is that
//      rank (node->Thread == PIC::ThisThread).  Ghost cells are intentionally not used
//      as sources, because the SWMF owner block is the authoritative copy.
//
//   4. MPI_Allreduce the dense B array and an integer ownership mask.  After the
//      reduction, every rank has a global cell-centered B cache for every interior
//      cell in every used AMR leaf.
//
//   5. Populate all allocated blocks, including their ghost cells.  Each ghost-cell
//      center is mapped back to the owning leaf/interior-cell index through the AMR
//      tree, then filled from the global cache.  This gives the standard AMPS
//      cell-centered interpolation routines the same data layout they normally see.
//
// The global B arrays exist only inside PrepareGlobalSWMFCoupledMagneticFieldForCutoff().
// The persistent state is the magnetic field written back into each cell's associated
// data buffer.  After that, Earth::Mode3D::RunCutoffRigidity() can use the existing
// mesh-based GetB path without knowing whether the field came from SWMF, an analytic
// debug override, or a standalone field initializer.
//======================================================================================

typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cAMRNode;

// Number of interior cell centers stored per AMR block.  The global dense cache
// below deliberately indexes only physical/interior cells.  Ghost cells are later
// reconstructed by mapping their center coordinates to an owning interior cell.
long int InteriorCellsPerBlock_() {
  return static_cast<long int>(_BLOCK_CELLS_X_) *
         static_cast<long int>(_BLOCK_CELLS_Y_) *
         static_cast<long int>(_BLOCK_CELLS_Z_);
}

// Convert a block-local interior index (i,j,k) to a compact row-major scalar index.
// The order must be consistent everywhere because it is used to build MPI-reduced
// dense arrays, not just local temporary arrays.
long int InteriorCellIndex_(int i,int j,int k) {
  return static_cast<long int>(i) +
         static_cast<long int>(_BLOCK_CELLS_X_) *
         (static_cast<long int>(j) +
          static_cast<long int>(_BLOCK_CELLS_Y_) * static_cast<long int>(k));
}

// Dense global index of an interior cell.  node->Temp_ID is assigned by
// AssignGlobalLeafTempIds_() immediately before packing and is valid only for this
// global-B materialization pass.
long int GlobalCellIndex_(cAMRNode *node,int i,int j,int k) {
  return node->Temp_ID * InteriorCellsPerBlock_() + InteriorCellIndex_(i,j,k);
}

// Assign deterministic dense IDs to all leaf blocks that participate in the
// calculation.  The SWMF/AMPS AMR tree is replicated across MPI ranks, therefore the
// same recursive traversal gives every process the same Temp_ID for the same leaf.
// Non-used leaves and non-leaf nodes are marked with Temp_ID=-1 so they cannot be
// accidentally inserted into the global cache.
void AssignGlobalLeafTempIds_(cAMRNode *node,long int& nUsedLeafBlocks) {
  if (node==NULL) return;

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->IsUsedInCalculationFlag==true) {
      node->Temp_ID = nUsedLeafBlocks++;
    }
    else {
      node->Temp_ID = -1;
    }

    return;
  }

  node->Temp_ID = -1;
  for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    AssignGlobalLeafTempIds_(node->downNode[nDownNode],nUsedLeafBlocks);
  }
}

// Collect the used leaves after IDs are assigned.  Keeping a vector avoids repeated
// full-tree traversals in the packing/allocation/population stages and keeps the
// subsequent code independent of the AMR-tree recursion details.
void CollectUsedLeafNodes_(cAMRNode *node,std::vector<cAMRNode*>& nodes) {
  if (node==NULL) return;

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((node->IsUsedInCalculationFlag==true) && (node->Temp_ID>=0)) nodes.push_back(node);
    return;
  }

  for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    CollectUsedLeafNodes_(node->downNode[nDownNode],nodes);
  }
}

// Ensure that every used AMR leaf has an allocated block on this MPI rank.
//
// In ordinary SWMF-coupled execution, node->block is allocated only for the local
// domain and boundary/ghost blocks.  The cutoff tracer is different: it may follow a
// particle through regions owned by other MPI ranks, so every rank needs block data
// for every used leaf.  Temporarily enabling AllowBlockAllocation lets us allocate
// those missing blocks while preserving the caller's original allocation policy.
long int AllocateMissingBlocks_(const std::vector<cAMRNode*>& nodes) {
  long int nAllocated=0;

  if (PIC::Mesh::mesh==NULL) return 0;

  // Respect the caller's original allocation policy.  We override it only inside
  // this helper because the global cutoff snapshot requires replicated blocks.
  const bool savedAllowBlockAllocation = PIC::Mesh::mesh->AllowBlockAllocation;
  PIC::Mesh::mesh->AllowBlockAllocation = true;

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode *node = *it;

    // Blocks owned by other ranks are normally absent here.  Allocate them so that
    // findTreeNode()+cell interpolation can work locally during backward tracing.
    if (node->block==NULL) {
      PIC::Mesh::mesh->AllocateBlock(node);
      if (node->block!=NULL) nAllocated++;
    }

    // Mirror the AMR-node global ID into the block for diagnostics and for any
    // future cell-level code that reaches the block before the node.
    if (node->block!=NULL) node->block->Temp_ID = node->Temp_ID;
  }

  PIC::Mesh::mesh->AllowBlockAllocation = savedAllowBlockAllocation;

  return nAllocated;
}

// MPI_Allreduce wrapper for large double arrays.
//
// MPI count arguments are int in the MPI-1/2 interface used here.  A large SWMF AMR
// mesh can exceed a safe single-call count, so the reduction is chunked.  The operation
// is SUM because only the owner rank contributes nonzero B for each interior cell.
void AllreduceDoubleVector_(const std::vector<double>& local,std::vector<double>& global) {
  global.assign(local.size(),0.0);

  const long int n = static_cast<long int>(local.size());
  const long int chunkMax = 100000000;

  for (long int offset=0;offset<n;offset+=chunkMax) {
    const long int chunk = std::min(chunkMax,n-offset);
    MPI_Allreduce((void*)(&local[offset]),&global[offset],static_cast<int>(chunk),
                  MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  }
}

// Same chunked reduction for the integer ownership mask.  The mask counts how many
// ranks supplied each global cell.  The normal value after reduction should be one;
// division by the mask makes the code robust to duplicate contributors, while the
// final completeness check still verifies that every interior cell was provided.
void AllreduceIntVector_(const std::vector<int>& local,std::vector<int>& global) {
  global.assign(local.size(),0);

  const long int n = static_cast<long int>(local.size());
  const long int chunkMax = 100000000;

  for (long int offset=0;offset<n;offset+=chunkMax) {
    const long int chunk = std::min(chunkMax,n-offset);
    MPI_Allreduce((void*)(&local[offset]),&global[offset],static_cast<int>(chunk),
                  MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  }
}

// Locate which interior cell of 'node' contains coordinate x.
//
// This is used mainly for ghost-cell reconstruction.  A ghost cell attached to one
// allocated block is not globally unique; its center lies inside a neighboring owner
// leaf.  We find that owner leaf first, then use this helper to convert the ghost-cell
// center into the owner's interior (i,j,k).  A small tolerance handles roundoff when a
// point lies exactly on a block boundary.
bool GetInteriorCellIndexForPoint_(const double *x,cAMRNode *node,int ijk[3]) {
  const int nCells[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  for (int idim=0;idim<3;idim++) {
    const double dx = (node->xmax[idim]-node->xmin[idim]) / static_cast<double>(nCells[idim]);

    if (!(dx>0.0)) return false;

    int idx = static_cast<int>(std::floor((x[idim]-node->xmin[idim])/dx));

    if (idx<0) {
      if (x[idim] < node->xmin[idim]-1.0e-10*dx) return false;
      idx = 0;
    }

    if (idx>=nCells[idim]) {
      if (x[idim] > node->xmax[idim]+1.0e-10*dx) return false;
      idx = nCells[idim]-1;
    }

    ijk[idim]=idx;
  }

  return true;
}

// Pack authoritative SWMF magnetic field values from this rank into a dense global
// array layout.
//
// Only owner blocks are packed (node->Thread == PIC::ThisThread).  Even if ghost or
// boundary copies exist on this rank, they are intentionally ignored as sources so that
// each global interior cell has a single authoritative contributor.  The B values are
// read from the same associated-data offset filled by PIC::CPLR::SWMF::RecieveCenterPointData().
long int PackOwnedInteriorMagneticField_(const std::vector<cAMRNode*>& nodes,
                                         std::vector<double>& localB,
                                         std::vector<int>& localMask) {
  long int nPacked=0;

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode *node = *it;

    // The SWMF coupler's authoritative interior values live on the block owner.
    // Do not pack ghost copies; otherwise the MPI sum would double count cells.
    if (node->Thread!=PIC::ThisThread) continue;
    if (node->block==NULL) continue;

    for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
          PIC::Mesh::cDataCenterNode *cell =
              node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell==NULL) continue;

          const long int c = GlobalCellIndex_(node,i,j,k);
          const long int b = 3*c;

          // Copy the three magnetic-field components from the SWMF-associated-data
          // buffer into the dense layout.  Units are unchanged: the coupler stores
          // B in the units expected by the existing AMPS field-evaluation path.
          memcpy(&localB[b],
                 cell->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset,
                 3*sizeof(double));

          // Mark this global cell as present.  After MPI_Allreduce, a value of 1
          // means exactly one owner provided the field; 0 means missing.
          localMask[c]=1;
          nPacked++;
        }
      }
    }
  }

  return nPacked;
}

// Write the MPI-reduced global magnetic field back into every allocated block.
//
// The loop includes both interior and ghost cells because AMPS interpolation stencils
// can request ghost-cell values near block boundaries.  For every cell center, we find
// the global owner leaf and copy the owner interior-cell B into the local cell's SWMF
// associated-data slot.  Cells whose centers cannot be mapped to a used owner leaf are
// set to zero and counted; these are normally outside the calculation domain and should
// not be touched by valid cutoff trajectories/interpolation stencils.
long int PopulateAllocatedBlocksFromGlobalMagneticField_(const std::vector<cAMRNode*>& nodes,
                                                        const std::vector<double>& globalB,
                                                        const std::vector<int>& globalMask) {
  long int nMissing=0;

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode *node = *it;

    if (node->block==NULL) continue;

    for (int i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (int j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++) {
        for (int k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          PIC::Mesh::cDataCenterNode *cell =
              node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell==NULL) continue;

          double *x = cell->GetX();

          // For interior cells, 'owner' is normally 'node'.  For ghost cells, this
          // maps the ghost-cell center to the neighboring AMR leaf that physically
          // owns that point.
          cAMRNode *owner = PIC::Mesh::mesh->findTreeNode(x,node);

          if ((owner==NULL) || (owner->Temp_ID<0) ||
              (owner->IsUsedInCalculationFlag==false)) {
            // The cell center is outside the used AMR domain.  Keep the data buffer
            // finite and count it for diagnostics; valid cutoff trajectories should
            // exit the domain before relying on such a value.
            double zero[3] = {0.0,0.0,0.0};
            memcpy(cell->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset,
                   zero,3*sizeof(double));
            nMissing++;
            continue;
          }

          int ijk[3];
          if (GetInteriorCellIndexForPoint_(x,owner,ijk)==false) {
            double zero[3] = {0.0,0.0,0.0};
            memcpy(cell->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset,
                   zero,3*sizeof(double));
            nMissing++;
            continue;
          }

          const long int c = GlobalCellIndex_(owner,ijk[0],ijk[1],ijk[2]);

          if ((c<0) || (c>=static_cast<long int>(globalMask.size())) ||
              (globalMask[c]<=0)) {
            double zero[3] = {0.0,0.0,0.0};
            memcpy(cell->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset,
                   zero,3*sizeof(double));
            nMissing++;
            continue;
          }

          double b[3];
          // Usually globalMask[c] == 1.  Averaging keeps the operation safe if the
          // same cell is contributed more than once by a future coupler variant.
          b[0]=globalB[3*c+0]/static_cast<double>(globalMask[c]);
          b[1]=globalB[3*c+1]/static_cast<double>(globalMask[c]);
          b[2]=globalB[3*c+2]/static_cast<double>(globalMask[c]);

          memcpy(cell->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset,
                 b,3*sizeof(double));
        }
      }
    }
  }

  return nMissing;
}

// Debug-field replacement for all currently allocated blocks.
//
// This is used by RedefineSWMFCoupledMagneticFieldToAnalyticDipole().  It is kept
// separate from PrepareGlobalSWMFCoupledMagneticFieldForCutoff(): the debug path writes
// B directly from an analytic callback, while the production path gathers the B that
// was imported from SWMF.  The callback signature matches PIC::CPLR::SWMF::Debug.
long int RedefineAllAllocatedMagneticField_(void (*f)(double*,double*)) {
  long int nCells=0;

  if (f==NULL) return 0;

  // Step 1: assign deterministic dense global IDs to used AMR leaves.
  long int nUsedLeafBlocks=0;
  AssignGlobalLeafTempIds_(PIC::Mesh::mesh->rootTree,nUsedLeafBlocks);

  // Step 2: collect the same ordered list of used leaves on every rank.
  std::vector<cAMRNode*> nodes;
  nodes.reserve(static_cast<size_t>(std::max(static_cast<long int>(0),nUsedLeafBlocks)));
  CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);

  double x[3],b[3];

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode *node = *it;
    if (node->block==NULL) continue;

    for (int i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (int j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++) {
        for (int k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          PIC::Mesh::cDataCenterNode *cell =
              node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell==NULL) continue;

          cell->GetX(x);
          f(x,b);
          memcpy(cell->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset,
                 b,3*sizeof(double));
          nCells++;
        }
      }
    }
  }

  return nCells;
}

} // anonymous namespace

//======================================================================================
// amps_pre_init  —  MUST be called before amps_init_mesh()
//======================================================================================
//
// Performs the subset of 3d_forward initialization that must precede the AMPS
// mesh-build phase.  main_lib.cpp::amps_init_mesh() calls this at its very start
// under a #if _PIC_COUPLER_MODE__SWMF_ guard.
//
// After this function returns:
//   • Earth::ModelMode == BoundaryInjectionMode
//     → amps_init_mesh() routes to MAIN_LIB_GEO::amps_init_mesh(), which calls
//       Earth::Init_AfterParser() (fixes Gap 3).
//   • cDensity3D::nEnergyBins/Emin_J/Emax_J reflect the input file
//     → initCellSamplingDataBuffer() allocates correctly-sized per-cell buffers
//       (fixes Gap 2).
//   • BoundingBoxInjection::prm is set
//     → main_lib.cpp::amps_init() can call InitDirectionIMF() correctly
//       (fixes Gaps 4 and 6).
//
void amps_pre_init() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  const char* inputFile = "AMPS_PARAM.in";

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForwardSWMF] amps_pre_init: reading '" << inputFile << "'.\n";

  // Parse once.  amps_init() will use s_prm directly — no second read.
  s_prm = EarthUtil::ParseAmpsParamFile(inputFile);

  // Mark the parameter object as SWMF-coupled: forces field access through
  // PIC::CPLR instead of Tsyganenko/DATAFILE in the shared 3-D helpers.
  s_prm.field.model                       = "SWMF";
  s_prm.calc.fieldEvalMethod              = "SWMF";
  s_prm.mode3d.forceAnalyticMagneticField = false;

  if (IsCutoffTarget_(s_prm)) {
    // Coupled cutoff-rigidity mode.  Keep the regular AMPS/SWMF mesh path, but
    // route per-coupling-step work to RunCutoffRigidity() instead of the
    // 3d_forward injector.  The cutoff tracer will sample the current SWMF
    // fields through Earth::Mode3D::EvaluateBackgroundMagneticFieldSI().
    Earth::ModelMode = Earth::CutoffRigidityMode;
    ApplyParsedDomainForCutoff_(s_prm);

    s_pre_initialized = true;

    if (PIC::ThisThread == 0)
      std::cout << "[Mode3DForwardSWMF] amps_pre_init complete:"
                << " ModelMode=CutoffRigidityMode"
                << ", field source=PIC::CPLR/SWMF"
                << ", output will be timestamped per amps_time_step().\n";
    return;
  }

  // Historical coupled 3d_forward mode.
  // ── Fix Gap 1 ────────────────────────────────────────────────────────────
  // Set Earth::ModelMode before amps_init_mesh() is called.
  // main_lib.cpp::amps_init_mesh() checks this at its entry:
  //   if (Earth::ModelMode == BoundaryInjectionMode)  →  MAIN_LIB_GEO path
  //   else                                            →  cutoff-rigidity path
  // The cutoff-rigidity path misses Earth::Init_AfterParser() and allocates
  // GeospaceFlag cell data that 3d_forward does not use.
  Earth::ModelMode = Earth::BoundaryInjectionMode;

  // ── Fix Gap 2 ────────────────────────────────────────────────────────────
  // Configure the cDensity3D energy grid before amps_init_mesh() calls
  // PIC::Mesh::initCellSamplingDataBuffer().  That function queries
  // cDensity3D::RequestSamplingData (already pushed into
  // PIC::IndividualModelSampling::RequestSamplingData by the first line of
  // amps_init_mesh()) to calculate each cell's sampling-buffer byte count:
  //     bytes = nEnergyBins * sizeof(double)
  // If nEnergyBins is still at its compile-time default here, every cell's
  // buffer will be under-allocated.
  Earth::Mode3DForward::cDensity3D::ConfigureEnergyGrid(s_prm);

  // ── Fix Gaps 4 & 6 ───────────────────────────────────────────────────────
  // Register prm with the bounding-box injector before amps_init_mesh() and
  // amps_init() run.  main_lib.cpp::amps_init() calls
  // Earth::BoundingBoxInjection::InitDirectionIMF() in the BoundaryInjectionMode
  // branch (now reachable because of the Gap-1 fix above), and InitDirectionIMF()
  // needs a valid prm to evaluate the background field at the domain boundary.
  // In the SWMF build ConfigureBackgroundFieldModel() is a no-op, so
  // InitDirectionIMF() falls through to the PIC::CPLR path or the safe default.
  Earth::BoundingBoxInjection::SetPrm(s_prm);

  s_pre_initialized = true;

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForwardSWMF] amps_pre_init complete:"
              << " ModelMode=BoundaryInjectionMode"
              << ", nEnergyBins=" << Earth::Mode3DForward::cDensity3D::nEnergyBins
              << ".\n";
#endif
}


//======================================================================================
// amps_init  —  called at the END of main_lib.cpp::amps_init()
//======================================================================================
//
// Completes the SWMF-coupled 3d_forward runtime after the full AMPS mesh and PIC
// core have been initialized.  Uses s_prm cached by amps_pre_init() — AMPS_PARAM.in
// is not re-read.
//
void amps_init() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  static bool initialized = false;
  if (initialized) {
    if (PIC::ThisThread == 0)
      std::cout << "[Mode3DForwardSWMF] Coupled 3d_forward runtime was already "
                   "initialized; skipping duplicate initialization.\n";
    return;
  }
  initialized = true;

  // Guard: amps_pre_init() must have run first.
  if (!s_pre_initialized)
    exit(__LINE__, __FILE__,
         "[Mode3DForwardSWMF] amps_init() called before amps_pre_init(). "
         "Ensure amps_pre_init() is called at the start of amps_init_mesh().");

  // Use the prm already parsed in amps_pre_init() — do not re-read the file.
  const EarthUtil::AmpsParam& prm = s_prm;

  if (IsCutoffTarget_(prm)) {
    // Coupled cutoff-rigidity mode needs only the generic AMPS/PIC
    // initialization already performed in main_lib.cpp before this hook is
    // called.  Do not register forward-injection callbacks.  The actual cutoff
    // calculation is launched from amps_cutoff_time_step() each time SWMF calls
    // amps_time_step(), after the current MHD fields have been imported.
    Earth::ModelMode = Earth::CutoffRigidityMode;
    ApplyParsedDomainForCutoff_(prm);

    for (int s=0; s<PIC::nTotalSpecies; ++s) {
      PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(s,1.0);
    }

    PIC::Mover::BackwardTimeIntegrationMode = _PIC_MODE_OFF_;

    if (PIC::ThisThread == 0)
      std::cout << "[Mode3DForwardSWMF] Initialized SWMF-coupled cutoff-rigidity runtime. "
                << "Each amps_time_step() call will write a timestamped cutoff snapshot.\n";
    return;
  }

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForwardSWMF] amps_init: completing coupled 3d_forward setup.\n";

  // Route AMPS particle motion through Earth::ParticleMover() →
  // Earth3DForward::MoverManager().  Earth::ModelMode was already set to
  // BoundaryInjectionMode in amps_pre_init(); confirm it is still correct.
  Earth::ModelMode = Earth::BoundaryInjectionMode;

  if (!Earth::Earth3DForward::SetMoverByName(prm.mode3dForward.particleMover)) {
    throw std::runtime_error(
        "Unknown 3d_forward_swmf particle mover '" +
        prm.mode3dForward.particleMover +
        "'. Valid movers are BORIS, RK4, GC/GC4, and HYBRID.");
  }

  Earth::Earth3DForward::ConfigureFieldEvaluation(prm);

  Earth::Mode3DForward::sInjectionEnergyDistribution =
      Earth::Mode3DForward::ParseInjectionEnergyDistribution(
          prm.mode3dForward.injectionEnergyDistribution);

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  if (Earth::Mode3DForward::sInjectionEnergyDistribution ==
      Earth::Mode3DForward::InjectionEnergyDistribution::LOG_UNIFORM) {
    throw std::runtime_error(
        "3d_forward_swmf LOG_UNIFORM injection requires AMPS to be compiled with "
        "_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_.");
  }
#endif

  Earth::Mode3DForward::sInitializeParticleTrajectories =
      prm.mode3dForward.initializeParticleTrajectories;
  Earth::Mode3DForward::sMaxParticleTrajectories = prm.mode3dForward.nParticleTrajectories;

  // Register the 3d_forward source-rate and particle-injection callbacks.
  // These override any BoundaryInjectionMode callbacks set by main_lib.cpp::amps_init()
  // (LocalBlockInjectionRate / InjectionProcessor) with the correct 3d_forward ones.
  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate =
      Earth::Mode3DForward::BoundaryInjectionSourceRate;
  PIC::BC::UserDefinedParticleInjectionFunction =
      Earth::Mode3DForward::InjectParticles;

  // In an SWMF-coupled build this call is a deliberate no-op: it leaves the field
  // source neutral so all field access goes through PIC::CPLR/SWMF.
  Earth::Mode3DForward::ConfigureBackgroundFieldModel(prm);

  // The 3d_forward injector does not use the IMF direction vector b to construct
  // particle velocities (it samples the full cosine-weighted inward hemisphere per
  // face).  Set a safe default so diagnostic messages that print b are well-defined.
  // Note: BoundingBoxInjection::SetPrm() was already called in amps_pre_init() and
  // InitDirectionIMF() was called by main_lib.cpp::amps_init() for BoundaryInjectionMode.
  Earth::BoundingBoxInjection::b[0] = 1.0;
  Earth::BoundingBoxInjection::b[1] = 0.0;
  Earth::BoundingBoxInjection::b[2] = 0.0;

#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  Earth::ParticleTracker::ConfigureRuntimeTrajectoryTracking(
      prm.mode3dForward.initializeParticleTrajectories,
      prm.mode3dForward.nParticleTrajectories,
      true);
#else
  if (Earth::Mode3DForward::sInitializeParticleTrajectories && PIC::ThisThread == 0) {
    std::cerr << "[Mode3DForwardSWMF] WARNING: trajectory initialization was "
              << "requested, but AMPS was compiled with "
              << "_PIC_PARTICLE_TRACKER_MODE_ != _PIC_MODE_ON_. No trajectory "
              << "records will be initialized.\n";
  }
#endif

  // cDensity3D::ConfigureEnergyGrid() was already called in amps_pre_init().
  // Calling it again here would be harmless but wasteful; the values are unchanged.

  // Enable the AMPS per-particle sampling dispatcher and register the 3d_forward
  // density sampler.  The duplicate-check handles the case where the generic
  // startup path already inserted the callback (cutoff-rigidity amps_init_mesh path).
  Earth::Sampling::ParticleData::SamplingMode = true;
  bool densitySamplerAlreadyRegistered = false;
  for (auto fn : Earth::Sampling::ParticleData::SampleParticleDataCallbacks) {
    if (fn == Mode3DForward::cDensity3D::SampleParticleData) {
      densitySamplerAlreadyRegistered = true;
      break;
    }
  }
  if (!densitySamplerAlreadyRegistered) {
    Earth::Sampling::ParticleData::SampleParticleDataCallbacks.push_back(
        Mode3DForward::cDensity3D::SampleParticleData);
  }

  // Initialize the spectrum before evaluating the source rate and particle weight.
  // 3d_forward convention: #SPECTRUM supplies the shape/normalization and
  // #DENSITY_3D supplies the actual simulated energy range.
  const auto forwardSpectrum =
      Earth::Mode3DForward::BuildForwardInjectionSpectrumMap(prm);
  InitGlobalSpectrumFromKeyValueMap(forwardSpectrum);

  Earth::Mode3DForward::sSpecies = 0;
  Earth::Mode3DForward::sDt      = Earth::Mode3DForward::EvaluateTimeStep(prm);
  PIC::ParticleWeightTimeStep::GlobalTimeStep[Earth::Mode3DForward::sSpecies] =
      Earth::Mode3DForward::sDt;

  Earth::Mode3DForward::sNParticlesPerIter = prm.mode3dForward.nParticlesPerIter;
  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber =
      Earth::Mode3DForward::sNParticlesPerIter;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(
      Earth::Mode3DForward::sSpecies);
  Earth::Mode3DForward::sParticleWeight =
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[Earth::Mode3DForward::sSpecies];

  Earth::Mode3DForward::InitAbsorptionSphere(prm);

  // ── Fix Gap 5 ────────────────────────────────────────────────────────────
  // The cutoff-rigidity path in main_lib.cpp::amps_init_mesh() sets
  //   PIC::Mover::ProcessOutsideDomainParticles =
  //       Earth::CutoffRigidity::ProcessOutsideDomainParticles
  // With the amps_pre_init() fix, amps_init_mesh() now takes the
  // MAIN_LIB_GEO path, which does NOT set this pointer (correct for
  // 3d_forward).  This reset is kept as an explicit guard in case another
  // code path sets it between amps_pre_init() and here.
  PIC::Mover::ProcessOutsideDomainParticles = nullptr;

  Mode3DForward::cDensity3D::Init(prm);
  Mode3DForward::cSphereFlux3D::Init(
      prm, Earth::Mode3DForward::sAbsorptionSphere, Earth::Mode3DForward::sDt);
  Earth::Mode3DForward::InitBoundaryInjectionTable();

  PIC::Mover::BackwardTimeIntegrationMode = _PIC_MODE_OFF_;
  PIC::SamplingMode = _RESTART_SAMPLING_MODE_;

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] Initialized SWMF-coupled 3d_forward runtime. "
              << "Field source=PIC::CPLR/SWMF"
              << ", mover="       << Earth::Earth3DForward::GetMoverName()
              << ", nParticlesPerIter=" << Earth::Mode3DForward::sNParticlesPerIter
              << ", particleWeight="    << Earth::Mode3DForward::sParticleWeight
              << ", dt="               << Earth::Mode3DForward::sDt
              << ", energySampling="
              << Earth::Mode3DForward::InjectionEnergyDistributionName(
                     Earth::Mode3DForward::sInjectionEnergyDistribution)
              << ", nEnergyBins=" << Mode3DForward::cDensity3D::nEnergyBins
              << "\n";
  }
#endif
}


bool IsForwardMode() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return false;
#else
  return s_pre_initialized && IsForwardTarget_(s_prm);
#endif
}

bool IsCutoffRigidityMode() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return false;
#else
  return s_pre_initialized && IsCutoffTarget_(s_prm);
#endif
}

// Public preparation hook called immediately before the SWMF-coupled cutoff solver.
//
// This routine intentionally performs all expensive globalisation work once per cutoff
// snapshot, not inside the particle mover.  After it returns, field access during the
// cutoff calculation is a local memory lookup/interpolation on every MPI rank.
void PrepareGlobalSWMFCoupledMagneticFieldForCutoff(bool verbose) {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  if (!s_pre_initialized) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before amps_pre_init().");
  }

  if (PIC::Mesh::mesh==NULL) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before PIC::Mesh::mesh is initialized.");
  }

  if (PIC::CPLR::SWMF::MagneticFieldOffset<0) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before the SWMF magnetic-field buffer is allocated.");
  }

  if (PIC::CPLR::SWMF::FirstCouplingOccured==false) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before the first SWMF coupling data receive.");
  }

  // Step 1: assign deterministic dense global IDs to used AMR leaves.
  long int nUsedLeafBlocks=0;
  AssignGlobalLeafTempIds_(PIC::Mesh::mesh->rootTree,nUsedLeafBlocks);

  // Step 2: collect the same ordered list of used leaves on every rank.
  std::vector<cAMRNode*> nodes;
  nodes.reserve(static_cast<size_t>(std::max(static_cast<long int>(0),nUsedLeafBlocks)));
  CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);

  // Dense global storage contains only interior cells.  Ghost-cell values are
  // regenerated from these interior cells after the MPI gather is complete.
  const long int nInteriorCells = nUsedLeafBlocks * InteriorCellsPerBlock_();

  std::vector<double> localB(static_cast<size_t>(3*nInteriorCells),0.0);
  std::vector<double> globalB;
  std::vector<int> localMask(static_cast<size_t>(nInteriorCells),0);
  std::vector<int> globalMask;

  // Step 3: allocate missing nonlocal blocks so every rank has the full AMR mesh
  // represented by actual cDataBlockAMR objects.
  const long int nNewlyAllocatedBlocks = AllocateMissingBlocks_(nodes);

  // Step 4: pack this rank's authoritative owner cells into the dense local cache.
  const long int nPackedLocal = PackOwnedInteriorMagneticField_(nodes,localB,localMask);

  // Step 5: make the dense cell-centered B cache globally available on every rank.
  AllreduceDoubleVector_(localB,globalB);
  AllreduceIntVector_(localMask,globalMask);

  // Step 6: write the global B values back into all allocated cells, including
  // ghost cells, so existing AMPS interpolation routines can be used unchanged.
  const long int nMissingLocal =
      PopulateAllocatedBlocksFromGlobalMagneticField_(nodes,globalB,globalMask);

  long int nPackedGlobal=0;
  MPI_Allreduce((void*)&nPackedLocal,&nPackedGlobal,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  long int nNewlyAllocatedMin=0,nNewlyAllocatedMax=0;
  MPI_Allreduce((void*)&nNewlyAllocatedBlocks,&nNewlyAllocatedMin,1,MPI_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&nNewlyAllocatedBlocks,&nNewlyAllocatedMax,1,MPI_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  long int nMissingMin=0,nMissingMax=0;
  MPI_Allreduce((void*)&nMissingLocal,&nMissingMin,1,MPI_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&nMissingLocal,&nMissingMax,1,MPI_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  // A complete coupled snapshot must supply every used leaf's interior-cell B.
  // If not, cutoff tracing would silently encounter uninitialized field values, so
  // fail early with an explicit diagnostic.
  if (nPackedGlobal!=nInteriorCells) {
    std::ostringstream msg;
    msg << "[Mode3DForwardSWMF] Global SWMF magnetic-field gather is incomplete: "
        << "received " << nPackedGlobal << " owner interior cells, expected "
        << nInteriorCells << ".";
    const std::string text = msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  // Print a compact summary from rank 0.  The min..max ranges show whether all ranks
  // allocated the same number of missing blocks and whether any ranks had ghost cells
  // outside the used AMR domain.
  if ((verbose==true) && (PIC::ThisThread==0)) {
    std::cout << "[Mode3DForwardSWMF] Prepared global SWMF B field for cutoff:"
              << " usedLeafBlocks=" << nUsedLeafBlocks
              << ", ownerInteriorCells=" << nPackedGlobal
              << ", newlyAllocatedBlocksPerRank=" << nNewlyAllocatedMin;

    if (nNewlyAllocatedMax!=nNewlyAllocatedMin) {
      std::cout << ".." << nNewlyAllocatedMax;
    }

    std::cout << ", missingGhostCellsPerRank=" << nMissingMin;

    if (nMissingMax!=nMissingMin) {
      std::cout << ".." << nMissingMax;
    }

    std::cout << ".\n";
    std::cout.flush();
  }
#endif
}

// Optional debug override.  It is deliberately not called by default in the coupled
// cutoff path: production SWMF-coupled cutoff now uses
// PrepareGlobalSWMFCoupledMagneticFieldForCutoff().  Call this manually when a known
// analytic dipole is desired for debugging the cutoff machinery independently of SWMF.
void RedefineSWMFCoupledMagneticFieldToAnalyticDipole() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  if (!s_pre_initialized) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] RedefineSWMFCoupledMagneticFieldToAnalyticDipole() called before amps_pre_init().");
  }

  if (PIC::Mesh::mesh==NULL) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] RedefineSWMFCoupledMagneticFieldToAnalyticDipole() called before PIC::Mesh::mesh is initialized.");
  }

  if (PIC::CPLR::SWMF::MagneticFieldOffset<0) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] RedefineSWMFCoupledMagneticFieldToAnalyticDipole() called before the SWMF magnetic-field buffer is allocated.");
  }

  Earth::GridlessMode::Dipole::SetMomentScale(s_prm.field.dipoleMoment_Me);
  Earth::GridlessMode::Dipole::SetTiltDeg(s_prm.field.dipoleTilt_deg);

  const long int nCells = RedefineAllAllocatedMagneticField_(AnalyticDipoleMagneticField_);

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] Replaced SWMF-coupled cell-centered B field "
              << "with analytic dipole:"
              << " DIPOLE_MOMENT=" << s_prm.field.dipoleMoment_Me
              << ", DIPOLE_TILT=" << s_prm.field.dipoleTilt_deg << " deg"
              << ", updatedCellsOnThisRank=" << nCells << ".\n";
    std::cout.flush();
  }
#endif
}

void amps_cutoff_time_step() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  if (!s_pre_initialized) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] amps_cutoff_time_step() called before amps_pre_init().");
  }

  if (!IsCutoffTarget_(s_prm)) return;

  const long int callIndex = s_cutoff_call_counter;
  const double   tSim_s    = s_cutoff_elapsed_time_s;
  const std::string suffix = FormatCutoffOutputSuffix_(callIndex,tSim_s);

  Earth::Mode3D::SetCutoffOutputFileSuffix(suffix);

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] cutoff snapshot " << callIndex
              << ": t_sim=" << tSim_s << " s, suffix='" << suffix << "'.\n";
    std::cout.flush();
  }

  try {
    // The cutoff solver can run trajectories across the entire AMR domain on each
    // MPI rank.  Materialize the distributed SWMF B field into replicated cell data
    // before starting any backward tracing.
    PrepareGlobalSWMFCoupledMagneticFieldForCutoff(true);

    // Run the Mode3D cutoff solver with progress reporting enabled for the coupled
    // snapshot.  RunCutoffRigidity() defaults to no progress bar for other callers.
    Earth::Mode3D::RunCutoffRigidity(s_prm,true);
  }
  catch (const std::exception& e) {
    std::ostringstream msg;
    msg << "[Mode3DForwardSWMF] SWMF-coupled cutoff-rigidity calculation failed: "
        << e.what();
    const std::string text = msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  ++s_cutoff_call_counter;
  s_cutoff_elapsed_time_s += CoupledCutoffOutputDtSeconds_(s_prm);
#endif
}

} // namespace Mode3DForwardSWMF
} // namespace Earth
