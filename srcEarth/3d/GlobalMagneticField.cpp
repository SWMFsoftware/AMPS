//======================================================================================
// GlobalMagneticField.cpp
//======================================================================================
//
// Generic materialization of a distributed AMPS cell-centered magnetic field into a
// replicated, read-only field snapshot for 3-D backward cutoff-rigidity tracing.
//
// The implementation is deliberately independent of the physical field source.  The
// caller passes the byte offset of the three B components in each cell-associated data
// buffer.  Therefore the same data-management machinery serves
//
//   * standalone 3-D cutoff, where Mode3D::InitMeshFields() writes DATAFILE B; and
//   * SWMF-coupled cutoff, where PIC::CPLR::SWMF writes coupled B.
//
// Only the interior cells of owner blocks are used as MPI sources.  Ghost cells are
// reconstructed after the gather by mapping each ghost-cell center to the owner AMR
// leaf that physically contains the point.  This avoids double counting ghost copies
// and gives AMPS' existing interpolation stencils the ghost-cell values they expect.
//======================================================================================

#include "GlobalMagneticField.h"

#include "pic.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

namespace Earth {
namespace Mode3D {
namespace GlobalMagneticField {
namespace {

// Short alias used throughout this file.  The helper works on the AMR nodes owned by
// the global PIC mesh instance; it does not create or modify the AMR topology itself.
typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cAMRNode;

// Return a safe diagnostic tag.  Avoiding a null pointer in log/error construction
// keeps the materialization helper usable from defensive/debug call sites.
std::string SafeTag_(const char* diagnosticTag) {
  return (diagnosticTag!=NULL && diagnosticTag[0]!='\0') ?
         std::string(diagnosticTag) : std::string("Mode3D::GlobalMagneticField");
}

// Number of physical/interior cell centers stored per AMR block.  The global dense
// cache indexes only these cells; ghost cells are filled later from the owner cells.
long int InteriorCellsPerBlock_() {
  return static_cast<long int>(_BLOCK_CELLS_X_) *
         static_cast<long int>(_BLOCK_CELLS_Y_) *
         static_cast<long int>(_BLOCK_CELLS_Z_);
}

// Convert a block-local interior index (i,j,k) into a compact row-major scalar index.
// The exact order is not physically important, but it must be deterministic and the
// same on every rank because the arrays are reduced with MPI_Allreduce.
long int InteriorCellIndex_(int i,int j,int k) {
  return static_cast<long int>(i) +
         static_cast<long int>(_BLOCK_CELLS_X_) *
         (static_cast<long int>(j) +
          static_cast<long int>(_BLOCK_CELLS_Y_) * static_cast<long int>(k));
}

// Dense global cell index.  node->Temp_ID is assigned by AssignGlobalLeafTempIds_()
// immediately before the gather; it is intentionally treated as temporary scratch
// state for this materialization pass.
long int GlobalCellIndex_(cAMRNode *node,int i,int j,int k) {
  return node->Temp_ID * InteriorCellsPerBlock_() + InteriorCellIndex_(i,j,k);
}

// Assign deterministic dense IDs to all used AMR leaf blocks.
//
// AMPS stores the AMR tree topology on every MPI rank even when the data blocks are
// distributed.  A recursive traversal in a fixed child order therefore assigns the
// same Temp_ID to the same physical leaf on every process.  These IDs become the block
// component of the dense global cache index.
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

  // Non-leaf nodes must never be used as dense-cache entries.
  node->Temp_ID = -1;

  for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    AssignGlobalLeafTempIds_(node->downNode[nDownNode],nUsedLeafBlocks);
  }
}

// Collect the used leaves after IDs are assigned.  Keeping an ordered vector avoids
// repeated full-tree traversals and guarantees all later loops use the same order on
// all ranks.
void CollectUsedLeafNodes_(cAMRNode *node,std::vector<cAMRNode*>& nodes) {
  if (node==NULL) return;

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((node->IsUsedInCalculationFlag==true) && (node->Temp_ID>=0)) {
      nodes.push_back(node);
    }
    return;
  }

  for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    CollectUsedLeafNodes_(node->downNode[nDownNode],nodes);
  }
}

// Ensure every used AMR leaf has an allocated cDataBlockAMR on this rank.
//
// In normal AMPS distributed-domain mode only local owner blocks and a limited set of
// neighbor/ghost blocks are allocated.  Backward cutoff tracing is nonlocal: one rank
// may trace a particle through any block while processing its assigned observation
// point.  Therefore we allocate missing blocks after mesh generation, just before the
// cutoff run, rather than forcing independent full-domain MPI initialization.
long int AllocateMissingBlocks_(const std::vector<cAMRNode*>& nodes) {
  long int nAllocated=0;

  if (PIC::Mesh::mesh==NULL) return 0;

  // Preserve the caller's allocation policy.  Some AMPS code intentionally disables
  // block allocation while manipulating only the tree; we override that policy only
  // within this helper because a replicated cutoff snapshot explicitly requires it.
  const bool savedAllowBlockAllocation = PIC::Mesh::mesh->AllowBlockAllocation;
  PIC::Mesh::mesh->AllowBlockAllocation = true;

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode *node = *it;

    if (node->block==NULL) {
      PIC::Mesh::mesh->AllocateBlock(node);
      if (node->block!=NULL) nAllocated++;
    }

    // Mirror the temporary global ID into the block.  The interpolation code uses the
    // node pointer, but this is useful for diagnostics and for future code that starts
    // from block-level data.
    if (node->block!=NULL) node->block->Temp_ID = node->Temp_ID;
  }

  PIC::Mesh::mesh->AllowBlockAllocation = savedAllowBlockAllocation;

  return nAllocated;
}

// Chunked MPI_Allreduce for large double vectors.  MPI count arguments are int in the
// MPI interface used by AMPS; chunking keeps very large AMR caches below that limit.
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

// Chunked MPI_Allreduce for integer masks.  A reduced mask value of 1 is the normal
// case: exactly one owner rank supplied that interior cell.  Values larger than one
// are averaged when writing B back; a zero triggers an explicit incomplete-gather
// failure before the cutoff calculation begins.
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

// Locate the interior cell of 'node' that contains coordinate x.
//
// This is needed for ghost-cell reconstruction.  A ghost cell belongs to the local
// allocated block object, but its center physically lies inside a neighboring owner
// leaf.  Once findTreeNode() identifies that owner leaf, this routine converts the
// coordinate into the owner's interior (i,j,k).  Small boundary tolerances absorb
// roundoff when x lies numerically on a cell or block face.
bool GetInteriorCellIndexForPoint_(const double *x,cAMRNode *node,int ijk[3]) {
  const int nCells[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  for (int idim=0;idim<3;idim++) {
    const double dx = (node->xmax[idim]-node->xmin[idim]) /
                      static_cast<double>(nCells[idim]);

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

// Pack authoritative magnetic-field values from owner blocks into the dense local
// contribution arrays.
//
// Only node->Thread == PIC::ThisThread contributes.  Non-owner blocks may exist on a
// rank either because AMPS created ghost/boundary blocks or because this helper has
// just allocated a full replicated mesh.  Packing them would double-count cells and
// could use stale/uninitialized ghost data, so they are deliberately ignored.
long int PackOwnedInteriorMagneticField_(const std::vector<cAMRNode*>& nodes,
                                         long int magneticFieldDataOffset,
                                         std::vector<double>& localB,
                                         std::vector<int>& localMask) {
  long int nPacked=0;

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode *node = *it;

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

          // The caller supplied the physical storage location.  In standalone Mode3D
          // this points at DATAFILE B; in SWMF-coupled cutoff it points at SWMF B.
          std::memcpy(&localB[b],
                      cell->GetAssociatedDataBufferPointer()+magneticFieldDataOffset,
                      3*sizeof(double));

          localMask[c]=1;
          nPacked++;
        }
      }
    }
  }

  return nPacked;
}

// Write the MPI-reduced magnetic field back into every allocated block, including
// ghost cells.  Ghost cells are reconstructed by asking the global AMR tree which
// owner leaf contains the ghost-cell center and copying the corresponding owner
// interior-cell value from the global cache.
long int PopulateAllocatedBlocksFromGlobalMagneticField_(
    const std::vector<cAMRNode*>& nodes,
    long int magneticFieldDataOffset,
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

          // For interior cells this usually returns 'node'.  For ghost cells it maps
          // the ghost-center coordinate to the neighboring owner leaf.
          cAMRNode *owner = PIC::Mesh::mesh->findTreeNode(x,node);

          if ((owner==NULL) || (owner->Temp_ID<0) ||
              (owner->IsUsedInCalculationFlag==false)) {
            // The cell center lies outside the used AMR domain.  Store a finite value
            // and count it.  Valid cutoff trajectories should leave the box before
            // depending on such a ghost value.
            double zero[3] = {0.0,0.0,0.0};
            std::memcpy(cell->GetAssociatedDataBufferPointer()+magneticFieldDataOffset,
                        zero,3*sizeof(double));
            nMissing++;
            continue;
          }

          int ijk[3];
          if (GetInteriorCellIndexForPoint_(x,owner,ijk)==false) {
            double zero[3] = {0.0,0.0,0.0};
            std::memcpy(cell->GetAssociatedDataBufferPointer()+magneticFieldDataOffset,
                        zero,3*sizeof(double));
            nMissing++;
            continue;
          }

          const long int c = GlobalCellIndex_(owner,ijk[0],ijk[1],ijk[2]);

          if ((c<0) || (c>=static_cast<long int>(globalMask.size())) ||
              (globalMask[c]<=0)) {
            double zero[3] = {0.0,0.0,0.0};
            std::memcpy(cell->GetAssociatedDataBufferPointer()+magneticFieldDataOffset,
                        zero,3*sizeof(double));
            nMissing++;
            continue;
          }

          double b[3];
          b[0]=globalB[3*c+0]/static_cast<double>(globalMask[c]);
          b[1]=globalB[3*c+1]/static_cast<double>(globalMask[c]);
          b[2]=globalB[3*c+2]/static_cast<double>(globalMask[c]);

          std::memcpy(cell->GetAssociatedDataBufferPointer()+magneticFieldDataOffset,
                      b,3*sizeof(double));
        }
      }
    }
  }

  return nMissing;
}

} // anonymous namespace

long int DataFileMagneticFieldDataOffset() {
  return PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin +
         PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset +
         PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
}

MaterializationStats MaterializeCellCenteredMagneticFieldForCutoff(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    bool verbose) {

  const std::string tag = SafeTag_(diagnosticTag);

  if (PIC::Mesh::mesh==NULL) {
    const std::string msg = "[" + tag + "] global magnetic-field materialization called before PIC::Mesh::mesh is initialized.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  if (PIC::Mesh::mesh->rootTree==NULL) {
    const std::string msg = "[" + tag + "] global magnetic-field materialization called before the AMR root tree exists.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  if (magneticFieldDataOffset<0) {
    const std::string msg = "[" + tag + "] global magnetic-field materialization received a negative B-field data offset.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  MaterializationStats stats;

  // Step 1: give every used AMR leaf the same dense global ID on every rank.
  long int nUsedLeafBlocks=0;
  AssignGlobalLeafTempIds_(PIC::Mesh::mesh->rootTree,nUsedLeafBlocks);
  stats.usedLeafBlocks = nUsedLeafBlocks;

  // Step 2: build the ordered leaf list used by all subsequent passes.
  std::vector<cAMRNode*> nodes;
  nodes.reserve(static_cast<size_t>(std::max(static_cast<long int>(0),nUsedLeafBlocks)));
  CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);

  // Dense storage for owner interior-cell values.  There are no ghost cells here;
  // ghost values are derived after the MPI reduction.
  const long int nInteriorCells = nUsedLeafBlocks * InteriorCellsPerBlock_();
  stats.expectedInteriorCells = nInteriorCells;

  std::vector<double> localB(static_cast<size_t>(3*nInteriorCells),0.0);
  std::vector<double> globalB;
  std::vector<int> localMask(static_cast<size_t>(nInteriorCells),0);
  std::vector<int> globalMask;

  // Step 3: allocate nonlocal blocks on this rank.  This is what replaces the old
  // independentDomainMode=true startup: normal MPI domain decomposition is preserved
  // until this explicit replicated cutoff snapshot is needed.
  const long int nNewlyAllocatedBlocks = AllocateMissingBlocks_(nodes);

  // Step 4: pack this rank's authoritative owner blocks.
  const long int nPackedLocal =
      PackOwnedInteriorMagneticField_(nodes,magneticFieldDataOffset,localB,localMask);

  // Step 5: replicate the owner-cell cache onto every rank.
  AllreduceDoubleVector_(localB,globalB);
  AllreduceIntVector_(localMask,globalMask);

  // Step 6: populate every allocated block, including ghost cells, from the cache.
  const long int nMissingLocal =
      PopulateAllocatedBlocksFromGlobalMagneticField_(nodes,magneticFieldDataOffset,
                                                     globalB,globalMask);

  long int nPackedGlobal=0;
  MPI_Allreduce((void*)&nPackedLocal,&nPackedGlobal,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  stats.ownerInteriorCells = nPackedGlobal;

  MPI_Allreduce((void*)&nNewlyAllocatedBlocks,&stats.newlyAllocatedBlocksMin,1,MPI_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&nNewlyAllocatedBlocks,&stats.newlyAllocatedBlocksMax,1,MPI_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&nMissingLocal,&stats.missingGhostCellsMin,1,MPI_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&nMissingLocal,&stats.missingGhostCellsMax,1,MPI_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  // Every used interior cell must have exactly one owner contribution.  Without this
  // check a missing local field would silently become zero after the Allreduce and the
  // cutoff classification could be wrong.
  if (nPackedGlobal!=nInteriorCells) {
    std::ostringstream msg;
    msg << "[" << tag << "] global magnetic-field gather is incomplete: received "
        << nPackedGlobal << " owner interior cells, expected " << nInteriorCells
        << ". Check that the normal AMPS mesh distribution allocated owner blocks "
        << "before the cutoff materialization step.";
    const std::string text = msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  if ((verbose==true) && (PIC::ThisThread==0)) {
    std::cout << "[" << tag << "] Prepared global cell-centered B field for cutoff:"
              << " usedLeafBlocks=" << stats.usedLeafBlocks
              << ", ownerInteriorCells=" << stats.ownerInteriorCells
              << ", newlyAllocatedBlocksPerRank=" << stats.newlyAllocatedBlocksMin;

    if (stats.newlyAllocatedBlocksMax!=stats.newlyAllocatedBlocksMin) {
      std::cout << ".." << stats.newlyAllocatedBlocksMax;
    }

    std::cout << ", missingGhostCellsPerRank=" << stats.missingGhostCellsMin;

    if (stats.missingGhostCellsMax!=stats.missingGhostCellsMin) {
      std::cout << ".." << stats.missingGhostCellsMax;
    }

    std::cout << ".\n";
    std::cout.flush();
  }

  return stats;
}

long int RedefineAllAllocatedMagneticField(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    void (*fieldCallback)(double*,double*),
    bool allocateMissingBlocks,
    bool verbose) {

  const std::string tag = SafeTag_(diagnosticTag);

  if (fieldCallback==NULL) return 0;

  if (PIC::Mesh::mesh==NULL || PIC::Mesh::mesh->rootTree==NULL) {
    const std::string msg = "[" + tag + "] magnetic-field redefinition called before the AMR mesh is initialized.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  if (magneticFieldDataOffset<0) {
    const std::string msg = "[" + tag + "] magnetic-field redefinition received a negative B-field data offset.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  long int nUsedLeafBlocks=0;
  AssignGlobalLeafTempIds_(PIC::Mesh::mesh->rootTree,nUsedLeafBlocks);

  std::vector<cAMRNode*> nodes;
  nodes.reserve(static_cast<size_t>(std::max(static_cast<long int>(0),nUsedLeafBlocks)));
  CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);

  long int nNewlyAllocatedBlocks=0;
  if (allocateMissingBlocks==true) {
    nNewlyAllocatedBlocks = AllocateMissingBlocks_(nodes);
  }

  long int nCells=0;
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
          fieldCallback(x,b);
          std::memcpy(cell->GetAssociatedDataBufferPointer()+magneticFieldDataOffset,
                      b,3*sizeof(double));
          nCells++;
        }
      }
    }
  }

  long int nCellsGlobal=0;
  long int nAllocatedMin=0,nAllocatedMax=0;
  MPI_Allreduce((void*)&nCells,&nCellsGlobal,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&nNewlyAllocatedBlocks,&nAllocatedMin,1,MPI_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&nNewlyAllocatedBlocks,&nAllocatedMax,1,MPI_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  if ((verbose==true) && (PIC::ThisThread==0)) {
    std::cout << "[" << tag << "] Redefined cell-centered B field on replicated cutoff mesh:"
              << " usedLeafBlocks=" << nUsedLeafBlocks
              << ", localCellsOnRank0=" << nCells
              << ", totalCellsOverRanks=" << nCellsGlobal
              << ", newlyAllocatedBlocksPerRank=" << nAllocatedMin;
    if (nAllocatedMax!=nAllocatedMin) std::cout << ".." << nAllocatedMax;
    std::cout << ".\n";
    std::cout.flush();
  }

  return nCells;
}

} // namespace GlobalMagneticField
} // namespace Mode3D
} // namespace Earth
