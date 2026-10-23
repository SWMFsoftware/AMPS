//$Id$
//interpolation routines

/*
 * pic_interpolation_routines.cpp
 *
 *  Created on: Jul 18, 2015
 *      Author: vtenishe
 */

/*
 * =====================================================================================
 * STENCIL-CELL CACHE OF THE CELL-CENTERED LINEAR INTERPOLATION (added 2026-10-22)
 * =====================================================================================
 *
 * 1. What is cached and why
 * -------------------------
 * CellCentered::Linear::InitStencil() builds the interpolation stencil of a point x in
 * block `node` in one of three ways:
 *   - trilinear on the 8 cell centers around x (GetTriliniarInterpolationStencil()),
 *     used inside a block and next to neighbors of the same resolution;
 *   - trilinear on the lattice of the coarser block, where each coarse support point
 *     that lies in a finer block is replaced by the average of its 2x2x2 fine cells
 *     (GetTriliniarInterpolationMutiBlockStencil()), used at coarse/fine interfaces;
 *   - a blend of the two within one cell of the interface (continuity).
 * It returns two representations of the stencil:
 *   - the pointer stencil (cStencil): cDataCenterNode pointers of the local block and
 *     its ghost cells (node->block->GetCenterNode()); used by all applications;
 *   - the optional row stencil (cRowStencil): decomposition-independent (leaf node,
 *     i,j,k) identifiers that do not need the block to be allocated on this rank;
 *     used by the Mode3D cutoff-rigidity tracing of srcEarth/srcGeospace to read a
 *     compact global field array.
 *
 * Building a stencil requires tree searches whose results depend only on the geometry
 * of the AMR tree, not on x within the cell, not on the field data and not on which
 * MPI rank owns a block:
 *   (a) row-stencil cells: every entry is given as cell (i,j,k) of a node with
 *       i,j,k in [-1,_BLOCK_CELLS_] (ghost indices included).  AddRowStencilCell()
 *       converts the cell center to the leaf that contains it with
 *       mesh->findTreeNode() and converts the coordinate to the cell index in that
 *       leaf (LocateRowStencilCell()).  A trilinear stencil has 8 such cells, a
 *       multiblock stencil up to 64, a blended stencil up to 72.
 *   (b) multiblock support points: GetTriliniarInterpolationMutiBlockStencil() calls
 *       mesh->findTreeNode() for each of its 8 support points (cell (i,j,k) of the
 *       coarse node) to find the block that covers it.
 *   (c) coarser-neighbor search: InitStencil() inspects up to 6 faces, 12 edges and 8
 *       corners of the block with node->GetNeibFace()/GetNeibEdge()/GetNeibCorner().
 *       With _MESH_GLOBAL_NODE_CONNECTION_INFO_MODE_==_OFF_AMR_MESH_ (the default)
 *       each of these is again an integer-index tree search.
 * A tree search climbs from the start node until the searched index lies inside the
 * node and then descends to the leaf, i.e. several dependent memory reads per level.
 * In Mode3D tracing every field evaluation builds a stencil, a trajectory makes
 * 1e3-1e5 field evaluations, and the same blocks are visited again and again: the
 * same searches are repeated millions of times with the same arguments.  A profile of
 * Mode3D MESH tracing (srcGeospace C1 test, port of this code) attributed about 55% of
 * the executed instructions to these searches and coordinate-to-index conversions,
 * while the arithmetic of the trajectory integrator was below 3%.
 *
 * The cache stores the results of (a), (b) and (c) for recently used nodes so that a
 * repeated request is answered by a table lookup.
 *
 * 2. Exactness
 * ------------
 * A cached value is computed by the original code at its first use and returned for
 * every later request with the same arguments:
 *   - (a) is filled only by AddRowStencilCell() and (b) only by
 *     GetTriliniarInterpolationMutiBlockStencil(); each caller computes the
 *     coordinate from (node,i,j,k) with its own expression, so a cached value is the
 *     result of exactly the computation it replaces;
 *   - state that is not tree geometry is NOT cached and is read at use time as before:
 *     IsUsedInCalculationFlag of the leaf found (re-tested on every use),
 *     RefinmentLevel and the coordinates of neighbors, node->block and the cell data
 *     (the pointer stencil still reads node->block->GetCenterNode()), all weights.
 * The stencils, their weights and the order in which entries are added and merged are
 * therefore unchanged and the interpolated values are bitwise identical with and
 * without the cache (verified: srcEarth and srcGeospace Mode3D MESH tests, 1 and 8
 * threads, cache on/off and original code; srcGeospace engine-vs-legacy comparison of
 * 8,888 trajectories, 0 mismatches).
 *
 * 3. Organization
 * ---------------
 *   - Per thread: every thread (OpenMP or std::thread) owns its table (thread_local),
 *     so there are no locks and no shared writes.
 *   - Direct-mapped table of nStencilCacheSlots=256 entries; the slot of a node is the
 *     upper 8 bits of a 64-bit multiplicative hash of the node address.  A request for
 *     another node mapping to the same slot evicts the previous entry.
 *   - One entry per node holds, for each of the (_BLOCK_CELLS_X_+2)(_BLOCK_CELLS_Y_+2)
 *     (_BLOCK_CELLS_Z_+2) cells (512 for 6x6x6 blocks): the row-cell owner and index
 *     (a), the support-point node (b), and flags telling which of them are already
 *     known; plus the 26 neighbors (c).  Values are filled lazily, only for the cells
 *     and neighbors actually requested.
 *   - Memory: about 10.5 KB per entry for 6x6x6 blocks, i.e. about 2.7 MB per thread
 *     that interpolates; allocated at the first use by the thread, released when the
 *     thread exits.  Independent of the mesh size.
 *   - Host code only: the cache code is excluded from CUDA device code (__CUDA_ARCH__),
 *     where the original searches are used.
 *   - A pointer to an entry is never kept across a request for another node (the slot
 *     may be reassigned in between).
 *
 * 4. What triggers re-caching
 * ---------------------------
 * An entry is valid for (mesh, node, mesh->nMeshModificationCounter).  It is cleared
 * and refilled lazily when
 *   - the slot is requested for another node (eviction);
 *   - PIC::Mesh::mesh points to another mesh object;
 *   - mesh->nMeshModificationCounter changed.  cMeshAMRgeneric increments it in
 *     splitTreeNode() (refinement), deleteTreeNode() (coarsening), AllocateBlock()
 *     and DeallocateBlock(), setMeshName()/generateMeshName() (mesh creation), and in
 *     CreateNewParallelDistributionLists() (every load-balancing redistribution).
 *     The cached values depend on none of block allocation and domain decomposition
 *     (the tree is replicated on all ranks; allocation and redistribution change
 *     node->block and node->Thread only), so clearing them in those cases is
 *     conservative; the cost is the lazy refill of the entries that are used again.
 *     An application that allocates or deallocates blocks between successive
 *     interpolations loses part of the benefit but not the correctness.
 * Code that changes the tree geometry (node coordinates, upNode/downNode links)
 * without incrementing nMeshModificationCounter would leave stale entries; such code
 * must increment the counter or run with the cache disabled.  Changes of
 * IsUsedInCalculationFlag or of the field data need no action.
 *
 * 5. Control
 * ----------
 * Environment variable AMPS_INTERPOLATION_STENCIL_CACHE=0 disables the cache (read
 * once per process); any other value or no value enables it.
 *
 * 6. Expected benefit
 * -------------------
 *   - Row-stencil users (srcEarth, srcGeospace with AMPS_EARTH_ENGINE=LEGACY): Mode3D
 *     MESH trajectory tracing measured 1.65-1.68x faster on one thread and
 *     1.53-1.64x faster with 8 threads (C1 test, 14k-36k blocks, 6x6x6 cells).
 *     Mesh generation and field initialization are not affected.
 *   - Applications that use only the pointer stencil: (a) does not apply (that path
 *     reads the block's ghost cells without tree searches), so the saving is limited
 *     to (b) and (c), i.e. to interpolation points near coarse/fine interfaces.  The
 *     gain depends on the fraction of such points and has not been measured.
 *   - Uniform meshes and points deep inside blocks take the trilinear path without
 *     (b) and (c); the cache adds one table lookup per row-stencil cell there.
 *
 * =====================================================================================
 * ROW-ONLY CONSTRUCTION (added 2026-10-23)
 * =====================================================================================
 * InitStencil(x,node,Stencil,RowStencil,BuildPointerStencil) and the trilinear and
 * multiblock routines build the pointer stencil and, if RowStencil!=NULL, the row
 * stencil in the same pass.  Callers that need only the row (the overload
 * InitStencil(x,node,RowStencil), used by the srcEarth/srcGeospace Mode3D field
 * evaluation) previously received a pointer stencil that they discarded: it was
 * built whenever the block was allocated on the rank (cell pointers, the
 * out-of-domain test of every cell, pointer merging, the coarse/fine blending of the
 * pointer stencil, the constant fallbacks with an extra tree search).  With
 * BuildPointerStencil=false all pointer-stencil steps are skipped and Stencil is
 * returned empty; the row is built by the same code and is identical.  The parameter
 * defaults to true, so existing callers are unchanged.
 * Measured (srcEarth C1 Mode3D MESH tracing, with the stencil-cell cache): 1.33-1.35x
 * on one thread, 1.18x with 8 threads.
 *
 * =====================================================================================
 * ROW-STRUCTURE REUSE (row-only construction, added 2026-10-23)
 * =====================================================================================
 * A trilinear row is determined by (node,i0,j0,k0), the lower index of the 2x2x2
 * support, and a multiblock row by (node,ijkStencilMin): within such a region of x
 * only the weights change, not the cells, the order in which they are added, or which
 * additions merge into an existing entry.  For row-only construction the trilinear
 * and multiblock routines keep, per thread, a direct-mapped table of recorded
 * structures (nRowStructureSlots=1024 slots, about 1.5 KB each):
 *   - the first construction with a key runs the original code and records the
 *     outcome of every AddRowStencilCell() call (appended as entry t, merged into
 *     entry t, or nothing added) and IsUsedInCalculationFlag of every node whose flag
 *     the construction read;
 *   - later constructions with the same key compute the weights with the same code
 *     and replay the recorded operations in the same order (append, or Weight+=w into
 *     the recorded entry); normalization, blending and merging that follow are
 *     unchanged.  The row is bitwise identical to the original construction.
 *   - A structure is rebuilt when one of the recorded flags changed and is discarded
 *     with the same triggers as the stencil-cell cache (mesh or
 *     nMeshModificationCounter changed, slot taken by another key).  A multiblock
 *     construction with a missing support point (the constant fallback depends on x)
 *     is not replayed.
 * Measured (same tests): another 1.1-1.25x.  Stencil-cell cache, row-only
 * construction and reuse together: 2.3-2.9x of the tracing time of the original code.
 * AMPS_INTERPOLATION_STENCIL_CACHE=0 disables the reuse as well.
 *
 * =====================================================================================
 * NORMALIZATION (changed 2026-10-23)
 * =====================================================================================
 * The weights of a complete trilinear or multiblock stencil sum to 1.  A stencil is
 * renormalized when, and only when, a contribution with a positive weight was not
 * added ("dropped weight"):
 *   - row stencil: AddRowStencilCell() added nothing (cell center outside the global
 *     box, owning leaf not used in the calculation, index out of range);
 *   - pointer stencil: the cell is not available in the local block or its ghost
 *     cells (GetCenterNode()==NULL), or cStencilGeneric::AddCell() rejected it
 *     (center outside the global box of a non-periodic mesh; AddCell() returns false).
 * This applies to the trilinear and to the multiblock stencils.  Before, the trilinear
 * row was normalized when Length!=8, which also normalized complete rows whose cells
 * had merged (ghost cells resolving to one coarser cell); the trilinear pointer
 * stencil tested StencilTable->Length, the stencil of OpenMP thread 0 rather than the
 * one being built; and multiblock stencils were never normalized.  A multiblock
 * pointer stencil with an unavailable cell still uses the constant fallback.  The
 * srcGeospace trajectory engine (gpu/EngineMesh.h) applies the same rule to its port
 * of the row stencil.
 * =====================================================================================
 */

#include "pic.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>

PIC::InterpolationRoutines::CellCentered::cStencil _TARGET_DEVICE_ _CUDA_MANAGED_ *PIC::InterpolationRoutines::CellCentered::StencilTable=NULL;
PIC::InterpolationRoutines::CornerBased::cStencil _TARGET_DEVICE_ _CUDA_MANAGED_ *PIC::InterpolationRoutines::CornerBased::StencilTable=NULL;

int PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent=0;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::nBlockFoundMax];
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::last=NULL;

namespace {
  typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cInterpolationTreeNode;

#ifndef __CUDA_ARCH__
  /*
   * Stencil-cell cache (host code). See "STENCIL-CELL CACHE" at the top of this file
   * for what is cached, why, the exactness argument and what triggers re-caching.
   * Users: AddRowStencilCell() (row cells), FindStencilSupportNode() (multiblock
   * support points) and GetStencilNeib() (coarser-neighbor search).
   */
  const int nStencilCacheCells=(_BLOCK_CELLS_X_+2)*(_BLOCK_CELLS_Y_+2)*(_BLOCK_CELLS_Z_+2);
  const int nStencilCacheSlots=256;   //the slot index is the upper 8 bits of a 64-bit hash
  enum {StencilCacheRowKnown=1,StencilCacheSupportKnown=2};

  static_assert((_BLOCK_CELLS_X_<256)&&(_BLOCK_CELLS_Y_<256)&&(_BLOCK_CELLS_Z_<256),"cell indices are stored as unsigned char");

  struct cStencilCacheEntry {
    const void *mesh;
    unsigned long int counter;
    const cInterpolationTreeNode *node;
    unsigned char known[nStencilCacheCells];
    cInterpolationTreeNode *rowNode[nStencilCacheCells];   //NULL: no row cell is added
    unsigned char rowIjk[nStencilCacheCells][3];
    cInterpolationTreeNode *supportNode[nStencilCacheCells];
    bool neibKnown[26];
    cInterpolationTreeNode *neib[26];                      //faces 0-5, edges 6-17, corners 18-25
  };

  bool StencilCacheEnabled() {
    static const bool enabled=[] {
      const char *e=getenv("AMPS_INTERPOLATION_STENCIL_CACHE");
      return !((e!=NULL)&&(strcmp(e,"0")==0));
    }();

    return enabled;
  }

  /* The calling thread's cache entry of node. The returned pointer must not be kept
   * across another call of this function (the slot may be reassigned). */
  cStencilCacheEntry *GetStencilCacheEntry(const cInterpolationTreeNode *node) {
    if ((node==NULL)||(StencilCacheEnabled()==false)) return NULL;

    thread_local std::unique_ptr<cStencilCacheEntry[]> table;
    if (!table) table.reset(new cStencilCacheEntry[nStencilCacheSlots]());

    const uint64_t key=static_cast<uint64_t>(reinterpret_cast<uintptr_t>(node));
    cStencilCacheEntry &e=table[(key*0x9E3779B97F4A7C15ULL)>>56];
    const void *mesh=PIC::Mesh::mesh;
    const unsigned long int counter=PIC::Mesh::mesh->nMeshModificationCounter;

    if ((e.node!=node)||(e.mesh!=mesh)||(e.counter!=counter)) {
      e.node=node,e.mesh=mesh,e.counter=counter;
      memset(e.known,0,sizeof(e.known));
      memset(e.neibKnown,0,sizeof(e.neibKnown));
    }

    return &e;
  }

  /* Index of cell (i,j,k), i,j,k in [-1,_BLOCK_CELLS_], in a cache entry; -1 otherwise. */
  int StencilCacheCellIndex(int i,int j,int k) {
    if ((i<-1)||(i>_BLOCK_CELLS_X_)||(j<-1)||(j>_BLOCK_CELLS_Y_)||(k<-1)||(k>_BLOCK_CELLS_Z_)) return -1;
    return (i+1)+(_BLOCK_CELLS_X_+2)*((j+1)+(_BLOCK_CELLS_Y_+2)*(k+1));
  }
#endif

  /*
   * Add a row-stencil entry identified by a physical cell-center coordinate.
   *
   * The AMR tree is replicated by AMPS, while node->block is distributed. This
   * helper intentionally uses only the tree geometry. It resolves a cell center
   * to the leaf tree node that owns that center and converts the coordinate to
   * an interior (i,j,k) index in that node. The resulting identifier is therefore
   * independent of which MPI rank owns the corresponding cDataBlockAMR.
   */
  _TARGET_HOST_ _TARGET_DEVICE_
  bool LocateRowStencilCell(const double *xIn,cInterpolationTreeNode *StartNode,
      cInterpolationTreeNode *&CellNode,int *ijk);

  _TARGET_HOST_ _TARGET_DEVICE_
  bool AddRowStencilCellByCoordinate(PIC::InterpolationRoutines::cRowStencil *RowStencil,double Weight,
      const double *xIn,cInterpolationTreeNode *StartNode) {
    if (RowStencil==NULL) return false;

    cInterpolationTreeNode *CellNode;
    int ijk[3];

    if (LocateRowStencilCell(xIn,StartNode,CellNode,ijk)==false) return false;
    if (CellNode->IsUsedInCalculationFlag==false) return false;

    RowStencil->AddCell(Weight,CellNode,ijk[0],ijk[1],ijk[2]);
    return true;
  }

  /*
   * The leaf and the interior (i,j,k) index that AddRowStencilCellByCoordinate()
   * adds for xIn, without its IsUsedInCalculationFlag test (which depends on the
   * state of the leaf, not on the tree geometry). Returns false when no cell is
   * added for any value of that flag.
   */
  _TARGET_HOST_ _TARGET_DEVICE_
  bool LocateRowStencilCell(const double *xIn,cInterpolationTreeNode *StartNode,
      cInterpolationTreeNode *&CellNode,int *ijk) {
    double x[3]={xIn[0],xIn[1],xIn[2]};

    /*
     * Non-periodic ghost centers outside the global domain do not correspond to
     * physical cells and are omitted, matching cStencilGeneric::AddCell(). For a
     * periodic mesh, map the center back into the primary global box so the row
     * references the physical periodic image rather than a ghost index.
     */
    for (int d=0;d<3;d++) {
      const double xmin=PIC::Mesh::mesh->xGlobalMin[d];
      const double xmax=PIC::Mesh::mesh->xGlobalMax[d];

      if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
        const double length=xmax-xmin;

        if (length>0.0) {
          while (x[d]<xmin) x[d]+=length;
          while (x[d]>xmax) x[d]-=length;
        }
      }
      else if ((x[d]<xmin)||(x[d]>xmax)) return false;
    }

    CellNode=PIC::Mesh::mesh->findTreeNode(x,StartNode);

    if (CellNode==NULL) return false;

    const int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

    for (int d=0;d<3;d++) {
      const double blockLength=CellNode->xmax[d]-CellNode->xmin[d];
      const double xLoc=(x[d]-CellNode->xmin[d])/blockLength*nCells[d];

      ijk[d]=(int)floor(xLoc);

      /* Protect against roundoff when a coordinate is infinitesimally outside a
       * block face. A genuine out-of-range index is rejected. */
      if ((ijk[d]<0)&&(xLoc>-1.0e-10)) ijk[d]=0;
      if ((ijk[d]>=nCells[d])&&(xLoc<nCells[d]+1.0e-10)) ijk[d]=nCells[d]-1;
      if ((ijk[d]<0)||(ijk[d]>=nCells[d])) return false;
    }

    return true;
  }

  /*
   * Add a row-stencil entry from indices measured in SourceNode. The indices may
   * be ghost indices such as -1 or _BLOCK_CELLS_X_. They are first converted to
   * a physical center coordinate and then canonicalized to the owning leaf node
   * by AddRowStencilCellByCoordinate().
   */
  _TARGET_HOST_ _TARGET_DEVICE_
  bool AddRowStencilCell(PIC::InterpolationRoutines::cRowStencil *RowStencil,double Weight,
      cInterpolationTreeNode *SourceNode,int i,int j,int k) {
    if ((RowStencil==NULL)||(SourceNode==NULL)) return false;

    const int ijk[3]={i,j,k};
    const int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    double x[3];

    #ifndef __CUDA_ARCH__
    cStencilCacheEntry *e=GetStencilCacheEntry(SourceNode);
    const int c=(e==NULL) ? -1 : StencilCacheCellIndex(i,j,k);

    if (c>=0) {
      if ((e->known[c]&StencilCacheRowKnown)==0) {
        cInterpolationTreeNode *CellNode;
        int ijkCell[3];

        for (int d=0;d<3;d++) {
          const double dx=(SourceNode->xmax[d]-SourceNode->xmin[d])/nCells[d];
          x[d]=SourceNode->xmin[d]+(ijk[d]+0.5)*dx;
        }

        if (LocateRowStencilCell(x,SourceNode,CellNode,ijkCell)==false) CellNode=NULL;

        e->rowNode[c]=CellNode;
        if (CellNode!=NULL) for (int d=0;d<3;d++) e->rowIjk[c][d]=static_cast<unsigned char>(ijkCell[d]);
        e->known[c]|=StencilCacheRowKnown;
      }

      cInterpolationTreeNode *CellNode=e->rowNode[c];

      if ((CellNode==NULL)||(CellNode->IsUsedInCalculationFlag==false)) return false;

      RowStencil->AddCell(Weight,CellNode,e->rowIjk[c][0],e->rowIjk[c][1],e->rowIjk[c][2]);
      return true;
    }
    #endif

    for (int d=0;d<3;d++) {
      const double dx=(SourceNode->xmax[d]-SourceNode->xmin[d])/nCells[d];
      x[d]=SourceNode->xmin[d]+(ijk[d]+0.5)*dx;
    }

    return AddRowStencilCellByCoordinate(RowStencil,Weight,x,SourceNode);
  }

  /*
   * findTreeNode(xLogical,node) for the support point (i,j,k) of node in
   * GetTriliniarInterpolationMutiBlockStencil(); the caller computes xLogical from
   * (node,i,j,k), so the cached result is that of the same search.
   */
  _TARGET_HOST_ _TARGET_DEVICE_
  cInterpolationTreeNode *FindStencilSupportNode(double *xLogical,cInterpolationTreeNode *node,const int *ijk) {
    #ifndef __CUDA_ARCH__
    cStencilCacheEntry *e=GetStencilCacheEntry(node);
    const int c=(e==NULL) ? -1 : StencilCacheCellIndex(ijk[0],ijk[1],ijk[2]);

    if (c>=0) {
      if ((e->known[c]&StencilCacheSupportKnown)==0) {
        e->supportNode[c]=PIC::Mesh::mesh->findTreeNode(xLogical,node);
        e->known[c]|=StencilCacheSupportKnown;
      }

      return e->supportNode[c];
    }
    #endif

    return PIC::Mesh::mesh->findTreeNode(xLogical,node);
  }

  /* node->GetNeibFace(n,0,0) (n=0-5), GetNeibEdge(n-6,0) (n=6-17) and
   * GetNeibCorner(n-18) (n=18-25) as used by InitStencil(). */
  _TARGET_HOST_ _TARGET_DEVICE_
  cInterpolationTreeNode *GetStencilNeib(cInterpolationTreeNode *node,int n) {
    #ifndef __CUDA_ARCH__
    cStencilCacheEntry *e=GetStencilCacheEntry(node);

    if ((e!=NULL)&&(e->neibKnown[n]==true)) return e->neib[n];
    #endif

    cInterpolationTreeNode *res;

    if (n<6) res=node->GetNeibFace(n,0,0,PIC::Mesh::mesh);
    else if (n<18) res=node->GetNeibEdge(n-6,0,PIC::Mesh::mesh);
    else res=node->GetNeibCorner(n-18,PIC::Mesh::mesh);

    #ifndef __CUDA_ARCH__
    if (e!=NULL) e->neib[n]=res,e->neibKnown[n]=true;
    #endif

    return res;
  }

#ifndef __CUDA_ARCH__
  /*
   * Row-structure reuse (row-only construction, host code). See "ROW-STRUCTURE
   * REUSE" at the top of this file.
   *
   * A trilinear row is determined by (node,i0,j0,k0) and a multiblock row by
   * (node,ijkStencilMin): the cells, the order in which they are added and which
   * additions merge into an existing entry do not depend on x within that region.
   * The first construction records, for every AddRowStencilCell() call, whether it
   * appended an entry, merged into entry t, or added nothing; later constructions
   * with the same key compute the weights with the same code and replay these
   * operations. A structure also records the IsUsedInCalculationFlag of every node
   * whose flag the construction read, and is rebuilt when one of them changed.
   */
  const int nRowStructureSlots=1024;   //the slot index is the upper 10 bits of a 64-bit hash
  const int nMaxRowContributions=64;   //multiblock: 8 support points x 8 fine cells
  const int nMaxRowChecks=nMaxRowContributions+8;
  enum {RowStructureTrilinear=1,RowStructureMultiBlock=2};
  enum {RowStructureEmpty=0,RowStructureReusable=1,RowStructureOriginal=2};
  enum {RowModeOriginal=0,RowModeReplay=1,RowModeRecord=2};

  struct cRowStructure {
    const void *mesh;
    unsigned long int counter;
    const cInterpolationTreeNode *node;
    int kind,i0,j0,k0;
    int state;                       //RowStructureEmpty/Reusable/Original (run the original code)
    bool recordOk;
    int nContrib;
    signed char target[nMaxRowContributions];   //entry index, -1: nothing added
    int length;
    cInterpolationTreeNode *entryNode[nMaxRowContributions];
    int entryIjk[nMaxRowContributions][3];
    unsigned char fineSupport[8];    //multiblock: the support point is represented by 8 fine cells
    int nCheck;
    const cInterpolationTreeNode *checkNode[nMaxRowChecks];
    bool checkUsed[nMaxRowChecks];
  };

  /* The calling thread's structure slot of (kind,node,i0,j0,k0); NULL when the cache
   * is disabled. The slot is reset when it holds another key or mesh state. */
  cRowStructure *GetRowStructure(int kind,const cInterpolationTreeNode *node,int i0,int j0,int k0) {
    if ((node==NULL)||(StencilCacheEnabled()==false)) return NULL;

    thread_local std::unique_ptr<cRowStructure[]> table;
    if (!table) table.reset(new cRowStructure[nRowStructureSlots]());

    const uint64_t key=static_cast<uint64_t>(reinterpret_cast<uintptr_t>(node))*0x9E3779B97F4A7C15ULL^
        (static_cast<uint64_t>(((kind*1024+(i0+1))*1024+(j0+1))*1024+(k0+1))*0xC2B2AE3D27D4EB4FULL);
    cRowStructure &e=table[(key*0x9E3779B97F4A7C15ULL)>>54];
    const void *mesh=PIC::Mesh::mesh;
    const unsigned long int counter=PIC::Mesh::mesh->nMeshModificationCounter;

    if ((e.node!=node)||(e.kind!=kind)||(e.i0!=i0)||(e.j0!=j0)||(e.k0!=k0)||(e.mesh!=mesh)||(e.counter!=counter)) {
      e.node=node,e.kind=kind,e.i0=i0,e.j0=j0,e.k0=k0,e.mesh=mesh,e.counter=counter;
      e.state=RowStructureEmpty;
    }

    return &e;
  }

  /* RowModeReplay/RowModeOriginal for a recorded structure whose node flags are
   * unchanged; otherwise start recording (RowModeRecord). */
  int StartRowStructure(cRowStructure *rs) {
    if (rs==NULL) return RowModeOriginal;

    if (rs->state!=RowStructureEmpty) {
      bool valid=true;

      for (int n=0;n<rs->nCheck;n++) if (rs->checkNode[n]->IsUsedInCalculationFlag!=rs->checkUsed[n]) {
        valid=false;
        break;
      }

      if (valid==true) return (rs->state==RowStructureReusable) ? RowModeReplay : RowModeOriginal;
    }

    rs->state=RowStructureEmpty,rs->recordOk=true,rs->nContrib=0,rs->length=0,rs->nCheck=0;
    return RowModeRecord;
  }

  /* Record that the construction read IsUsedInCalculationFlag of n. */
  void AddRowStructureCheck(cRowStructure *rs,const cInterpolationTreeNode *n) {
    if (n==NULL) return;
    for (int i=0;i<rs->nCheck;i++) if (rs->checkNode[i]==n) return;

    if (rs->nCheck==nMaxRowChecks) {
      rs->recordOk=false;
      return;
    }

    rs->checkNode[rs->nCheck]=n,rs->checkUsed[rs->nCheck]=n->IsUsedInCalculationFlag;
    rs->nCheck++;
  }

  /* Record the outcome of AddRowStencilCell(row,w,SourceNode,i,j,k), which returned
   * `added` and found the row with LengthBefore entries. */
  void RecordRowContribution(cRowStructure *rs,const PIC::InterpolationRoutines::cRowStencil *row,
      int LengthBefore,bool added,const cInterpolationTreeNode *SourceNode,int i,int j,int k) {
    if (rs->nContrib==nMaxRowContributions) {
      rs->recordOk=false;
      return;
    }

    //the resolution of the cell is in the stencil-cell cache (filled by AddRowStencilCell)
    cStencilCacheEntry *e=GetStencilCacheEntry(SourceNode);
    const int c=(e==NULL) ? -1 : StencilCacheCellIndex(i,j,k);

    if ((c<0)||((e->known[c]&StencilCacheRowKnown)==0)) {
      rs->recordOk=false;
      return;
    }

    cInterpolationTreeNode *owner=e->rowNode[c];
    int t=-1;

    AddRowStructureCheck(rs,owner);

    if (added==true) {
      const int ijk[3]={e->rowIjk[c][0],e->rowIjk[c][1],e->rowIjk[c][2]};

      if (row->Length>LengthBefore) t=LengthBefore;
      else for (int n=0;n<row->Length;n++) {
        const PIC::InterpolationRoutines::cRowStencil::cElement &el=row->Element[n];

        if ((el.node==owner)&&(el.i==ijk[0])&&(el.j==ijk[1])&&(el.k==ijk[2])) {
          t=n;
          break;
        }
      }

      if ((t<0)||(t>=nMaxRowContributions)||(row->Element[t].node!=owner)||(row->Element[t].i!=ijk[0])||
          (row->Element[t].j!=ijk[1])||(row->Element[t].k!=ijk[2])) {
        rs->recordOk=false;
        return;
      }

      if (t==LengthBefore) {
        rs->entryNode[t]=owner;
        for (int d=0;d<3;d++) rs->entryIjk[t][d]=ijk[d];
        rs->length=row->Length;
      }
    }

    rs->target[rs->nContrib++]=static_cast<signed char>(t);
  }

  void EndRowStructure(cRowStructure *rs,bool reusable) {
    rs->state=((rs->recordOk==true)&&(reusable==true)) ? RowStructureReusable : RowStructureOriginal;
  }

  /* Replay contribution c with weight w: the operation AddRowStencilCell() performed
   * when the structure was recorded (append, merge into an entry, or nothing).
   * Returns true when nothing is added and w>0 (a dropped weight). */
  inline bool ApplyRowContribution(const cRowStructure *rs,int c,double w,
      PIC::InterpolationRoutines::cRowStencil *row) {
    const int t=rs->target[c];

    if (t<0) return (w>0.0);   //true: a positive weight is dropped
    if (w<0.0) exit(__LINE__,__FILE__,"Error: found a negative interpolation weight");

    if (t==row->Length) {
      PIC::InterpolationRoutines::cRowStencil::cElement &el=row->Element[t];

      el.node=rs->entryNode[t];
      el.i=rs->entryIjk[t][0],el.j=rs->entryIjk[t][1],el.k=rs->entryIjk[t][2];
      el.Weight=w;
      row->Length++;
    }
    else if (t<row->Length) row->Element[t].Weight+=w;
    else exit(__LINE__,__FILE__,"Error: inconsistent row structure");

    return false;
  }
#endif

  /* Build the one-cell row used by the legacy constant-interpolation fallback. */
  _TARGET_HOST_ _TARGET_DEVICE_
  void BuildConstantRowStencil(double *x,cInterpolationTreeNode *node,
      PIC::InterpolationRoutines::cRowStencil *RowStencil) {
    if (RowStencil==NULL) return;

    RowStencil->flush();
    AddRowStencilCellByCoordinate(RowStencil,1.0,x,node);
  }

  /*
   * Add a center-node pointer while merging duplicate pointers. The legacy
   * multiblock routine obtains this behavior through cStencilGeneric::Add(); the
   * helper preserves it when the stencil is assembled directly element by
   * element alongside the row stencil.
   */
  _TARGET_HOST_ _TARGET_DEVICE_
  bool AddPhysicalStencilCell(PIC::InterpolationRoutines::CellCentered::cStencil &Stencil,
      double Weight,PIC::Mesh::cDataCenterNode *cell,int LocalCellID) {
    for (int iElement=0;iElement<Stencil.Length;iElement++) {
      if (Stencil.cell[iElement]==cell) {
        Stencil.Weight[iElement]+=Weight;
        return true;
      }
    }

    //false: the cell is outside the global domain and its weight is dropped
    return Stencil.AddCell(Weight,cell,LocalCellID);
  }
}

//the locally ordered interpolation coeffcients for the corner based interpolation procedure
//thread_local double PIC::InterpolationRoutines::CornerBased::InterpolationCoefficientTable_LocalNodeOrder[8]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

// macro switch is needed in the case some other interpolation is used
// and interface function is not compiled
#if _PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_
//definition of the interface functions that are used to access the FORTRAN written linear cell cenetered interpolation library
extern "C"{
  void interface__cell_centered_linear_interpolation__init_stencil_(int* nDim, double* XyzIn_D, int* nIndexes, int* nCell_D,  int* nGridOut, double* Weight_I, int* iIndexes_II, int* IsSecondOrder, int* UseGhostCell);
}
#endif//_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_


//initialize the interpolation module
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::Init() {

  //init the stencil table
//  CellCentered::StencilTable=new CellCentered::cStencil[PIC::nTotalThreadsOpenMP];
//  CornerBased::StencilTable=new CornerBased::cStencil[PIC::nTotalThreadsOpenMP];

  #ifndef __CUDA_ARCH__
  if (CellCentered::StencilTable==NULL) {
    amps_new_managed(CellCentered::StencilTable,PIC::nTotalThreadsOpenMP); 
    amps_new_managed(CornerBased::StencilTable,PIC::nTotalThreadsOpenMP);
  }
  #endif
}

//determine stencil for the cell centered piecewise constant interpolation
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil) {
  int i,j,k;
  long int nd;
  PIC::Mesh::cDataCenterNode *cell;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  //flush the stencil
  Stencil.flush();

  //find the block
  if (node==NULL) {
    node=PIC::Mesh::mesh->findTreeNode(x);

    if (node==NULL) exit(__LINE__,__FILE__,"Error: the location is outside of the computational domain");
    if (node->block==NULL) {
      char msg[200];

      #ifndef __CUDA_ARCH__
      sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
      exit(__LINE__,__FILE__,msg);
      #else
      printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
      exit(__LINE__,__FILE__);
      #endif
    }
  }
  else if (node->block==NULL) {
    char msg[200];

    #ifndef __CUDA_ARCH__
    sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
    exit(__LINE__,__FILE__,msg);
    #else 
    printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
    exit(__LINE__,__FILE__,msg); 
    #endif
  }

  //find cell
  nd = PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false);
  cell=node->block->GetCenterNode(nd);//getCenterNodeLocalNumber(i,j,k));

  //add the cell to the stencil
  if (cell!=NULL) {
    Stencil.AddCell(1.0,cell,nd);
  }
  else exit(__LINE__,__FILE__,"Error: cell is not initialized");
}

//determine the stencil for the cell centered linear interpolation using interpolation library from ../share/Library/src/
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(double *XyzIn_D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil,PIC::InterpolationRoutines::cRowStencil *RowStencil,bool BuildPointerStencil) {

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();

  if  (_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__SWMF_) {
    exit(__LINE__,__FILE__,"Error: the procedure is not adapted fro using with OpenMP: INTERFACE::BlockFound... need to be converted into an array as PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP] ");
  }
  #else
  int ThreadOpenMP=0;
  #endif

  Stencil.flush();
  if (RowStencil!=NULL) RowStencil->flush();

  // macro switch is needed in the case some other interpolation is used
  // and interface function is not compiled
  if (_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
    //find the block if needed
    if (node==NULL) {
      node=PIC::Mesh::mesh->findTreeNode(XyzIn_D);

      if (node==NULL) exit(__LINE__,__FILE__,"Error: the location is outside of the computational domain");
      if ((node->block==NULL)&&(RowStencil==NULL)) {
        char msg[200];

        #ifndef __CUDA_ARCH__
        sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
        exit(__LINE__,__FILE__,msg);
        #else 
        printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
        exit(__LINE__,__FILE__);
        #endif
      }
    }
    else if ((node->block==NULL)&&(RowStencil==NULL)) {
      char msg[200];

      #ifndef __CUDA_ARCH__
      sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
      exit(__LINE__,__FILE__,msg);
      #else 
      printf("Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::GPU::ThisThread,node->Thread);
      exit(__LINE__,__FILE__,msg);
      #endif
    }

    //check whether the point is located deep in the block -> use three linear interpolation
    double iLoc,jLoc,kLoc;
    double *xmin,*xmax;

    xmin=node->xmin;
    xmax=node->xmax;

    iLoc=(XyzIn_D[0]-xmin[0])/(xmax[0]-xmin[0])*_BLOCK_CELLS_X_;
    jLoc=(XyzIn_D[1]-xmin[1])/(xmax[1]-xmin[1])*_BLOCK_CELLS_Y_;
    kLoc=(XyzIn_D[2]-xmin[2])/(xmax[2]-xmin[2])*_BLOCK_CELLS_Z_;

    #if _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
    //in case the mesh is uniform (_AMR_MESH_TYPE_ == _AMR_MESH_TYPE__UNIFORM_) -> use simple trilinear interpolation
    if (_AMR_MESH_TYPE_ == _AMR_MESH_TYPE__UNIFORM_) {
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil,RowStencil,BuildPointerStencil);
      return;
    }

    //if the point of interest is deep inside the block or all neighbors of the block has the same resolution level -> use simple trilinear interpolation
    if ((node->RefinmentLevel==node->minNeibRefinmentLevel)&&(node->RefinmentLevel==node->maxNeibRefinmentLevel)) {
      //all blocks around has the same level of the resolution
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil,RowStencil,BuildPointerStencil);
      return;
    }
    else if ((1.0<iLoc)&&(iLoc<_BLOCK_CELLS_X_-1) && (1.0<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-1) && (1.0<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-1)) {
      //the point is deep inside the block
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil,RowStencil,BuildPointerStencil);
      return;
    }
    else if (node->RefinmentLevel==node->minNeibRefinmentLevel) {
      if  ((0.5<iLoc)&&(iLoc<_BLOCK_CELLS_X_-0.5) && (0.5<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-0.5) && (0.5<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-0.5)) {
        //1. the point is deeper than a half cell size insode the block
        //2. the block is geometrically largest among the neibours
        //3. -> it is only the intermal points of the block that will be used in the stencil
        GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil,RowStencil,BuildPointerStencil);
        return;
      }
      else  {
        //the block is largest between the neibours
        //getermine the size limit for the interpolation stencil
        double xStencilMin[3],xStencilMax[3];
        double dxCell[3];

        dxCell[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
        dxCell[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
        dxCell[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

        #pragma ivdep
        for (int idim=0;idim<3;idim++) {
          int iInterval;

          iInterval=(int)((XyzIn_D[idim]-(xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
          xStencilMin[idim]=(xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
          xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
        }

        GetTriliniarInterpolationMutiBlockStencil(XyzIn_D,xStencilMin,xStencilMax,node,Stencil,RowStencil,BuildPointerStencil);
        return;
      }
    }
    else {
      //1. find a coarser block that is close to the point, and can be used for interpolation
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *CoarserBlock=NULL;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *NeibNode;
      int idim,iFace;
      double dxCell[3];

      dxCell[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
      dxCell[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
      dxCell[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

      int nBlockCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
      double dmin=10.0*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;

      bool CornerTestFlagTable[8]={false,false,false,false,  false,false,false,false};
      bool EdgeTestFlagTable[12]={false,false,false,false, false,false,false,false, false,false,false,false};

      for (idim=0;idim<3;idim++) {
        switch (idim) {
        case 0:
          if (iLoc<=1.0) iFace=0;
          else if (iLoc>=_BLOCK_CELLS_X_-1.0) iFace=1;
          else continue;

          break;
        case 1:
          if (jLoc<=1.0) iFace=2;
          else if (jLoc>=_BLOCK_CELLS_Y_-1.0) iFace=3;
          else continue;

          break;
        case 2:
          if (kLoc<=1.0) iFace=4;
          else if (kLoc>=_BLOCK_CELLS_Z_-1.0) iFace=5;
          else continue;
        }

        //check blocks connected through the face
        NeibNode=GetStencilNeib(node,iFace);

        double xLoc[3]={iLoc,jLoc,kLoc};

        if (NeibNode!=NULL) if ((NeibNode->RefinmentLevel<node->RefinmentLevel)&&(NeibNode->IsUsedInCalculationFlag==true)) {
          //found a coarser block
          int cnt=0;

          for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

          if (cnt==3) {
            //the block can be used for interpolation
            switch(iFace) {
            case 0:
              if ((xLoc[0]<1.0)&&(xLoc[0]<dmin)) dmin=xLoc[0],CoarserBlock=NeibNode;
              break;

            case 1:
              if ((xLoc[0]>nBlockCells[0]-1)&&(nBlockCells[0]-xLoc[0]<dmin)) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
              break;

            case 2:
              if ((xLoc[1]<1.0)&&(xLoc[1]<dmin)) dmin=xLoc[1],CoarserBlock=NeibNode;
              break;

            case 3:
              if ((xLoc[1]>nBlockCells[1]-1)&&(nBlockCells[1]-xLoc[1]<dmin)) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
              break;

            case 4:
              if ((xLoc[2]<1.0)&&(xLoc[2]<dmin)) dmin=xLoc[2],CoarserBlock=NeibNode;
              break;

            case 5:
              if ((xLoc[2]>nBlockCells[2]-1)&&(nBlockCells[2]-xLoc[2]<dmin)) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
              break;

            default:
              exit(__LINE__,__FILE__,"error: something is wrong");
            }
          }
        }

        //check blocks connected though the edges
        static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};


        for (int iEdge=0;iEdge<4;iEdge++) if (EdgeTestFlagTable[faceEdges[iFace][iEdge]]==false) {
          EdgeTestFlagTable[faceEdges[iFace][iEdge]]=true;
          NeibNode=GetStencilNeib(node,6+faceEdges[iFace][iEdge]);

          if (NeibNode!=NULL) if ((NeibNode->RefinmentLevel<node->RefinmentLevel)&&(NeibNode->IsUsedInCalculationFlag==true)) {
            //found a coarser block
            int cnt=0;

            for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

            if (cnt==3) {
              //the block can be used for interpolation
              switch (faceEdges[iFace][iEdge]) {
              case 0:
                if ((xLoc[1]<1.0)&&(xLoc[2]<1.0)) {
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 1:
                if ((xLoc[1]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 2:
                if ((xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 3:
                if ((xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;


              case 4:
                if ((xLoc[0]<1.0)&&(xLoc[2]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 5:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 6:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 7:
                if ((xLoc[0]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 8:
                if ((xLoc[0]<1.0)&&(xLoc[1]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                }

                break;

              case 9:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                }

                break;

              case 10:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]>nBlockCells[1]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                }

                break;

              case 11:
                if ((xLoc[0]<1.0)&&(xLoc[1]>nBlockCells[1]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                }


                break;

              default:
                exit(__LINE__,__FILE__,"error: something is wrong");
              }
            }
          }
        }

        //check blocks connected through the corners
        static const int FaceNodeMap[6][4]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};

        for (int iCorner=0;iCorner<4;iCorner++) if (CornerTestFlagTable[FaceNodeMap[iFace][iCorner]]==false) {
          CornerTestFlagTable[FaceNodeMap[iFace][iCorner]]=true;
          NeibNode=GetStencilNeib(node,18+FaceNodeMap[iFace][iCorner]);

          if (NeibNode!=NULL) if ((NeibNode->RefinmentLevel<node->RefinmentLevel)&&(NeibNode->IsUsedInCalculationFlag==true)) {
            //found a coarser block
            int cnt=0;

            for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

            if (cnt==3) {
              //the block can be used for interpolation
              switch (FaceNodeMap[iFace][iCorner]) {
              case 0:
                if ((xLoc[0]<1.0)&&(xLoc[1]<1.0)&&(xLoc[2]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 1:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]<1.0)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 2:
                if ((xLoc[0]<1.0)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]<1.0)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 3:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]<1.0)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (xLoc[2]<dmin) dmin=xLoc[2],CoarserBlock=NeibNode;
                }

                break;


              case 4:
                if ((xLoc[0]<1.0)&&(xLoc[1]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 5:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]<1.0)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (xLoc[1]<dmin) dmin=xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 6:
                if ((xLoc[0]<1.0)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (xLoc[0]<dmin) dmin=xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              case 7:
                if ((xLoc[0]>nBlockCells[0]-1)&&(xLoc[1]>nBlockCells[1]-1)&&(xLoc[2]>nBlockCells[2]-1)) {
                  if (nBlockCells[0]-xLoc[0]<dmin) dmin=nBlockCells[0]-xLoc[0],CoarserBlock=NeibNode;
                  if (nBlockCells[1]-xLoc[1]<dmin) dmin=nBlockCells[1]-xLoc[1],CoarserBlock=NeibNode;
                  if (nBlockCells[2]-xLoc[2]<dmin) dmin=nBlockCells[2]-xLoc[2],CoarserBlock=NeibNode;
                }

                break;

              default:
                exit(__LINE__,__FILE__,"error: something is wrong");
              }
            }
          }
        }
      }

      if (CoarserBlock!=NULL) {
        //getermine the size limit for the interpolation stencil
        double xStencilMin[3],xStencilMax[3];
        double dxCell[3];

        dxCell[0]=(CoarserBlock->xmax[0]-CoarserBlock->xmin[0])/_BLOCK_CELLS_X_;
        dxCell[1]=(CoarserBlock->xmax[1]-CoarserBlock->xmin[1])/_BLOCK_CELLS_Y_;
        dxCell[2]=(CoarserBlock->xmax[2]-CoarserBlock->xmin[2])/_BLOCK_CELLS_Z_;

        #pragma ivdep
        for (idim=0;idim<3;idim++) {
          int iInterval;

          iInterval=(int)((XyzIn_D[idim]-(CoarserBlock->xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
          xStencilMin[idim]=(CoarserBlock->xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
          xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
        }

        GetTriliniarInterpolationMutiBlockStencil(XyzIn_D,xStencilMin,xStencilMax,CoarserBlock,Stencil,RowStencil,BuildPointerStencil);


        /* the following is used to get continuous interpolsation
        1. if the point of interpolation is in the "coarse" block -> use GetTriliniarInterpolationMutiBlockStencil()
        2. if the point of interpolation is in the "fine" block ->
           a) def "dmin" is the min length to the boundary of the block
           b) if dmin<0.5 -> "coarse" stencil interpolation
           c) if 0.5<dmin<1-> combination of "fine" and "coarse" stencils:
             interpolation=(dmin-0.5)/0.5*fine stencil+[1-(dmin-0.5)/0.5)]*corse stencil
           d) is 1<dmin -> "fine" stencil interpolation
        */

        if ((0.5<dmin)&&(dmin<=1.0)) {
          PIC::InterpolationRoutines::CellCentered::cStencil FineStencil;
          PIC::InterpolationRoutines::cRowStencil FineRowStencil;
          PIC::InterpolationRoutines::cRowStencil *FineRowStencilPtr=(RowStencil==NULL) ? NULL : &FineRowStencil;

          GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,FineStencil,FineRowStencilPtr,BuildPointerStencil);

          if (BuildPointerStencil==true) {
            Stencil.MultiplyScalar(1.0-(dmin-0.5)/0.5);
            FineStencil.MultiplyScalar((dmin-0.5)/0.5);

            Stencil.Add(&FineStencil);
          }

          /* Apply exactly the same continuous coarse/fine blending to the
           * decomposition-independent row. */
          if (RowStencil!=NULL) {
            RowStencil->MultiplyScalar(1.0-(dmin-0.5)/0.5);
            FineRowStencil.MultiplyScalar((dmin-0.5)/0.5);
            RowStencil->Add(&FineRowStencil);
          }
        }

        return;
      }

      //2.threre is no a coarser block that can be used for interpolation
      //The current block is the coarstest that can be used to construct the interpolation stencil
      double xStencilMin[3],xStencilMax[3];

      dxCell[0]=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
      dxCell[1]=(node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_;
      dxCell[2]=(node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_;

      #pragma ivdep
      for (idim=0;idim<3;idim++) {
        int iInterval;

        iInterval=(int)((XyzIn_D[idim]-(node->xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
        xStencilMin[idim]=(node->xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
        xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
      }

      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil,RowStencil,BuildPointerStencil);
      return;
    }

    #elif _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__SWMF_
    if ( ((node->RefinmentLevel==node->minNeibRefinmentLevel)&&(node->RefinmentLevel==node->maxNeibRefinmentLevel)) ||
        (0.5<iLoc)&&(iLoc<_BLOCK_CELLS_X_-0.5) && (0.5<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-0.5) && (0.5<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-0.5) ) {
      GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node,Stencil,RowStencil,BuildPointerStencil);
      return;
    }

    // if the point of interest is very close to the cell center
    // then truncated stencil to the single point
    if( fabs(iLoc-0.5-(long)iLoc) < PrecisionCellCenter &&
        fabs(jLoc-0.5-(long)jLoc) < PrecisionCellCenter &&
        fabs(kLoc-0.5-(long)kLoc) < PrecisionCellCenter   ){
          BuildConstantRowStencil(XyzIn_D,node,RowStencil);
          if ((BuildPointerStencil==true)&&(node->block!=NULL)) PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node,Stencil);
          else Stencil.flush();
          return;
    }

    //re-init variables in the INTERFACE and flush the Stencil
    INTERFACE::iBlockFoundCurrent=0;
    for(int iBlock = 0; iBlock < INTERFACE::nBlockFoundMax; iBlock++ ) INTERFACE::BlockFound[iBlock] = NULL;

    Stencil.flush();
    INTERFACE::last=node;

    //prepare the call of FORTRAN interpolate_amr subroutine
    int nDim       = DIM;
    // number of indices to identify cell ON A GIVEN PROCESSOR:
    // block id + DIM indices of cell in a block
    int nIndexes   = DIM + 1;
    int nCell_D[3] = {_BLOCK_CELLS_X_, _BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

    // IMPORTANT NOTE: iIndexes_II has an extra index per cell which is
    // a processor ID => need one more index in addition to nIndexes
    int iIndexes_II[(DIM+1+1)*PIC::InterpolationRoutines::nMaxStencilLength];
    double WeightStencil[PIC::InterpolationRoutines::nMaxStencilLength];
    int IsSecondOrder;
    int UseGhostCell=1;
    int nCellStencil;

    // call the interpolation subroutine
    interface__cell_centered_linear_interpolation__init_stencil_(&nDim,XyzIn_D,&nIndexes, nCell_D, &nCellStencil, WeightStencil, iIndexes_II, &IsSecondOrder, &UseGhostCell);

    // size of the stencil and weights are known
    // need to identify cells in the stencil
    // NOTE: FORTRAN is column-major for 2D arrays like iIndexes_II
    bool PhysicalStencilAvailable=true;
    bool TopologicalStencilAvailable=true;

    for (int iCellStencil = 0; iCellStencil < nCellStencil; iCellStencil++) {
      int ind[3]={0,0,0};
      PIC::Mesh::cDataBlockAMR  *block;
      PIC::Mesh::cDataCenterNode *cell;

  //    int iThread = iIndexes_II[0    +iCellStencil*(nIndexes+1)];
      for(int i = 0; i < DIM; i++)
        //cell indices are 1-based in FORTRAN
        ind[i]    = iIndexes_II[1+i  +iCellStencil*(nIndexes+1)]-1;
      int iBlock  = iIndexes_II[1+DIM+iCellStencil*(nIndexes+1)];

      /* The FORTRAN library has already returned a topological cell
       * identifier. Store it before checking local allocation so the optional
       * row remains complete on ranks that do not own the block. */
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *StencilNode=NULL;

      if ((iBlock>=0)&&(iBlock<INTERFACE::nBlockFoundMax)) StencilNode=INTERFACE::BlockFound[iBlock];

      if (StencilNode!=NULL) {
        if (RowStencil!=NULL) RowStencil->AddCell(WeightStencil[iCellStencil],StencilNode,ind[0],ind[1],ind[2]);
      }
      else {
        TopologicalStencilAvailable=false;
      }

      //retrieve the pointer to the current cell
      if ((BuildPointerStencil==false)||(StencilNode==NULL)||((block=StencilNode->block)==NULL)) {
        PhysicalStencilAvailable=false;
        continue;
      }

      cell=block->GetCenterNode(_getCenterNodeLocalNumber(ind[0],ind[1],ind[2]));

      if (cell==NULL) {
        //there is not enough local information to build the pointer stencil
        PhysicalStencilAvailable=false;
        continue;
      }

      int nd=_getCenterNodeLocalNumber(ind[0],ind[1],ind[2]);
      Stencil.AddCell(WeightStencil[iCellStencil],StencilNode->block->GetCenterNode(nd),nd);
    }

    if ((RowStencil!=NULL)&&(TopologicalStencilAvailable==false)) BuildConstantRowStencil(XyzIn_D,node,RowStencil);

    if (PhysicalStencilAvailable==false) {
      /* Preserve the legacy fallback whenever the containing block is local. If
       * the caller explicitly requested a row for a remote containing block,
       * there is no valid local pointer stencil; return it empty while keeping
       * the complete row. */
      if ((BuildPointerStencil==true)&&(node->block!=NULL)) PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node,Stencil);
      else Stencil.flush();

      return;
    }
    #else  //_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
    exit(__LINE__,__FILE__,"Error: the option is unknown");
    #endif //_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
  }
  else exit(__LINE__,__FILE__,"ERROR: cell centered linear interpolation is currently available only through interface, add corresponding block to the input file!");
}

//triliniar interpolation used inside blocks
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Linear::GetTriliniarInterpolationStencil(double iLoc,double jLoc,double kLoc,double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil &Stencil,PIC::InterpolationRoutines::cRowStencil *RowStencil,bool BuildPointerStencil) {
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;

  //row-only construction (BuildPointerStencil==false): treat the block as not allocated
  block=(BuildPointerStencil==true) ? node->block : NULL;
  if ((block==NULL)&&(RowStencil==NULL)) exit(__LINE__,__FILE__,"Error: the block is not allocated");

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  //flush the stencil
  Stencil.flush();
  if (RowStencil!=NULL) RowStencil->flush();

  //determine the aray of the cell's pointer that will be used in the interpolation stencil
  double w[3],InterpolationWeight;
  int i,j,k,i0,j0,k0,nd;

  i0=(iLoc<0.5) ? -1 : (int)(iLoc-0.50);
  j0=(jLoc<0.5) ? -1 : (int)(jLoc-0.50);
  k0=(kLoc<0.5) ? -1 : (int)(kLoc-0.50);

  //get coefficients used in determening of the interpolation weights
  w[0]=iLoc-(i0+0.5);
  w[1]=jLoc-(j0+0.5);
  w[2]=kLoc-(k0+0.5);

  //dropped-weight normalization criterion (see "NORMALIZATION" at the top of this file)
  bool RowDropped=false,PointerDropped=false;

  //row-structure reuse (row-only construction): the row is determined by (node,i0,j0,k0)
  int RowMode=0,iContribution=0;

  #ifndef __CUDA_ARCH__
  cRowStructure *RowStructure=NULL;

  if ((BuildPointerStencil==false)&&(RowStencil!=NULL)) {
    RowStructure=GetRowStructure(RowStructureTrilinear,node,i0,j0,k0);
    RowMode=StartRowStructure(RowStructure);
  }
  #endif

  for (i=0;i<2;i++) for (j=0;j<2;j++) for (k=0;k<2;k++) {
    nd=_getCenterNodeLocalNumber(i0+i,j0+j,k0+k);
    switch (i+2*j+4*k) {
      case 0+0*2+0*4:
        InterpolationWeight=(1.0-w[0])*(1.0-w[1])*(1.0-w[2]);
        break;
      case 1+0*2+0*4:
        InterpolationWeight=w[0]*(1.0-w[1])*(1.0-w[2]);
        break;
      case 0+1*2+0*4:
        InterpolationWeight=(1.0-w[0])*w[1]*(1.0-w[2]);
        break;
      case 1+1*2+0*4:
        InterpolationWeight=w[0]*w[1]*(1.0-w[2]);
        break;

      case 0+0*2+1*4:
        InterpolationWeight=(1.0-w[0])*(1.0-w[1])*w[2];
        break;
      case 1+0*2+1*4:
        InterpolationWeight=w[0]*(1.0-w[1])*w[2];
        break;
      case 0+1*2+1*4:
        InterpolationWeight=(1.0-w[0])*w[1]*w[2];
        break;
      case 1+1*2+1*4:
        InterpolationWeight=w[0]*w[1]*w[2];
        break;

      default:
        exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    /* Populate the topological row even if this rank has no cDataCenterNode for
     * the cell. Ghost indices are canonicalized to the owning leaf node. */
    #ifndef __CUDA_ARCH__
    if (RowMode==RowModeReplay) {
      if (ApplyRowContribution(RowStructure,iContribution,InterpolationWeight,RowStencil)==true) RowDropped=true;
    }
    else
    #endif
    {
      const int LengthBefore=(RowStencil!=NULL) ? RowStencil->Length : 0;
      const bool added=AddRowStencilCell(RowStencil,InterpolationWeight,node,i0+i,j0+j,k0+k);

      if ((RowStencil!=NULL)&&(added==false)&&(InterpolationWeight>0.0)) RowDropped=true;

      #ifndef __CUDA_ARCH__
      if (RowMode==RowModeRecord) RecordRowContribution(RowStructure,RowStencil,LengthBefore,added,node,i0+i,j0+j,k0+k);
      #else
      (void)LengthBefore;
      #endif
    }

    iContribution++;

    cell=(block==NULL) ? NULL : block->GetCenterNode(nd);

    if (cell!=NULL) {
      if ((Stencil.AddCell(InterpolationWeight,cell,nd)==false)&&(InterpolationWeight>0.0)) PointerDropped=true;
    }
    else if (InterpolationWeight>0.0) PointerDropped=true;
  }

  #ifndef __CUDA_ARCH__
  if (RowMode==RowModeRecord) EndRowStructure(RowStructure,true);
  #endif
  (void)iContribution;

  /* Normalize the row when a positive weight was dropped (support point outside
   * the global box or in a block not used in the calculation). Remote allocation
   * never affects this normalization. */
  if ((RowStencil!=NULL)&&(RowStencil->Length>0)&&(RowDropped==true)) RowStencil->Normalize();

  if (Stencil.Length==0) {
    //no local cell has been found -> retain the legacy constant fallback when possible
    if (block!=NULL) PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node,Stencil);
    return;
  }
  else if (PointerDropped==true) {
    //a positive weight was dropped (cell not available locally or outside the domain) -> renormalize
    Stencil.Normalize();
  }
}



_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CellCentered::Linear::GetTriliniarInterpolationMutiBlockStencil(double *x,double *xStencilMin,double *xStencilMax,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CellCentered::cStencil &Stencil,PIC::InterpolationRoutines::cRowStencil *RowStencil,bool BuildPointerStencil) {
  long int nd;
  PIC::Mesh::cDataCenterNode *cell;

  Stencil.flush();
  if (RowStencil!=NULL) RowStencil->flush();

  /*
   * xStencilMin and xStencilMax are retained in the public interface for
   * backward compatibility. The existing AMPS algorithm reconstructs the same
   * limits from x and node, so these arguments remain intentionally unused.
   */
  (void)xStencilMin;
  (void)xStencilMax;

  const int nCells[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
  double dxCell[3];
  int idim;

  for (idim=0;idim<3;idim++) dxCell[idim]=(node->xmax[idim]-node->xmin[idim])/nCells[idim];

  //determine the lower logical coarse-center index and local coordinates
  int ijkStencilMin[3];
  double xStencilLower[3],xLoc[3];

  for (idim=0;idim<3;idim++) {
    ijkStencilMin[idim]=(x[idim]-node->xmin[idim]<0.5*dxCell[idim]) ? -1 : (int)((x[idim]-node->xmin[idim]-0.5*dxCell[idim])/dxCell[idim]);
    xStencilLower[idim]=node->xmin[idim]+(ijkStencilMin[idim]+0.5)*dxCell[idim];

    xLoc[idim]=(x[idim]-xStencilLower[idim])/dxCell[idim];
    if (xLoc[idim]<0.0) xLoc[idim]=0.0;
    if (xLoc[idim]>1.0) xLoc[idim]=1.0;
  }

  bool GeometryAvailable=true;
  bool PhysicalStencilAvailable=true;

  //dropped-weight normalization criterion (see "NORMALIZATION" at the top of this file)
  bool RowDropped=false,PointerDropped=false;

  //row-structure reuse (row-only construction): the row is determined by (node,ijkStencilMin)
  #ifndef __CUDA_ARCH__
  int RowMode=RowModeOriginal,iContribution=0;
  cRowStructure *RowStructure=NULL;

  if ((BuildPointerStencil==false)&&(RowStencil!=NULL)) {
    RowStructure=GetRowStructure(RowStructureMultiBlock,node,ijkStencilMin[0],ijkStencilMin[1],ijkStencilMin[2]);
    RowMode=StartRowStructure(RowStructure);
  }
  #endif

  /*
   * Visit the eight logical coarse-lattice support points in the same order as
   * the legacy implementation. The final trilinear coefficient is known here,
   * so row entries can be appended directly without allocating a second array of
   * full-size row stencils on the stack.
   */
  for (int di=0;di<2;di++) for (int dj=0;dj<2;dj++) for (int dk=0;dk<2;dk++) {
    int ijk[3]={ijkStencilMin[0]+di,ijkStencilMin[1]+dj,ijkStencilMin[2]+dk};
    double xLogical[3];

    for (idim=0;idim<3;idim++) xLogical[idim]=node->xmin[idim]+(ijk[idim]+0.5)*dxCell[idim];

    double StencilElementWeight;

    switch (di+2*dj+4*dk) {
    case 0+0*2+0*4:
      StencilElementWeight=(1.0-xLoc[0])*(1.0-xLoc[1])*(1.0-xLoc[2]);
      break;
    case 1+0*2+0*4:
      StencilElementWeight=xLoc[0]*(1.0-xLoc[1])*(1.0-xLoc[2]);
      break;
    case 0+1*2+0*4:
      StencilElementWeight=(1.0-xLoc[0])*xLoc[1]*(1.0-xLoc[2]);
      break;
    case 1+1*2+0*4:
      StencilElementWeight=xLoc[0]*xLoc[1]*(1.0-xLoc[2]);
      break;

    case 0+0*2+1*4:
      StencilElementWeight=(1.0-xLoc[0])*(1.0-xLoc[1])*xLoc[2];
      break;
    case 1+0*2+1*4:
      StencilElementWeight=xLoc[0]*(1.0-xLoc[1])*xLoc[2];
      break;
    case 0+1*2+1*4:
      StencilElementWeight=(1.0-xLoc[0])*xLoc[1]*xLoc[2];
      break;
    case 1+1*2+1*4:
      StencilElementWeight=xLoc[0]*xLoc[1]*xLoc[2];
      break;

    default:
      exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    #ifndef __CUDA_ARCH__
    const int iSupportPoint=4*di+2*dj+dk;

    if (RowMode==RowModeReplay) {
      if (RowStructure->fineSupport[iSupportPoint]==0) {
        if (ApplyRowContribution(RowStructure,iContribution++,StencilElementWeight,RowStencil)==true) RowDropped=true;
      }
      else {
        const double FineCellWeight=(1.0/8.0)*StencilElementWeight;

        for (int q=0;q<8;q++) if (ApplyRowContribution(RowStructure,iContribution++,FineCellWeight,RowStencil)==true) RowDropped=true;
      }

      continue;
    }
    #endif

    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *StencilNode=FindStencilSupportNode(xLogical,node,ijk);

    #ifndef __CUDA_ARCH__
    if (RowMode==RowModeRecord) AddRowStructureCheck(RowStructure,StencilNode);
    #endif

    if ((StencilNode==NULL)||(StencilNode->IsUsedInCalculationFlag==false)) {
      GeometryAvailable=false;
      continue;
    }

    if (StencilNode->RefinmentLevel==node->RefinmentLevel) {
      /*
       * Same-level data are read by the legacy implementation through the
       * reference block and its ghost cells. Keep that pointer path unchanged.
       * The row path canonicalizes the possibly-ghost index to the owning node.
       */
      {
        const int LengthBefore=(RowStencil!=NULL) ? RowStencil->Length : 0;
        const bool added=AddRowStencilCell(RowStencil,StencilElementWeight,node,ijk[0],ijk[1],ijk[2]);

        if ((RowStencil!=NULL)&&(added==false)&&(StencilElementWeight>0.0)) RowDropped=true;

        #ifndef __CUDA_ARCH__
        if (RowMode==RowModeRecord) {
          RowStructure->fineSupport[iSupportPoint]=0;
          RecordRowContribution(RowStructure,RowStencil,LengthBefore,added,node,ijk[0],ijk[1],ijk[2]);
        }
        #else
        (void)LengthBefore,(void)added;
        #endif
      }

      nd=_getCenterNodeLocalNumber(ijk[0],ijk[1],ijk[2]);
      cell=((BuildPointerStencil==false)||(node->block==NULL)) ? NULL : node->block->GetCenterNode(nd);

      if (cell!=NULL) {
        if ((AddPhysicalStencilCell(Stencil,StencilElementWeight,cell,nd)==false)&&(StencilElementWeight>0.0)) PointerDropped=true;
      }
      else PhysicalStencilAvailable=false;
    }
    else {
      /*
       * A logical coarse center covered by a finer block is represented by the
       * average of its 2x2x2 fine cells. This is the existing AMPS 2:1 AMR rule.
       * The row entries are created even if StencilNode->block is remote.
       */
      int iNeib[3];

      for (idim=0;idim<3;idim++) {
        iNeib[idim]=2*((int)((xLogical[idim]-StencilNode->xmin[idim])/dxCell[idim]));
      }

      for (int ii=0;ii<2;ii++) for (int jj=0;jj<2;jj++) for (int kk=0;kk<2;kk++) {
        const int iFine=iNeib[0]+ii;
        const int jFine=iNeib[1]+jj;
        const int kFine=iNeib[2]+kk;
        const double FineCellWeight=(1.0/8.0)*StencilElementWeight;

        {
          const int LengthBefore=(RowStencil!=NULL) ? RowStencil->Length : 0;
          const bool added=AddRowStencilCell(RowStencil,FineCellWeight,StencilNode,iFine,jFine,kFine);

          if ((RowStencil!=NULL)&&(added==false)&&(FineCellWeight>0.0)) RowDropped=true;

          #ifndef __CUDA_ARCH__
          if (RowMode==RowModeRecord) {
            RowStructure->fineSupport[iSupportPoint]=1;
            RecordRowContribution(RowStructure,RowStencil,LengthBefore,added,StencilNode,iFine,jFine,kFine);
          }
          #else
          (void)LengthBefore,(void)added;
          #endif
        }

        nd=_getCenterNodeLocalNumber(iFine,jFine,kFine);
        cell=((BuildPointerStencil==false)||(StencilNode->block==NULL)) ? NULL : StencilNode->block->GetCenterNode(nd);

        if (cell!=NULL) {
          if ((AddPhysicalStencilCell(Stencil,FineCellWeight,cell,nd)==false)&&(FineCellWeight>0.0)) PointerDropped=true;
        }
        else PhysicalStencilAvailable=false;
      }
    }
  }

  #ifndef __CUDA_ARCH__
  //a structure with a missing support point is not reused: the constant fallback depends on x
  if (RowMode==RowModeRecord) EndRowStructure(RowStructure,GeometryAvailable);
  #endif

  if (GeometryAvailable==false) {
    /* The old implementation falls back to constant interpolation when any
     * logical support point is outside the active mesh. Mirror that choice for
     * both representations. */
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *InterpolationNode=PIC::Mesh::mesh->findTreeNode(x,node);

    BuildConstantRowStencil(x,InterpolationNode,RowStencil);

    if ((BuildPointerStencil==true)&&(InterpolationNode!=NULL)&&(InterpolationNode->block!=NULL)) {
      PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,InterpolationNode,Stencil);
    }
    else Stencil.flush();

    return;
  }

  //a positive weight was dropped (a cell outside the global domain) -> renormalize
  if ((RowStencil!=NULL)&&(RowStencil->Length>0)&&(RowDropped==true)) RowStencil->Normalize();

  if ((PhysicalStencilAvailable==false)&&(BuildPointerStencil==true)) {
    /* Preserve the legacy constant fallback for local interpolation while the
     * optional row keeps the complete decomposition-independent linear stencil. */
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *InterpolationNode=PIC::Mesh::mesh->findTreeNode(x,node);

    if ((InterpolationNode!=NULL)&&(InterpolationNode->block!=NULL)) {
      PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,InterpolationNode,Stencil);
    }
    else Stencil.flush();
  }
  else if ((BuildPointerStencil==true)&&(Stencil.Length>0)&&(PointerDropped==true)) Stencil.Normalize();
}

//init stencil for the corner based interpolation
_TARGET_HOST_ _TARGET_DEVICE_
void PIC::InterpolationRoutines::CornerBased::InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,PIC::InterpolationRoutines::CornerBased::cStencil &Stencil,double *InterpolationCoefficientTable) {
  int iStencil,jStencil,kStencil,iX[3],nd,idim;
  double w,xLoc[3],dx[3],*xMinNode,*xMaxNode;
  PIC::Mesh::cDataCornerNode* CornerNode;
  PIC::Mesh::cDataBlockAMR *block;

  if (node==NULL) node=PIC::Mesh::mesh->findTreeNode(x);

  xMinNode=node->xmin;
  xMaxNode=node->xmax;

  //the size of the cells
  dx[0]=(xMaxNode[0]-xMinNode[0])/_BLOCK_CELLS_X_;
  dx[1]=(xMaxNode[1]-xMinNode[1])/_BLOCK_CELLS_Y_;
  dx[2]=(xMaxNode[2]-xMinNode[2])/_BLOCK_CELLS_Z_;

  //get the local coordinate for the interpolation point location
  #pragma ivdep
  for (idim=0;idim<3;idim++) {
    if ((x[idim]<xMinNode[idim])||(x[idim]>xMaxNode[idim])) exit(__LINE__,__FILE__,"Error: the point is out of block");
    if (fabs(x[idim]-xMaxNode[idim])<1e-10*dx[idim]) x[idim]= xMaxNode[idim]-1e-10*dx[idim];    
    xLoc[idim]=(x[idim]-xMinNode[idim])/dx[idim];
    iX[idim]=(int)(xLoc[idim]);
    xLoc[idim]-=iX[idim];
  }
  
  //build interpolation stencil
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  Stencil.flush();

  //return an empty stencil if the block is not allocated
  if ((block=node->block)==NULL) return; // PIC::InterpolationRoutines::CornerBased::StencilTable+ThreadOpenMP;

  for (iStencil=0;iStencil<2;iStencil++) for (jStencil=0;jStencil<2;jStencil++) for (kStencil=0;kStencil<2;kStencil++) {
    //get the local ID of the 'corner node'
    nd=_getCornerNodeLocalNumber(iStencil+iX[0],jStencil+iX[1],kStencil+iX[2]);
    CornerNode=block->GetCornerNode(nd);

    if (CornerNode!=NULL) {
      switch (iStencil+2*jStencil+4*kStencil) {
      case 0+0*2+0*4:
        w=(1.0-xLoc[0])*(1.0-xLoc[1])*(1.0-xLoc[2]);
        InterpolationCoefficientTable[0]=w;
        break;
      case 1+0*2+0*4:
        w=xLoc[0]*(1.0-xLoc[1])*(1.0-xLoc[2]);
        InterpolationCoefficientTable[1]=w;
        break;
      case 0+1*2+0*4:
        w=(1.0-xLoc[0])*xLoc[1]*(1.0-xLoc[2]);
        InterpolationCoefficientTable[3]=w;
        break;
      case 1+1*2+0*4:
        w=xLoc[0]*xLoc[1]*(1.0-xLoc[2]);
        InterpolationCoefficientTable[2]=w;
        break;

      case 0+0*2+1*4:
        w=(1.0-xLoc[0])*(1.0-xLoc[1])*xLoc[2];
        InterpolationCoefficientTable[4]=w;
        break;
      case 1+0*2+1*4:
        w=xLoc[0]*(1.0-xLoc[1])*xLoc[2];
        InterpolationCoefficientTable[5]=w;
        break;
      case 0+1*2+1*4:
        w=(1.0-xLoc[0])*xLoc[1]*xLoc[2];
        InterpolationCoefficientTable[7]=w;
        break;
      case 1+1*2+1*4:
        w=xLoc[0]*xLoc[1]*xLoc[2];
        InterpolationCoefficientTable[6]=w;
        break;

      default:
        exit(__LINE__,__FILE__,"Error: the option is not defined");
      }

      Stencil.AddCell(w,CornerNode,nd);
    }
    else {
      switch (iStencil+2*jStencil+4*kStencil) {
      case 0+0*2+0*4:
        InterpolationCoefficientTable[0]=0.0;
        break;
      case 1+0*2+0*4:
        InterpolationCoefficientTable[1]=0.0;
        break;
      case 0+1*2+0*4:
        InterpolationCoefficientTable[3]=0.0;
        break;
      case 1+1*2+0*4:
        InterpolationCoefficientTable[2]=0.0;
        break;

      case 0+0*2+1*4:
        InterpolationCoefficientTable[4]=0.0;
        break;
      case 1+0*2+1*4:
        InterpolationCoefficientTable[5]=0.0;
        break;
      case 0+1*2+1*4:
        InterpolationCoefficientTable[7]=0.0;
        break;
      case 1+1*2+1*4:
        InterpolationCoefficientTable[6]=0.0;
        break;

      default:
        exit(__LINE__,__FILE__,"Error: the option is not defined");
      }
    }
  }

  Stencil.Normalize();
}
