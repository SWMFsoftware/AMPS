#ifndef _SRC_GEOSPACE_GPU_ENGINE_MESH_H_
#define _SRC_GEOSPACE_GPU_ENGINE_MESH_H_

//======================================================================================
// EngineMesh.h -- device port of the Mode3D compact-field interpolation.
//
// The legacy Mode3D code evaluates B with
//
//   node = PIC::Mesh::mesh->findTreeNode(x,lastNode)
//   PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,cRowStencil)
//   B = sum_s weight_s * GlobalMagneticField_[3*(Temp_ID*Nc + i + Nx*(j + Ny*k))]
//
// The AMPS tree is a pointer structure on the host.  Here it is flattened into index
// arrays (EgMeshView) and the following AMPS routines are ported line by line:
//
//   cMeshAMRgeneric::findTreeNode(double*,start), findTreeNode(int*,start)
//   cTreeNodeAMR::neibNodeFace/neibNodeEdge/neibNodeCorner (integer-index form used
//       with _MESH_GLOBAL_NODE_CONNECTION_INFO_MODE_=_OFF_AMR_MESH_)
//   pic_interpolation_routines.cpp:
//       AddRowStencilCellByCoordinate, AddRowStencilCell, BuildConstantRowStencil,
//       CellCentered::Linear::InitStencil (_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_
//       ROUTINE__AMPS_, non-uniform AMR mesh, non-periodic),
//       GetTriliniarInterpolationStencil, GetTriliniarInterpolationMutiBlockStencil
//       (row-stencil parts only; the pointer stencil is not used by Mode3D tracing)
//   cRowStencil::AddCell/MultiplyScalar/Add/Normalize
//   GlobalMagneticField.cpp: CheckedGlobalCellIndex_, InterpolateField_
//
// The same operations in the same order are used, so the interpolated field is bitwise
// identical to the legacy tracer in double precision.  This includes the normalization
// rule of the row stencil (core "NORMALIZATION", 2026-10-23): a trilinear or multiblock
// row is renormalized exactly when a contribution with a positive weight was not added;
// Trilinear() and MultiBlock() must be kept in step with the core.  The tree search result does not depend
// on the start-node hint (the integer search climbs until the index box contains the
// point and then descends), so the hint is used only for speed.
//
// Block size: _BLOCK_CELLS_X/Y/Z_ of the AMPS build are passed in at run time through
// EgMeshView::nCells and checked against the compile-time constants on the host.
//======================================================================================

#include "EngineTypes.h"
#include <algorithm>
#include <climits>
#include <thread>
#include <vector>

namespace Geospace {
namespace GPU {

// Upper bound of a Mode3D row: multiblock (8 support points x 8 fine cells) blended
// with one trilinear fine stencil (8) is at most 72 distinct cells.  AMPS allows 128.
enum : int { kMaxRowStencil=80 };

// Stencil-cell cache (optional, see BuildStencilCache below).
//
// InitStencil resolves every stencil cell, given as (node,i,j,k) with i,j,k in
// [-1,nCells], to the leaf that owns the cell centre: AddCell() computes the centre and
// AddCellByCoordinate() runs a tree search and converts the centre to the owner's cell
// index; MultiBlock() runs the same tree search for its support points, and the coarse
// neighbour search calls NeibFace/NeibEdge/NeibCorner.  For a fixed mesh these results
// depend only on (node,i,j,k) (resp. node and the face/edge/corner), so they are
// computed once per used leaf with the same functions and looked up afterwards.  The
// stencil, its weights and the order of their accumulation are unchanged, so the
// interpolated field is bitwise identical with and without the cache.
//
//   interior cells (0<=i<nCells): resolve to the node itself; verified per leaf when
//                                 the cache is built (a leaf that fails is not cached)
//   ghost cells (one index -1 or nCells): one EgGhostRef per cell
struct EgGhostRef {
  int node;           // Find(centre,node): -1 when outside the tree
  signed char i,j,k;  // owner cell index (valid when ok)
  signed char ok;     // AddCellByCoordinate(centre) would add (node,i,j,k)
};
enum : int { kNeibPerLeaf=26 };   // 6 faces, 12 edges, 8 corners

template <class R>
struct EgMeshView {
  int nNodes;
  int root;
  int nCells[3];                 // _BLOCK_CELLS_X_, _Y_, _Z_
  R xGlobalMin[3];
  R xGlobalMax[3];
  R dxMaxRefinement[3];

  const R* xmin;                 // [3*nNodes]
  const R* xmax;                 // [3*nNodes]
  const int* ixMin;              // xMinGlobalIndex [3*nNodes]
  const int* sizeIndex;          // NodeGeometricSizeIndex [nNodes]
  const int* up;                 // parent index or -1
  const int* down;               // [8*nNodes], -1 when absent
  const int* level;              // RefinmentLevel
  const int* minNeibLevel;       // minNeibRefinmentLevel
  const int* maxNeibLevel;       // maxNeibRefinmentLevel
  const unsigned char* used;     // IsUsedInCalculationFlag
  const long long* tempId;       // Temp_ID (compact leaf index or -1)

  long long usedLeafBlocks;      // GlobalUsedLeafBlocks_
  long long interiorCells;       // GlobalInteriorCellCount_
  const R* B;                    // [3*interiorCells] internal units
  const unsigned char* presence; // GlobalCellPresence_

  // Stencil-cell cache; all nullptr when absent (uncached path).
  const int* cacheSlot;          // [nNodes] leaf slot or -1
  const short* ghostSlot;        // [(nx+2)(ny+2)(nz+2)] ghost index or -1 (interior)
  const EgGhostRef* ghost;       // [nGhost*nSlots]
  const int* neib;               // [kNeibPerLeaf*nSlots] NeibFace 0-5, NeibEdge 6-17, NeibCorner 18-25
  int nGhost;
};

template <class R>
struct EgRowStencil {
  int Length;
  int node[kMaxRowStencil];
  int i[kMaxRowStencil];
  int j[kMaxRowStencil];
  int k[kMaxRowStencil];
  R w[kMaxRowStencil];

  EG_HD void flush() { Length=0; }

  // cRowStencil::AddCell: merge an identical (node,i,j,k) entry, else append.
  EG_HD void AddCell(R weight,int n,int ii,int jj,int kk,int& error) {
    if (weight<R(0)) { error=kErrNegativeWeight; return; }
    for (int e=0;e<Length;e++) {
      if (node[e]==n && i[e]==ii && j[e]==jj && k[e]==kk) {
        w[e]+=weight;
        return;
      }
    }
    if (Length==kMaxRowStencil) { error=kErrRowStencilOverflow; return; }
    node[Length]=n; i[Length]=ii; j[Length]=jj; k[Length]=kk; w[Length]=weight;
    Length++;
  }

  EG_HD void MultiplyScalar(R a) {
    for (int e=0;e<Length;e++) w[e]*=a;
  }

  EG_HD void Add(const EgRowStencil& t,int& error) {
    for (int e=0;e<t.Length;e++) AddCell(t.w[e],t.node[e],t.i[e],t.j[e],t.k[e],error);
  }

  EG_HD R GetWeightSum() const {
    R res=R(0);
    for (int e=0;e<Length;e++) res+=w[e];
    return res;
  }

  EG_HD void Normalize() {
    const R norm=GetWeightSum();
    if (norm>R(0)) {
      for (int e=0;e<Length;e++) w[e]/=norm;
    }
  }
};

// x86-64 conversion semantics of (int)floor(v) used by AMPS: values outside the int
// range (including NaN) produce INT_MIN ("integer indefinite").  Making this explicit
// keeps the device result identical to the host for degenerate positions.
template <class R> EG_HD int FloorToIntX86(R v) {
  const R f=std::floor(v);
  if (!(f>=R(-2147483648.0) && f<R(2147483648.0))) return INT_MIN;
  return static_cast<int>(f);
}

template <class R>
class EgMesh {
 public:
  EG_HD explicit EgMesh(const EgMeshView<R>& m) : m_(m) {}

  // cMeshAMRgeneric::findTreeNode(int *ix,startNode)
  EG_HD int FindByIndex(const int ix[3],int start) const {
    int n=(start<0) ? m_.root : start;
    for (;;) {
      bool inblock=true;
      for (int d=0;d<3;d++) {
        const int lo=m_.ixMin[3*n+d];
        if ((ix[d]<lo)||(ix[d]>=lo+m_.sizeIndex[n])) { inblock=false; break; }
      }
      if (inblock) {
        const int half=m_.sizeIndex[n]/2;
        const int i=(ix[0]-m_.ixMin[3*n+0]<half) ? 0 : 1;
        const int j=(ix[1]-m_.ixMin[3*n+1]<half) ? 0 : 1;
        const int k=(ix[2]-m_.ixMin[3*n+2]<half) ? 0 : 1;
        const int t=m_.down[8*n+i+2*(j+2*k)];
        if (t>=0) { n=t; continue; }
        return n;
      }
      else {
        if (m_.up[n]>=0) { n=m_.up[n]; continue; }
        return -1;
      }
    }
  }

  // cMeshAMRgeneric::findTreeNode(double *x,startNode)
  EG_HD int Find(const R x[3],int start) const {
    int ix[3];
    for (int d=0;d<3;d++)
      ix[d]=FloorToIntX86<R>((x[d]-m_.xGlobalMin[d])/m_.dxMaxRefinement[d]);
    int res=FindByIndex(ix,start);
    bool flag=false;
    if (res>=0) for (int d=0;d<3;d++) {
      if (x[d]<m_.xmin[3*res+d]) ix[d]--,flag=true;
      if (x[d]>=m_.xmax[3*res+d]) ix[d]++,flag=true;
    }
    if (flag) res=FindByIndex(ix,res);
    return res;
  }

  // cTreeNodeAMR::neibNodeFace(iFace+2*(jFace+2*nface)) with iFace=jFace=0, i.e.
  // GetNeibFace(nface,0,0,mesh) as called by InitStencil.
  EG_HD int NeibFace(int n,int nface) const {
    int ix[3]={m_.ixMin[3*n+0],m_.ixMin[3*n+1],m_.ixMin[3*n+2]};
    const int s=m_.sizeIndex[n];
    switch (nface) {
      case 0: ix[0]-=1; break;
      case 1: ix[0]-=1; ix[0]+=1+s; break;
      case 2: ix[1]-=1; break;
      case 3: ix[1]-=1; ix[1]+=1+s; break;
      case 4: ix[2]-=1; break;
      case 5: ix[2]-=1; ix[2]+=1+s; break;
    }
    return FindByIndex(ix,n);
  }

  // cTreeNodeAMR::neibNodeEdge(iEdge+2*nedge) with iEdge=0 (GetNeibEdge(nedge,0,mesh)).
  EG_HD int NeibEdge(int n,int nedge) const {
    int inc[3];
    EdgeIncrement(nedge,inc);
    int ix[3];
    const int s=m_.sizeIndex[n];
    for (int d=0;d<3;d++) {
      ix[d]=m_.ixMin[3*n+d];
      if (inc[d]==-1) ix[d]-=1;
      else if (inc[d]==1) ix[d]+=s;
    }
    return FindByIndex(ix,n);
  }

  // cTreeNodeAMR::neibNodeCorner(i)
  EG_HD int NeibCorner(int n,int c) const {
    int ix[3];
    const int s=m_.sizeIndex[n];
    for (int d=0;d<3;d++) {
      const int bit=(c>>d)&1;
      ix[d]=bit ? m_.ixMin[3*n+d]+s : m_.ixMin[3*n+d]-1;
    }
    return FindByIndex(ix,n);
  }

  // AddRowStencilCellByCoordinate()
  EG_HD bool AddCellByCoordinate(EgRowStencil<R>& row,R weight,const R xIn[3],int start,
                                 int& error) const {
    int cellNode,ijk[3];
    if (!ResolveCellByCoordinate(xIn,start,cellNode,ijk)) return false;
    row.AddCell(weight,cellNode,ijk[0],ijk[1],ijk[2],error);
    return true;
  }

  // The owner leaf and cell index that AddRowStencilCellByCoordinate() adds for xIn;
  // false when it adds nothing.
  EG_HD bool ResolveCellByCoordinate(const R xIn[3],int start,int& cellNode,
                                     int ijk[3]) const {
    R x[3]={xIn[0],xIn[1],xIn[2]};
    for (int d=0;d<3;d++) {
      const R xmn=m_.xGlobalMin[d];
      const R xmx=m_.xGlobalMax[d];
      if ((x[d]<xmn)||(x[d]>xmx)) return false;
    }
    cellNode=Find(x,start);
    if ((cellNode<0)||(m_.used[cellNode]==0)) return false;

    // AMPS protects against round-off at a block face with a 1e-10 cell tolerance.  In
    // double this is kept literally (bitwise identity).  In float, 6+1e-10 rounds to
    // 6.0 and a coordinate that rounds onto the upper face would lose its stencil cell,
    // so the tolerance is scaled to the float resolution of the cell coordinate.
    const R faceTol=(sizeof(R)==sizeof(double)) ? R(1.0e-10)
        : R(64)*std::numeric_limits<R>::epsilon()*R(m_.nCells[0]);
    for (int d=0;d<3;d++) {
      const R blockLength=m_.xmax[3*cellNode+d]-m_.xmin[3*cellNode+d];
      const R xLoc=(x[d]-m_.xmin[3*cellNode+d])/blockLength*m_.nCells[d];
      ijk[d]=(int)std::floor(xLoc);
      if ((ijk[d]<0)&&(xLoc>-faceTol)) ijk[d]=0;
      if ((ijk[d]>=m_.nCells[d])&&(xLoc<m_.nCells[d]+faceTol)) ijk[d]=m_.nCells[d]-1;
      if ((ijk[d]<0)||(ijk[d]>=m_.nCells[d])) return false;
    }
    return true;
  }

  // Centre of cell (i,j,k) of sourceNode (AddRowStencilCell; the same expressions give
  // the support points xLogical of MultiBlock()).
  EG_HD void CellCentre(int sourceNode,int i,int j,int k,R x[3]) const {
    const int ijk[3]={i,j,k};
    for (int d=0;d<3;d++) {
      const R dx=(m_.xmax[3*sourceNode+d]-m_.xmin[3*sourceNode+d])/m_.nCells[d];
      x[d]=m_.xmin[3*sourceNode+d]+(ijk[d]+R(0.5))*dx;
    }
  }

  // Cache entry of cell (i,j,k) of `node`: nullptr when the node is not cached or the
  // index is outside [-1,nCells]; *interior=true (and nullptr) for an interior cell,
  // which resolves to (node,i,j,k) itself.
  EG_HD const EgGhostRef* CachedCell(int node,int i,int j,int k,bool& interior) const {
    interior=false;
    if (m_.cacheSlot==nullptr || node<0) return nullptr;
    const int slot=m_.cacheSlot[node];
    if (slot<0) return nullptr;
    const int nx=m_.nCells[0],ny=m_.nCells[1],nz=m_.nCells[2];
    if ((i<-1)||(i>nx)||(j<-1)||(j>ny)||(k<-1)||(k>nz)) return nullptr;
    const int g=m_.ghostSlot[(i+1)+(nx+2)*((j+1)+(ny+2)*(k+1))];
    if (g<0) { interior=true; return nullptr; }
    return m_.ghost+(long long)slot*m_.nGhost+g;
  }

  // AddRowStencilCell()
  EG_HD bool AddCell(EgRowStencil<R>& row,R weight,int sourceNode,int i,int j,int k,
                     int& error) const {
    if (sourceNode<0) return false;
    bool interior;
    const EgGhostRef* c=CachedCell(sourceNode,i,j,k,interior);
    if (interior) { row.AddCell(weight,sourceNode,i,j,k,error); return true; }
    if (c!=nullptr) {
      if (!c->ok) return false;
      row.AddCell(weight,c->node,c->i,c->j,c->k,error);
      return true;
    }
    R x[3];
    CellCentre(sourceNode,i,j,k,x);
    return AddCellByCoordinate(row,weight,x,sourceNode,error);
  }

  // Find(xLogical,node) for the centre xLogical of cell (i,j,k) of node (MultiBlock).
  EG_HD int FindCellCentre(int node,int i,int j,int k,const R xLogical[3]) const {
    bool interior;
    const EgGhostRef* c=CachedCell(node,i,j,k,interior);
    if (interior) return node;
    if (c!=nullptr) return c->node;
    return Find(xLogical,node);
  }

  // NeibFace/NeibEdge/NeibCorner through the cache (entries 0-5, 6-17, 18-25).
  EG_HD int CachedNeib(int n,int entry) const {
    if (m_.cacheSlot!=nullptr && n>=0) {
      const int slot=m_.cacheSlot[n];
      if (slot>=0) return m_.neib[(long long)slot*kNeibPerLeaf+entry];
    }
    if (entry<6) return NeibFace(n,entry);
    if (entry<18) return NeibEdge(n,entry-6);
    return NeibCorner(n,entry-18);
  }

  // BuildConstantRowStencil()
  EG_HD void BuildConstant(const R x[3],int node,EgRowStencil<R>& row,int& error) const {
    row.flush();
    AddCellByCoordinate(row,R(1.0),x,node,error);
  }

  // GetTriliniarInterpolationStencil() -- row part
  EG_HD void Trilinear(R iLoc,R jLoc,R kLoc,int node,EgRowStencil<R>& row,
                       int& error) const {
    row.flush();
    int i0,j0,k0;
    R w[3];
    i0=(iLoc<R(0.5)) ? -1 : (int)(iLoc-R(0.50));
    j0=(jLoc<R(0.5)) ? -1 : (int)(jLoc-R(0.50));
    k0=(kLoc<R(0.5)) ? -1 : (int)(kLoc-R(0.50));
    w[0]=iLoc-(i0+R(0.5));
    w[1]=jLoc-(j0+R(0.5));
    w[2]=kLoc-(k0+R(0.5));

    // Normalization -- the rule of the AMPS core (src/pic/pic_interpolation_routines.cpp,
    // "NORMALIZATION").  The 8 weights of a complete stencil sum to 1, so the row is
    // renormalized when, and only when, a contribution with a positive weight was not
    // added (AddCell() returned false: cell centre outside the global box, owning leaf
    // not used in the calculation, index out of range).  Entries that merged (ghost
    // cells resolving to the same coarser cell) do not trigger it.  The engine must
    // apply exactly the core's rule: it reproduces the legacy row stencil bitwise
    // (AMPS_EARTH_ENGINE_VERIFY, the Mode3D mesh self-check).
    bool dropped=false;
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) for (int k=0;k<2;k++) {
      R weight;
      switch (i+2*j+4*k) {
        case 0+0*2+0*4: weight=(R(1.0)-w[0])*(R(1.0)-w[1])*(R(1.0)-w[2]); break;
        case 1+0*2+0*4: weight=w[0]*(R(1.0)-w[1])*(R(1.0)-w[2]); break;
        case 0+1*2+0*4: weight=(R(1.0)-w[0])*w[1]*(R(1.0)-w[2]); break;
        case 1+1*2+0*4: weight=w[0]*w[1]*(R(1.0)-w[2]); break;
        case 0+0*2+1*4: weight=(R(1.0)-w[0])*(R(1.0)-w[1])*w[2]; break;
        case 1+0*2+1*4: weight=w[0]*(R(1.0)-w[1])*w[2]; break;
        case 0+1*2+1*4: weight=(R(1.0)-w[0])*w[1]*w[2]; break;
        default:        weight=w[0]*w[1]*w[2]; break;
      }
      if (!AddCell(row,weight,node,i0+i,j0+j,k0+k,error) && weight>R(0)) dropped=true;
    }
    if ((row.Length>0)&&dropped) row.Normalize();
  }

  // GetTriliniarInterpolationMutiBlockStencil() -- row part
  EG_HD void MultiBlock(const R x[3],int node,EgRowStencil<R>& row,int& error) const {
    row.flush();
    R dxCell[3];
    for (int d=0;d<3;d++) dxCell[d]=(m_.xmax[3*node+d]-m_.xmin[3*node+d])/m_.nCells[d];

    int ijkStencilMin[3];
    R xStencilLower[3],xLoc[3];
    for (int d=0;d<3;d++) {
      const R xmn=m_.xmin[3*node+d];
      ijkStencilMin[d]=(x[d]-xmn<R(0.5)*dxCell[d]) ? -1 :
          (int)((x[d]-xmn-R(0.5)*dxCell[d])/dxCell[d]);
      xStencilLower[d]=xmn+(ijkStencilMin[d]+R(0.5))*dxCell[d];
      xLoc[d]=(x[d]-xStencilLower[d])/dxCell[d];
      if (xLoc[d]<R(0.0)) xLoc[d]=R(0.0);
      if (xLoc[d]>R(1.0)) xLoc[d]=R(1.0);
    }

    bool geometryAvailable=true;
    // Normalization -- the rule of the AMPS core (see Trilinear()): the row is
    // renormalized when a contribution with a positive weight was not added.  A
    // missing support point (geometryAvailable==false) selects the constant fallback
    // instead, as in the core.
    bool dropped=false;
    for (int di=0;di<2;di++) for (int dj=0;dj<2;dj++) for (int dk=0;dk<2;dk++) {
      const int ijk[3]={ijkStencilMin[0]+di,ijkStencilMin[1]+dj,ijkStencilMin[2]+dk};
      R xLogical[3];
      for (int d=0;d<3;d++) xLogical[d]=m_.xmin[3*node+d]+(ijk[d]+R(0.5))*dxCell[d];

      R weight;
      switch (di+2*dj+4*dk) {
        case 0+0*2+0*4: weight=(R(1.0)-xLoc[0])*(R(1.0)-xLoc[1])*(R(1.0)-xLoc[2]); break;
        case 1+0*2+0*4: weight=xLoc[0]*(R(1.0)-xLoc[1])*(R(1.0)-xLoc[2]); break;
        case 0+1*2+0*4: weight=(R(1.0)-xLoc[0])*xLoc[1]*(R(1.0)-xLoc[2]); break;
        case 1+1*2+0*4: weight=xLoc[0]*xLoc[1]*(R(1.0)-xLoc[2]); break;
        case 0+0*2+1*4: weight=(R(1.0)-xLoc[0])*(R(1.0)-xLoc[1])*xLoc[2]; break;
        case 1+0*2+1*4: weight=xLoc[0]*(R(1.0)-xLoc[1])*xLoc[2]; break;
        case 0+1*2+1*4: weight=(R(1.0)-xLoc[0])*xLoc[1]*xLoc[2]; break;
        default:        weight=xLoc[0]*xLoc[1]*xLoc[2]; break;
      }

      // xLogical is CellCentre(node,ijk): the same expressions, hence the same bits.
      const int stencilNode=FindCellCentre(node,ijk[0],ijk[1],ijk[2],xLogical);
      if ((stencilNode<0)||(m_.used[stencilNode]==0)) {
        geometryAvailable=false;
        continue;
      }

      if (m_.level[stencilNode]==m_.level[node]) {
        if (!AddCell(row,weight,node,ijk[0],ijk[1],ijk[2],error) && weight>R(0)) dropped=true;
      }
      else {
        int iNeib[3];
        for (int d=0;d<3;d++)
          iNeib[d]=2*((int)((xLogical[d]-m_.xmin[3*stencilNode+d])/dxCell[d]));
        for (int ii=0;ii<2;ii++) for (int jj=0;jj<2;jj++) for (int kk=0;kk<2;kk++) {
          const R fineCellWeight=(R(1.0)/R(8.0))*weight;
          if (!AddCell(row,fineCellWeight,stencilNode,iNeib[0]+ii,iNeib[1]+jj,iNeib[2]+kk,error) &&
              fineCellWeight>R(0)) dropped=true;
        }
      }
    }

    if (!geometryAvailable) {
      const int interpolationNode=Find(x,node);
      BuildConstant(x,interpolationNode,row,error);
    }
    else if ((row.Length>0)&&dropped) row.Normalize();
  }

  // CellCentered::Linear::InitStencil() -- AMPS routine, row part, node!=NULL.
  EG_HD void InitStencil(const R x[3],int node,EgRowStencil<R>& row,
                         EgRowStencil<R>& fineRow,int& error) const {
    row.flush();
    const R* xmn=m_.xmin+3*node;
    const R* xmx=m_.xmax+3*node;
    const int nx=m_.nCells[0],ny=m_.nCells[1],nz=m_.nCells[2];

    const R iLoc=(x[0]-xmn[0])/(xmx[0]-xmn[0])*nx;
    const R jLoc=(x[1]-xmn[1])/(xmx[1]-xmn[1])*ny;
    const R kLoc=(x[2]-xmn[2])/(xmx[2]-xmn[2])*nz;

    const int lev=m_.level[node];
    if ((lev==m_.minNeibLevel[node])&&(lev==m_.maxNeibLevel[node])) {
      Trilinear(iLoc,jLoc,kLoc,node,row,error);
      return;
    }
    else if ((R(1.0)<iLoc)&&(iLoc<nx-1) && (R(1.0)<jLoc)&&(jLoc<ny-1) &&
             (R(1.0)<kLoc)&&(kLoc<nz-1)) {
      Trilinear(iLoc,jLoc,kLoc,node,row,error);
      return;
    }
    else if (lev==m_.minNeibLevel[node]) {
      if ((R(0.5)<iLoc)&&(iLoc<nx-R(0.5)) && (R(0.5)<jLoc)&&(jLoc<ny-R(0.5)) &&
          (R(0.5)<kLoc)&&(kLoc<nz-R(0.5))) {
        Trilinear(iLoc,jLoc,kLoc,node,row,error);
        return;
      }
      else {
        MultiBlock(x,node,row,error);
        return;
      }
    }
    else {
      int coarserBlock=-1;
      R dxCell[3];
      dxCell[0]=(xmx[0]-xmn[0])/nx;
      dxCell[1]=(xmx[1]-xmn[1])/ny;
      dxCell[2]=(xmx[2]-xmn[2])/nz;

      const int nBlockCells[3]={nx,ny,nz};
      R dmin=R(10.0*nx*ny*nz);
      bool cornerTested[8]={false,false,false,false,false,false,false,false};
      bool edgeTested[12]={false,false,false,false,false,false,
                           false,false,false,false,false,false};
      const R xLoc[3]={iLoc,jLoc,kLoc};

      for (int idim=0;idim<3;idim++) {
        int iFace;
        switch (idim) {
          case 0:
            if (iLoc<=R(1.0)) iFace=0;
            else if (iLoc>=nx-R(1.0)) iFace=1;
            else continue;
            break;
          case 1:
            if (jLoc<=R(1.0)) iFace=2;
            else if (jLoc>=ny-R(1.0)) iFace=3;
            else continue;
            break;
          default:
            if (kLoc<=R(1.0)) iFace=4;
            else if (kLoc>=nz-R(1.0)) iFace=5;
            else continue;
            break;
        }

        // blocks connected through the face
        int neib=CachedNeib(node,iFace);
        if (neib>=0) if ((m_.level[neib]<lev)&&(m_.used[neib]!=0)) {
          if (Covers(neib,x,dxCell)) {
            switch (iFace) {
              case 0:
                if ((xLoc[0]<R(1.0))&&(xLoc[0]<dmin)) dmin=xLoc[0],coarserBlock=neib;
                break;
              case 1:
                if ((xLoc[0]>nBlockCells[0]-1)&&(nBlockCells[0]-xLoc[0]<dmin))
                  dmin=nBlockCells[0]-xLoc[0],coarserBlock=neib;
                break;
              case 2:
                if ((xLoc[1]<R(1.0))&&(xLoc[1]<dmin)) dmin=xLoc[1],coarserBlock=neib;
                break;
              case 3:
                if ((xLoc[1]>nBlockCells[1]-1)&&(nBlockCells[1]-xLoc[1]<dmin))
                  dmin=nBlockCells[1]-xLoc[1],coarserBlock=neib;
                break;
              case 4:
                if ((xLoc[2]<R(1.0))&&(xLoc[2]<dmin)) dmin=xLoc[2],coarserBlock=neib;
                break;
              case 5:
                if ((xLoc[2]>nBlockCells[2]-1)&&(nBlockCells[2]-xLoc[2]<dmin))
                  dmin=nBlockCells[2]-xLoc[2],coarserBlock=neib;
                break;
            }
          }
        }

        // blocks connected through the edges
        for (int iEdge=0;iEdge<4;iEdge++) {
          const int edge=FaceEdge(iFace,iEdge);
          if (edgeTested[edge]) continue;
          edgeTested[edge]=true;
          neib=CachedNeib(node,6+edge);
          if (neib>=0) if ((m_.level[neib]<lev)&&(m_.used[neib]!=0)) {
            if (Covers(neib,x,dxCell)) EdgeCandidate(edge,xLoc,nBlockCells,neib,dmin,coarserBlock);
          }
        }

        // blocks connected through the corners
        for (int iCorner=0;iCorner<4;iCorner++) {
          const int corner=FaceCorner(iFace,iCorner);
          if (cornerTested[corner]) continue;
          cornerTested[corner]=true;
          neib=CachedNeib(node,18+corner);
          if (neib>=0) if ((m_.level[neib]<lev)&&(m_.used[neib]!=0)) {
            if (Covers(neib,x,dxCell)) CornerCandidate(corner,xLoc,nBlockCells,neib,dmin,coarserBlock);
          }
        }
      }

      if (coarserBlock>=0) {
        MultiBlock(x,coarserBlock,row,error);
        if ((R(0.5)<dmin)&&(dmin<=R(1.0))) {
          Trilinear(iLoc,jLoc,kLoc,node,fineRow,error);
          row.MultiplyScalar(R(1.0)-(dmin-R(0.5))/R(0.5));
          fineRow.MultiplyScalar((dmin-R(0.5))/R(0.5));
          row.Add(fineRow,error);
        }
        return;
      }

      Trilinear(iLoc,jLoc,kLoc,node,row,error);
      return;
    }
  }

  // CheckedGlobalCellIndex_() + InterpolateField_() accumulation.
  // Returns false when the row is empty (legacy: InterpolateMagneticField()==false).
  // Fast path of Interpolate() for the common case: InitStencil() selects the
  // trilinear stencil of `node` (its first three branches, not the coarser-neighbour
  // search) and all 8 cells lie inside `node`, a used leaf.  Then every AddCell() of the
  // general path finds `node` itself (the cell centre is inside it) and the cell's own
  // index, the 8 entries are distinct and added in loop order, and the stencil is not
  // normalized.  B is therefore accumulated here with the same weights in the same order
  // and with the same consistency checks, but without the 8 tree searches and without
  // building the stencil in (GPU local) memory.  Returns false when the case does not
  // apply (or a weight is negative); the caller then runs the general path.
  EG_HD bool InterpolateInteriorTrilinear(const R x[3],int node,R B[3],int& error) const {
    if (m_.used[node]==0) return false;
    const R* xmn=m_.xmin+3*node;
    const R* xmx=m_.xmax+3*node;
    const int nx=m_.nCells[0],ny=m_.nCells[1],nz=m_.nCells[2];
    const R iLoc=(x[0]-xmn[0])/(xmx[0]-xmn[0])*nx;
    const R jLoc=(x[1]-xmn[1])/(xmx[1]-xmn[1])*ny;
    const R kLoc=(x[2]-xmn[2])/(xmx[2]-xmn[2])*nz;
    const int lev=m_.level[node];
    const bool trilinear=
        ((lev==m_.minNeibLevel[node])&&(lev==m_.maxNeibLevel[node])) ||
        ((R(1.0)<iLoc)&&(iLoc<nx-1) && (R(1.0)<jLoc)&&(jLoc<ny-1) &&
         (R(1.0)<kLoc)&&(kLoc<nz-1)) ||
        ((lev==m_.minNeibLevel[node]) &&
         (R(0.5)<iLoc)&&(iLoc<nx-R(0.5)) && (R(0.5)<jLoc)&&(jLoc<ny-R(0.5)) &&
         (R(0.5)<kLoc)&&(kLoc<nz-R(0.5)));
    if (!trilinear) return false;
    // Trilinear(): lower cell indices and weights, same expressions.
    const int i0=(iLoc<R(0.5)) ? -1 : (int)(iLoc-R(0.50));
    const int j0=(jLoc<R(0.5)) ? -1 : (int)(jLoc-R(0.50));
    const int k0=(kLoc<R(0.5)) ? -1 : (int)(kLoc-R(0.50));
    if ((i0<0)||(i0+1>=nx)||(j0<0)||(j0+1>=ny)||(k0<0)||(k0+1>=nz)) return false;
    R w[3];
    w[0]=iLoc-(i0+R(0.5));
    w[1]=jLoc-(j0+R(0.5));
    w[2]=kLoc-(k0+R(0.5));
    R weight[8];
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) for (int k=0;k<2;k++) {
      const int e=i+2*j+4*k;
      R wt;
      switch (e) {
        case 0+0*2+0*4: wt=(R(1.0)-w[0])*(R(1.0)-w[1])*(R(1.0)-w[2]); break;
        case 1+0*2+0*4: wt=w[0]*(R(1.0)-w[1])*(R(1.0)-w[2]); break;
        case 0+1*2+0*4: wt=(R(1.0)-w[0])*w[1]*(R(1.0)-w[2]); break;
        case 1+1*2+0*4: wt=w[0]*w[1]*(R(1.0)-w[2]); break;
        case 0+0*2+1*4: wt=(R(1.0)-w[0])*(R(1.0)-w[1])*w[2]; break;
        case 1+0*2+1*4: wt=w[0]*(R(1.0)-w[1])*w[2]; break;
        case 0+1*2+1*4: wt=(R(1.0)-w[0])*w[1]*w[2]; break;
        default:        wt=w[0]*w[1]*w[2]; break;
      }
      if (wt<R(0)) return false;   // the general path reports kErrNegativeWeight
      weight[e]=wt;
    }

    const long long tid=m_.tempId[node];
    if ((tid<0)||(tid>=m_.usedLeafBlocks)) { error=kErrMeshConsistency; return true; }
    const long long nc=(long long)nx*ny*nz;
    B[0]=B[1]=B[2]=R(0);
    // Entry order of the general path: the i/j/k loops of Trilinear().
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) for (int k=0;k<2;k++) {
      const R wt=weight[i+2*j+4*k];
      const long long cell=tid*nc+(long long)(i0+i)+(long long)nx*
                           ((long long)(j0+j)+(long long)ny*(long long)(k0+k));
      if ((cell<0)||(cell>=m_.interiorCells)||(m_.presence[cell]!=1)) {
        error=kErrMeshConsistency; return true;
      }
      const long long v=3*cell;
      B[0]+=wt*m_.B[v+0];
      B[1]+=wt*m_.B[v+1];
      B[2]+=wt*m_.B[v+2];
    }
    return true;
  }

  EG_HD bool Interpolate(const R x[3],int node,R B[3],EgRowStencil<R>& row,
                         EgRowStencil<R>& fineRow,int& error) const {
    if (InterpolateInteriorTrilinear(x,node,B,error)) return error==kErrNone;
    InitStencil(x,node,row,fineRow,error);
    if (row.Length<=0) return false;
    B[0]=B[1]=B[2]=R(0);
    const long long nc=(long long)m_.nCells[0]*m_.nCells[1]*m_.nCells[2];
    for (int s=0;s<row.Length;s++) {
      const int n=row.node[s];
      const int i=row.i[s],j=row.j[s],k=row.k[s];
      if ((i<0)||(i>=m_.nCells[0])||(j<0)||(j>=m_.nCells[1])||(k<0)||(k>=m_.nCells[2])) {
        error=kErrMeshConsistency; return false;
      }
      const long long tid=m_.tempId[n];
      if ((tid<0)||(tid>=m_.usedLeafBlocks)) { error=kErrMeshConsistency; return false; }
      const long long cell=tid*nc+(long long)i+(long long)m_.nCells[0]*
                           ((long long)j+(long long)m_.nCells[1]*(long long)k);
      if ((cell<0)||(cell>=m_.interiorCells)||(m_.presence[cell]!=1)) {
        error=kErrMeshConsistency; return false;
      }
      const long long v=3*cell;
      B[0]+=row.w[s]*m_.B[v+0];
      B[1]+=row.w[s]*m_.B[v+1];
      B[2]+=row.w[s]*m_.B[v+2];
    }
    return true;
  }

 private:
  EG_HD bool Covers(int neib,const R x[3],const R dxCell[3]) const {
    int cnt=0;
    for (int ii=0;ii<3;ii++)
      if ((m_.xmin[3*neib+ii]-dxCell[ii]<=x[ii]) && (m_.xmax[3*neib+ii]+dxCell[ii]>=x[ii])) cnt++;
    return cnt==3;
  }

  static EG_HD void EdgeIncrement(int edge,int inc[3]) {
    // Increment[12][3] of cTreeNodeAMR::neibNodeEdge
    switch (edge) {
      case 0:  inc[0]=0;  inc[1]=-1; inc[2]=-1; break;
      case 1:  inc[0]=0;  inc[1]=1;  inc[2]=-1; break;
      case 2:  inc[0]=0;  inc[1]=1;  inc[2]=1;  break;
      case 3:  inc[0]=0;  inc[1]=-1; inc[2]=1;  break;
      case 4:  inc[0]=-1; inc[1]=0;  inc[2]=-1; break;
      case 5:  inc[0]=1;  inc[1]=0;  inc[2]=-1; break;
      case 6:  inc[0]=1;  inc[1]=0;  inc[2]=1;  break;
      case 7:  inc[0]=-1; inc[1]=0;  inc[2]=1;  break;
      case 8:  inc[0]=-1; inc[1]=-1; inc[2]=0;  break;
      case 9:  inc[0]=1;  inc[1]=-1; inc[2]=0;  break;
      case 10: inc[0]=1;  inc[1]=1;  inc[2]=0;  break;
      default: inc[0]=-1; inc[1]=1;  inc[2]=0;  break;
    }
  }

  // faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}}
  static EG_HD int FaceEdge(int face,int n) {
    const int t[24]={4,11,7,8, 5,10,6,9, 0,9,3,8, 1,10,2,11, 0,5,1,4, 3,6,2,7};
    return t[4*face+n];
  }

  // FaceNodeMap[6][4]={{0,2,4,6},{1,3,5,7},{0,1,4,5},{2,3,6,7},{0,1,2,3},{4,5,6,7}}
  static EG_HD int FaceCorner(int face,int n) {
    const int t[24]={0,2,4,6, 1,3,5,7, 0,1,4,5, 2,3,6,7, 0,1,2,3, 4,5,6,7};
    return t[4*face+n];
  }

  static EG_HD void Low(R v,R& dmin,int& cb,int neib) { if (v<dmin) dmin=v,cb=neib; }

  static EG_HD void EdgeCandidate(int edge,const R xLoc[3],const int nb[3],int neib,
                                  R& dmin,int& cb) {
    const R x0=xLoc[0],x1=xLoc[1],x2=xLoc[2];
    const R h0=nb[0]-x0,h1=nb[1]-x1,h2=nb[2]-x2;
    switch (edge) {
      case 0: if ((x1<R(1.0))&&(x2<R(1.0))) { Low(x1,dmin,cb,neib); Low(x2,dmin,cb,neib); } break;
      case 1: if ((x1<R(1.0))&&(x2>nb[2]-1)) { Low(x1,dmin,cb,neib); Low(h2,dmin,cb,neib); } break;
      case 2: if ((x1>nb[1]-1)&&(x2>nb[2]-1)) { Low(h1,dmin,cb,neib); Low(h2,dmin,cb,neib); } break;
      case 3: if ((x1>nb[1]-1)&&(x2<R(1.0))) { Low(h1,dmin,cb,neib); Low(x2,dmin,cb,neib); } break;
      case 4: if ((x0<R(1.0))&&(x2<R(1.0))) { Low(x0,dmin,cb,neib); Low(x2,dmin,cb,neib); } break;
      case 5: if ((x0>nb[0]-1)&&(x2<R(1.0))) { Low(h0,dmin,cb,neib); Low(x2,dmin,cb,neib); } break;
      case 6: if ((x0>nb[0]-1)&&(x2>nb[2]-1)) { Low(h0,dmin,cb,neib); Low(h2,dmin,cb,neib); } break;
      case 7: if ((x0<R(1.0))&&(x2>nb[2]-1)) { Low(x0,dmin,cb,neib); Low(h2,dmin,cb,neib); } break;
      case 8: if ((x0<R(1.0))&&(x1<R(1.0))) { Low(x0,dmin,cb,neib); Low(x1,dmin,cb,neib); } break;
      case 9: if ((x0>nb[0]-1)&&(x1<R(1.0))) { Low(h0,dmin,cb,neib); Low(x1,dmin,cb,neib); } break;
      case 10: if ((x0>nb[0]-1)&&(x1>nb[1]-1)) { Low(h0,dmin,cb,neib); Low(h1,dmin,cb,neib); } break;
      default: if ((x0<R(1.0))&&(x1>nb[1]-1)) { Low(x0,dmin,cb,neib); Low(h1,dmin,cb,neib); } break;
    }
  }

  static EG_HD void CornerCandidate(int corner,const R xLoc[3],const int nb[3],int neib,
                                    R& dmin,int& cb) {
    // Corner c: bit d set -> upper side in dimension d.
    bool ok=true;
    R v[3];
    for (int d=0;d<3;d++) {
      if ((corner>>d)&1) { ok=ok && (xLoc[d]>nb[d]-1); v[d]=nb[d]-xLoc[d]; }
      else { ok=ok && (xLoc[d]<R(1.0)); v[d]=xLoc[d]; }
    }
    if (ok) { Low(v[0],dmin,cb,neib); Low(v[1],dmin,cb,neib); Low(v[2],dmin,cb,neib); }
  }

  EgMeshView<R> m_;
};

// Host storage of the stencil-cell cache.
struct EgStencilCacheHost {
  std::vector<int> cacheSlot;      // [nNodes]
  std::vector<short> ghostSlot;    // [(nx+2)(ny+2)(nz+2)]
  std::vector<EgGhostRef> ghost;   // [nGhost*nSlots]
  std::vector<int> neib;           // [kNeibPerLeaf*nSlots]
  int nGhost=0;
  long long nSlots=0;              // cached leaves
  long long nRejected=0;           // used leaves left uncached (interior check failed)
};

// Build the stencil-cell cache of every used leaf of `v` (host).  Every entry is
// computed by the uncached EgMesh functions on the same arrays, so a lookup returns
// exactly what the uncached path computes.  A leaf is cached only when every interior
// cell resolves to the leaf itself with the same index (the interior shortcut).
template <class R>
void BuildStencilCache(EgMeshView<R> v,EgStencilCacheHost& c,int nThreads) {
  c=EgStencilCacheHost{};
  v.cacheSlot=nullptr; v.ghostSlot=nullptr; v.ghost=nullptr; v.neib=nullptr; v.nGhost=0;
  const int nx=v.nCells[0],ny=v.nCells[1],nz=v.nCells[2];
  if (v.nNodes<=0 || nx<1 || ny<1 || nz<1 || nx>127 || ny>127 || nz>127) return;
  const EgMesh<R> mesh(v);

  c.ghostSlot.assign((std::size_t)(nx+2)*(ny+2)*(nz+2),-1);
  std::vector<int> gi,gj,gk;
  for (int k=-1;k<=nz;++k) for (int j=-1;j<=ny;++j) for (int i=-1;i<=nx;++i) {
    const bool interior=(i>=0 && i<nx && j>=0 && j<ny && k>=0 && k<nz);
    if (interior) continue;
    c.ghostSlot[(i+1)+(nx+2)*((j+1)+(ny+2)*(k+1))]=static_cast<short>(gi.size());
    gi.push_back(i); gj.push_back(j); gk.push_back(k);
  }
  if (gi.size()>32767) { c=EgStencilCacheHost{}; return; }
  c.nGhost=static_cast<int>(gi.size());

  // Candidate leaves in node order.
  c.cacheSlot.assign((std::size_t)v.nNodes,-1);
  std::vector<int> nodeOf;
  for (int n=0;n<v.nNodes;++n) {
    if (v.used[n]==0) continue;
    bool leaf=true;
    for (int q=0;q<8;++q) if (v.down[8*n+q]>=0) { leaf=false; break; }
    if (!leaf) continue;
    c.cacheSlot[(std::size_t)n]=static_cast<int>(nodeOf.size());
    nodeOf.push_back(n);
  }
  c.nSlots=static_cast<long long>(nodeOf.size());
  c.ghost.assign((std::size_t)c.nSlots*c.nGhost,EgGhostRef{-1,0,0,0,0});
  c.neib.assign((std::size_t)c.nSlots*kNeibPerLeaf,-1);
  std::vector<char> reject((std::size_t)c.nSlots,0);

  auto work=[&](long long s0,long long s1) {
    for (long long s=s0;s<s1;++s) {
      const int n=nodeOf[(std::size_t)s];
      bool ok=true;
      for (int k=0;k<nz && ok;++k) for (int j=0;j<ny && ok;++j) for (int i=0;i<nx && ok;++i) {
        R x[3];
        int cn,ijk[3];
        mesh.CellCentre(n,i,j,k,x);
        ok=mesh.ResolveCellByCoordinate(x,n,cn,ijk) && cn==n &&
           ijk[0]==i && ijk[1]==j && ijk[2]==k;
      }
      EgGhostRef* g=c.ghost.data()+s*c.nGhost;
      for (int e=0;e<c.nGhost && ok;++e) {
        R x[3];
        int cn,ijk[3];
        mesh.CellCentre(n,gi[e],gj[e],gk[e],x);
        g[e].node=mesh.Find(x,n);
        const bool added=mesh.ResolveCellByCoordinate(x,n,cn,ijk);
        if (added) {
          // ResolveCellByCoordinate runs the same Find(x,n).
          if (cn!=g[e].node) { ok=false; break; }
          g[e].i=static_cast<signed char>(ijk[0]);
          g[e].j=static_cast<signed char>(ijk[1]);
          g[e].k=static_cast<signed char>(ijk[2]);
        }
        g[e].ok=added ? 1 : 0;
      }
      int* nb=c.neib.data()+s*kNeibPerLeaf;
      for (int f=0;f<6;++f) nb[f]=mesh.NeibFace(n,f);
      for (int e=0;e<12;++e) nb[6+e]=mesh.NeibEdge(n,e);
      for (int q=0;q<8;++q) nb[18+q]=mesh.NeibCorner(n,q);
      if (!ok) reject[(std::size_t)s]=1;
    }
  };
  nThreads=std::max(1,std::min<int>(nThreads,64));
  if (nThreads==1 || c.nSlots<1024) work(0,c.nSlots);
  else {
    std::vector<std::thread> th;
    const long long chunk=(c.nSlots+nThreads-1)/nThreads;
    for (int t=0;t<nThreads;++t) {
      const long long s0=t*chunk,s1=std::min(c.nSlots,s0+chunk);
      if (s0<s1) th.emplace_back(work,s0,s1);
    }
    for (auto& t:th) t.join();
  }
  for (long long s=0;s<c.nSlots;++s)
    if (reject[(std::size_t)s]) { c.cacheSlot[(std::size_t)nodeOf[(std::size_t)s]]=-1; ++c.nRejected; }
}

} // namespace GPU
} // namespace Geospace

#endif
