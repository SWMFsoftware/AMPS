#include "pic.h"   
#include <cassert>
#include <array>

namespace PIC {
namespace Mesh {

//=========================== Face table (1D/2D/3D) ============================
//
// Compile-time dimensionality via _MESH_DIMENSION_ (1,2,3).
// Faces are numbered 0..2*_MESH_DIMENSION_-1, matching mesh_boundary_blocks.cpp style:
//
//  face=0 → X-   (normal axis n=0, sign -)
//  face=1 → X+   (n=0, +)
//  face=2 → Y-   (n=1, -) if _MESH_DIMENSION_>=2
//  face=3 → Y+   (n=1, +) if _MESH_DIMENSION_>=2
//  face=4 → Z-   (n=2, -) if _MESH_DIMENSION_==3
//  face=5 → Z+   (n=2, +) if _MESH_DIMENSION_==3
//
// Each face also carries its two in-plane axes (t1, t2) for building (iFace,jFace).

struct FaceDesc {
  int n;    // normal axis (0,1,2)
  int t1;   // first in-plane axis
  int t2;   // second in-plane axis (=-1 in 1D)
  int bit;  // BoundaryFaceMask bit for this face
};

static inline void build_face_table(std::array<FaceDesc,6>& F, int& nFaces) {
  for (auto& e : F) { e = {-1,-1,-1, Face_None}; }

  // X- , X+
  F[0] = {0, (_MESH_DIMENSION_>=2?1:-1), (_MESH_DIMENSION_==3?2:-1), Face_XMin};
  F[1] = {0, (_MESH_DIMENSION_>=2?1:-1), (_MESH_DIMENSION_==3?2:-1), Face_XMax};

#if _MESH_DIMENSION_ >= 2
  // Y- , Y+
  F[2] = {1, 0, (_MESH_DIMENSION_==3?2:-1), Face_YMin};
  F[3] = {1, 0, (_MESH_DIMENSION_==3?2:-1), Face_YMax};
#endif

#if _MESH_DIMENSION_ == 3
  // Z- , Z+
  F[4] = {2, 0, 1, Face_ZMin};
  F[5] = {2, 0, 1, Face_ZMax};
#endif

  nFaces = 2 * _MESH_DIMENSION_;
}

//=========================== Utilities =======================================

static inline int clamp_to_face_index(int cornerIndex, int N) {
  // cornerIndex ∈ {0,N}; face-local valid range is [0, max(N-1,0)]
  return (cornerIndex == 0) ? 0 : (N > 0 ? (N - 1) : 0);
}

static inline std::array<int,3> cell_counts() {
  return { _BLOCK_CELLS_X_, _BLOCK_CELLS_Y_, _BLOCK_CELLS_Z_ };
}

static inline std::array<int,3> make_corner(int i, int j, int k) {
  return { i, j, k };
}

// Compute face-local (iFace,jFace) from a corner for a given FaceDesc.
static inline void corner_to_face_local(const FaceDesc& fd,
                                        const std::array<int,3>& N,
                                        const std::array<int,3>& idx,
                                        int& iFace, int& jFace)
{
  if (fd.t1 < 0) { iFace = 0; jFace = 0; return; }  // 1D: no in-plane axes

  const int iAxis = fd.t1;
  iFace = clamp_to_face_index(idx[iAxis], N[iAxis]);

  if (fd.t2 >= 0) {
    const int jAxis = fd.t2;
    jFace = clamp_to_face_index(idx[jAxis], N[jAxis]);
  } else {
    jFace = 0;
  }
}

// Compute which global faces the (i,j,k) corner of node 'nd' lies on.
// A face contributes if the corner sits on that block face AND the face has no neighbor (nullptr).
static inline int GetCornerBoundaryMask(cTreeNodeAMR<cDataBlockAMR>* nd, int i, int j, int k) {
  int mask = Face_None;

  std::array<FaceDesc,6> FT; int nFaces=0;
  build_face_table(FT, nFaces);

  const auto N   = cell_counts();
  const auto idx = make_corner(i,j,k);

  for (int f = 0; f < nFaces; ++f) {
    const auto& fd = FT[f];

    const int nAxis = fd.n;
    const int Nn    = N[nAxis];
    const int idxN  = idx[nAxis];

    const bool isMinFace = (f % 2 == 0); // even=min, odd=max
    const bool onFace =
      (isMinFace && idxN == 0) ||
      (!isMinFace && idxN == Nn);

    if (!onFace) continue;

    int iFace=0, jFace=0;
    corner_to_face_local(fd, N, idx, iFace, jFace);
    auto* neib = nd->GetNeibFace(f, iFace, jFace, PIC::Mesh::mesh);

    if (neib == nullptr) {
      mask |= fd.bit;
    }
  }

  return mask;
}

// Coordinates of corner (i,j,k); robust for 1D/2D/3D.
static inline void GetCornerPosition(cTreeNodeAMR<cDataBlockAMR>* nd, int i, int j, int k, double x[3]) {
  const double* xmin = nd->xmin;
  const double* xmax = nd->xmax;
  const int Nx = _BLOCK_CELLS_X_;
  const int Ny = _BLOCK_CELLS_Y_;
  const int Nz = _BLOCK_CELLS_Z_;

  const auto safe_div = [](double num, int den) -> double {
    return (den > 0) ? (num / den) : 0.0;
  };

  const double dx = safe_div(xmax[0] - xmin[0], Nx);
  const double dy = safe_div(xmax[1] - xmin[1], Ny);
  const double dz = safe_div(xmax[2] - xmin[2], Nz);

  x[0] = xmin[0] + i * dx;
  x[1] = xmin[1] + j * dy;
  x[2] = xmin[2] + k * dz;
}

//=========================== Public API ======================================

bool InitBoundaryCornerVector(std::vector<cBoundaryCornerInfo>& out) {
  out.clear();

  // Iterate *this thread's* nodes, as requested.
  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != nullptr; node = node->nextNodeThisThread)
  {
    // Safety: require payload block and leaf (parents can appear in distribution lists)
    if (node->block == nullptr) continue;
    if (node->lastBranchFlag() != _BOTTOM_BRANCH_TREE_) continue; // keep leaves only

    const int Nx = _BLOCK_CELLS_X_;
    const int Ny = _BLOCK_CELLS_Y_;
    const int Nz = _BLOCK_CELLS_Z_;

    const int iSet[2] = {0, Nx};
    const int jSet[2] = {0, Ny};
    const int kSet[2] = {0, Nz};

#if _MESH_DIMENSION_ == 1
    // 1D: (i in {0,Nx}), j=k=0
    for (int ii = 0; ii < 2; ++ii) {
      const int i = iSet[ii], j = 0, k = 0;
      const int faceMask = GetCornerBoundaryMask(node, i, j, k);
      if (faceMask == Face_None) continue;

      cBoundaryCornerInfo rec{};
      rec.node = node; rec.i = i; rec.j = j; rec.k = k; rec.faceMask = faceMask;
      GetCornerPosition(node, i, j, k, rec.x);
      out.push_back(rec);
    }

#elif _MESH_DIMENSION_ == 2
    // 2D: (i in {0,Nx}) × (j in {0,Ny}), k=0
    for (int ii = 0; ii < 2; ++ii) {
      for (int jj = 0; jj < 2; ++jj) {
        const int i = iSet[ii], j = jSet[jj], k = 0;
        const int faceMask = GetCornerBoundaryMask(node, i, j, k);
        if (faceMask == Face_None) continue;

        cBoundaryCornerInfo rec{};
        rec.node = node; rec.i = i; rec.j = j; rec.k = k; rec.faceMask = faceMask;
        GetCornerPosition(node, i, j, k, rec.x);
        out.push_back(rec);
      }
    }

#else // _MESH_DIMENSION_ == 3
    // 3D: (i in {0,Nx}) × (j in {0,Ny}) × (k in {0,Nz})
    for (int ii = 0; ii < 2; ++ii) {
      for (int jj = 0; jj < 2; ++jj) {
        for (int kk = 0; kk < 2; ++kk) {
          const int i = iSet[ii], j = jSet[jj], k = kSet[kk];
          const int faceMask = GetCornerBoundaryMask(node, i, j, k);
          if (faceMask == Face_None) continue;

          cBoundaryCornerInfo rec{};
          rec.node = node; rec.i = i; rec.j = j; rec.k = k; rec.faceMask = faceMask;
          GetCornerPosition(node, i, j, k, rec.x);
          out.push_back(rec);
        }
      }
    }
#endif
  }

  return true;
}

std::string FaceMaskToString(int mask) {
  std::string s;
  auto add = [&](const char* t){ if(!s.empty()) s += "|"; s += t; };
  if (mask & Face_XMin) add("X-");
  if (mask & Face_XMax) add("X+");
#if _MESH_DIMENSION_ >= 2
  if (mask & Face_YMin) add("Y-");
  if (mask & Face_YMax) add("Y+");
#endif
#if _MESH_DIMENSION_ == 3
  if (mask & Face_ZMin) add("Z-");
  if (mask & Face_ZMax) add("Z+");
#endif
  if (s.empty()) s = "None";
  return s;
}

bool IsOnFace(int faceMask, BoundaryFaceMask face) {
  return (faceMask & static_cast<int>(face)) != 0;
}

void DecodeBoundaryFaces(int faceMask, int& nBoundaryFaces, int BoundaryFaceList[3]) {
  nBoundaryFaces = 0;
  auto push = [&](int bit) {
    if ((faceMask & bit) && nBoundaryFaces < 3) BoundaryFaceList[nBoundaryFaces++] = bit;
  };

  // Deterministic order, pruned to compile-time dimension
  push(Face_XMin); push(Face_XMax);
#if _MESH_DIMENSION_ >= 2
  push(Face_YMin); push(Face_YMax);
#endif
#if _MESH_DIMENSION_ == 3
  push(Face_ZMin); push(Face_ZMax);
#endif
  // assert(nBoundaryFaces <= 3);
}

} // namespace Mesh
} // namespace PIC

