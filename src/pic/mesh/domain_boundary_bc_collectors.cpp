/***************************************************************************************
*  domain_boundary_bc_collectors.cpp  (C++14-safe)
*
*  Collect and mark domain-boundary face-layer and corner cells (Dirichlet & Neumann).
*  - Robust domain check via node->xmin/xmax vs mesh->xGlobalMin/Max with eps.
*  - Works in 1D/2D/3D, leaf blocks only.
*  - Appends records to std::vector<cBoundaryCellInfo>.
*  - Uses your BC APIs:
*      Dirichlet: cell->SetBCTypeDirichlet(PIC::Mesh::BCTypeDirichlet);
*      Neumann  : cell->SetBCTypeNeumann (PIC::Mesh::BCTypeNeumann, i_in, j_in, k_in);
*
*  PURPOSE
*  -------
*  Provide utilities to collect and mark (1) *face-layer* cells and (2) *corner* cells
*  that lie on the **PHYSICAL DOMAIN BOUNDARY** of an AMR block assigned to the current
*  MPI rank. Two BC flavors are supported:
*
*    • Dirichlet: set a fixed field value at boundary cells.
*      -> cell->SetBCTypeDirichlet(PIC::Mesh::BCTypeDirichlet);
*
*    • Neumann: set a normal-derivative (flux) BC at boundary cells, referencing a
*      cell **inside** the domain to define ∂f/∂n.
*      -> cell->SetBCTypeNeumann(BCTypeDirichlet,i_in, j_in, k_in);
*
*  For every selected cell we append a record {i,j,k, node*, cell*} into a user-supplied
*  vector. Selection works in 1D/2D/3D and is restricted to **requested domain faces**.
*
*  WHY “DOMAIN” (NOT JUST BLOCK EDGE)?
*  -----------------------------------
*  A block face may abut another block, ghost layer, or periodic image. We therefore
*  test domain boundary status by comparing the block’s physical extents
*  (node->xmin[], node->xmax[]) to global extents (mesh->xGlobalMin[], xGlobalMax[])
*  with a tiny tolerance (eps) to absorb roundoff. This is robust to AMR and MPI layout.
*
*  FACE INDEXING (matches block faces)
*  -----------------------------------
*     0 : -X (left)     1 : +X (right)
*     2 : -Y (bottom)   3 : +Y (top)        [if _MESH_DIMENSION_ >= 2]
*     4 : -Z (back)     5 : +Z (front)      [if _MESH_DIMENSION_ == 3]
*
*  FUNCTIONS
*  ---------
*  (Dirichlet)
*    A1. CollectAndMarkDomainBoundaryDirichletCellsOnFaces(faces, nFaces, out)
*        Select the **face-layer** of center cells on requested domain faces and mark BC.
*
*    A2. CollectAndMarkDomainBoundaryDirichletCornerCells(faces, nFaces, out)
*        Select **corner** center cells formed *only* by requested faces (e.g., {0,3,5})
*        and mark BC.
*
*  (Neumann)
*    B1. CollectAndMarkDomainBoundaryNeumannCellsOnFaces(faces, nFaces, out)
*        Same selection as A1, but set **Neumann** BC. The inside reference index
*        (i_in, j_in, k_in) is one cell inward along the face normal:
*          -X: (1,  j,  k)   +X: (NX-2, j,   k)
*          -Y: (i,  1,  k)   +Y: (i,   NY-2, k)
*          -Z: (i,  j,  1)   +Z: (i,   j,    NZ-2)
*        (Guards ensure NX,NY,NZ >= 2 in the respective directions.)
*
*    B2. CollectAndMarkDomainBoundaryNeumannCornerCells(faces, nFaces, out)
*        Corner selection as in A2, but compute a single **inside** reference index
*        that is one cell inward along *each* boundary direction participating in the
*        corner:
*          i_in = ( corner uses -X ? 1 : corner uses +X ? NX-2 : i_b )
*          j_in = ( corner uses -Y ? 1 : corner uses +Y ? NY-2 : j_b )
*          k_in = ( corner uses -Z ? 1 : corner uses +Z ? NZ-2 : k_b )
*        where (i_b, j_b, k_b) is the boundary corner cell index
*        (e.g., (-X,+Y,+Z) → (0, NY-1, NZ-1)).
*
*  NUMERICAL TOLERANCE
*  -------------------
*  epsX = 1e-12 * (1 + |xGlobalMax[0] - xGlobalMin[0]|), similarly for Y,Z.
*
*  COMPLEXITY
*  ----------
*  Faces: O(B * face_area_cells). Corners: O(B). Per-block dedup uses a byte map.
*
*  EXAMPLE
*  -------
*    std::vector<cBoundaryCellInfo> faceD, cornerD, faceN, cornerN;
*    int faces[] = {0, 5}; // -X and +Z
*
*    PIC::Mesh::CollectAndMarkDomainBoundaryDirichletCellsOnFaces(faces, 2, faceD);
*    PIC::Mesh::CollectAndMarkDomainBoundaryDirichletCornerCells(faces, 2, cornerD);
*
*    PIC::Mesh::CollectAndMarkDomainBoundaryNeumannCellsOnFaces(faces, 2, faceN);
*    PIC::Mesh::CollectAndMarkDomainBoundaryNeumannCornerCells(faces, 2, cornerN);
*
*    // After this:
*    //  • faceD/faceN contain all face-layer boundary cells on -X and +Z.
*    //  • cornerD/cornerN contain cells at (-X,+Z) (2D) or (-X,+Z) edges/corners (3D)
*    //    that are made only from the requested faces; Neumann entries also store the
*    //    interior index used via SetBCTypeNeumann(BCTypeDirichlet,i_in, j_in, k_in);.
*
****************************************************************************************/

#include <vector>
#include <cmath>
#include <cstdio>

#include "../pic.h"                     // AMR mesh types and PIC::Mesh::mesh

namespace PIC { namespace Mesh {

// Pull exact pointer types from your existing record type; this prevents mismatches.
using Record = PIC::Mesh::cBoundaryCellInfo;
using NodePtr = decltype(Record{}.node);
using CellPtr = decltype(Record{}.cell);

namespace { // helpers
  inline bool approx_le(double a, double b, double eps) { return (a - b) <= eps; }
  inline bool approx_ge(double a, double b, double eps) { return (b - a) <= eps; }

  struct DomainFaceFlags { bool Xm=false, Xp=false, Ym=false, Yp=false, Zm=false, Zp=false; };

  // Accept any ::cTreeNodeAMR<BlockT>* (C++14-friendly).
  template<class BlockT>
  inline DomainFaceFlags computeDomainFaceFlags(const ::cTreeNodeAMR<BlockT>* node) {
    DomainFaceFlags f;

    const double Lx  = PIC::Mesh::mesh->xGlobalMax[0] - PIC::Mesh::mesh->xGlobalMin[0];
    const double epsX = 1e-12 * (1.0 + std::fabs(Lx));
    f.Xm = approx_le(node->xmin[0], PIC::Mesh::mesh->xGlobalMin[0], epsX);
    f.Xp = approx_ge(node->xmax[0], PIC::Mesh::mesh->xGlobalMax[0], epsX);

#if _MESH_DIMENSION_ >= 2
    const double Ly  = PIC::Mesh::mesh->xGlobalMax[1] - PIC::Mesh::mesh->xGlobalMin[1];
    const double epsY = 1e-12 * (1.0 + std::fabs(Ly));
    f.Ym = approx_le(node->xmin[1], PIC::Mesh::mesh->xGlobalMin[1], epsY);
    f.Yp = approx_ge(node->xmax[1], PIC::Mesh::mesh->xGlobalMax[1], epsY);
#endif
#if _MESH_DIMENSION_ == 3
    const double Lz  = PIC::Mesh::mesh->xGlobalMax[2] - PIC::Mesh::mesh->xGlobalMin[2];
    const double epsZ = 1e-12 * (1.0 + std::fabs(Lz));
    f.Zm = approx_le(node->xmin[2], PIC::Mesh::mesh->xGlobalMin[2], epsZ);
    f.Zp = approx_ge(node->xmax[2], PIC::Mesh::mesh->xGlobalMax[2], epsZ);
#endif
    return f;
  }

  inline bool hasFace(const int* faces, int nFaces, int f) {
    for (int i = 0; i < nFaces; ++i) if (faces[i] == f) return true;
    return false;
  }
} // anon

// ============================================================================
// A1) Dirichlet on requested DOMAIN faces: face-layer cells
// ============================================================================
bool CollectAndMarkDomainBoundaryDirichletCellsOnFaces(const int* faces, int nFaces,
                                                       std::vector<Record>& out)
{
  out.clear();
  if (!faces || nFaces <= 0) { std::printf("$PREFIX: DirichletFaces: no faces\n"); return false; }
  if (!PIC::Mesh::mesh || !PIC::Mesh::mesh->ParallelNodesDistributionList) {
    std::printf("$PREFIX: DirichletFaces: mesh not initialized\n"); return false; }

  constexpr int NX = _BLOCK_CELLS_X_;
  constexpr int NY = _BLOCK_CELLS_Y_;
  constexpr int NZ = _BLOCK_CELLS_Z_;
  std::size_t added = 0;

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != nullptr; node = node->nextNodeThisThread)
  {
    if (!node->block) continue;
    if (node->lastBranchFlag() != _BOTTOM_BRANCH_TREE_) continue;

    const auto f = computeDomainFaceFlags(node);

    bool touches = false;
    for (int i = 0; i < nFaces && !touches; ++i) {
      switch (faces[i]) {
        case 0: touches |= f.Xm; break; case 1: touches |= f.Xp; break;
#if _MESH_DIMENSION_ >= 2
        case 2: touches |= f.Ym; break; case 3: touches |= f.Yp; break;
#endif
#if _MESH_DIMENSION_ == 3
        case 4: touches |= f.Zm; break; case 5: touches |= f.Zp; break;
#endif
        default: break;
      }
    }
    if (!touches) continue;

    std::vector<unsigned char> seen(static_cast<size_t>(NX) * NY * NZ, 0);
    auto push_once = [&](int i, int j, int k, CellPtr cell) {
      const size_t idx = static_cast<size_t>(i) + NX * (static_cast<size_t>(j) + static_cast<size_t>(NY) * static_cast<size_t>(k));
      if (seen[idx]) return; seen[idx] = 1;

      cell->SetBCTypeDirichlet(PIC::Mesh::BCTypeDirichlet);

      Record rec;
      rec.i = i; rec.j = j; rec.k = k;
      rec.node = static_cast<NodePtr>(node);
      rec.cell = cell;
      out.push_back(rec);
      ++added;
    };

    for (int fi = 0; fi < nFaces; ++fi) {
      switch (faces[fi]) {
        case 0: // -X
          if (!f.Xm) break;
          for (int k=0;k<NZ;++k) for (int j=0;j<NY;++j)
            if (auto* c = node->block->GetCenterNode(0,j,k)) push_once(0,j,k, static_cast<CellPtr>(c));
          break;

        case 1: // +X
          if (!f.Xp) break;
          for (int k=0;k<NZ;++k) for (int j=0;j<NY;++j) {
            const int i = NX-1;
            if (auto* c = node->block->GetCenterNode(i,j,k)) push_once(i,j,k, static_cast<CellPtr>(c));
          }
          break;

#if _MESH_DIMENSION_ >= 2
        case 2: // -Y
          if (!f.Ym) break;
          for (int k=0;k<NZ;++k) for (int i=0;i<NX;++i)
            if (auto* c = node->block->GetCenterNode(i,0,k)) push_once(i,0,k, static_cast<CellPtr>(c));
          break;

        case 3: // +Y
          if (!f.Yp) break;
          for (int k=0;k<NZ;++k) for (int i=0;i<NX;++i) {
            const int j = NY-1;
            if (auto* c = node->block->GetCenterNode(i,j,k)) push_once(i,j,k, static_cast<CellPtr>(c));
          }
          break;
#endif

#if _MESH_DIMENSION_ == 3
        case 4: // -Z
          if (!f.Zm) break;
          for (int j=0;j<NY;++j) for (int i=0;i<NX;++i)
            if (auto* c = node->block->GetCenterNode(i,j,0)) push_once(i,j,0, static_cast<CellPtr>(c));
          break;

        case 5: // +Z
          if (!f.Zp) break;
          for (int j=0;j<NY;++j) for (int i=0;i<NX;++i) {
            const int k = NZ-1;
            if (auto* c = node->block->GetCenterNode(i,j,k)) push_once(i,j,k, static_cast<CellPtr>(c));
          }
          break;
#endif
        default: break;
      }
    }
  }

  std::printf("$PREFIX: DirichletFaces: marked %zu cells on rank %d\n", added, PIC::ThisThread);
  return true;
}

// ============================================================================
// A2) Dirichlet on requested DOMAIN corners: single cell per corner
// ============================================================================
bool CollectAndMarkDomainBoundaryDirichletCornerNodesOnFaces(const int* faces, int nFaces,
                                                      std::vector<Record>& out)
{
  out.clear();
  if (!faces || nFaces <= 0) { std::printf("$PREFIX: DirichletCorners: no faces\n"); return false; }
  if (!PIC::Mesh::mesh || !PIC::Mesh::mesh->ParallelNodesDistributionList) {
    std::printf("$PREFIX: DirichletCorners: mesh not initialized\n"); return false; }

  constexpr int NX = _BLOCK_CELLS_X_;
  constexpr int NY = _BLOCK_CELLS_Y_;
  constexpr int NZ = _BLOCK_CELLS_Z_;
  std::size_t added = 0;

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != nullptr; node = node->nextNodeThisThread)
  {
    if (!node->block) continue;
    if (node->lastBranchFlag() != _BOTTOM_BRANCH_TREE_) continue;

    const auto f = computeDomainFaceFlags(node);

    std::vector<unsigned char> seen(static_cast<size_t>(NX) * NY * NZ, 0);
    auto push_corner_once = [&](int i, int j, int k) {
      const size_t idx = static_cast<size_t>(i) + NX * (static_cast<size_t>(j) + static_cast<size_t>(NY) * static_cast<size_t>(k));
      if (seen[idx]) return; seen[idx] = 1;
      if (auto* c = node->block->GetCenterNode(i,j,k)) {
        auto* cell = static_cast<CellPtr>(c);
        cell->SetBCTypeDirichlet(PIC::Mesh::BCTypeDirichlet);

        Record rec;
        rec.i = i; rec.j = j; rec.k = k;
        rec.node = static_cast<NodePtr>(node);
        rec.cell = cell;
        out.push_back(rec);
        ++added;
      }
    };

#if _MESH_DIMENSION_ == 1
    if (f.Xm && hasFace(faces,nFaces,0)) push_corner_once(0,0,0);
    if (f.Xp && hasFace(faces,nFaces,1)) push_corner_once(NX-1,0,0);
#endif

#if _MESH_DIMENSION_ == 2
    if (f.Xm&&f.Ym&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,2)) push_corner_once(0,0,0);
    if (f.Xm&&f.Yp&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,3)) push_corner_once(0,NY-1,0);
    if (f.Xp&&f.Ym&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,2)) push_corner_once(NX-1,0,0);
    if (f.Xp&&f.Yp&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,3)) push_corner_once(NX-1,NY-1,0);
#endif

#if _MESH_DIMENSION_ == 3
    if (f.Xm&&f.Ym&&f.Zm&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,4)) push_corner_once(0,0,0);
    if (f.Xm&&f.Ym&&f.Zp&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,5)) push_corner_once(0,0,NZ-1);
    if (f.Xm&&f.Yp&&f.Zm&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,4)) push_corner_once(0,NY-1,0);
    if (f.Xm&&f.Yp&&f.Zp&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,5)) push_corner_once(0,NY-1,NZ-1);
    if (f.Xp&&f.Ym&&f.Zm&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,4)) push_corner_once(NX-1,0,0);
    if (f.Xp&&f.Ym&&f.Zp&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,5)) push_corner_once(NX-1,0,NZ-1);
    if (f.Xp&&f.Yp&&f.Zm&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,4)) push_corner_once(NX-1,NY-1,0);
    if (f.Xp&&f.Yp&&f.Zp&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,5)) push_corner_once(NX-1,NY-1,NZ-1);
#endif
  }

  std::printf("$PREFIX: DirichletCorners: marked %zu cells on rank %d\n", added, PIC::ThisThread);
  return true;
}

// ============================================================================
// B1) Neumann on requested DOMAIN faces: face-layer cells with inside ref idx
// ============================================================================
bool CollectAndMarkDomainBoundaryNeumannCellsOnFaces(const int* faces, int nFaces,
                                                     std::vector<Record>& out)
{
  out.clear();
  if (!faces || nFaces <= 0) { std::printf("$PREFIX: NeumannFaces: no faces\n"); return false; }
  if (!PIC::Mesh::mesh || !PIC::Mesh::mesh->ParallelNodesDistributionList) {
    std::printf("$PREFIX: NeumannFaces: mesh not initialized\n"); return false; }

  constexpr int NX = _BLOCK_CELLS_X_;
  constexpr int NY = _BLOCK_CELLS_Y_;
  constexpr int NZ = _BLOCK_CELLS_Z_;

  if (NX < 2 || (_MESH_DIMENSION_ >= 2 && NY < 2) || (_MESH_DIMENSION_ == 3 && NZ < 2)) {
    std::printf("$PREFIX: NeumannFaces: block too thin to choose inside indices\n");
    return false;
  }

  std::size_t added = 0;

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != nullptr; node = node->nextNodeThisThread)
  {
    if (!node->block) continue;
    if (node->lastBranchFlag() != _BOTTOM_BRANCH_TREE_) continue;

    const auto f = computeDomainFaceFlags(node);

    bool touches = false;
    for (int i = 0; i < nFaces && !touches; ++i) {
      switch (faces[i]) {
        case 0: touches |= f.Xm; break; case 1: touches |= f.Xp; break;
#if _MESH_DIMENSION_ >= 2
        case 2: touches |= f.Ym; break; case 3: touches |= f.Yp; break;
#endif
#if _MESH_DIMENSION_ == 3
        case 4: touches |= f.Zm; break; case 5: touches |= f.Zp; break;
#endif
        default: break;
      }
    }
    if (!touches) continue;

    std::vector<unsigned char> seen(static_cast<size_t>(NX) * NY * NZ, 0);
    auto push_once = [&](int i, int j, int k, int ii, int jj, int kk, CellPtr cell) {
      const size_t idx = static_cast<size_t>(i) + NX * (static_cast<size_t>(j) + static_cast<size_t>(NY) * static_cast<size_t>(k));
      if (seen[idx]) return; seen[idx] = 1;

      cell->SetBCTypeNeumann(PIC::Mesh::BCTypeNeumann, ii, jj, kk);

      Record rec;
      rec.i = i; rec.j = j; rec.k = k;
      rec.node = static_cast<NodePtr>(node);
      rec.cell = cell;
      out.push_back(rec);
      ++added;
    };

    for (int fi = 0; fi < nFaces; ++fi) {
      switch (faces[fi]) {
        case 0: // -X: boundary i=0, inside i=1
          if (!f.Xm) break;
          for (int k=0;k<NZ;++k) for (int j=0;j<NY;++j)
            if (auto* c = node->block->GetCenterNode(0,j,k))
              push_once(0,j,k, 1, j, k, static_cast<CellPtr>(c));
          break;

        case 1: // +X: boundary i=NX-1, inside i=NX-2
          if (!f.Xp) break;
          for (int k=0;k<NZ;++k) for (int j=0;j<NY;++j) {
            const int i = NX-1, ii = NX-2;
            if (auto* c = node->block->GetCenterNode(i,j,k))
              push_once(i,j,k, ii, j, k, static_cast<CellPtr>(c));
          }
          break;

#if _MESH_DIMENSION_ >= 2
        case 2: // -Y: boundary j=0, inside j=1
          if (!f.Ym) break;
          for (int k=0;k<NZ;++k) for (int i=0;i<NX;++i)
            if (auto* c = node->block->GetCenterNode(i,0,k))
              push_once(i,0,k, i, 1, k, static_cast<CellPtr>(c));
          break;

        case 3: // +Y: boundary j=NY-1, inside j=NY-2
          if (!f.Yp) break;
          for (int k=0;k<NZ;++k) for (int i=0;i<NX;++i) {
            const int j = NY-1, jj = NY-2;
            if (auto* c = node->block->GetCenterNode(i,j,k))
              push_once(i,j,k, i, jj, k, static_cast<CellPtr>(c));
          }
          break;
#endif

#if _MESH_DIMENSION_ == 3
        case 4: // -Z: boundary k=0, inside k=1
          if (!f.Zm) break;
          for (int j=0;j<NY;++j) for (int i=0;i<NX;++i)
            if (auto* c = node->block->GetCenterNode(i,j,0))
              push_once(i,j,0, i, j, 1, static_cast<CellPtr>(c));
          break;

        case 5: // +Z: boundary k=NZ-1, inside k=NZ-2
          if (!f.Zp) break;
          for (int j=0;j<NY;++j) for (int i=0;i<NX;++i) {
            const int k = NZ-1, kk = NZ-2;
            if (auto* c = node->block->GetCenterNode(i,j,k))
              push_once(i,j,k, i, j, kk, static_cast<CellPtr>(c));
          }
          break;
#endif
        default: break;
      }
    }
  }

  std::printf("$PREFIX: NeumannFaces: marked %zu cells on rank %d\n", added, PIC::ThisThread);
  return true;
}

// ============================================================================
// B2) Neumann on requested DOMAIN corners: one inside ref index per corner
// ============================================================================
bool CollectAndMarkDomainBoundaryNeumannCornerNodesOnFaces(const int* faces, int nFaces,
                                                    std::vector<Record>& out)
{
  out.clear();
  if (!faces || nFaces <= 0) { std::printf("$PREFIX: NeumannCorners: no faces\n"); return false; }
  if (!PIC::Mesh::mesh || !PIC::Mesh::mesh->ParallelNodesDistributionList) {
    std::printf("$PREFIX: NeumannCorners: mesh not initialized\n"); return false; }

  constexpr int NX = _BLOCK_CELLS_X_;
  constexpr int NY = _BLOCK_CELLS_Y_;
  constexpr int NZ = _BLOCK_CELLS_Z_;

  if (NX < 2 || (_MESH_DIMENSION_ >= 2 && NY < 2) || (_MESH_DIMENSION_ == 3 && NZ < 2)) {
    std::printf("$PREFIX: NeumannCorners: block too thin to choose inside indices\n");
    return false;
  }

  std::size_t added = 0;

  auto inward_i = [&](bool usesXm, bool usesXp, int i_b) {
    if (usesXm) return 1; if (usesXp) return NX-2; return i_b;
  };
#if _MESH_DIMENSION_ >= 2
  auto inward_j = [&](bool usesYm, bool usesYp, int j_b) {
    if (usesYm) return 1; if (usesYp) return NY-2; return j_b;
  };
#endif
#if _MESH_DIMENSION_ == 3
  auto inward_k = [&](bool usesZm, bool usesZp, int k_b) {
    if (usesZm) return 1; if (usesZp) return NZ-2; return k_b;
  };
#endif

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != nullptr; node = node->nextNodeThisThread)
  {
    if (!node->block) continue;
    if (node->lastBranchFlag() != _BOTTOM_BRANCH_TREE_) continue;

    const auto f = computeDomainFaceFlags(node);

    std::vector<unsigned char> seen(static_cast<size_t>(NX) * NY * NZ, 0);
    auto push_corner_once = [&](int i_b, int j_b, int k_b,
                                bool usesXm, bool usesXp,
                                bool usesYm, bool usesYp,
                                bool usesZm, bool usesZp)
    {
      const size_t idx = static_cast<size_t>(i_b) + NX * (static_cast<size_t>(j_b) + static_cast<size_t>(NY) * static_cast<size_t>(k_b));
      if (seen[idx]) return; seen[idx] = 1;
      if (auto* c = node->block->GetCenterNode(i_b, j_b, k_b)) {
        auto* cell = static_cast<CellPtr>(c);
        const int ii = inward_i(usesXm, usesXp, i_b);
#if _MESH_DIMENSION_ >= 2
        const int jj = inward_j(usesYm, usesYp, j_b);
#else
        const int jj = 0;
#endif
#if _MESH_DIMENSION_ == 3
        const int kk = inward_k(usesZm, usesZp, k_b);
#else
        const int kk = 0;
#endif
        cell->SetBCTypeNeumann(PIC::Mesh::BCTypeNeumann, ii, jj, kk);

        Record rec;
        rec.i = i_b; rec.j = j_b; rec.k = k_b;
        rec.node = static_cast<NodePtr>(node);
        rec.cell = cell;
        out.push_back(rec);
        ++added;
      }
    };

#if _MESH_DIMENSION_ == 1
    if (f.Xm && hasFace(faces,nFaces,0)) push_corner_once(0,0,0, /*Xm*/true,false, false,false, false,false);
    if (f.Xp && hasFace(faces,nFaces,1)) push_corner_once(NX-1,0,0, false,true,  false,false, false,false);
#endif

#if _MESH_DIMENSION_ == 2
    if (f.Xm&&f.Ym&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,2))
      push_corner_once(0,0,0, true,false, true,false, false,false);
    if (f.Xm&&f.Yp&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,3))
      push_corner_once(0,NY-1,0, true,false, false,true, false,false);
    if (f.Xp&&f.Ym&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,2))
      push_corner_once(NX-1,0,0, false,true, true,false, false,false);
    if (f.Xp&&f.Yp&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,3))
      push_corner_once(NX-1,NY-1,0, false,true, false,true, false,false);
#endif

#if _MESH_DIMENSION_ == 3
    if (f.Xm&&f.Ym&&f.Zm&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,4))
      push_corner_once(0,0,0, true,false,  true,false,  true,false);
    if (f.Xm&&f.Ym&&f.Zp&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,5))
      push_corner_once(0,0,NZ-1, true,false,  true,false,  false,true);
    if (f.Xm&&f.Yp&&f.Zm&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,4))
      push_corner_once(0,NY-1,0, true,false,  false,true,  true,false);
    if (f.Xm&&f.Yp&&f.Zp&&hasFace(faces,nFaces,0)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,5))
      push_corner_once(0,NY-1,NZ-1, true,false,  false,true,  false,true);
    if (f.Xp&&f.Ym&&f.Zm&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,4))
      push_corner_once(NX-1,0,0, false,true,  true,false,  true,false);
    if (f.Xp&&f.Ym&&f.Zp&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,2)&&hasFace(faces,nFaces,5))
      push_corner_once(NX-1,0,NZ-1, false,true,  true,false,  false,true);
    if (f.Xp&&f.Yp&&f.Zm&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,4))
      push_corner_once(NX-1,NY-1,0, false,true,  false,true,  true,false);
    if (f.Xp&&f.Yp&&f.Zp&&hasFace(faces,nFaces,1)&&hasFace(faces,nFaces,3)&&hasFace(faces,nFaces,5))
      push_corner_once(NX-1,NY-1,NZ-1, false,true,  false,true,  false,true);
#endif
  }

  std::printf("$PREFIX: NeumannCorners: marked %zu cells on rank %d\n", added, PIC::ThisThread);
  return true;
}

}} // namespace PIC::Mesh

