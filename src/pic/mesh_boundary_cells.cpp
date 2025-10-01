#include "pic.h"

namespace PIC {
namespace Mesh {

bool InitBoundaryCellVector(std::vector<cBoundaryCellInfo>* BoundaryCellVector) {
  if (BoundaryCellVector != nullptr) {
    BoundaryCellVector->clear();
  }

  if (PIC::Mesh::mesh == nullptr) {
    std::printf("$PREFIX: InitBoundaryCellVector: ERROR: mesh is null\n");
    return false;
  }
  if (PIC::Mesh::mesh->ParallelNodesDistributionList == nullptr) {
    std::printf("$PREFIX: InitBoundaryCellVector: ERROR: ParallelNodesDistributionList is null\n");
    return false;
  }

  // Faces & block cell counts by dimension
#if _MESH_DIMENSION_ == 1
  constexpr int nFaces = 2;  // -X, +X
  const int NX = _BLOCK_CELLS_X_;
  const int NY = 1;
  const int NZ = 1;
#elif _MESH_DIMENSION_ == 2
  constexpr int nFaces = 4;  // -X, +X, -Y, +Y
  const int NX = _BLOCK_CELLS_X_;
  const int NY = _BLOCK_CELLS_Y_;
  const int NZ = 1;
#else
  constexpr int nFaces = 6;  // -X, +X, -Y, +Y, -Z, +Z
  const int NX = _BLOCK_CELLS_X_;
  const int NY = _BLOCK_CELLS_Y_;
  const int NZ = _BLOCK_CELLS_Z_;
#endif

  const int STRIDE_J = NX;
  const int STRIDE_K = NX * NY;

  // Helper: is there any neighbor across face f for this node?
  auto faceHasNeighbor = [](cBoundaryCellInfo::Node_t* node, int f) -> bool {
#if _MESH_DIMENSION_ == 1
    return node->GetNeibFace(f, 0, 0, PIC::Mesh::mesh) != nullptr;
#elif _MESH_DIMENSION_ == 2
    for (int ii = 0; ii < 2; ++ii)
      if (node->GetNeibFace(f, ii, 0, PIC::Mesh::mesh) != nullptr) return true;
    return false;
#else
    for (int jj = 0; jj < 2; ++jj)
      for (int ii = 0; ii < 2; ++ii)
        if (node->GetNeibFace(f, ii, jj, PIC::Mesh::mesh) != nullptr) return true;
    return false;
#endif
  };

  // First pass: reset all boundary flags to false
  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != nullptr; node = node->nextNodeThisThread)
  {
    if (node->block == nullptr) continue;

    for (int k = 0; k < NZ; ++k)
    for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i) {
      auto* c = node->block->GetCenterNode(i, j, k);
      if (c) c->SetBoundaryFlag(false);
    }
  }

  // Second pass: identify and mark boundary cells
  std::size_t added = 0;
  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != nullptr; node = node->nextNodeThisThread)
  {
    if (node->block == nullptr) continue; // skip non-leaves or unallocated

    // Per-block dedup mask: avoids double-adding edge/corner cells
    std::vector<unsigned char> seen(static_cast<size_t>(NX) * NY * NZ, 0);

    auto add_cell = [&](int i, int j, int k) {
      const int idx = i + STRIDE_J * (j + NY * k);
      if (seen[static_cast<size_t>(idx)]) return;
      seen[static_cast<size_t>(idx)] = 1;

      auto* c = node->block->GetCenterNode(i, j, k);
      if (!c) return; // guard against absent/ghost allocations
      
      // Set boundary flag to true
      c->SetBoundaryFlag(true);
      
      if (BoundaryCellVector != nullptr) {
        cBoundaryCellInfo rec;
        rec.node = node; rec.i = i; rec.j = j; rec.k = k; rec.cell = c;
        BoundaryCellVector->push_back(rec);
      }
      ++added;
    };

    // Face IDs: 0=-X, 1=+X, 2=-Y, 3=+Y, 4=-Z, 5=+Z
    // Append ONLY cells that lie on faces with NO neighbor (true domain boundary)

    // -X: i=0
    if (!faceHasNeighbor(node, 0)) {
      const int i = 0;
      for (int k = 0; k < NZ; ++k)
      for (int j = 0; j < NY; ++j)
        add_cell(i, j, k);
    }

    // +X: i=NX-1
    if (!faceHasNeighbor(node, 1)) {
      const int i = NX - 1;
      for (int k = 0; k < NZ; ++k)
      for (int j = 0; j < NY; ++j)
        add_cell(i, j, k);
    }

#if _MESH_DIMENSION_ >= 2
    // -Y: j=0
    if (!faceHasNeighbor(node, 2)) {
      const int j = 0;
      for (int k = 0; k < NZ; ++k)
      for (int i = 0; i < NX; ++i)
        add_cell(i, j, k);
    }

    // +Y: j=NY-1
    if (!faceHasNeighbor(node, 3)) {
      const int j = NY - 1;
      for (int k = 0; k < NZ; ++k)
      for (int i = 0; i < NX; ++i)
        add_cell(i, j, k);
    }
#endif

#if _MESH_DIMENSION_ == 3
    // -Z: k=0
    if (!faceHasNeighbor(node, 4)) {
      const int k = 0;
      for (int j = 0; j < NY; ++j)
      for (int i = 0; i < NX; ++i)
        add_cell(i, j, k);
    }

    // +Z: k=NZ-1
    if (!faceHasNeighbor(node, 5)) {
      const int k = NZ - 1;
      for (int j = 0; j < NY; ++j)
      for (int i = 0; i < NX; ++i)
        add_cell(i, j, k);
    }
#endif
  }

  std::printf("$PREFIX: InitBoundaryCellVector: collected %zu boundary cells on rank %d\n",
              added, PIC::ThisThread);
  return true;
}

} // namespace Mesh
} // namespace PIC
