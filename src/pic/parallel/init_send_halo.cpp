

#include "../pic.h"

namespace PIC {
namespace Parallel {

// Dispatcher
void InitSendHaloLayer(cHalo& SendHalo, bool communicate_entire_block) {
#if _MESH_DIMENSION_ == 1
  InitSendHaloLayer_1D(SendHalo, communicate_entire_block);
#elif _MESH_DIMENSION_ == 2
  InitSendHaloLayer_2D(SendHalo, communicate_entire_block);
#elif _MESH_DIMENSION_ == 3
  InitSendHaloLayer_3D(SendHalo, communicate_entire_block);
#else
  SendHalo.clear();
  if (PIC::Mesh::mesh) SendHalo.nMeshModificationCounter = PIC::Mesh::mesh->nMeshModificationCounter;
#endif
}

// -------------------------
// 1D implementation
// -------------------------
void InitSendHaloLayer_1D(cHalo& SendHalo, bool communicate_entire_block) {
  SendHalo.clear();

  exit(__LINE__,__FILE__,"errror: the functino was never tested or ran -- check if it works");

  auto* mesh = PIC::Mesh::mesh;
  if (mesh == nullptr) return;

  const int ThisThread    = PIC::ThisThread;
  const int nTotalThreads = mesh->nTotalThreads;
  if (nTotalThreads <= 1) { SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter; return; }

  const int NG = (int)_GHOST_CELLS_X_;
  if (!communicate_entire_block && NG <= 0) { SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter; return; }

  const int Nx = _BLOCK_CELLS_X_;

  // ---- helpers (all lambdas) ----
  auto apply_side_slab = [&](int N, int sign, int& iMin, int& iMax) {
    if (communicate_entire_block) { iMin = 0; iMax = N - 1; return; }
    const int ng = std::min(N, NG);
    iMin = 0; iMax = N - 1;
    if (sign < 0) { iMin = 0;      iMax = ng - 1; }
    if (sign > 0) { iMin = N - ng; iMax = N - 1;  }
  };

  auto derive_corner_from_cell = [&](const int cellMin[3], const int cellMax[3],
                                     int cornerMin[3], int cornerMax[3]) {
    // active: x; inactive: y,z -> 0..0
    cornerMin[0] = std::max(0, std::min(cellMin[0], Nx));
    cornerMax[0] = std::max(0, std::min(cellMax[0] + 1, Nx));
    cornerMin[1] = 0; cornerMax[1] = 0;
    cornerMin[2] = 0; cornerMax[2] = 0;
  };

  auto add_or_merge = [&](cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode, int ToThread,
                          const int cellMin[3], const int cellMax[3]) {
    if (!startNode) return;
    if (ToThread < 0 || ToThread >= nTotalThreads) return;
    if (ToThread == ThisThread) return;
    if (mesh->ParallelSendRecvMap[ThisThread][ToThread] == false) return;

    int cornerMin[3], cornerMax[3];
    derive_corner_from_cell(cellMin, cellMax, cornerMin, cornerMax);

    for (auto& e : SendHalo.list) {
      if (e.startNode == startNode && e.ToThread == ToThread) {
        e.i_cell_min[0]   = std::min(e.i_cell_min[0],   cellMin[0]);
        e.i_cell_max[0]   = std::max(e.i_cell_max[0],   cellMax[0]);
        e.i_corner_min[0] = std::min(e.i_corner_min[0], cornerMin[0]);
        e.i_corner_max[0] = std::max(e.i_corner_max[0], cornerMax[0]);
        return;
      }
    }

    cHaloEntry e;
    e.startNode = startNode;
    e.ToThread  = ToThread;

    e.i_cell_min[0] = cellMin[0]; e.i_cell_max[0] = cellMax[0];
    e.i_cell_min[1] = 0;          e.i_cell_max[1] = 0;
    e.i_cell_min[2] = 0;          e.i_cell_max[2] = 0;

    e.i_corner_min[0] = cornerMin[0]; e.i_corner_max[0] = cornerMax[0];
    e.i_corner_min[1] = 0;            e.i_corner_max[1] = 0;
    e.i_corner_min[2] = 0;            e.i_corner_max[2] = 0;

    SendHalo.list.push_back(e);
  };

  // Face indexing: 2 faces, 1 segment
  const int nFaces = 2;
  const int nSeg   = 1;
  const int nFaceIdx = nFaces * nSeg;

  for (auto* node = mesh->ParallelNodesDistributionList[ThisThread];
       node != nullptr; node = node->nextNodeThisThread) {

    if (!node->block) continue;
    if (!node->IsUsedInCalculationFlag) continue;
    if (node->Thread != ThisThread) continue;

    // Face neighbors: 0=x-, 1=x+
    for (int iface = 0; iface < nFaceIdx; ++iface) {
      auto* neib = node->neibNodeFace(iface, mesh);
      if (!neib) continue;

      const int ToThread = neib->Thread;
      if (ToThread == ThisThread) continue;

      const int signX = (iface == 0 ? -1 : +1);

      int cellMin[3] = {0,0,0}, cellMax[3] = {0,0,0};
      apply_side_slab(Nx, signX, cellMin[0], cellMax[0]);

      add_or_merge(node, ToThread, cellMin, cellMax);
    }

    // Corner neighbors: 2 corners (equivalent to faces in 1D)
    for (int ic = 0; ic < 2; ++ic) {
      auto* neib = node->neibNodeCorner(ic, mesh);
      if (!neib) continue;

      const int ToThread = neib->Thread;
      if (ToThread == ThisThread) continue;

      const int signX = (ic & 1) ? +1 : -1;

      int cellMin[3] = {0,0,0}, cellMax[3] = {0,0,0};
      apply_side_slab(Nx, signX, cellMin[0], cellMax[0]);

      add_or_merge(node, ToThread, cellMin, cellMax);
    }
  }

  // Defensive cleanup
  SendHalo.list.erase(
    std::remove_if(SendHalo.list.begin(), SendHalo.list.end(),
      [&](const cHaloEntry& e) {
        if (!e.startNode) return true;
        if (e.ToThread < 0 || e.ToThread >= nTotalThreads) return true;
        if (e.i_cell_max[0] < e.i_cell_min[0]) return true;
        return false;
      }),
    SendHalo.list.end()
  );

  // Snapshot mesh modification counter at exit (per request)
  SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter;
}

// -------------------------
// 2D implementation
// -------------------------
void InitSendHaloLayer_2D(cHalo& SendHalo, bool communicate_entire_block) {
  SendHalo.clear();

  exit(__LINE__,__FILE__,"errror: the functino was never tested or ran -- check if it works");

  auto* mesh = PIC::Mesh::mesh;
  if (mesh == nullptr) return;

  const int ThisThread    = PIC::ThisThread;
  const int nTotalThreads = mesh->nTotalThreads;
  if (nTotalThreads <= 1) { SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter; return; }

  const int NG = std::max((int)_GHOST_CELLS_X_, (int)_GHOST_CELLS_Y_);
  if (!communicate_entire_block && NG <= 0) { SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter; return; }

  const int Nx = _BLOCK_CELLS_X_;
  const int Ny = _BLOCK_CELLS_Y_;

  // ---- helpers (all lambdas) ----
  auto apply_side_slab = [&](int N, int sign, int& iMin, int& iMax) {
    if (communicate_entire_block) { iMin = 0; iMax = N - 1; return; }
    const int ng = std::min(N, NG);
    iMin = 0; iMax = N - 1;
    if (sign < 0) { iMin = 0;      iMax = ng - 1; }
    if (sign > 0) { iMin = N - ng; iMax = N - 1;  }
  };

  auto derive_corner_from_cell = [&](const int cellMin[3], const int cellMax[3],
                                     int cornerMin[3], int cornerMax[3]) {
    cornerMin[0] = std::max(0, std::min(cellMin[0], Nx));
    cornerMax[0] = std::max(0, std::min(cellMax[0] + 1, Nx));
    cornerMin[1] = std::max(0, std::min(cellMin[1], Ny));
    cornerMax[1] = std::max(0, std::min(cellMax[1] + 1, Ny));
    cornerMin[2] = 0; cornerMax[2] = 0;
  };

  auto add_or_merge = [&](cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode, int ToThread,
                          const int cellMin[3], const int cellMax[3]) {
    if (!startNode) return;
    if (ToThread < 0 || ToThread >= nTotalThreads) return;
    if (ToThread == ThisThread) return;
    if (mesh->ParallelSendRecvMap[ThisThread][ToThread] == false) return;

    int cornerMin[3], cornerMax[3];
    derive_corner_from_cell(cellMin, cellMax, cornerMin, cornerMax);

    for (auto& e : SendHalo.list) {
      if (e.startNode == startNode && e.ToThread == ToThread) {
        e.i_cell_min[0]   = std::min(e.i_cell_min[0],   cellMin[0]);
        e.i_cell_max[0]   = std::max(e.i_cell_max[0],   cellMax[0]);
        e.i_cell_min[1]   = std::min(e.i_cell_min[1],   cellMin[1]);
        e.i_cell_max[1]   = std::max(e.i_cell_max[1],   cellMax[1]);

        e.i_corner_min[0] = std::min(e.i_corner_min[0], cornerMin[0]);
        e.i_corner_max[0] = std::max(e.i_corner_max[0], cornerMax[0]);
        e.i_corner_min[1] = std::min(e.i_corner_min[1], cornerMin[1]);
        e.i_corner_max[1] = std::max(e.i_corner_max[1], cornerMax[1]);
        return;
      }
    }

    cHaloEntry e;
    e.startNode = startNode;
    e.ToThread  = ToThread;

    e.i_cell_min[0] = cellMin[0]; e.i_cell_max[0] = cellMax[0];
    e.i_cell_min[1] = cellMin[1]; e.i_cell_max[1] = cellMax[1];
    e.i_cell_min[2] = 0;          e.i_cell_max[2] = 0;

    e.i_corner_min[0] = cornerMin[0]; e.i_corner_max[0] = cornerMax[0];
    e.i_corner_min[1] = cornerMin[1]; e.i_corner_max[1] = cornerMax[1];
    e.i_corner_min[2] = 0;            e.i_corner_max[2] = 0;

    SendHalo.list.push_back(e);
  };

  // Face enumeration in D dims:
  //   nFaces = 2*D
  //   nSeg   = 2^(D-1)
  const int nFaces = 4;
  const int nSeg   = 2;            // 2^(2-1)
  const int nFaceIdx = nFaces * nSeg;

  auto faceDir_to_sign = [&](int faceDir, int sign[3]) {
    sign[0]=0; sign[1]=0; sign[2]=0;
    if      (faceDir == 0) sign[0] = -1; // x-
    else if (faceDir == 1) sign[0] = +1; // x+
    else if (faceDir == 2) sign[1] = -1; // y-
    else if (faceDir == 3) sign[1] = +1; // y+
  };

  auto build_cell_range_from_sign = [&](const int sign[3], int cellMin[3], int cellMax[3]) {
    apply_side_slab(Nx, sign[0], cellMin[0], cellMax[0]);
    apply_side_slab(Ny, sign[1], cellMin[1], cellMax[1]);
    cellMin[2] = 0; cellMax[2] = 0;
  };

  for (auto* node = mesh->ParallelNodesDistributionList[ThisThread];
       node != nullptr; node = node->nextNodeThisThread) {

    if (!node->block) continue;
    if (!node->IsUsedInCalculationFlag) continue;
    if (node->Thread != ThisThread) continue;

    // Face neighbors: faceDir = iface / nSeg
    for (int iface = 0; iface < nFaceIdx; ++iface) {
      auto* neib = node->neibNodeFace(iface, mesh);
      if (!neib) continue;

      const int ToThread = neib->Thread;
      if (ToThread == ThisThread) continue;

      const int faceDir = iface / nSeg;
      int sign[3];
      faceDir_to_sign(faceDir, sign);

      int cellMin[3], cellMax[3];
      build_cell_range_from_sign(sign, cellMin, cellMax);

      add_or_merge(node, ToThread, cellMin, cellMax);
    }

    // Corner neighbors: 2^2 = 4 corners
    for (int ic = 0; ic < 4; ++ic) {
      auto* neib = node->neibNodeCorner(ic, mesh);
      if (!neib) continue;

      const int ToThread = neib->Thread;
      if (ToThread == ThisThread) continue;

      int sign[3] = {0,0,0};
      sign[0] = (ic & 1) ? +1 : -1;
      sign[1] = (ic & 2) ? +1 : -1;

      int cellMin[3], cellMax[3];
      build_cell_range_from_sign(sign, cellMin, cellMax);

      add_or_merge(node, ToThread, cellMin, cellMax);
    }
  }

  SendHalo.list.erase(
    std::remove_if(SendHalo.list.begin(), SendHalo.list.end(),
      [&](const cHaloEntry& e) {
        if (!e.startNode) return true;
        if (e.ToThread < 0 || e.ToThread >= nTotalThreads) return true;
        if (e.i_cell_max[0] < e.i_cell_min[0]) return true;
        if (e.i_cell_max[1] < e.i_cell_min[1]) return true;
        return false;
      }),
    SendHalo.list.end()
  );

  SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter;
}

// -------------------------
// 3D implementation
// -------------------------
void InitSendHaloLayer_3D(cHalo& SendHalo, bool communicate_entire_block) {
  SendHalo.clear();

  auto* mesh = PIC::Mesh::mesh;
  if (mesh == nullptr) return;

  const int ThisThread    = PIC::ThisThread;
  const int nTotalThreads = mesh->nTotalThreads;
  if (nTotalThreads <= 1) { SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter; return; }

  const int NG = std::max((int)_GHOST_CELLS_X_,
                 std::max((int)_GHOST_CELLS_Y_, (int)_GHOST_CELLS_Z_));
  if (!communicate_entire_block && NG <= 0) { SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter; return; }

  const int Nx = _BLOCK_CELLS_X_;
  const int Ny = _BLOCK_CELLS_Y_;
  const int Nz = _BLOCK_CELLS_Z_;

  // ---- helpers (all lambdas) ----
  auto apply_side_slab = [&](int N, int sign, int& iMin, int& iMax) {
    if (communicate_entire_block) { iMin = 0; iMax = N - 1; return; }
    const int ng = std::min(N, NG);
    iMin = 0; iMax = N - 1;
    if (sign < 0) { iMin = 0;      iMax = ng - 1; }
    if (sign > 0) { iMin = N - ng; iMax = N - 1;  }
  };

  auto derive_corner_from_cell = [&](const int cellMin[3], const int cellMax[3],
                                     int cornerMin[3], int cornerMax[3]) {
    cornerMin[0] = std::max(0, std::min(cellMin[0], Nx));
    cornerMax[0] = std::max(0, std::min(cellMax[0] + 1, Nx));
    cornerMin[1] = std::max(0, std::min(cellMin[1], Ny));
    cornerMax[1] = std::max(0, std::min(cellMax[1] + 1, Ny));
    cornerMin[2] = std::max(0, std::min(cellMin[2], Nz));
    cornerMax[2] = std::max(0, std::min(cellMax[2] + 1, Nz));
  };

  auto add_or_merge = [&](cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode, int ToThread,
                          const int cellMin[3], const int cellMax[3]) {
    if (!startNode) return;
    if (ToThread < 0 || ToThread >= nTotalThreads) return;
    if (ToThread == ThisThread) return;
    if (mesh->ParallelSendRecvMap[ThisThread][ToThread] == false) return;

    int cornerMin[3], cornerMax[3];
    derive_corner_from_cell(cellMin, cellMax, cornerMin, cornerMax);

    for (auto& e : SendHalo.list) {
      if (e.startNode == startNode && e.ToThread == ToThread) {
        for (int d=0; d<3; ++d) {
          e.i_cell_min[d]   = std::min(e.i_cell_min[d],   cellMin[d]);
          e.i_cell_max[d]   = std::max(e.i_cell_max[d],   cellMax[d]);
          e.i_corner_min[d] = std::min(e.i_corner_min[d], cornerMin[d]);
          e.i_corner_max[d] = std::max(e.i_corner_max[d], cornerMax[d]);
        }
        return;
      }
    }

    cHaloEntry e;
    e.startNode = startNode;
    e.ToThread  = ToThread;
    for (int d=0; d<3; ++d) {
      e.i_cell_min[d]   = cellMin[d];
      e.i_cell_max[d]   = cellMax[d];
      e.i_corner_min[d] = cornerMin[d];
      e.i_corner_max[d] = cornerMax[d];
    }
    SendHalo.list.push_back(e);
  };

  // Face enumeration in D dims: nFaces=2*D, nSeg=2^(D-1)
  const int nFaces = 6;
  const int nSeg   = 4;     // 2^(3-1)
  const int nFaceIdx = nFaces * nSeg;

  auto faceDir_to_sign = [&](int faceDir, int sign[3]) {
    sign[0]=0; sign[1]=0; sign[2]=0;
    if      (faceDir == 0) sign[0] = -1; // x-
    else if (faceDir == 1) sign[0] = +1; // x+
    else if (faceDir == 2) sign[1] = -1; // y-
    else if (faceDir == 3) sign[1] = +1; // y+
    else if (faceDir == 4) sign[2] = -1; // z-
    else if (faceDir == 5) sign[2] = +1; // z+
  };

  auto build_cell_range_from_sign = [&](const int sign[3], int cellMin[3], int cellMax[3]) {
    apply_side_slab(Nx, sign[0], cellMin[0], cellMax[0]);
    apply_side_slab(Ny, sign[1], cellMin[1], cellMax[1]);
    apply_side_slab(Nz, sign[2], cellMin[2], cellMax[2]);
  };

  // EdgeIncrement mapping used to convert edge-id to sign vector (two offset dims, one along-edge dim)
  const int EdgeIncrement[12][3] = {
    { 0,-1,-1}, { 0, 1,-1}, { 0, 1, 1}, { 0,-1, 1},
    {-1, 0,-1}, { 1, 0,-1}, { 1, 0, 1}, {-1, 0, 1},
    {-1,-1, 0}, { 1,-1, 0}, { 1, 1, 0}, {-1, 1, 0}
  };

  for (auto* node = mesh->ParallelNodesDistributionList[ThisThread];
       node != nullptr; node = node->nextNodeThisThread) {

    if (!node->block) continue;
    if (!node->IsUsedInCalculationFlag) continue;
    if (node->Thread != ThisThread) continue;

    // Face neighbors: faceDir = iface / nSeg
    for (int iface = 0; iface < nFaceIdx; ++iface) {
      auto* neib = node->neibNodeFace(iface, mesh);
      if (!neib) continue;

      const int ToThread = neib->Thread;
      if (ToThread == ThisThread) continue;

      const int faceDir = iface / nSeg;
      int sign[3];
      faceDir_to_sign(faceDir, sign);

      int cellMin[3], cellMax[3];
      build_cell_range_from_sign(sign, cellMin, cellMax);

      add_or_merge(node, ToThread, cellMin, cellMax);
    }

    // Edge neighbors: 12 edges * 2 segments
    for (int iedgeSeg = 0; iedgeSeg < 12 * 2; ++iedgeSeg) {
      auto* neib = node->neibNodeEdge(iedgeSeg, mesh);
      if (!neib) continue;

      const int ToThread = neib->Thread;
      if (ToThread == ThisThread) continue;

      const int iedge = iedgeSeg / 2;

      int sign[3] = {0,0,0};
      for (int d=0; d<3; ++d) {
        if (EdgeIncrement[iedge][d] < 0) sign[d] = -1;
        else if (EdgeIncrement[iedge][d] > 0) sign[d] = +1;
        else sign[d] = 0;
      }

      int cellMin[3], cellMax[3];
      build_cell_range_from_sign(sign, cellMin, cellMax);

      add_or_merge(node, ToThread, cellMin, cellMax);
    }

    // Corner neighbors: 8 corners
    for (int ic = 0; ic < 8; ++ic) {
      auto* neib = node->neibNodeCorner(ic, mesh);
      if (!neib) continue;

      const int ToThread = neib->Thread;
      if (ToThread == ThisThread) continue;

      int sign[3];
      sign[0] = (ic & 1) ? +1 : -1;
      sign[1] = (ic & 2) ? +1 : -1;
      sign[2] = (ic & 4) ? +1 : -1;

      int cellMin[3], cellMax[3];
      build_cell_range_from_sign(sign, cellMin, cellMax);

      add_or_merge(node, ToThread, cellMin, cellMax);
    }
  }

  // Defensive cleanup
  SendHalo.list.erase(
    std::remove_if(SendHalo.list.begin(), SendHalo.list.end(),
      [&](const cHaloEntry& e) {
        if (!e.startNode) return true;
        if (e.ToThread < 0 || e.ToThread >= nTotalThreads) return true;
        for (int d=0; d<3; ++d) {
          if (e.i_cell_max[d]   < e.i_cell_min[d])   return true;
          if (e.i_corner_max[d] < e.i_corner_min[d]) return true;
        }
        return false;
      }),
    SendHalo.list.end()
  );

  SendHalo.nMeshModificationCounter = mesh->nMeshModificationCounter;
}

} // namespace Parallel
} // namespace PIC

