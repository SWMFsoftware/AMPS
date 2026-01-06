#include <vector>
#include <cstring>
#include <algorithm>

#include "../pic.h"

namespace PIC {
namespace Parallel {

// ============================================================================
// NEW overload: uses PIC::Parallel::cHalo SendHalo (SEND-side descriptors)
//
// Key behavior change vs old implementation:
//   - Old: SEND selection relied on mesh->DomainBoundaryLayerNodesList[To] and sent
//          full-node arrays for each block.
//   - New: SEND selection uses SendHalo.list entries (owned block, ToThread, ranges),
//          and sends ONLY the sub-range defined by ghost thickness (or full block if
//          SendHalo was built with communicate_entire_block=true).
//
// Mesh-change robustness:
//   - If SendHalo.nMeshModificationCounter != mesh->nMeshModificationCounter,
//     the halo list is rebuilt via InitSendHaloLayer(SendHalo).
//
// Wire format per neighbor rank 'To':
//   [int nBlocks] +
//   repeat nBlocks times:
//     [cAMRnodeID blockId] +
//     [int iMin] [int iMax] [int jMin] [int jMax] [int kMin] [int kMax]  // range for chosen node type
//     [payload for all nodes in that range in deterministic i/j/k order]
//       where each node payload is manager.nDoubles doubles
//
// Receiver applies ONLY into ghost blocks owned by 'From' (node->Thread == From).
// ============================================================================
void SyncNodeHalo_DomainBoundaryLayer(const cNodeHaloSyncManager& manager,
                                      cHalo& SendHalo,
                                      int tagBase /*=31000*/) {
#if _AMR_PARALLEL_DATA_EXCHANGE_MODE_ == _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_

  // ===========================================================================
  // 0) Basic sanity / early exits
  // ===========================================================================
  auto* mesh = PIC::Mesh::mesh;
  if (mesh == nullptr) return;

  const int ThisThread    = PIC::ThisThread;
  const int nTotalThreads = mesh->nTotalThreads;

  if (nTotalThreads <= 1) return;
  if (manager.nDoubles <= 0) return;

  const bool isCornerNodes = (manager.NodeType == eNodeType::Corner);

  // ===========================================================================
  // 1) Ensure SendHalo is initialized for the current mesh topology
  //   Rebuild if:
  //     (A) mesh topology changed
  //     (B) communicate_entire_block requested by manager differs from halo mode
  // --------------------------------------------------------------------------
  const bool meshChanged =
      (SendHalo.nMeshModificationCounter != mesh->nMeshModificationCounter);

  const bool modeChanged =
      (SendHalo.communicate_entire_block != manager.communicate_entire_block);

  if (meshChanged || modeChanged) {
    // Build halo in the exact mode requested by manager
    InitSendHaloLayer(SendHalo, manager.communicate_entire_block);
  }

  // Two tags per exchange: size + payload
  const int tagSize = tagBase + 1;
  const int tagData = tagBase + 2;

  const int bytesPerNodePayload = manager.nDoubles * (int)sizeof(double);

  // Each block record now includes: cAMRnodeID + 6 ints (range) + payload bytes
  const int bytesPerRangeHeader = 6 * (int)sizeof(int);

  // A small zero payload to keep stream aligned when a node pointer is null (rare)
  std::vector<char> zeroPayload((size_t)bytesPerNodePayload, (char)0);

  // ===========================================================================
  // 2) Helpers (ALL lambdas)
  // ===========================================================================

  auto pack_int = [](std::vector<char>& buf, int v) {
    const char* p = reinterpret_cast<const char*>(&v);
    buf.insert(buf.end(), p, p + (int)sizeof(int));
  };

  auto unpack_int = [](const char*& p) -> int {
    int v;
    std::memcpy(&v, p, sizeof(int));
    p += sizeof(int);
    return v;
  };

  auto pack_node_payload = [&](std::vector<char>& buf, char* associatedBufferPtr) {
    const char* src =
        reinterpret_cast<const char*>(associatedBufferPtr + manager.RelativeOffsetBytesFromAssociated);
    buf.insert(buf.end(), src, src + bytesPerNodePayload);
  };

  auto pack_node_zeros = [&](std::vector<char>& buf) {
    buf.insert(buf.end(), zeroPayload.begin(), zeroPayload.end());
  };

  auto apply_node_payload = [&](char* associatedBufferPtr, const char* payloadBytes) {
    double* dst =
        reinterpret_cast<double*>(associatedBufferPtr + manager.RelativeOffsetBytesFromAssociated);
    const double* src = reinterpret_cast<const double*>(payloadBytes);

    if (manager.Op == eHaloOp::Replace) {
      std::memcpy(dst, src, (size_t)bytesPerNodePayload);
    }
    else {
      for (int q = 0; q < manager.nDoubles; ++q) dst[q] += src[q];
    }
  };

  auto pack_amr_id = [&](std::vector<char>& buf, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    cAMRnodeID nodeid;
    mesh->GetAMRnodeID(nodeid, node);
    buf.insert(buf.end(),
               reinterpret_cast<const char*>(&nodeid),
               reinterpret_cast<const char*>(&nodeid) + (int)sizeof(cAMRnodeID));
  };

  auto unpack_amr_id = [&](const char*& p, const char* pend, cAMRnodeID& out) -> bool {
    if (p + (int)sizeof(cAMRnodeID) > pend) return false;
    std::memcpy(&out, p, sizeof(cAMRnodeID));
    p += sizeof(cAMRnodeID);
    return true;
  };

  // Compute node-count in the sent range (inclusive bounds)
  auto nNodesInRange = [&](int iMin, int iMax, int jMin, int jMax, int kMin, int kMax) -> size_t {
    if (iMax < iMin || jMax < jMin || kMax < kMin) return 0;
    return (size_t)(iMax - iMin + 1) *
           (size_t)(jMax - jMin + 1) *
           (size_t)(kMax - kMin + 1);
  };

  // Pack one halo entry record: [id] + [range(6 ints)] + [payloads for nodes in range]
  auto pack_halo_entry = [&](std::vector<char>& buf, const cHaloEntry& e) {
    auto* node = e.startNode;
    if (node == nullptr || node->block == nullptr) return;

    // 1) block id
    pack_amr_id(buf, node);

    // 2) range header (for selected node type)
    int iMin, iMax, jMin, jMax, kMin, kMax;
    if (isCornerNodes) {
      iMin = e.i_corner_min[0]; iMax = e.i_corner_max[0];
      jMin = e.i_corner_min[1]; jMax = e.i_corner_max[1];
      kMin = e.i_corner_min[2]; kMax = e.i_corner_max[2];
    }
    else {
      iMin = e.i_cell_min[0]; iMax = e.i_cell_max[0];
      jMin = e.i_cell_min[1]; jMax = e.i_cell_max[1];
      kMin = e.i_cell_min[2]; kMax = e.i_cell_max[2];
    }

    pack_int(buf, iMin); pack_int(buf, iMax);
    pack_int(buf, jMin); pack_int(buf, jMax);
    pack_int(buf, kMin); pack_int(buf, kMax);

    // 3) payload
    if (isCornerNodes) {
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i) {

            auto* cn = node->block->GetCornerNode(mesh->getCornerNodeLocalNumber(i, j, k));
            if (cn == nullptr) { pack_node_zeros(buf); continue; }
            pack_node_payload(buf, cn->GetAssociatedDataBufferPointer());
          }
    }
    else {
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i) {

            auto* cn = node->block->GetCenterNode(mesh->getCenterNodeLocalNumber(i, j, k));
            if (cn == nullptr) { pack_node_zeros(buf); continue; }
            pack_node_payload(buf, cn->GetAssociatedDataBufferPointer());
          }
    }
  };

  // Unpack/apply one received buffer from rank 'From'
  auto unpack_buffer_from = [&](const std::vector<char>& rbuf, int From) {
    const char* p    = rbuf.data();
    const char* pend = rbuf.data() + rbuf.size();

    if (p + (int)sizeof(int) > pend) return;
    const int nBlocks = unpack_int(p);

    for (int b = 0; b < nBlocks; ++b) {
      // ---- read id ----
      cAMRnodeID nodeid;
      if (!unpack_amr_id(p, pend, nodeid)) return;

      // ---- read range ----
      if (p + bytesPerRangeHeader > pend) return;

      const int iMin = unpack_int(p);
      const int iMax = unpack_int(p);
      const int jMin = unpack_int(p);
      const int jMax = unpack_int(p);
      const int kMin = unpack_int(p);
      const int kMax = unpack_int(p);

      const size_t nNodes = nNodesInRange(iMin,iMax,jMin,jMax,kMin,kMax);
      const size_t payloadBytes = nNodes * (size_t)bytesPerNodePayload;

      // ---- locate local block ----
      auto* node = mesh->findAMRnodeWithID(nodeid);

      // Apply ONLY to ghost blocks owned by sender 'From'
      if (node == nullptr || node->block == nullptr || node->Thread != From) {
        // Skip payload bytes
        if (p + (ptrdiff_t)payloadBytes > pend) return;
        p += payloadBytes;
        continue;
      }

      // ---- apply node payloads in same i/j/k order ----
      if (isCornerNodes) {
        for (int k = kMin; k <= kMax; ++k)
          for (int j = jMin; j <= jMax; ++j)
            for (int i = iMin; i <= iMax; ++i) {

              if (p + bytesPerNodePayload > pend) return;

              auto* cn = node->block->GetCornerNode(mesh->getCornerNodeLocalNumber(i, j, k));
              if (cn != nullptr) apply_node_payload(cn->GetAssociatedDataBufferPointer(), p);

              p += bytesPerNodePayload;
            }
      }
      else {
        for (int k = kMin; k <= kMax; ++k)
          for (int j = jMin; j <= jMax; ++j)
            for (int i = iMin; i <= iMax; ++i) {

              if (p + bytesPerNodePayload > pend) return;

              auto* cn = node->block->GetCenterNode(mesh->getCenterNodeLocalNumber(i, j, k));
              if (cn != nullptr) apply_node_payload(cn->GetAssociatedDataBufferPointer(), p);

              p += bytesPerNodePayload;
            }
      }
    }
  };

  // ===========================================================================
  // 3) Build send buffers per neighbor based on SendHalo.list and exchange sizes
  // ===========================================================================
  std::vector<int> recvBytes(nTotalThreads, 0);
  std::vector<std::vector<char>> sendBuf(nTotalThreads);
  std::vector<std::vector<char>> recvBuf(nTotalThreads);

  for (int To = 0; To < nTotalThreads; ++To) {
    if (To == ThisThread) continue;
    if (mesh->ParallelSendRecvMap[ThisThread][To] == false) continue;

    // Count how many halo entries go to 'To'
    int nBlocksToSend = 0;
    for (const auto& e : SendHalo.list) {
      if (e.ToThread != To) continue;
      if (e.startNode == nullptr || e.startNode->block == nullptr) continue;
      if (e.startNode->IsUsedInCalculationFlag == false) continue;

      // SEND list should already be owned by ThisThread, but keep it defensive:
      if (e.startNode->Thread != ThisThread) continue;

      ++nBlocksToSend;
    }

    // Pack header: number of block records
    pack_int(sendBuf[To], nBlocksToSend);

    // Pack each block record
    for (const auto& e : SendHalo.list) {
      if (e.ToThread != To) continue;
      if (e.startNode == nullptr || e.startNode->block == nullptr) continue;
      if (e.startNode->IsUsedInCalculationFlag == false) continue;
      if (e.startNode->Thread != ThisThread) continue;

      pack_halo_entry(sendBuf[To], e);
    }

    // Exchange sizes (blocking Sendrecv)
    int sbytes = (int)sendBuf[To].size();
    int rbytes = 0;

    MPI_Sendrecv(&sbytes, 1, MPI_INT, To, tagSize,
                 &rbytes, 1, MPI_INT, To, tagSize,
                 MPI_GLOBAL_COMMUNICATOR, MPI_STATUS_IGNORE);

    recvBytes[To] = rbytes;
    recvBuf[To].resize((size_t)rbytes);
  }

  // ===========================================================================
  // 4) Exchange payload buffers (nonblocking) and wait
  // ===========================================================================
  std::vector<MPI_Request> req;
  req.reserve(2 * (size_t)nTotalThreads);

  for (int To = 0; To < nTotalThreads; ++To) {
    if (To == ThisThread) continue;
    if (mesh->ParallelSendRecvMap[ThisThread][To] == false) continue;

    if (recvBytes[To] > 0) {
      MPI_Request r;
      MPI_Irecv(recvBuf[To].data(), recvBytes[To], MPI_BYTE, To, tagData,
                MPI_GLOBAL_COMMUNICATOR, &r);
      req.push_back(r);
    }

    if (!sendBuf[To].empty()) {
      MPI_Request s;
      MPI_Isend(sendBuf[To].data(), (int)sendBuf[To].size(), MPI_BYTE, To, tagData,
                MPI_GLOBAL_COMMUNICATOR, &s);
      req.push_back(s);
    }
  }

  if (!req.empty()) MPI_Waitall((int)req.size(), req.data(), MPI_STATUSES_IGNORE);

  // ===========================================================================
  // 5) Unpack/apply received buffers
  // ===========================================================================
  for (int From = 0; From < nTotalThreads; ++From) {
    if (From == ThisThread) continue;
    if (mesh->ParallelSendRecvMap[ThisThread][From] == false) continue;
    if (recvBuf[From].empty()) continue;

    unpack_buffer_from(recvBuf[From], From);
  }

#endif
}

// ============================================================================
// Backward-compatible wrapper: keeps old signature intact.
// Internally uses a static SendHalo that is rebuilt when mesh changes.
// ============================================================================
void SyncNodeHalo_DomainBoundaryLayer(const cNodeHaloSyncManager& manager,
                                      int tagBase /*=31000*/) {
  static cHalo _SendHalo;
  SyncNodeHalo_DomainBoundaryLayer(manager, _SendHalo, tagBase);
}

} // namespace Parallel
} // namespace PIC

