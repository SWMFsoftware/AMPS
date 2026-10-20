// Read-only views of private data of the unchanged LiouvilleCOPS classes (see access_shim.cpp).
#ifndef COPS_CUDA_ACCESS_SHIM_H
#define COPS_CUDA_ACCESS_SHIM_H

#include <vector>

class cBVH;
class cGaugeModel;

struct ShimNode {
  double bmin[3], bmax[3];
  int left, right, first, count;
};

// Copy of cBVH::Nodes (flattened) and cBVH::FaceOrder.
void ShimBVH(const cBVH& bvh, std::vector<ShimNode>& nodes, std::vector<int>& faceOrder);

// cGaugeModel::ShieldTable (empty unless NGAcceptance=srccg_shield).
const std::vector<double>& ShimShieldTable(const cGaugeModel& gauges);

#endif
