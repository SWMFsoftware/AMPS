// Read-only access to two private members of the UNCHANGED LiouvilleCOPS classes (ROADMAP C11).
//
// The GPU implementation must use exactly the objects main() built:
//   * the BVH of cBVH (Nodes, FaceOrder): the device traversal then visits the same
//     nodes and leaves in the same order as cBVH::FirstHit/AnyHit;
//   * the NG shield table of cGaugeModel (built by a Monte Carlo with mt19937_64 in
//     cGaugeModel::BuildShieldTable).
// Both are private.  ../LiouvilleCOPS must not change, so this one translation unit
// re-declares access for reading only.  Under the Itanium C++ ABI (GCC) access
// specifiers do not change the layout of these classes.  The values extracted here
// are cross-checked at run time against the classes' PUBLIC methods
// (cBVH::FirstHit, cGaugeModel::NGWeight) by gpu_liouville.cu before any use.
#define private public
#include "BVH.h"
#include "Gauges.h"
#undef private

#include "access_shim.h"

void ShimBVH(const cBVH& bvh, std::vector<ShimNode>& nodes, std::vector<int>& faceOrder) {
  nodes.resize(bvh.Nodes.size());
  for (size_t i = 0; i < bvh.Nodes.size(); i++) {
    const auto& n = bvh.Nodes[i];
    ShimNode& s = nodes[i];
    s.bmin[0] = n.bmin.x; s.bmin[1] = n.bmin.y; s.bmin[2] = n.bmin.z;
    s.bmax[0] = n.bmax.x; s.bmax[1] = n.bmax.y; s.bmax[2] = n.bmax.z;
    s.left = n.left; s.right = n.right; s.first = n.first; s.count = n.count;
  }
  faceOrder = bvh.FaceOrder;
}

const std::vector<double>& ShimShieldTable(const cGaugeModel& gauges) {
  return gauges.ShieldTable;
}
