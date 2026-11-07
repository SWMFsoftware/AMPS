#include "mesh_test_fixture.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace AMPS {
namespace CoreMeshTests {
namespace {

typedef cAmpsMesh<PIC::Mesh::cDataCornerNode, PIC::Mesh::cDataCenterNode,
                  PIC::Mesh::cDataBlockAMR> Mesh;

int CoarsenDataOffset = -1;

double ReadBlockDatum(const Node* node) {
  AMPS_EXPECT_TRUE(node != NULL);
  AMPS_EXPECT_TRUE(node->block != NULL);
  AMPS_EXPECT_TRUE(CoarsenDataOffset >= 0);
  double value = 0.0;
  std::memcpy(&value,
              node->block->GetAssociatedDataBufferPointer()+CoarsenDataOffset,
              sizeof(value));
  return value;
}

void WriteBlockDatum(Node* node,double value) {
  AMPS_EXPECT_TRUE(node != NULL);
  AMPS_EXPECT_TRUE(node->block != NULL);
  AMPS_EXPECT_TRUE(CoarsenDataOffset >= 0);
  std::memcpy(node->block->GetAssociatedDataBufferPointer()+CoarsenDataOffset,
              &value,sizeof(value));
}

double CoarsenResolution(double*) { return 1.0; }

void InitTopologyOnly(Mesh& mesh) {
  double xmin[3] = {0.0, 0.0, 0.0};
  double xmax[3] = {1.0, 1.0, 1.0};
  mesh.AllowBlockAllocation = false;
  mesh.init(xmin, xmax, CoarsenResolution);
}

void InitAllocated(Mesh& mesh) {
  double xmin[3] = {0.0, 0.0, 0.0};
  double xmax[3] = {1.0, 1.0, 1.0};
  mesh.AllowBlockAllocation = true;
  mesh.init(xmin, xmax, CoarsenResolution);
}

void CollectLeaves(Node* node, std::vector<Node*>& leaves) {
  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    leaves.push_back(node);
    return;
  }
  for (int child = 0; child < (1 << _MESH_DIMENSION_); ++child) {
    AMPS_EXPECT_TRUE(node->downNode[child] != NULL);
    CollectLeaves(node->downNode[child], leaves);
  }
}

double NodeVolume(const Node* node) {
  double result = 1.0;
  for (int d = 0; d < _MESH_DIMENSION_; ++d)
    result *= node->xmax[d] - node->xmin[d];
  return result;
}

double LeafVolume(Mesh& mesh) {
  std::vector<Node*> leaves;
  CollectLeaves(mesh.rootTree, leaves);
  double result = 0.0;
  for (std::size_t n = 0; n < leaves.size(); ++n)
    result += NodeVolume(leaves[n]);
  return result;
}

std::uint64_t RootReceipt(const Node* root) {
  // This checksum is deliberately independent of allocation addresses and of
  // the production mesh-signature routine.  A complete split/coarsen cycle
  // must restore these exact dyadic identity fields.
  std::uint64_t result = UINT64_C(1469598103934665603);
  const std::uint64_t prime = UINT64_C(1099511628211);
  result = (result ^ static_cast<unsigned char>(root->RefinmentLevel)) * prime;
  result = (result ^ static_cast<std::uint64_t>(root->NodeGeometricSizeIndex)) *
      prime;
  for (int d = 0; d < _MESH_DIMENSION_; ++d)
    result = (result ^ static_cast<std::uint64_t>(root->xMinGlobalIndex[d])) *
        prime;
  return result;
}

struct RestrictionReceipt {
  int prolong_calls;
  int calls;
  int observed_children;
  double child_volume;
  double parent_value;
};

void ProlongLoadMeasure(Node* parent, Node* const* children,
                        int nChildren, void* raw) {
  RestrictionReceipt* receipt =
      static_cast<RestrictionReceipt*>(raw);
  ++receipt->prolong_calls;
  receipt->observed_children = nChildren;
  receipt->parent_value = ReadBlockDatum(parent);

  // Constant preservation is the zero-gradient limit of every consistent
  // prolongation scheme.  The application callback, not the generic mesh,
  // owns the decision that this particular datum is intensive.
  for (int child = 0; child < nChildren; ++child) {
    AMPS_EXPECT_TRUE(children[child] != NULL);
    AMPS_EXPECT_TRUE(children[child]->upNode == parent);
    WriteBlockDatum(children[child],receipt->parent_value);
  }
}

void RestrictLoadMeasure(Node* parent, Node* const* children,
                         int nChildren, void* raw) {
  RestrictionReceipt* receipt =
      static_cast<RestrictionReceipt*>(raw);
  ++receipt->calls;
  receipt->observed_children = nChildren;

  /* The reserved block slot is real opaque AMPS associated storage.  Treat it
   * here as an intensive scalar q and apply the independent finite-volume
   * restriction q_p=sum(q_c V_c)/sum(V_c).  The unit cube makes V dimensioned
   * in model-length^D; no AMPS interpolation routine participates. */
  double weighted = 0.0;
  double volume = 0.0;
  for (int child = 0; child < nChildren; ++child) {
    AMPS_EXPECT_TRUE(children[child] != NULL);
    AMPS_EXPECT_TRUE(children[child]->upNode == parent);
    const double child_volume = NodeVolume(children[child]);
    weighted += ReadBlockDatum(children[child]) * child_volume;
    volume += child_volume;
  }
  AMPS_EXPECT_TRUE(volume > 0.0);
  WriteBlockDatum(parent,weighted / volume);
  receipt->child_volume = volume;
}

}  // namespace

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-COARSEN"),
               "balanced leaf siblings collapse without a partition hole") {
  const int branching = 1 << _MESH_DIMENSION_;
  const double tolerance =
      32.0 * std::numeric_limits<double>::epsilon();

  // First prove that a collapse which would create a two-level face jump is
  // rejected before the topology or generation counter changes.
  Mesh rejected;
  InitTopologyOnly(rejected);
  AMPS_EXPECT_TRUE(rejected.splitTreeNode(rejected.rootTree));
  AMPS_EXPECT_TRUE(rejected.splitTreeNode(rejected.rootTree->downNode[0]));
  AMPS_EXPECT_TRUE(rejected.splitTreeNode(rejected.rootTree->downNode[1]));
  AMPS_EXPECT_TRUE(
      rejected.splitTreeNode(rejected.rootTree->downNode[1]->downNode[0]));
  const long rejected_generation = rejected.GetMeshID();
  Node* rejected_parent = rejected.rootTree->downNode[0];
  AMPS_EXPECT_FALSE(rejected.coarsenTreeNode(rejected_parent));
  AMPS_EXPECT_EQ(rejected.GetMeshID(), rejected_generation);
  AMPS_EXPECT_FALSE(rejected_parent->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);

  // A balanced non-root family collapses to its level-one parent, restores
  // the uniform root partition, and retains reciprocal face search.
  Mesh mesh;
  InitTopologyOnly(mesh);
  AMPS_EXPECT_TRUE(mesh.splitTreeNode(mesh.rootTree));
  AMPS_EXPECT_TRUE(mesh.splitTreeNode(mesh.rootTree->downNode[0]));
  Node* lower_parent = mesh.rootTree->downNode[0];
  AMPS_EXPECT_TRUE(mesh.coarsenTreeNode(lower_parent));
  AMPS_EXPECT_TRUE(lower_parent->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);

  std::vector<Node*> leaves;
  CollectLeaves(mesh.rootTree, leaves);
  AMPS_EXPECT_EQ(static_cast<int>(leaves.size()), branching);
  AMPS_EXPECT_NEAR(LeafVolume(mesh), 1.0, tolerance);
  for (int d = 0; d < _MESH_DIMENSION_; ++d) {
    Node* expected = mesh.rootTree->downNode[1 << d];
    AMPS_EXPECT_TRUE(
        lower_parent->GetNeibFace(2 * d + 1, 0, 0, &mesh) == expected);
    AMPS_EXPECT_TRUE(expected->GetNeibFace(2 * d, 0, 0, &mesh) ==
                     lower_parent);
  }

  // The legacy public entry takes any child of the sibling family.  It now
  // routes through the same validated coarsener instead of the historical
  // unconditional fatal path.
  AMPS_EXPECT_TRUE(mesh.deleteTreeNode(mesh.rootTree->downNode[0]));
  AMPS_EXPECT_TRUE(mesh.rootTree->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
  AMPS_EXPECT_NEAR(NodeVolume(mesh.rootTree), 1.0, tolerance);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-DATA"),
               "block data prolongation and restriction use explicit laws") {
  Mesh mesh;
  InitAllocated(mesh);
  WriteBlockDatum(mesh.rootTree,7.25);
  RestrictionReceipt receipt = {0, 0, 0, 0.0, 0.0};
  AMPS_EXPECT_TRUE(mesh.splitTreeNode(
      mesh.rootTree, ProlongLoadMeasure, &receipt));

  const int children = 1 << _MESH_DIMENSION_;
  AMPS_EXPECT_EQ(receipt.prolong_calls, 1);
  AMPS_EXPECT_EQ(receipt.observed_children, children);
  AMPS_EXPECT_NEAR(receipt.parent_value, 7.25, 0.0);
  for (int child = 0; child < children; ++child)
    AMPS_EXPECT_NEAR(
        ReadBlockDatum(mesh.rootTree->downNode[child]), 7.25, 0.0);

  // Allocated payload may never take the topology-only shortcut.  With no
  // restriction law the operation is rejected before changing generation,
  // connectivity, blocks, or their data bytes.
  const long before_missing_policy = mesh.GetMeshID();
  AMPS_EXPECT_FALSE(mesh.coarsenTreeNode(mesh.rootTree));
  AMPS_EXPECT_EQ(mesh.GetMeshID(), before_missing_policy);
  AMPS_EXPECT_FALSE(mesh.rootTree->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
  for (int child = 0; child < children; ++child) {
    AMPS_EXPECT_TRUE(mesh.rootTree->downNode[child] != NULL);
    AMPS_EXPECT_TRUE(mesh.rootTree->downNode[child]->block != NULL);
    AMPS_EXPECT_NEAR(
        ReadBlockDatum(mesh.rootTree->downNode[child]), 7.25, 0.0);
  }

  double independent_weighted_sum = 0.0;
  double independent_volume = 0.0;
  for (int child = 0; child < children; ++child) {
    Node* node = mesh.rootTree->downNode[child];
    // A non-symmetric sequence prevents an erroneous first-child copy from
    // matching the finite-volume oracle by accident.
    WriteBlockDatum(node,2.0 + 0.75 * child);
    independent_weighted_sum +=
        ReadBlockDatum(node) * NodeVolume(node);
    independent_volume += NodeVolume(node);
  }
  const double expected = independent_weighted_sum / independent_volume;
  AMPS_EXPECT_TRUE(mesh.coarsenTreeNode(
      mesh.rootTree, RestrictLoadMeasure, &receipt));
  AMPS_EXPECT_EQ(receipt.prolong_calls, 1);
  AMPS_EXPECT_EQ(receipt.calls, 1);
  AMPS_EXPECT_EQ(receipt.observed_children, children);
  AMPS_EXPECT_NEAR(receipt.child_volume, 1.0,
                   32.0 * std::numeric_limits<double>::epsilon());
  AMPS_EXPECT_NEAR(ReadBlockDatum(mesh.rootTree), expected,
                   32.0 * std::numeric_limits<double>::epsilon());
  AMPS_EXPECT_TRUE(mesh.rootTree->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-CYCLE"),
               "repeated refine/coarsen cycles restore topology and volume") {
  Mesh mesh;
  InitTopologyOnly(mesh);
  Node* const root = mesh.rootTree;
  const std::uint64_t expected_receipt = RootReceipt(root);
  const double tolerance =
      32.0 * std::numeric_limits<double>::epsilon();

  for (int cycle = 0; cycle < 32; ++cycle) {
    AMPS_EXPECT_TRUE(mesh.splitTreeNode(root));
    AMPS_EXPECT_NEAR(LeafVolume(mesh), 1.0, tolerance);
    AMPS_EXPECT_TRUE(mesh.coarsenTreeNode(root));
    AMPS_EXPECT_TRUE(mesh.rootTree == root);
    AMPS_EXPECT_TRUE(root->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
    AMPS_EXPECT_EQ(RootReceipt(root), expected_receipt);
    AMPS_EXPECT_NEAR(NodeVolume(root), 1.0, tolerance);
    AMPS_EXPECT_EQ(mesh.GetMeshID(), 2L * (cycle + 1));
  }
}

void LinkMeshCoarsenDataTests() {}

void RegisterMeshCoarsenData() {
  if (CoarsenDataOffset<0)
    CoarsenDataOffset=
        PIC::Mesh::cDataBlockAMR::RequestInternalBlockData(sizeof(double));
}

}  // namespace CoreMeshTests
}  // namespace AMPS
