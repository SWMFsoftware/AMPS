#include "mesh_test_fixture.h"

#include <cmath>
#include <limits>
#include <set>

namespace AMPS {
namespace CoreMeshTests {

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-REFINE"),
               "two production splits create the independently counted tree") {
  const std::size_t branching = static_cast<std::size_t>(1 << _MESH_DIMENSION_);
  std::vector<Node*> all = AllNodes();
  std::vector<Node*> leaves = Leaves();
  AMPS_EXPECT_EQ(all.size(), 1 + 2 * branching);
  AMPS_EXPECT_EQ(leaves.size(), 2 * branching - 1);
  AMPS_EXPECT_FALSE(Root()->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
  AMPS_EXPECT_FALSE(Root()->downNode[0]->lastBranchFlag() ==
                    _BOTTOM_BRANCH_TREE_);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-VOLUME"),
               "leaf D-volumes equal dyadic theory and partition the unit domain") {
  std::vector<Node*> leaves = Leaves();
  double total = 0.0;
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    const double expected = std::ldexp(
        1.0, -_MESH_DIMENSION_ * static_cast<int>(leaves[i]->RefinmentLevel));
    const double actual = DyadicLeafVolume(leaves[i]);
    AMPS_EXPECT_NEAR(actual, expected, 8.0 * std::numeric_limits<double>::epsilon());
    total += actual;
  }
  AMPS_EXPECT_NEAR(total, 1.0, 32.0 * std::numeric_limits<double>::epsilon());
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-COUNTERS"),
               "mesh generation increments exactly once for each accepted split") {
  EnsureTopologyMesh();
  // init() resets the generation to zero.  The fixture accepts exactly the
  // root and lower-child splits and performs no distribution or allocation.
  AMPS_EXPECT_EQ(PIC::Mesh::mesh->GetMeshID(), 2L);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-LOAD"),
               "bottom-list membership is unique before load values are assigned") {
  EnsureTopologyMesh();

  // Load balancing consumes BranchBottomNodeList rather than recursively
  // walking the tree.  Prove that this production list is a one-to-one image
  // of the independent recursive leaf traversal before checking callback
  // values; otherwise a duplicate and an omission could preserve aggregate
  // load and escape detection.
  const std::vector<Node*> leaves = Leaves();
  std::set<Node*> expected(leaves.begin(), leaves.end());
  std::set<Node*> listed;
  std::size_t list_length = 0;
  for (Node* node = PIC::Mesh::mesh->BranchBottomNodeList;
       node != NULL; node = node->nextBranchBottomNode) {
    ++list_length;
    AMPS_EXPECT_TRUE(node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
    AMPS_EXPECT_TRUE(expected.find(node) != expected.end());
    AMPS_EXPECT_TRUE(listed.insert(node).second);
  }
  AMPS_EXPECT_EQ(list_length, leaves.size());
  AMPS_EXPECT_EQ(listed.size(), expected.size());
  for (std::set<Node*>::const_iterator it = expected.begin();
       it != expected.end(); ++it)
    AMPS_EXPECT_TRUE(listed.find(*it) != listed.end());

  PIC::Mesh::mesh->SetParallelLoadMeasure(LevelLoad);
  for (std::size_t i = 0; i < leaves.size(); ++i)
    AMPS_EXPECT_NEAR(leaves[i]->ParallelLoadMeasure,
                     1.0 + leaves[i]->RefinmentLevel, 0.0);

  // Upper branches are explicitly assigned zero by SetParallelLoadMeasure;
  // stale work estimates must not enter a later space-filling-curve balance.
  AMPS_EXPECT_NEAR(Root()->ParallelLoadMeasure, 0.0, 0.0);
  AMPS_EXPECT_NEAR(Root()->downNode[0]->ParallelLoadMeasure, 0.0, 0.0);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("ADAPT-MAX-LEVEL"),
               "maximum-level refinement is rejected without changing generation") {
  EnsureTopologyMesh();
  Node synthetic;
  synthetic.RefinmentLevel = _MAX_REFINMENT_LEVEL_;
  const long before = PIC::Mesh::mesh->GetMeshID();
  AMPS_EXPECT_FALSE(PIC::Mesh::mesh->splitTreeNode(&synthetic));
  AMPS_EXPECT_EQ(PIC::Mesh::mesh->GetMeshID(), before);
  AMPS_EXPECT_TRUE(synthetic.lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
}

void LinkMeshAdaptationTests() {}

}  // namespace CoreMeshTests
}  // namespace AMPS
