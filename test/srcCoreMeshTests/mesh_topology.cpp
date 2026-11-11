#include "mesh_test_fixture.h"

#include <cmath>
#include <cstdio>
#include <set>

namespace AMPS {
namespace CoreMeshTests {
namespace {

typedef cAmpsMesh<PIC::Mesh::cDataCornerNode, PIC::Mesh::cDataCenterNode,
                  PIC::Mesh::cDataBlockAMR> Mesh;

// One model-length unit is larger than the root cell's D-dimensional
// characteristic size.  buildMesh() must therefore perform its mandatory root
// split and then terminate without an application-specific refinement.  This
// bounded law exercises the production construction workflow independently of
// the direct-split fixture used by the remaining topology cases.
double BuildWorkflowResolution(double*) { return 1.0; }

}  // namespace

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("TOPOLOGY-UNIFORM"),
               "root split is the exact uniform dyadic partition") {
  Node* root = Root();
  const int children = 1 << _MESH_DIMENSION_;
  for (int child = 0; child < children; ++child) {
    Node* node = root->downNode[child];
    AMPS_EXPECT_TRUE(node != NULL);
    AMPS_EXPECT_EQ(static_cast<int>(node->RefinmentLevel), 1);
    for (int d = 0; d < _MESH_DIMENSION_; ++d) {
      const int upper = (child >> d) & 1;
      AMPS_EXPECT_NEAR(node->xmin[d], upper ? 0.5 : 0.0, 0.0);
      AMPS_EXPECT_NEAR(node->xmax[d], upper ? 1.0 : 0.5, 0.0);
    }
  }

  // Exercise the standard application path as well as the fixture's direct
  // split.  The private authority prevents this smoke from changing the
  // process-global PIC mesh used by subsequent exact-case invocations.
  Mesh generated;
  double xmin[3] = {0.0, 0.0, 0.0};
  double xmax[3] = {1.0, 1.0, 1.0};
  generated.AllowBlockAllocation = false;
  generated.init(xmin, xmax, BuildWorkflowResolution);
  generated.buildMesh();

  AMPS_EXPECT_FALSE(generated.rootTree->lastBranchFlag() ==
                    _BOTTOM_BRANCH_TREE_);
  AMPS_EXPECT_EQ(generated.GetMeshID(), 1L);
  AMPS_EXPECT_EQ(generated.getMeshMaxRefinmentLevel(), 1);

  std::set<Node*> generated_leaves;
  std::size_t generated_list_length = 0;
  for (Node* leaf = generated.BranchBottomNodeList; leaf != NULL;
       leaf = leaf->nextBranchBottomNode) {
    ++generated_list_length;
    AMPS_EXPECT_TRUE(leaf->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
    AMPS_EXPECT_TRUE(generated_leaves.insert(leaf).second);
    AMPS_EXPECT_EQ(static_cast<int>(leaf->RefinmentLevel), 1);
  }
  AMPS_EXPECT_EQ(generated_list_length,
                 static_cast<std::size_t>(1 << _MESH_DIMENSION_));
  AMPS_EXPECT_EQ(generated_leaves.size(), generated_list_length);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("TOPOLOGY-NONUNIFORM"),
               "locally refined leaves obey a strict 2:1 dyadic hierarchy") {
  std::vector<Node*> leaves = Leaves();
  const std::size_t n = static_cast<std::size_t>(1 << _MESH_DIMENSION_);
  AMPS_EXPECT_EQ(leaves.size(), 2 * n - 1);

  std::size_t level_one = 0;
  std::size_t level_two = 0;
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    if (leaves[i]->RefinmentLevel == 1) ++level_one;
    if (leaves[i]->RefinmentLevel == 2) ++level_two;
  }
  AMPS_EXPECT_EQ(level_one, n - 1);
  AMPS_EXPECT_EQ(level_two, n);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("SEARCH"),
               "search obeys half-open ownership at interiors, interfaces, and limits") {
  std::vector<Node*> leaves = Leaves();

  // First exercise ordinary interior points.  Midpoints are deliberately
  // dyadic and therefore cannot move across a leaf boundary by roundoff.
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    double midpoint[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < _MESH_DIMENSION_; ++d)
      midpoint[d] = 0.5 * (leaves[i]->xmin[d] + leaves[i]->xmax[d]);
    AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(midpoint) == leaves[i]);
  }

  // AMPS owns Cartesian boxes as [xmin,xmax): an exact internal plane belongs
  // to the child on its positive side.  Check that convention independently
  // at every face of both split levels.  Other coordinates stay strictly in
  // the all-lower branch so the expected child bit is unambiguous.
  for (int normal = 0; normal < _MESH_DIMENSION_; ++normal) {
    double fine_face[3] = {0.125, 0.125, 0.125};
    fine_face[normal] = 0.25;
    AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(fine_face) ==
                     Root()->downNode[0]->downNode[1 << normal]);

    double coarse_face[3] = {0.25, 0.25, 0.25};
    coarse_face[normal] = 0.5;
    AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(coarse_face) ==
                     Root()->downNode[1 << normal]);
  }

#if _MESH_DIMENSION_ >= 2
  // Intersections of two root planes are a corner in 2-D and an edge in 3-D.
  // The positive-side convention must be applied in both directions at once.
  for (int first = 0; first < _MESH_DIMENSION_; ++first)
    for (int second = first + 1; second < _MESH_DIMENSION_; ++second) {
      double interface[3] = {0.25, 0.25, 0.25};
      interface[first] = interface[second] = 0.5;
      AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(interface) ==
                       Root()->downNode[(1 << first) | (1 << second)]);
    }
#endif

#if _MESH_DIMENSION_ == 3
  double central_corner[3] = {0.5, 0.5, 0.5};
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(central_corner) ==
                   Root()->downNode[7]);
#endif

  // The lower domain limit is included and the upper one is excluded.  Test
  // one representable value on each side as well as the limits themselves;
  // this catches accidental inclusive-upper tests and integer-index clamping.
  double lower[3] = {0.0, 0.0, 0.0};
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(lower) ==
                   Root()->downNode[0]->downNode[0]);
  for (int d = 0; d < _MESH_DIMENSION_; ++d) {
    double outside_lower[3] = {0.125, 0.125, 0.125};
    outside_lower[d] = std::nextafter(0.0, -1.0);
    AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(outside_lower) == NULL);

    double upper_limit[3] = {0.75, 0.75, 0.75};
    upper_limit[d] = 1.0;
    AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(upper_limit) == NULL);
  }

  double inside_upper[3] = {0.0, 0.0, 0.0};
  for (int d = 0; d < _MESH_DIMENSION_; ++d)
    inside_upper[d] = std::nextafter(1.0, 0.0);
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(inside_upper) ==
                   Root()->downNode[(1 << _MESH_DIMENSION_) - 1]);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("PARENT-CHILD"),
               "parent links, levels, and integer dyadic coordinates agree") {
  std::vector<Node*> all = AllNodes();
  Node* root = Root();
  AMPS_EXPECT_TRUE(root->upNode == NULL);
  for (std::size_t i = 0; i < all.size(); ++i) {
    Node* parent = all[i];
    if (parent->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) continue;
    for (int child = 0; child < (1 << _MESH_DIMENSION_); ++child) {
      Node* node = parent->downNode[child];
      AMPS_EXPECT_TRUE(node->upNode == parent);
      AMPS_EXPECT_EQ(static_cast<int>(node->RefinmentLevel),
                     static_cast<int>(parent->RefinmentLevel) + 1);
      AMPS_EXPECT_EQ(node->NodeGeometricSizeIndex,
                     parent->NodeGeometricSizeIndex / 2);
      for (int d = 0; d < _MESH_DIMENSION_; ++d) {
        const int expected = parent->xMinGlobalIndex[d] +
            (((child >> d) & 1) ? node->NodeGeometricSizeIndex : 0);
        AMPS_EXPECT_EQ(node->xMinGlobalIndex[d], expected);
      }
    }
  }
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("NEIGHBORS"),
               "face, coarse-fine, edge, and corner slots name geometric neighbors") {
  Node* refined_parent = Root()->downNode[0];

  // Same-level face reciprocity inside the refined parent.
  for (int d = 0; d < _MESH_DIMENSION_; ++d) {
    Node* lower = refined_parent->downNode[0];
    Node* upper = refined_parent->downNode[1 << d];
    Node* forward = lower->GetNeibFace(2 * d + 1, 0, 0, PIC::Mesh::mesh);
    Node* backward = upper->GetNeibFace(2 * d, 0, 0, PIC::Mesh::mesh);
    AMPS_EXPECT_TRUE(forward == upper);
    AMPS_EXPECT_TRUE(backward == lower);
  }

  // Across a 2:1 interface the coarse face exposes one slot per fine leaf:
  // one slot in 1-D, two in 2-D, and four in 3-D.  Tangent bit zero is the
  // lower half and bit one the upper half.  Conversely every sub-slot of a
  // fine face resolves to the same coarse leaf.
  for (int normal = 0; normal < _MESH_DIMENSION_; ++normal) {
    int tangent[2] = {-1, -1};
    int n_tangent = 0;
    for (int d = 0; d < _MESH_DIMENSION_; ++d)
      if (d != normal) tangent[n_tangent++] = d;

    Node* coarse = Root()->downNode[1 << normal];
    const int i_max = (_MESH_DIMENSION_ >= 2) ? 2 : 1;
    const int j_max = (_MESH_DIMENSION_ == 3) ? 2 : 1;
    for (int j = 0; j < j_max; ++j)
      for (int i = 0; i < i_max; ++i) {
        int fine_code = 1 << normal;
        if (n_tangent > 0 && i != 0) fine_code |= 1 << tangent[0];
        if (n_tangent > 1 && j != 0) fine_code |= 1 << tangent[1];
        Node* fine = refined_parent->downNode[fine_code];

        AMPS_EXPECT_TRUE(coarse->GetNeibFace(2 * normal, i, j,
                                             PIC::Mesh::mesh) == fine);
        AMPS_EXPECT_TRUE(fine->GetNeibFace(2 * normal + 1, i, j,
                                           PIC::Mesh::mesh) == coarse);
      }
  }

#if _MESH_DIMENSION_ >= 2
  // Corner number is the sign-bit code x+2*y(+4*z).  From child c, the
  // direction toward the parent's center and the diagonally opposite child
  // are both the D-bit complement of c.
  const int child_mask = (1 << _MESH_DIMENSION_) - 1;
  for (int child = 0; child <= child_mask; ++child) {
    const int opposite = child_mask ^ child;
    AMPS_EXPECT_TRUE(refined_parent->downNode[child]->GetNeibCorner(
                         opposite, PIC::Mesh::mesh) ==
                     refined_parent->downNode[opposite]);
  }
#endif

#if _MESH_DIMENSION_ == 3
  // The three internal edges incident on child 0 follow AMPS' documented
  // edge ordering: edge 2 crosses +y,+z; 6 crosses +x,+z; 10 crosses +x,+y.
  // Both edge segments must identify the same equal-level neighbor here.
  const int edge_number[3] = {2, 6, 10};
  const int edge_neighbor[3] = {6, 5, 3};
  for (int edge = 0; edge < 3; ++edge)
    for (int segment = 0; segment < 2; ++segment)
      AMPS_EXPECT_TRUE(refined_parent->downNode[0]->GetNeibEdge(
                           edge_number[edge], segment, PIC::Mesh::mesh) ==
                       refined_parent->downNode[edge_neighbor[edge]]);
#endif
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("CHECKSUM"),
               "canonical topology and native mesh files preserve unique identities") {
  const std::uint64_t first = CanonicalTopologyChecksum();
  const std::uint64_t second = CanonicalTopologyChecksum();
  AMPS_EXPECT_TRUE(first != 0);
  AMPS_EXPECT_EQ(first, second);

  // AMR IDs are the production serialization/search identity.  Pairwise
  // comparison is independent of pointer allocation order and catches a child
  // whose bit path was not included in its identifier.
  std::vector<Node*> all = AllNodes();
  for (std::size_t i = 0; i < all.size(); ++i)
    for (std::size_t j = i + 1; j < all.size(); ++j)
      AMPS_EXPECT_TRUE(all[i]->AMRnodeID != all[j]->AMRnodeID);

  // Exercise the real topology-file writer/reader, not only the address-free
  // checksum above. The comet integration exposed an uninitialized face
  // cache in compact DIM1/DIM2 builds: saveTreeStructure dereferenced garbage
  // even though geometric neighbor queries were correct. A separate mesh
  // keeps this persistence regression independent of other cases' allocation
  // state and preserves the process-global PIC mesh authority.
  Mesh original, restored;
  double xmin[3] = {0.0, 0.0, 0.0};
  double xmax[3] = {1.0, 1.0, 1.0};
  original.AllowBlockAllocation = restored.AllowBlockAllocation = false;
  original.init(xmin, xmax, BuildWorkflowResolution);
  original.buildMesh();
  restored.init(xmin, xmax, BuildWorkflowResolution);
  const char* filename = "core-mesh-native-roundtrip.msh";
  original.saveMeshFile(filename);
  restored.readMeshFile(filename);

  // The oracle is the exact root split of a unit D-cube. Child bit d chooses
  // [0,1/2] or [1/2,1] on axis d. Verify identities, parent links, midpoint
  // ownership and face neighbors against that geometry after deserialization.
  // Comparing only a second production-generated hash could hide a shared
  // serialization error, and comparing pointers across meshes would be wrong.
  const int child_count = 1 << _MESH_DIMENSION_;
  for (int child = 0; child < child_count; ++child) {
    Node* node = restored.rootTree->downNode[child];
    AMPS_EXPECT_TRUE(node != NULL);
    AMPS_EXPECT_TRUE(node->upNode == restored.rootTree);
    AMPS_EXPECT_EQ(static_cast<int>(node->RefinmentLevel), 1);
    AMPS_EXPECT_TRUE(node->AMRnodeID == original.rootTree->downNode[child]->AMRnodeID);
    double midpoint[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < _MESH_DIMENSION_; ++d) {
      const bool upper = ((child >> d) & 1) != 0;
      AMPS_EXPECT_NEAR(node->xmin[d], upper ? 0.5 : 0.0, 0.0);
      AMPS_EXPECT_NEAR(node->xmax[d], upper ? 1.0 : 0.5, 0.0);
      midpoint[d] = upper ? 0.75 : 0.25;
      const int face = 2*d + (upper ? 0 : 1);
      AMPS_EXPECT_TRUE(node->GetNeibFace(face, 0, 0, &restored) ==
          restored.rootTree->downNode[child ^ (1 << d)]);
      // The native lower-dimensional file actually stores these face slots.
      // Check them as well as the geometric accessor: zero-initialization
      // alone would prevent the crash while silently losing setter updates.
#if _MESH_DIMENSION_ == 1
      AMPS_EXPECT_TRUE(node->neibNodeFace[face] ==
          restored.rootTree->downNode[child ^ (1 << d)]);
#elif _MESH_DIMENSION_ == 2
      AMPS_EXPECT_TRUE(node->neibNodeFace[2*face] ==
          restored.rootTree->downNode[child ^ (1 << d)]);
#endif
    }
    AMPS_EXPECT_TRUE(restored.findTreeNode(midpoint) == node);
  }
  AMPS_EXPECT_EQ(std::remove(filename), 0);
}

void LinkMeshTopologyTests() {}

}  // namespace CoreMeshTests
}  // namespace AMPS
