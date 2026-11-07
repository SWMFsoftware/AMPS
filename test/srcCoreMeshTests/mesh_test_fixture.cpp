#include "mesh_test_fixture.h"

#include <cmath>
#include <cstring>

namespace AMPS {
namespace CoreMeshTests {
namespace {

bool topology_ready = false;
bool allocation_ready = false;
bool interpolation_ready = false;

// init() stores this callback for later adaptation.  The direct fixture splits
// are deliberate, so the value merely states the coarsest admissible physical
// cell length in model length units and does not drive construction.
double FixtureResolution(double*) { return 1.0; }

void Collect(Node* node, std::vector<Node*>& all, std::vector<Node*>& leaves) {
  AMPS_EXPECT_TRUE(node != NULL);
  all.push_back(node);
  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    leaves.push_back(node);
    return;
  }

  for (int child = 0; child < (1 << _MESH_DIMENSION_); ++child) {
    AMPS_EXPECT_TRUE(node->downNode[child] != NULL);
    Collect(node->downNode[child], all, leaves);
  }
}

std::uint64_t Fnv1a(std::uint64_t hash, std::uint64_t value) {
  // FNV-1a is used only as a deterministic receipt, not as a physical oracle.
  // Feed bytes explicitly so results do not depend on host endianness.
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

}  // namespace

void EnsureTopologyMesh() {
  if (topology_ready) return;

  AMPS_EXPECT_TRUE(PIC::Mesh::mesh != NULL);
  double xmin[3] = {0.0, 0.0, 0.0};
  double xmax[3] = {1.0, 1.0, 1.0};

  // Keep blocks unallocated while topology is changing.  This separates tree
  // connectivity from center/corner storage and avoids assigning any
  // unqualified prolongation semantics to associated application data.
  PIC::Mesh::mesh->AllowBlockAllocation = false;
  PIC::Mesh::mesh->init(xmin, xmax, FixtureResolution);
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree));
  AMPS_EXPECT_TRUE(
      PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree->downNode[0]));
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();
  PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
  topology_ready = true;
}

void EnsureAllocatedMesh() {
  EnsureTopologyMesh();
  if (allocation_ready) return;

  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation = true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();
  EnsureInterpolationRoutines();
  allocation_ready = true;
}

void EnsureInterpolationRoutines() {
  // Init() creates dimension-dependent interpolation lookup storage.  Keep it
  // process-global just as the production application does: exhaustive mask
  // tests may replace the active mesh, but must not repeatedly initialize the
  // interpolation subsystem or depend on another case having run first.
  if (interpolation_ready) return;
  PIC::InterpolationRoutines::Init();
  interpolation_ready = true;
}

Node* Root() {
  EnsureTopologyMesh();
  return PIC::Mesh::mesh->rootTree;
}

std::vector<Node*> Leaves() {
  EnsureTopologyMesh();
  std::vector<Node*> all;
  std::vector<Node*> leaves;
  Collect(PIC::Mesh::mesh->rootTree, all, leaves);
  return leaves;
}

std::vector<Node*> AllNodes() {
  EnsureTopologyMesh();
  std::vector<Node*> all;
  std::vector<Node*> leaves;
  Collect(PIC::Mesh::mesh->rootTree, all, leaves);
  return all;
}

double DyadicLeafVolume(const Node* node) {
  double volume = 1.0;
  for (int dimension = 0; dimension < _MESH_DIMENSION_; ++dimension)
    volume *= node->xmax[dimension] - node->xmin[dimension];
  return volume;
}

std::uint64_t CanonicalTopologyChecksum() {
  std::vector<Node*> all = AllNodes();
  std::vector<std::uint64_t> records;
  records.reserve(all.size());
  for (std::size_t n = 0; n < all.size(); ++n) {
    std::uint64_t record = static_cast<unsigned char>(all[n]->RefinmentLevel);
    record = record * UINT64_C(1315423911) +
             static_cast<std::uint64_t>(all[n]->NodeGeometricSizeIndex);
    for (int d = 0; d < _MESH_DIMENSION_; ++d)
      record = record * UINT64_C(1315423911) +
               static_cast<std::uint64_t>(all[n]->xMinGlobalIndex[d]);
    records.push_back(record);
  }
  std::sort(records.begin(), records.end());

  std::uint64_t checksum = UINT64_C(1469598103934665603);
  for (std::size_t n = 0; n < records.size(); ++n)
    checksum = Fnv1a(checksum, records[n]);
  return checksum;
}

double UnitLoad(Node*) { return 1.0; }
double LevelLoad(Node* node) { return 1.0 + node->RefinmentLevel; }

}  // namespace CoreMeshTests
}  // namespace AMPS
