#include "pic.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// T10 is a collective test executable: unlike process-local unit fixtures,
// every case is entered by every MPI rank.  Local assertions are accumulated
// instead of throwing so a failed rank cannot strand its peers in a later
// collective.  The final MAX reduction turns any local failure into a global
// failure while preserving a rank-by-rank receipt for the Python orchestrator.

namespace {

typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> Node;

int center_data_offset = -1;
int corner_data_offset = -1;

int RequestCenterData(int offset) {
  center_data_offset = offset;
  return 2 * static_cast<int>(sizeof(double));
}

int RequestCornerData(int offset) {
  corner_data_offset = offset;
  return 2 * static_cast<int>(sizeof(double));
}

struct CheckContext {
  bool ok;
  std::vector<std::string> errors;
  std::string receipt;

  CheckContext() : ok(true) {}

  void Check(bool condition, const char* format, ...) {
    if (condition) return;
    ok = false;
    char buffer[512];
    va_list arguments;
    va_start(arguments, format);
    std::vsnprintf(buffer, sizeof(buffer), format, arguments);
    va_end(arguments);
    if (errors.size() < 8) errors.push_back(buffer);
  }

  std::string Message() const {
    if (ok) return receipt.empty() ? "ok" : receipt;
    std::ostringstream out;
    for (std::size_t i = 0; i < errors.size(); ++i) {
      if (i != 0) out << "; ";
      out << errors[i];
    }
    return out.str();
  }
};

double UnitResolution(double*) { return 1.0; }

double PeriodicResolution(double*) {
  // The user domain is two metres wide and blocks contain two cells.  A
  // requested cell length of L/16 gives the periodic real domain and its
  // production ghost shell the same dyadic refinement level.
  return 0.125;
}

double UnitLoad(Node*) { return 1.0; }

double SkewLoad(Node* node) {
  // A deterministic, strongly nonuniform manufactured cost.  The function is
  // evaluated only for the current owner, uses replicated geometry, and has
  // no dependence on rank.  Consequently all decompositions target the same
  // weighted Morton partition while actually moving blocks at nRank > 1.
  const double x = 0.5 * (node->xmin[0] + node->xmax[0]);
  const double y = 0.5 * (node->xmin[1] + node->xmax[1]);
  return 1.0 + (x < -0.25 ? 23.0 : 0.0) + (y > 0.25 ? 7.0 : 0.0);
}

void CollectLeaves(Node* node, std::vector<Node*>& leaves) {
  if (node == NULL) return;
  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    if (node->IsUsedInCalculationFlag) leaves.push_back(node);
    return;
  }
  for (int child = 0; child < (1 << _MESH_DIMENSION_); ++child)
    CollectLeaves(node->downNode[child], leaves);
}

std::vector<Node*> Leaves() {
  std::vector<Node*> leaves;
  CollectLeaves(PIC::Mesh::mesh->rootTree, leaves);
  return leaves;
}

std::uint64_t Mix(std::uint64_t hash, std::uint64_t value) {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

std::uint64_t NodeKey(const Node* node) {
  // This identity is derived only from dyadic AMR geometry.  It deliberately
  // does not use pointer values, list order, Temp_ID, or owning rank.
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = Mix(hash, static_cast<unsigned char>(node->RefinmentLevel));
  hash = Mix(hash, static_cast<std::uint64_t>(node->NodeGeometricSizeIndex));
  for (int d = 0; d < 3; ++d)
    hash = Mix(hash, static_cast<std::uint64_t>(node->xMinGlobalIndex[d]));
  return hash;
}

double BlockCode(const Node* node) {
  // Keep the code below 2^53 so integer-valued doubles compare exactly after
  // MPI packing, addition, block migration, and restart serialization.
  std::uint64_t key = NodeKey(node);
  return static_cast<double>(1 + (key % UINT64_C(1000000000)));
}

void Store2(char* buffer, int offset, double a, double b) {
  double value[2] = {a, b};
  std::memcpy(buffer + offset, value, sizeof(value));
}

void Load2(const char* buffer, int offset, double& a, double& b) {
  double value[2];
  std::memcpy(value, buffer + offset, sizeof(value));
  a = value[0];
  b = value[1];
}

double CenterValue(const Node* node, int i, int j, int k) {
  return BlockCode(node) + i + 4 * j + 16 * k;
}

double CornerValue(const Node* node, int i, int j, int k) {
  return BlockCode(node) + 64 + i + 4 * j + 16 * k;
}

void SetBlockFields(Node* node, double center_shift, double corner_shift) {
  if (node == NULL || node->block == NULL) return;
  for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
    for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
      for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
        PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
            PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
        const double value = CenterValue(node, i, j, k);
        Store2(cell->GetAssociatedDataBufferPointer(), center_data_offset,
               value + center_shift, 2.0 * value + center_shift);
      }
  for (int k = 0; k <= _BLOCK_CELLS_Z_; ++k)
    for (int j = 0; j <= _BLOCK_CELLS_Y_; ++j)
      for (int i = 0; i <= _BLOCK_CELLS_X_; ++i) {
        PIC::Mesh::cDataCornerNode* corner = node->block->GetCornerNode(
            PIC::Mesh::mesh->getCornerNodeLocalNumber(i, j, k));
        if (corner == NULL) continue;
        const double value = CornerValue(node, i, j, k);
        Store2(corner->GetAssociatedDataBufferPointer(), corner_data_offset,
               value + corner_shift, 2.0 * value + corner_shift);
      }
}

void SetOwnerFields(double center_shift = 0.0, double corner_shift = 0.0) {
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread)
    SetBlockFields(node, center_shift, corner_shift);
}

void SetGhostCenterFields(double a, double b) {
  if (PIC::Mesh::mesh->DomainBoundaryLayerNodesList == NULL) return;
  for (int owner = 0; owner < PIC::nTotalThreads; ++owner) {
    if (owner == PIC::ThisThread) continue;
    for (Node* node = PIC::Mesh::mesh->DomainBoundaryLayerNodesList[owner];
         node != NULL; node = node->nextNodeThisThread) {
      if (node->block == NULL) continue;
      for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
        for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
          for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
            PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
                PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
            Store2(cell->GetAssociatedDataBufferPointer(), center_data_offset,
                   a, b);
          }
    }
  }
}

void CheckOwnerFields(CheckContext& check) {
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread) {
    check.Check(node->block != NULL, "owner block is not allocated");
    if (node->block == NULL) continue;
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          double a, b;
          PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
              PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
          Load2(cell->GetAssociatedDataBufferPointer(), center_data_offset, a, b);
          const double expected = CenterValue(node, i, j, k);
          check.Check(a == expected && b == 2.0 * expected,
                      "migrated center datum differs from owner oracle");
        }
  }
}

std::vector<std::pair<std::uint64_t, int> > OwnerMap(CheckContext& check) {
  std::vector<std::pair<std::uint64_t, int> > map;
  // Thread is the replicated ownership authority.  Only the current rank's
  // ParallelNodesDistributionList is a traversable load list after the
  // domain-boundary layer is installed: remote leaves reuse their intrusive
  // next/prev links for this rank's received-ghost lists.  Traversing remote
  // owner lists would therefore test an invalid representation rather than an
  // AMPS contract.  Build the global map from the replicated tree, then prove
  // that the one local load list contains every and only locally owned leaf.
  const std::vector<Node*> leaves = Leaves();
  std::set<Node*> expected_local;
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    check.Check(leaves[index]->Thread >= 0 &&
                    leaves[index]->Thread < PIC::nTotalThreads,
                "leaf has an out-of-range owning rank");
    map.push_back(std::make_pair(NodeKey(leaves[index]),
                                 leaves[index]->Thread));
    if (leaves[index]->Thread == PIC::ThisThread)
      expected_local.insert(leaves[index]);
  }
  std::set<Node*> observed_local;
  Node* previous = NULL;
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread) {
    check.Check(node->Thread == PIC::ThisThread,
                "local load list contains a remotely owned leaf");
    check.Check(node->prevNodeThisThread == previous,
                "local distribution-list backward link is inconsistent");
    check.Check(observed_local.insert(node).second,
                "local load list contains a duplicate leaf");
    previous = node;
  }
  check.Check(observed_local == expected_local,
              "local load list does not equal the replicated owner map");
  std::sort(map.begin(), map.end());
  check.Check(map.size() == leaves.size(),
              "replicated owner map does not contain every used leaf");
  return map;
}

std::uint64_t OwnerMapChecksum(
    const std::vector<std::pair<std::uint64_t, int> >& map) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  for (std::size_t i = 0; i < map.size(); ++i) {
    hash = Mix(hash, map[i].first);
    hash = Mix(hash, static_cast<std::uint64_t>(map[i].second));
  }
  return hash;
}

void CheckReplicatedChecksum(CheckContext& check, std::uint64_t local) {
  std::vector<unsigned long long> all(PIC::nTotalThreads);
  const unsigned long long send = static_cast<unsigned long long>(local);
  MPI_Allgather(&send, 1, MPI_UNSIGNED_LONG_LONG, &all[0], 1,
                MPI_UNSIGNED_LONG_LONG, MPI_GLOBAL_COMMUNICATOR);
  for (int rank = 0; rank < PIC::nTotalThreads; ++rank)
    check.Check(all[rank] == send, "replicated owner map differs on rank %d", rank);
}

bool TouchesRank(Node* node, int rank, int* same_level, int* coarse_fine,
                 int* face_count = NULL, int* edge_count = NULL,
                 int* corner_count = NULL) {
  bool found = false;
  for (int face = 0; face < 6; ++face)
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) {
        Node* neighbor = node->GetNeibFace(face, a, b, PIC::Mesh::mesh);
        if (neighbor != NULL && neighbor->Thread == rank) {
          found = true;
          if (face_count != NULL) ++*face_count;
          if (neighbor->RefinmentLevel == node->RefinmentLevel)
            ++*same_level;
          else
            ++*coarse_fine;
        }
      }
  for (int edge = 0; edge < 12; ++edge)
    for (int segment = 0; segment < 2; ++segment) {
      Node* neighbor = node->GetNeibEdge(edge, segment, PIC::Mesh::mesh);
      if (neighbor != NULL && neighbor->Thread == rank) {
        found = true;
        if (edge_count != NULL) ++*edge_count;
        if (neighbor->RefinmentLevel == node->RefinmentLevel)
          ++*same_level;
        else
          ++*coarse_fine;
      }
    }
  for (int corner = 0; corner < 8; ++corner) {
    Node* neighbor = node->GetNeibCorner(corner, PIC::Mesh::mesh);
    if (neighbor != NULL && neighbor->Thread == rank) {
      found = true;
      if (corner_count != NULL) ++*corner_count;
      if (neighbor->RefinmentLevel == node->RefinmentLevel)
        ++*same_level;
      else
        ++*coarse_fine;
    }
  }
  return found;
}

long int& CellHead(Node* node, int local = 0) {
  return node->block->FirstCellParticleTable[local];
}

void ParticlePosition(const Node* node, double* x) {
  for (int d = 0; d < 3; ++d) x[d] = 0.5 * (node->xmin[d] + node->xmax[d]);
}

void CellCenterPosition(const Node* node, int i, int j, int k, double* x) {
  // AMR node bounds are exact dyadic coordinates.  The cell-centre oracle is
  // x_min + (index+1/2) Delta x and is intentionally independent of any PIC
  // interpolation or center-node storage routine under test.
  const int index[3] = {i, j, k};
  const int cells[3] = {
      _BLOCK_CELLS_X_, _BLOCK_CELLS_Y_, _BLOCK_CELLS_Z_};
  for (int d = 0; d < 3; ++d)
    x[d] = node->xmin[d] + (index[d] + 0.5) *
        (node->xmax[d] - node->xmin[d]) / cells[d];
}

long int AddParticle(Node* node, double tag) {
  long int& head = CellHead(node);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(head, true);
  double x[3], v[3] = {tag, 0.5 * tag, -0.25 * tag};
  ParticlePosition(node, x);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
  PIC::ParticleBuffer::SetX(x, particle);
  PIC::ParticleBuffer::SetV(v, particle);
  return particle;
}

struct ParticleDigest {
  long long count;
  unsigned long long hash;
  double velocity_sum;
  ParticleDigest() : count(0), hash(0), velocity_sum(0.0) {}
};

void AccumulateList(long int head, ParticleDigest& digest) {
  for (long int particle = head; particle != -1;
       particle = PIC::ParticleBuffer::GetNext(particle)) {
    ++digest.count;
    const double* v = PIC::ParticleBuffer::GetV(particle);
    digest.velocity_sum += v[0];
    std::uint64_t bits = 0;
    std::memcpy(&bits, v, sizeof(bits));
    digest.hash ^= static_cast<unsigned long long>(Mix(NodeKey(
        PIC::Mesh::mesh->findTreeNode(PIC::ParticleBuffer::GetX(particle))), bits));
  }
}

ParticleDigest LocalOwnerParticleDigest() {
  ParticleDigest result;
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread) {
    if (node->block == NULL) continue;
    const int cells = _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_;
    for (int cell = 0; cell < cells; ++cell)
      AccumulateList(node->block->FirstCellParticleTable[cell], result);
  }
  return result;
}

ParticleDigest GlobalDigest() {
  const ParticleDigest local = LocalOwnerParticleDigest();
  ParticleDigest global;
  MPI_Allreduce(&local.count, &global.count, 1, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.hash, &global.hash, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local.velocity_sum, &global.velocity_sum, 1, MPI_DOUBLE,
                MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  return global;
}

void SeedOwnerParticles() {
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread) {
    const double tag = static_cast<double>(1 + NodeKey(node) % UINT64_C(1000000));
    AddParticle(node, tag);
  }
}

void DeleteList(long int& head) {
  long int particle = head;
  while (particle != -1) {
    const long int next = PIC::ParticleBuffer::GetNext(particle);
    PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(particle, true);
    particle = next;
  }
  head = -1;
}

void ClearOwnerParticles() {
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread) {
    if (node->block == NULL) continue;
    const int cells = _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_;
    for (int cell = 0; cell < cells; ++cell)
      DeleteList(node->block->FirstCellParticleTable[cell]);
  }
}

void InitializeMesh(CheckContext& check) {
  check.Check(center_data_offset >= 0 && corner_data_offset >= 0,
              "application associated-data offsets were not registered");
  double xmin[3] = {-1.0, -1.0, -1.0};
  double xmax[3] = {1.0, 1.0, 1.0};
  PIC::Mesh::mesh->AllowBlockAllocation = false;

#if _PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_ON_
  PIC::BC::ExternalBoundary::Periodic::Init(xmin, xmax, PeriodicResolution);
  PIC::Mesh::mesh->buildMesh();
#else
  PIC::Mesh::mesh->init(xmin, xmax, UnitResolution);
  check.Check(PIC::Mesh::mesh->splitTreeNode(PIC::Mesh::mesh->rootTree),
              "could not split root AMR node");
  check.Check(PIC::Mesh::mesh->splitTreeNode(
                  PIC::Mesh::mesh->rootTree->downNode[0]),
              "could not create the coarse/fine interface");
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();
  PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
#endif

  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation = true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();
  PIC::DomainBlockDecomposition::UpdateBlockTable();
#if _PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_ON_
  PIC::BC::ExternalBoundary::Periodic::InitBlockPairTable(true);
#endif
}

void CaseOwnerMap(CheckContext& check) {
  const std::vector<std::pair<std::uint64_t, int> > map = OwnerMap(check);
  const std::uint64_t checksum = OwnerMapChecksum(map);
  CheckReplicatedChecksum(check, checksum);
  std::vector<int> counts(PIC::nTotalThreads, 0);
  for (std::size_t i = 0; i < map.size(); ++i) ++counts[map[i].second];
  for (int rank = 0; rank < PIC::nTotalThreads; ++rank)
    check.Check(counts[rank] > 0, "rank %d owns no used block", rank);
  std::ostringstream receipt;
  receipt << "leaves=" << map.size() << " owner_checksum=" << checksum;
  check.receipt = receipt.str();
}

void CaseOwnerGhost(CheckContext& check) {
  OwnerMap(check);
  long long local_ghosts = 0;
  std::set<Node*> unique;
  if (PIC::Mesh::mesh->DomainBoundaryLayerNodesList != NULL) {
    for (int owner = 0; owner < PIC::nTotalThreads; ++owner) {
      for (Node* node = PIC::Mesh::mesh->DomainBoundaryLayerNodesList[owner];
           node != NULL; node = node->nextNodeThisThread) {
        ++local_ghosts;
        check.Check(owner != PIC::ThisThread,
                    "local owner appears in received-ghost list");
        check.Check(node->Thread == owner, "ghost list owner is incorrect");
        check.Check(node->block != NULL, "received ghost block is unallocated");
        check.Check(unique.insert(node).second,
                    "received ghost block is duplicated locally");
        int same = 0, coarse = 0;
        check.Check(TouchesRank(node, PIC::ThisThread, &same, &coarse),
                    "received ghost does not touch a local owner block");
      }
    }
  }
  long long global_ghosts = 0;
  MPI_Allreduce(&local_ghosts, &global_ghosts, 1, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  check.Check((PIC::nTotalThreads == 1 && global_ghosts == 0) ||
                  (PIC::nTotalThreads > 1 && global_ghosts > 0),
              "ghost inventory is inconsistent with rank count");
  std::ostringstream receipt;
  receipt << "global_received_ghosts=" << global_ghosts;
  check.receipt = receipt.str();
}

void RepartitionAndCheck(CheckContext& check, bool require_particle_invariance) {
  const std::vector<std::pair<std::uint64_t, int> > before = OwnerMap(check);
  const ParticleDigest particle_before = GlobalDigest();
  SetOwnerFields();
  PIC::Mesh::mesh->SetParallelLoadMeasure(SkewLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  const std::vector<std::pair<std::uint64_t, int> > after = OwnerMap(check);
  CheckReplicatedChecksum(check, OwnerMapChecksum(after));
  int changed = 0;
  for (std::size_t i = 0; i < before.size() && i < after.size(); ++i)
    if (before[i].second != after[i].second) ++changed;
  check.Check(PIC::nTotalThreads == 1 || changed > 0,
              "weighted repartition moved no ownership at nRank > 1");
  CheckOwnerFields(check);
  if (require_particle_invariance) {
    const ParticleDigest particle_after = GlobalDigest();
    check.Check(particle_after.count == particle_before.count,
                "repartition changed global particle count");
    check.Check(particle_after.hash == particle_before.hash,
                "repartition changed particle identity checksum");
    check.Check(particle_after.velocity_sum == particle_before.velocity_sum,
                "repartition changed the velocity receipt");
  }
  std::ostringstream receipt;
  receipt << "changed_owners=" << changed
          << " owner_checksum=" << OwnerMapChecksum(after);
  check.receipt = receipt.str();
}

void CaseOwnerRepartition(CheckContext& check) {
  RepartitionAndCheck(check, false);
}

void CaseCutCell(CheckContext& check) {
  CutCell::cTriangleFace triangle;
  double a[3] = {0.0, 0.0, 0.0};
  double b[3] = {1.0, 0.0, 0.0};
  double c[3] = {0.0, 1.0, 0.0};
  triangle.SetFaceNodes(a, b, c);
  unsigned long long local_mask = 0;
  double local_t_sum = 0.0;
  unsigned long long expected_mask = 0;
  double expected_t_sum = 0.0;
  for (int sample = 0; sample < 24; ++sample) {
    const double x = (sample % 6) * 0.2 - 0.05;
    const double y = (sample / 6) * 0.2 + 0.05;
    const bool expected = x >= 0.0 && y >= 0.0 && x + y <= 1.0;
    if (expected) {
      expected_mask |= (UINT64_C(1) << sample);
      expected_t_sum += 0.5;
    }
    if (sample % PIC::nTotalThreads != PIC::ThisThread) continue;
    double origin[3] = {x, y, 1.0};
    double direction[3] = {0.0, 0.0, -2.0};
    double t = 0.0;
    const bool hit = triangle.RayIntersection(origin, direction, t, 0.0);
    check.Check(hit == expected, "distributed ray classification differs at %d", sample);
    if (hit) {
      local_mask |= (UINT64_C(1) << sample);
      local_t_sum += t;
    }
  }
  unsigned long long global_mask = 0;
  double global_t_sum = 0.0;
  MPI_Allreduce(&local_mask, &global_mask, 1, MPI_UNSIGNED_LONG_LONG, MPI_BOR,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_t_sum, &global_t_sum, 1, MPI_DOUBLE, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  check.Check(global_mask == expected_mask,
              "distributed cut-cell hit mask differs from analytic triangle");
  check.Check(global_t_sum == expected_t_sum,
              "distributed cut-cell intersection distances are not invariant");
  std::ostringstream receipt;
  receipt << "hit_mask=" << global_mask << " t_sum=" << global_t_sum;
  check.receipt = receipt.str();
}

void CheckGhostHalo(CheckContext& check, double shift_a, double shift_b) {
  long long local_checked = 0;
  if (PIC::Mesh::mesh->DomainBoundaryLayerNodesList != NULL) {
    for (int owner = 0; owner < PIC::nTotalThreads; ++owner) {
      if (owner == PIC::ThisThread) continue;
      for (Node* node = PIC::Mesh::mesh->DomainBoundaryLayerNodesList[owner];
           node != NULL; node = node->nextNodeThisThread) {
        if (node->block == NULL) continue;
        for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
          for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
            for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
              double first, second;
              PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
                  PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
              Load2(cell->GetAssociatedDataBufferPointer(), center_data_offset,
                    first, second);
              const double value = CenterValue(node, i, j, k);
              check.Check(first == value + shift_a &&
                              second == 2.0 * value + shift_b,
                          "ghost halo value differs from owner field");
              ++local_checked;
            }
      }
    }
  }
  long long global_checked = 0;
  MPI_Allreduce(&local_checked, &global_checked, 1, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  check.Check(PIC::nTotalThreads == 1 || global_checked > 0,
              "no remote ghost values were checked");
  std::ostringstream receipt;
  receipt << "checked_values=" << global_checked;
  check.receipt = receipt.str();
}

void CaseHaloReplace(CheckContext& check) {
  SetOwnerFields();
  SetGhostCenterFields(-701.0, -907.0);
  PIC::Parallel::cNodeHaloSyncManager manager;
  manager.NodeType = PIC::Parallel::eNodeType::Center;
  manager.RelativeOffsetBytesFromAssociated = center_data_offset;
  manager.nDoubles = 2;
  manager.Op = PIC::Parallel::eHaloOp::Replace;
  manager.communicate_entire_block = true;
  PIC::Parallel::cHalo halo;
  PIC::Parallel::SyncNodeHalo_DomainBoundaryLayer(manager, halo, 37100);
  CheckGhostHalo(check, 0.0, 0.0);
}

void CaseHaloAdd(CheckContext& check) {
  SetOwnerFields();
  const double add_a = 7.0, add_b = 11.0;
  SetGhostCenterFields(add_a, add_b);
  PIC::Parallel::cNodeHaloSyncManager manager;
  manager.NodeType = PIC::Parallel::eNodeType::Center;
  manager.RelativeOffsetBytesFromAssociated = center_data_offset;
  manager.nDoubles = 2;
  manager.Op = PIC::Parallel::eHaloOp::Add;
  manager.communicate_entire_block = true;
  PIC::Parallel::cHalo halo;
  PIC::Parallel::SyncNodeHalo_DomainBoundaryLayer(manager, halo, 37200);
  CheckGhostHalo(check, add_a, add_b);
}

void CaseSendMask(CheckContext& check) {
  int local_same = 0, local_coarse = 0, local_relations = 0;
  int local_center_bits = 0, local_corner_bits = 0;
  int local_faces = 0, local_edges = 0, local_corners = 0;
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread) {
    for (int destination = 0; destination < PIC::nTotalThreads; ++destination) {
      if (destination == PIC::ThisThread) continue;
      int same = 0, coarse = 0, faces = 0, edges = 0, corners = 0;
      const bool touches = TouchesRank(node, destination, &same, &coarse,
                                       &faces, &edges, &corners);
      // InitLayerBlock is a send-table callback.  AMR calls it only for an
      // actual (node,destination) communication relation; its historical
      // no-neighbor fallback intentionally creates a basic full-block mask
      // and therefore is not a predicate for discovering communication.
      // Discover contacts independently above, and invoke the production
      // callback only within its valid precondition.
      if (!touches) continue;
      std::vector<unsigned char> center(
          PIC::Mesh::BlockElementSendMask::CenterNode::GetSize());
      std::vector<unsigned char> corner(
          PIC::Mesh::BlockElementSendMask::CornerNode::GetSize());
      PIC::Mesh::BlockElementSendMask::InitLayerBlock(
          node, destination, &center[0], &corner[0]);
      int center_bits = 0, corner_bits = 0;
      for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
        for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
          for (int i = 0; i < _BLOCK_CELLS_X_; ++i)
            if (PIC::Mesh::BlockElementSendMask::CenterNode::Test(
                    i, j, k, &center[0])) ++center_bits;
      for (int k = 0; k <= _BLOCK_CELLS_Z_; ++k)
        for (int j = 0; j <= _BLOCK_CELLS_Y_; ++j)
          for (int i = 0; i <= _BLOCK_CELLS_X_; ++i)
            if (PIC::Mesh::BlockElementSendMask::CornerNode::Test(
                    i, j, k, &corner[0])) ++corner_bits;
      check.Check(center_bits > 0 && corner_bits > 0,
                  "communicating block has an empty production send mask");
      ++local_relations;
      local_same += same;
      local_coarse += coarse;
      local_center_bits += center_bits;
      local_corner_bits += corner_bits;
      local_faces += faces;
      local_edges += edges;
      local_corners += corners;
    }
  }

  PIC::Parallel::cHalo halo;
  PIC::Parallel::InitSendHaloLayer(halo, false);
  for (std::size_t index = 0; index < halo.list.size(); ++index) {
    const PIC::Parallel::cHaloEntry& entry = halo.list[index];
    check.Check(entry.startNode != NULL &&
                    entry.startNode->Thread == PIC::ThisThread,
                "send-halo entry is not locally owned");
    int same = 0, coarse = 0;
    check.Check(TouchesRank(entry.startNode, entry.ToThread, &same, &coarse),
                "send-halo destination has no geometric contact");
    for (int d = 0; d < 3; ++d) {
      const int cells = d == 0 ? _BLOCK_CELLS_X_ :
                        (d == 1 ? _BLOCK_CELLS_Y_ : _BLOCK_CELLS_Z_);
      check.Check(entry.i_cell_min[d] >= 0 &&
                      entry.i_cell_max[d] < cells &&
                      entry.i_cell_min[d] <= entry.i_cell_max[d],
                  "send-halo cell range is outside the block");
      check.Check(entry.i_corner_min[d] >= 0 &&
                      entry.i_corner_max[d] <= cells &&
                      entry.i_corner_min[d] <= entry.i_corner_max[d],
                  "send-halo corner range is outside the block");
    }
  }

  int values[8] = {local_same, local_coarse, local_relations,
                   local_center_bits, local_corner_bits, local_faces,
                   local_edges, local_corners};
  int global[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  MPI_Allreduce(values, global, 8, MPI_INT, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  if (PIC::nTotalThreads > 1) {
    check.Check(global[1] > 0, "no coarse/fine inter-rank send mask was exercised");
    check.Check(global[2] > 0 && global[3] > 0 && global[4] > 0,
                "inter-rank send masks contain no communicated nodes");
    if (PIC::nTotalThreads >= 4)
      check.Check(global[0] > 0,
                  "four-rank matrix exercised no same-level send mask");
    if (PIC::nTotalThreads >= 4)
      check.Check(global[5] > 0 && global[6] > 0 && global[7] > 0,
                  "four-rank matrix did not cover face/edge/corner contacts");
  }
  std::ostringstream receipt;
  receipt << "relations=" << global[2] << " same=" << global[0]
          << " coarse_fine=" << global[1] << " faces=" << global[5]
          << " edges=" << global[6] << " corners=" << global[7];
  check.receipt = receipt.str();
}

void CaseBlockPack(CheckContext& check) {
  Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
  check.Check(node != NULL && node->block != NULL,
              "rank has no owner block for pack/unpack");
  if (node == NULL || node->block == NULL) return;
  // Initialize every owner so the post-unpack whole-owner check distinguishes
  // corruption of the selected record from uninitialized unrelated blocks.
  SetOwnerFields();
  AddParticle(node, 101.0 + PIC::ThisThread);
  AddParticle(node, 201.0 + PIC::ThisThread);
  Node* table[1] = {node};
  unsigned long int length[1] = {0};
  PIC::Mesh::MoveBlock::GetBlockDataSize(table, 1, length);
  check.Check(length[0] > 0, "move-block packed length is zero");
  std::vector<char> buffer(length[0]);
  const int packed = PIC::Mesh::MoveBlock::PackBlockData(table, 1, &buffer[0]);
  check.Check(static_cast<unsigned long int>(packed) == length[0],
              "move-block size oracle and packer disagree");
  check.Check(CellHead(node) == -1,
              "move-block packer did not detach serialized particles");
  SetBlockFields(node, -3000.0, -4000.0);
  const int unpacked = PIC::Mesh::MoveBlock::UnpackBlockData(table, 1, &buffer[0]);
  check.Check(unpacked == packed, "move-block unpack consumed a different size");
  ParticleDigest digest;
  AccumulateList(CellHead(node), digest);
  check.Check(digest.count == 2, "move-block unpack did not restore two particles");
  CheckOwnerFields(check);
  long long local_bytes = packed, global_bytes = 0;
  MPI_Allreduce(&local_bytes, &global_bytes, 1, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  std::ostringstream receipt;
  receipt << "packed_bytes=" << global_bytes;
  check.receipt = receipt.str();
}

long long CountGhostParticles() {
  long long count = 0;
  if (PIC::Mesh::mesh->DomainBoundaryLayerNodesList == NULL) return count;
  for (int owner = 0; owner < PIC::nTotalThreads; ++owner) {
    if (owner == PIC::ThisThread) continue;
    for (Node* node = PIC::Mesh::mesh->DomainBoundaryLayerNodesList[owner];
         node != NULL; node = node->nextNodeThisThread) {
      if (node->block == NULL) continue;
      const int cells = _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_;
      for (int cell = 0; cell < cells; ++cell)
        for (long int particle = node->block->FirstCellParticleTable[cell];
             particle != -1; particle = PIC::ParticleBuffer::GetNext(particle))
          ++count;
    }
  }
  return count;
}

void CaseParticleExchange(CheckContext& check) {
  long long local_seeded = 0;
  if (PIC::Mesh::mesh->DomainBoundaryLayerNodesList != NULL) {
    for (int owner = 0; owner < PIC::nTotalThreads; ++owner) {
      if (owner == PIC::ThisThread) continue;
      Node* ghost = PIC::Mesh::mesh->DomainBoundaryLayerNodesList[owner];
      if (ghost == NULL || ghost->block == NULL) continue;
      AddParticle(ghost, 1000.0 + 100.0 * PIC::ThisThread + owner);
      ++local_seeded;
    }
  }
  if (PIC::nTotalThreads == 1) {
    Node* owner = PIC::Mesh::mesh->ParallelNodesDistributionList[0];
    AddParticle(owner, 1000.0);
    local_seeded = 1;
  }
  long long global_seeded = 0;
  MPI_Allreduce(&local_seeded, &global_seeded, 1, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  PIC::Parallel::ExchangeParticleData();
  const ParticleDigest after = GlobalDigest();
  long long local_ghost_particles = CountGhostParticles();
  long long global_ghost_particles = 0;
  MPI_Allreduce(&local_ghost_particles, &global_ghost_particles, 1,
                MPI_LONG_LONG_INT, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  check.Check(after.count == global_seeded,
              "rank-boundary exchange duplicated or lost particles");
  check.Check(global_ghost_particles == 0,
              "sent particles remain in received-ghost lists");
  std::ostringstream receipt;
  receipt << "migrated_particles=" << global_seeded
          << " identity_hash=" << after.hash;
  check.receipt = receipt.str();
}

void CaseParticleRepartition(CheckContext& check) {
  SeedOwnerParticles();
  RepartitionAndCheck(check, true);
}

void CaseOuterLoss(CheckContext& check) {
  double x[3] = {0.95, 0.70, 0.70};
  Node* node = PIC::Mesh::mesh->findTreeNode(x);
  check.Check(node != NULL, "outer-loss start has no AMR leaf");
  int local_created = 0, local_deleted = 0;
  if (node != NULL && node->Thread == PIC::ThisThread) {
    check.Check(node->block != NULL, "outer-loss owner block is unallocated");
    if (node->block != NULL) {
      const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
      double v[3] = {0.20, 0.0, 0.0};
      PIC::ParticleBuffer::SetI(_H2O_SPEC_, particle);
      PIC::ParticleBuffer::SetX(x, particle);
      PIC::ParticleBuffer::SetV(v, particle);
      ++local_created;
      const int result = PIC::Mover::Simple(particle, 1.0, node);
      check.Check(result == _PARTICLE_LEFT_THE_DOMAIN_,
                  "outward mover did not report domain loss");
      if (!PIC::ParticleBuffer::IsParticleAllocated(particle)) ++local_deleted;
    }
  }
  int created = 0, deleted = 0;
  MPI_Allreduce(&local_created, &created, 1, MPI_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_deleted, &deleted, 1, MPI_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  check.Check(created == 1 && deleted == 1,
              "outer boundary did not delete exactly one global particle");
  check.receipt = "created=1 deleted=1";
}

void CasePeriodicExchange(CheckContext& check) {
#if _PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_ON_
  PIC::BC::ExternalBoundary::Periodic::cBlockPairTable* selected = NULL;
  for (int index = 0;
       index < PIC::BC::ExternalBoundary::Periodic::BlockPairTableLength;
       ++index) {
    PIC::BC::ExternalBoundary::Periodic::cBlockPairTable& pair =
        PIC::BC::ExternalBoundary::Periodic::BlockPairTable[index];
    if (pair.GhostBlock != NULL && pair.RealBlock != NULL &&
        pair.GhostBlock->Thread == PIC::ThisThread &&
        pair.GhostBlock->block != NULL && pair.RealBlock->block != NULL) {
      selected = &pair;
      break;
    }
  }
  check.Check(selected != NULL, "rank has no allocated periodic block pair");
  int local_seeded = 0;
  if (selected != NULL) {
    AddParticle(selected->GhostBlock, 5000.0 + PIC::ThisThread);
    local_seeded = 1;
  }
  int global_seeded = 0;
  MPI_Allreduce(&local_seeded, &global_seeded, 1, MPI_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  PIC::BC::ExternalBoundary::Periodic::ExchangeParticles();
  const long long local_allocated = PIC::ParticleBuffer::GetAllPartNum();
  long long global_allocated = 0;
  MPI_Allreduce(&local_allocated, &global_allocated, 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  check.Check(global_allocated == global_seeded,
              "periodic exchange duplicated or lost particles");
  if (selected != NULL)
    check.Check(CellHead(selected->GhostBlock) == -1,
                "periodic source ghost still contains its particle");
  std::ostringstream receipt;
  receipt << "periodic_particles=" << global_allocated
          << " pairs="
          << PIC::BC::ExternalBoundary::Periodic::BlockPairTableLength;
  check.receipt = receipt.str();
#else
  check.Check(false, "periodic case requires the periodic compile profile");
#endif
}

void CaseAffinityDiagnostics(CheckContext& check) {
  // Do not mutate process affinity in a portable test: schedulers and cgroups
  // own that policy.  Validate the production parser and per-rank discovery
  // used to construct an affinity request, including a noncontiguous set.
  check.Check(PIC::Parallel::CountCpusInList("4") == 1,
              "single-CPU affinity list count is incorrect");
  check.Check(PIC::Parallel::CountCpusInList("0-7") == 8,
              "contiguous affinity range count is incorrect");
  check.Check(PIC::Parallel::CountCpusInList("0-3,8-11") == 8,
              "noncontiguous affinity list count is incorrect");
  check.Check(PIC::Parallel::CountCpusInList("7-3") == 0,
              "reversed affinity range was not rejected by the counter");
  const std::string online = PIC::Parallel::GetNodeOnlineCpuList();
  const int local_count = PIC::Parallel::CountCpusInList(online);
  check.Check(!online.empty() && local_count > 0,
              "online CPU discovery returned no usable CPU");
  int minimum = 0, maximum = 0;
  MPI_Allreduce(&local_count, &minimum, 1, MPI_INT, MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_count, &maximum, 1, MPI_INT, MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  check.Check(minimum > 0 && maximum >= minimum,
              "MPI ranks returned invalid online CPU counts");
  std::ostringstream receipt;
  receipt << "online_cpu_count_min=" << minimum
          << " online_cpu_count_max=" << maximum;
  check.receipt = receipt.str();
}

void SaveRestart(const std::string& path) {
  if (PIC::ThisThread == 0) std::remove(path.c_str());
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  PIC::Restart::SaveParticleData(path.c_str());
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

void FinishRestart(const std::string& path) {
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread == 0) std::remove(path.c_str());
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

void CaseRestart(CheckContext& check, bool repartition) {
  SeedOwnerParticles();
  const ParticleDigest expected = GlobalDigest();
  std::ostringstream name;
  name << "t10-restart-" << (repartition ? "changed" : "same")
       << "-n" << PIC::nTotalThreads << ".amps";
  SaveRestart(name.str());
  int changed = 0;
  if (repartition) {
    const std::vector<std::pair<std::uint64_t, int> > before = OwnerMap(check);
    PIC::Mesh::mesh->SetParallelLoadMeasure(SkewLoad);
    PIC::Mesh::mesh->CreateNewParallelDistributionLists();
    PIC::DomainBlockDecomposition::UpdateBlockTable();
    const std::vector<std::pair<std::uint64_t, int> > after = OwnerMap(check);
    for (std::size_t i = 0; i < before.size(); ++i)
      if (before[i].second != after[i].second) ++changed;
    check.Check(PIC::nTotalThreads == 1 || changed > 0,
                "restart changed-decomposition fixture did not change owners");
  }
  ClearOwnerParticles();
  check.Check(GlobalDigest().count == 0, "particle clear before restart failed");
  PIC::Restart::ReadParticleData(name.str().c_str());
  const ParticleDigest restored = GlobalDigest();
  check.Check(restored.count == expected.count,
              "restart changed global particle count");
  check.Check(restored.hash == expected.hash,
              "restart changed particle identity checksum");
  check.Check(restored.velocity_sum == expected.velocity_sum,
              "restart changed global velocity receipt");
  FinishRestart(name.str());
  std::ostringstream receipt;
  receipt << "particles=" << restored.count << " hash=" << restored.hash
          << " changed_owners=" << changed;
  check.receipt = receipt.str();
}

void CaseRestartSame(CheckContext& check) { CaseRestart(check, false); }
void CaseRestartChanged(CheckContext& check) { CaseRestart(check, true); }

void CaseRankInvariant(CheckContext& check) {
  SeedOwnerParticles();
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_] = 2.0;
  const double mass = PIC::MolecularData::GetMass(_H2O_SPEC_);
  long long local_count = 0;
  double local_momentum[3] = {0.0, 0.0, 0.0};
  double local_energy = 0.0;
  double local_weight = 0.0, local_mass = 0.0;
  unsigned long long local_hash = 0;
  double local_rate = 0.0;
  for (Node* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node != NULL; node = node->nextNodeThisThread) {
    for (long int particle = CellHead(node); particle != -1;
         particle = PIC::ParticleBuffer::GetNext(particle)) {
      ++local_count;
      local_weight += 2.0;
      local_mass += 2.0 * mass;
      const double* velocity = PIC::ParticleBuffer::GetV(particle);
      double speed2 = 0.0;
      for (int d = 0; d < 3; ++d) {
        local_momentum[d] += 2.0 * mass * velocity[d];
        speed2 += velocity[d] * velocity[d];
      }
      local_energy += mass * speed2;  // 0.5 * weight(=2) * m * v^2
      local_hash ^= static_cast<unsigned long long>(NodeKey(node));
    }
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          double x[3];
          CellCenterPosition(node, i, j, k, x);
          local_rate += 1.0 + x[0] + 2.0 * x[1] + 4.0 * x[2];
        }
  }
  long long count = 0;
  double momentum[3] = {0.0, 0.0, 0.0}, energy = 0.0, rate = 0.0;
  double weight = 0.0, total_mass = 0.0;
  unsigned long long hash = 0;
  MPI_Allreduce(&local_count, &count, 1, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local_momentum, momentum, 3, MPI_DOUBLE, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_energy, &energy, 1, MPI_DOUBLE, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_weight, &weight, 1, MPI_DOUBLE, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_rate, &rate, 1, MPI_DOUBLE, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&local_hash, &hash, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                MPI_GLOBAL_COMMUNICATOR);

  const std::vector<Node*> leaves = Leaves();
  long long expected_count = static_cast<long long>(leaves.size());
  double expected_momentum[3] = {0.0, 0.0, 0.0};
  double expected_energy = 0.0, expected_rate = 0.0;
  unsigned long long expected_hash = 0;
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    const double tag = static_cast<double>(1 + NodeKey(leaves[index]) % UINT64_C(1000000));
    const double velocity[3] = {tag, 0.5 * tag, -0.25 * tag};
    double speed2 = 0.0;
    for (int d = 0; d < 3; ++d) {
      expected_momentum[d] += 2.0 * mass * velocity[d];
      speed2 += velocity[d] * velocity[d];
    }
    expected_energy += mass * speed2;
    expected_hash ^= static_cast<unsigned long long>(NodeKey(leaves[index]));
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          double x[3];
          CellCenterPosition(leaves[index], i, j, k, x);
          expected_rate += 1.0 + x[0] + 2.0 * x[1] + 4.0 * x[2];
        }
  }
  const double relative_tolerance = 2.0e-14;
  check.Check(count == expected_count, "rank-invariant global count differs");
  check.Check(weight == 2.0 * expected_count,
              "rank-invariant statistical weight differs");
  check.Check(std::fabs(total_mass - 2.0 * expected_count * mass) <=
                  relative_tolerance * std::max(
                      std::fabs(2.0 * expected_count * mass), 1.0e-300),
              "rank-invariant physical mass differs");
  check.Check(hash == expected_hash, "rank-invariant XOR checksum differs");
  for (int d = 0; d < 3; ++d)
    check.Check(std::fabs(momentum[d] - expected_momentum[d]) <=
                    relative_tolerance * std::max(
                        std::fabs(expected_momentum[d]), 1.0e-300),
                "rank-invariant momentum differs in component %d", d);
  check.Check(std::fabs(energy - expected_energy) <=
                  relative_tolerance * std::max(
                      std::fabs(expected_energy), 1.0e-300),
              "rank-invariant kinetic energy differs");
  check.Check(std::fabs(rate - expected_rate) <=
                  relative_tolerance * std::max(std::fabs(expected_rate), 1.0),
              "rank-invariant sampled rate differs");
  std::ostringstream receipt;
  // Twelve significant digits are sufficient to expose physical changes but
  // intentionally suppress harmless last-bit differences from MPI reduction
  // tree ordering across 1/2/4 ranks.
  receipt.setf(std::ios::scientific);
  receipt.precision(12);
  receipt << "count=" << count << " hash=" << hash
          << " weight=" << weight << " mass=" << total_mass
          << " px=" << momentum[0] << " py=" << momentum[1]
          << " pz=" << momentum[2] << " energy=" << energy
          << " sampled_rate=" << rate;
  check.receipt = receipt.str();
}

typedef void (*CaseFunction)(CheckContext&);

struct CaseDescriptor {
  const char* id;
  CaseFunction function;
  bool periodic;
};

const CaseDescriptor cases[] = {
    {"MESH-MPI-OWNER-MAP-001", CaseOwnerMap, false},
    {"MESH-MPI-OWNER-GHOST-001", CaseOwnerGhost, false},
    {"MESH-MPI-LOAD-REPARTITION-001", CaseOwnerRepartition, false},
    {"MESH-MPI-CUTCELL-CONSISTENCY-001", CaseCutCell, false},
    {"PIC-MPI-HALO-REPLACE-001", CaseHaloReplace, false},
    {"PIC-MPI-HALO-ADD-001", CaseHaloAdd, false},
    {"PIC-MPI-SEND-MASK-001", CaseSendMask, false},
    {"PIC-MPI-BLOCK-PACK-001", CaseBlockPack, false},
    {"PIC-MPI-PARTICLE-EXCHANGE-001", CaseParticleExchange, false},
    {"PIC-MPI-PARTICLE-REPARTITION-001", CaseParticleRepartition, false},
    {"PIC-MPI-OUTER-LOSS-001", CaseOuterLoss, false},
    {"PIC-MPI-PERIODIC-EXCHANGE-001", CasePeriodicExchange, true},
    {"PIC-MPI-AFFINITY-DIAGNOSTICS-001", CaseAffinityDiagnostics, false},
    {"PIC-MPI-RESTART-SAME-001", CaseRestartSame, false},
    {"PIC-MPI-RESTART-REPARTITION-001", CaseRestartChanged, false},
    {"PIC-MPI-RANK-INVARIANCE-001", CaseRankInvariant, false},
};

const CaseDescriptor* FindCase(const std::string& id) {
  for (std::size_t index = 0; index < sizeof(cases) / sizeof(cases[0]); ++index)
    if (id == cases[index].id) return cases + index;
  return NULL;
}

}  // namespace

int main(int argc, char** argv) {
  std::string case_id;
  bool list_only = false;
  for (int index = 1; index < argc; ++index) {
    if (std::strcmp(argv[index], "--list") == 0) {
      list_only = true;
    } else if (std::strcmp(argv[index], "--case") == 0 && index + 1 < argc) {
      case_id = argv[++index];
    } else {
      std::fprintf(stderr, "usage: %s [--list | --case ID]\n", argv[0]);
      return 2;
    }
  }

  if (list_only) {
    for (std::size_t index = 0; index < sizeof(cases) / sizeof(cases[0]); ++index)
      std::printf("%s\n", cases[index].id);
    return 0;
  }
  const CaseDescriptor* selected = FindCase(case_id);
  if (selected == NULL) {
    std::fprintf(stderr, "unknown or missing T10 case: %s\n", case_id.c_str());
    return 2;
  }

  PIC::Init_BeforeParser();
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestCenterData);
  PIC::IndividualModelSampling::RequestStaticCellCornerData->push_back(
      RequestCornerData);
  PIC::Mesh::initCellSamplingDataBuffer();
  PIC::ParticleBuffer::Init(200000);

  CheckContext check;
#if _PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_ON_
  check.Check(selected->periodic,
              "nonperiodic case selected in periodic compile profile");
#else
  check.Check(!selected->periodic,
              "periodic case selected in nonperiodic compile profile");
#endif
  InitializeMesh(check);
  selected->function(check);

  int local_failure = check.ok ? 0 : 1;
  int global_failure = 0;
  MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  std::printf("[T10][case=%s][rank=%d/%d] %s %s\n", selected->id,
              PIC::ThisThread, PIC::nTotalThreads,
              check.ok ? "PASS" : "FAIL", check.Message().c_str());
  std::fflush(stdout);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread == 0)
    std::printf("[%s] %s %s\n", selected->id,
                global_failure == 0 ? "PASS" : "FAIL",
                check.receipt.empty() ? "collective case complete" :
                                        check.receipt.c_str());
  MPI_Finalize();
  return global_failure;
}
