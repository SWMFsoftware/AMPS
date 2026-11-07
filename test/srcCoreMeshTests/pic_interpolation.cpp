#include "mesh_test_fixture.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <tuple>
#include <vector>

namespace AMPS {
namespace CoreMeshTests {

#if _MESH_DIMENSION_ == 3
namespace {

typedef cAmpsMesh<PIC::Mesh::cDataCornerNode, PIC::Mesh::cDataCenterNode,
                  PIC::Mesh::cDataBlockAMR> Mesh;

double Affine(const double* x) {
  return 1.25 + 2.0 * x[0] - 3.0 * x[1] + 0.5 * x[2];
}

double EvaluateCornerAffine(
    const PIC::InterpolationRoutines::CornerBased::cStencil& stencil) {
  double value = 0.0;
  for (int i = 0; i < stencil.Length; ++i)
    value += stencil.Weight[i] * Affine(stencil.cell[i]->GetX());
  return value;
}

double WeightSum(
    const PIC::InterpolationRoutines::CornerBased::cStencil& stencil) {
  double sum = 0.0;
  for (int i = 0; i < stencil.Length; ++i) sum += stencil.Weight[i];
  return sum;
}

struct RowKey {
  Node* node;
  int i;
  int j;
  int k;
  double weight;
};

bool RowKeyLess(const RowKey& a, const RowKey& b) {
  if (a.node != b.node) return a.node < b.node;
  if (a.i != b.i) return a.i < b.i;
  if (a.j != b.j) return a.j < b.j;
  return a.k < b.k;
}

std::vector<RowKey> CanonicalRow(
    const PIC::InterpolationRoutines::cRowStencil& row) {
  std::vector<RowKey> result;
  for (int n = 0; n < row.Length; ++n) {
    RowKey key = {row.Element[n].node, row.Element[n].i, row.Element[n].j,
                  row.Element[n].k, row.Element[n].Weight};
    result.push_back(key);
  }
  std::sort(result.begin(), result.end(), RowKeyLess);
  return result;
}

double LocalResolution(double*) { return 1.0; }

double RowCoordinate(const PIC::InterpolationRoutines::cRowStencil& row,
                     int dimension) {
  const int cells[3] = {
      _BLOCK_CELLS_X_, _BLOCK_CELLS_Y_, _BLOCK_CELLS_Z_};
  double result = 0.0;
  for (int n = 0; n < row.Length; ++n) {
    const Node* node = row.Element[n].node;
    const int index[3] = {
        row.Element[n].i, row.Element[n].j, row.Element[n].k};
    const double dx =
        (node->xmax[dimension] - node->xmin[dimension]) / cells[dimension];
    result += row.Element[n].Weight *
        (node->xmin[dimension] + (index[dimension] + 0.5) * dx);
  }
  return result;
}

bool SameCanonicalRow(const std::vector<RowKey>& a,
                      const std::vector<RowKey>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t n = 0; n < a.size(); ++n)
    if (a[n].node != b[n].node || a[n].i != b[n].i ||
        a[n].j != b[n].j || a[n].k != b[n].k ||
        a[n].weight != b[n].weight)
      return false;
  return true;
}

// Production interpolation reaches the active AMR tree through the global
// PIC mesh pointer.  This narrow scope lets the 256-mask oracle exercise a
// private topology without changing ownership of the application fixture.
class ScopedActiveMesh {
 public:
  explicit ScopedActiveMesh(Mesh* active) : saved_(PIC::Mesh::mesh) {
    PIC::Mesh::mesh = active;
  }
  ~ScopedActiveMesh() { PIC::Mesh::mesh = saved_; }

 private:
  Mesh* saved_;
};

}  // namespace

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-UNIFORM"),
               "trilinear corner stencil exactly reproduces an affine field") {
  EnsureAllocatedMesh();
  Node* node = Root()->downNode[0]->downNode[0];
  double x[3];
  for (int d = 0; d < 3; ++d)
    x[d] = node->xmin[d] + 0.37 * (node->xmax[d] - node->xmin[d]);

  PIC::InterpolationRoutines::CornerBased::cStencil stencil;
  PIC::InterpolationRoutines::CornerBased::InitStencil(x, node, stencil);
  AMPS_EXPECT_EQ(stencil.Length, 8);
  AMPS_EXPECT_NEAR(WeightSum(stencil), 1.0, 32.0 * DBL_EPSILON);
  AMPS_EXPECT_NEAR(EvaluateCornerAffine(stencil), Affine(x), 256.0 * DBL_EPSILON);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-BOUNDARY"),
               "upper-face stencil is normalized and remains in the final cell") {
  EnsureAllocatedMesh();
  Node* node = Root()->downNode[0]->downNode[0];
  double x[3] = {node->xmax[0], node->xmax[1], node->xmax[2]};
  PIC::InterpolationRoutines::CornerBased::cStencil stencil;
  PIC::InterpolationRoutines::CornerBased::InitStencil(x, node, stencil);
  AMPS_EXPECT_EQ(stencil.Length, 8);
  AMPS_EXPECT_NEAR(WeightSum(stencil), 1.0, 64.0 * DBL_EPSILON);
  AMPS_EXPECT_NEAR(EvaluateCornerAffine(stencil), Affine(x), 2.0e-9);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-AMR"),
               "row stencil crosses a 2:1 interface with conservative weights") {
  EnsureAllocatedMesh();
  Node* fine = Root()->downNode[0]->downNode[1];
  const double dx = (fine->xmax[0] - fine->xmin[0]) / _BLOCK_CELLS_X_;
  double x[3] = {fine->xmax[0] - 0.05 * dx,
                 0.5 * (fine->xmin[1] + fine->xmax[1]),
                 0.5 * (fine->xmin[2] + fine->xmax[2])};
  PIC::InterpolationRoutines::cRowStencil row;
  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x, fine, row);
  AMPS_EXPECT_TRUE(row.Length > 1);
  AMPS_EXPECT_NEAR(row.GetWeightSum(), 1.0, 128.0 * DBL_EPSILON);

  bool coarse_support = false;
  bool fine_support = false;
  for (int n = 0; n < row.Length; ++n) {
    AMPS_EXPECT_TRUE(row.Element[n].Weight >= 0.0);
    coarse_support = coarse_support || row.Element[n].node->RefinmentLevel == 1;
    fine_support = fine_support || row.Element[n].node->RefinmentLevel == 2;
  }
  AMPS_EXPECT_TRUE(coarse_support);
  AMPS_EXPECT_TRUE(fine_support);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-LEGACY-MASK"),
               "all 256 legacy masks reproduce coordinates across AMR interfaces") {
  EnsureInterpolationRoutines();
  const int octants = 8;

  // Keep every private mesh alive until the sweep completes.  The production
  // interpolation cache keys geometry by mesh address, node address, and mesh
  // modification counter.  Constructing one automatic Mesh at the same stack
  // address on every iteration can make two different C++ object lifetimes
  // appear identical when their masks have the same number of splits.  That is
  // not a production mesh lifecycle, and it would turn this geometry oracle
  // into an allocator-address test.  Distinct live objects give each mask an
  // unambiguous mesh authority while retaining the normal cache-enabled path.
  std::vector<std::unique_ptr<Mesh> > mask_meshes;
  mask_meshes.reserve(256);

  for (unsigned mask = 0; mask < 256; ++mask) {
    mask_meshes.push_back(std::unique_ptr<Mesh>(new Mesh));
    Mesh& mesh = *mask_meshes.back();
    ScopedActiveMesh active(&mesh);
    double xmin[3] = {0.0, 0.0, 0.0};
    double xmax[3] = {1.0, 1.0, 1.0};

    // Row-only interpolation consumes logical (leaf,i,j,k) cells and does not
    // require block storage.  Keeping this topology-only makes the exhaustive
    // mask sweep bounded while still running the production nonuniform
    // interpolation algorithm, unlike the former leaf-count proxy.
    mesh.AllowBlockAllocation = false;
    mesh.init(xmin, xmax, LocalResolution);
    AMPS_EXPECT_TRUE(mesh.splitTreeNode(mesh.rootTree));
    for (int octant = 0; octant < octants; ++octant)
      if ((mask & (1u << octant)) != 0) {
        AMPS_EXPECT_TRUE(mesh.splitTreeNode(mesh.rootTree->downNode[octant]));
      }

    // InitStencil selects between same-level trilinear and coarse/fine
    // multiblock interpolation from each leaf's min/max neighboring refinement
    // levels.  buildMesh(), readMeshFile(), and the shared fixture all finalize
    // these derived values after topology changes; this topology-only private
    // mesh must do the same.  Without this call a coarse leaf next to a refined
    // octant retains the constructor sentinel/default and can be misclassified
    // as uniform, producing a biased ghost-cell coordinate rather than testing
    // the actual AMR interpolation closure.
    mesh.SetNodeNeibResolutionLevelLimit();

    // Adapt the legacy srcInterpolateAMRTest invariant: interpolating each
    // Cartesian coordinate must return that coordinate.  For every octant,
    // probe within one interpolation cell of each central face, each central
    // edge, and the central corner.  The source-side sign is selected by the
    // octant code; interface_mask selects the codimension (1, 2, or 3).
    const int cells[3] = {
        _BLOCK_CELLS_X_, _BLOCK_CELLS_Y_, _BLOCK_CELLS_Z_};
    for (int octant = 0; octant < octants; ++octant) {
      for (int interface_mask = 1; interface_mask < 8; ++interface_mask) {
        double x[3];
        for (int d = 0; d < 3; ++d) {
          const bool upper = (octant & (1 << d)) != 0;
          x[d] = upper ? 0.6875 : 0.3125;
          if ((interface_mask & (1 << d)) != 0) {
            const double coarse_cell = 0.5 / cells[d];
            x[d] = 0.5 + (upper ? 1.0 : -1.0) * 0.20 * coarse_cell;
          }
        }

        Node* source = mesh.findTreeNode(x);
        AMPS_EXPECT_TRUE(source != NULL);
        PIC::InterpolationRoutines::cRowStencil row;
        PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(
            x, source, row);
        AMPS_EXPECT_TRUE(row.Length > 0);
        AMPS_EXPECT_NEAR(row.GetWeightSum(), 1.0, 512.0 * DBL_EPSILON);
        for (int n = 0; n < row.Length; ++n) {
          AMPS_EXPECT_TRUE(row.Element[n].node != NULL);
          AMPS_EXPECT_TRUE(row.Element[n].node->lastBranchFlag() ==
                           _BOTTOM_BRANCH_TREE_);
          AMPS_EXPECT_TRUE(row.Element[n].Weight >= -64.0 * DBL_EPSILON);
        }
        for (int d = 0; d < 3; ++d) {
          const double reconstructed = RowCoordinate(row, d);
          if (std::fabs(reconstructed - x[d]) > 2.0e-10) {
            // An exhaustive mask failure must identify the exact AMR geometry;
            // a bare coordinate mismatch cannot distinguish a face, edge, or
            // corner closure, nor prove whether cache state changed the row.
            std::ostringstream message;
            message << "legacy-mask coordinate reproduction failed: mask="
                    << mask << " source_octant=" << octant
                    << " interface_mask=" << interface_mask
                    << " dimension=" << d << " source_level="
                    << static_cast<int>(source->RefinmentLevel)
                    << " neighbor_levels=["
                    << static_cast<int>(source->minNeibRefinmentLevel) << ","
                    << static_cast<int>(source->maxNeibRefinmentLevel) << "]"
                    << " x=" << x[d] << " reconstructed="
                    << reconstructed << " row_length=" << row.Length
                    << " row={";
            for (int element = 0; element < row.Length; ++element) {
              if (element != 0) message << ";";
              const Node* owner = row.Element[element].node;
              const int index[3] = {row.Element[element].i,
                                    row.Element[element].j,
                                    row.Element[element].k};
              message << "L" << static_cast<int>(owner->RefinmentLevel)
                      << " i=" << index[0] << "," << index[1] << ","
                      << index[2] << " w=" << row.Element[element].Weight
                      << " c=";
              for (int coordinate = 0; coordinate < 3; ++coordinate) {
                if (coordinate != 0) message << ",";
                const double dx =
                    (owner->xmax[coordinate] - owner->xmin[coordinate]) /
                    cells[coordinate];
                message << owner->xmin[coordinate] +
                    (index[coordinate] + 0.5) * dx;
              }
            }
            message << "}";
            throw AMPS::Test::Failure(message.str());
          }
        }
      }
    }
  }
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-DUPLICATE"),
               "row insertion merges duplicate logical cells and adds weights") {
  EnsureTopologyMesh();
  Node* node = Root()->downNode[1];
  PIC::InterpolationRoutines::cRowStencil row;
  row.AddCell(0.25, node, 0, 1, 0);
  row.AddCell(0.75, node, 0, 1, 0);
  AMPS_EXPECT_EQ(row.Length, 1);
  AMPS_EXPECT_NEAR(row.Element[0].Weight, 1.0, 0.0);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-PERMUTATION"),
               "row merge is invariant to source-stencil insertion order") {
  EnsureTopologyMesh();
  Node* a = Root()->downNode[1];
  Node* b = Root()->downNode[2];
  PIC::InterpolationRoutines::cRowStencil first;
  PIC::InterpolationRoutines::cRowStencil second;
  first.AddCell(0.2, a, 0, 0, 0);
  first.AddCell(0.3, b, 1, 0, 0);
  first.AddCell(0.5, a, 0, 0, 0);
  second.AddCell(0.3, b, 1, 0, 0);
  second.AddCell(0.5, a, 0, 0, 0);
  second.AddCell(0.2, a, 0, 0, 0);

  const std::vector<RowKey> x = CanonicalRow(first);
  const std::vector<RowKey> y = CanonicalRow(second);
  AMPS_EXPECT_EQ(x.size(), y.size());
  for (std::size_t n = 0; n < x.size(); ++n) {
    AMPS_EXPECT_TRUE(x[n].node == y[n].node);
    AMPS_EXPECT_EQ(x[n].i, y[n].i);
    AMPS_EXPECT_EQ(x[n].j, y[n].j);
    AMPS_EXPECT_EQ(x[n].k, y[n].k);
    AMPS_EXPECT_NEAR(x[n].weight, y[n].weight, DBL_EPSILON);
  }
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-CACHE"),
               "cached row replay invalidates after adjacent AMR refinement") {
  // The cache has no public hit/miss counter.  Force its documented enabled
  // setting before the first interpolation in this isolated test process,
  // then test behavior which a stale entry cannot satisfy.
  ::setenv("AMPS_INTERPOLATION_STENCIL_CACHE", "1", 1);
  EnsureTopologyMesh();
  EnsureInterpolationRoutines();

  Node* node = Root()->downNode[0]->downNode[1];
  Node* coarse = Root()->downNode[1];
  const double dx = (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_;
  double x[3] = {node->xmax[0] - 0.05 * dx,
                 0.5 * (node->xmin[1] + node->xmax[1]),
                 0.5 * (node->xmin[2] + node->xmax[2])};
  PIC::InterpolationRoutines::cRowStencil first;
  PIC::InterpolationRoutines::cRowStencil replay;
  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x, node, first);
  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x, node, replay);
  const std::vector<RowKey> a = CanonicalRow(first);
  const std::vector<RowKey> b = CanonicalRow(replay);
  AMPS_EXPECT_TRUE(SameCanonicalRow(a, b));

  bool first_uses_coarse = false;
  for (int n = 0; n < first.Length; ++n)
    first_uses_coarse = first_uses_coarse || first.Element[n].node == coarse;
  AMPS_EXPECT_TRUE(first_uses_coarse);

  // Keep the query and its fine source-node pointer unchanged while refining
  // only the adjacent coarse support.  The mesh generation must change, and
  // the cache key's generation component must cause a new row to be built.
  const long generation = PIC::Mesh::mesh->GetMeshID();
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->splitTreeNode(coarse));
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->GetMeshID() > generation);
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->findTreeNode(x) == node);

  PIC::InterpolationRoutines::cRowStencil rebuilt;
  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x, node,
                                                                 rebuilt);
  const std::vector<RowKey> c = CanonicalRow(rebuilt);
  bool new_fine_support = false;
  for (int n = 0; n < rebuilt.Length; ++n) {
    Node* support = rebuilt.Element[n].node;
    AMPS_EXPECT_TRUE(support != NULL);
    AMPS_EXPECT_TRUE(support->lastBranchFlag() == _BOTTOM_BRANCH_TREE_);
    AMPS_EXPECT_TRUE(support != coarse);
    new_fine_support = new_fine_support || support->upNode == coarse;
  }
  AMPS_EXPECT_TRUE(new_fine_support);
  AMPS_EXPECT_FALSE(SameCanonicalRow(a, c));
  AMPS_EXPECT_NEAR(rebuilt.GetWeightSum(), 1.0, 256.0 * DBL_EPSILON);
}

AMPS_TEST_CASE(AMPS_MESH_CASE_ID("INTERPOLATION-INPUT-IMMUTABLE"),
               "corner interpolation never modifies an upper-boundary input") {
  EnsureAllocatedMesh();
  Node* node = Root()->downNode[0]->downNode[0];
  double x[3] = {node->xmax[0], node->xmax[1], node->xmax[2]};
  const double original[3] = {x[0], x[1], x[2]};
  PIC::InterpolationRoutines::CornerBased::cStencil stencil;
  PIC::InterpolationRoutines::CornerBased::InitStencil(x, node, stencil);
  for (int d = 0; d < 3; ++d) AMPS_EXPECT_NEAR(x[d], original[d], 0.0);
}

#endif  // _MESH_DIMENSION_ == 3

void LinkPICInterpolationTests() {}

}  // namespace CoreMeshTests
}  // namespace AMPS
