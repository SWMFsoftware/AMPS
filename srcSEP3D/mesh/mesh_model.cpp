#include "mesh_model.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>

namespace SEP3D {
namespace Mesh {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

double Component(const Core::Vec3& value, int axis) {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

void SetComponent(Core::Vec3* value, int axis, double component) {
  if (axis == 0) value->x = component;
  else if (axis == 1) value->y = component;
  else value->z = component;
}

double WrapAngle(double value) {
  const double period = 2.0 * Core::Const::kPi;
  value = std::fmod(value, period);
  return value < 0.0 ? value + period : value;
}

double Clamp(double value, double lower, double upper) {
  return std::max(lower, std::min(value, upper));
}

double BlockSide(const LeafBlock& block) {
  return block.maximumM.x - block.minimumM.x;
}

Core::Vec3 BlockCenter(const LeafBlock& block) {
  return 0.5 * (block.minimumM + block.maximumM);
}

double RequestedInBlock(const LeafBlock& block,
                        const ResolutionConfiguration& configuration) {
  // The centre and eight corners are sufficient for this monotone radial law
  // and conservatively sample a curved tube crossing a block.  The result is
  // the smallest requested cell size, so an unsampled minimum can only delay
  // refinement by one level and is caught by the balance shoulder tests.
  double requested = RequestedCellSizeM(BlockCenter(block), configuration);
  for (int child = 0; child < 8; ++child) {
    Core::Vec3 point;
    SetComponent(&point, 0, (child & 1) ? block.maximumM.x : block.minimumM.x);
    SetComponent(&point, 1, (child & 2) ? block.maximumM.y : block.minimumM.y);
    SetComponent(&point, 2, (child & 4) ? block.maximumM.z : block.minimumM.z);
    requested = std::min(requested,
                         RequestedCellSizeM(point, configuration));
  }
  return requested;
}

void Split(const LeafBlock& parent, std::vector<LeafBlock>* children) {
  const Core::Vec3 middle = BlockCenter(parent);
  for (int child = 0; child < 8; ++child) {
    LeafBlock node;
    node.level = parent.level + 1;
    node.path = (parent.path << 3) | static_cast<std::uint64_t>(child);
    for (int axis = 0; axis < 3; ++axis) {
      const bool upper = (child & (1 << axis)) != 0;
      SetComponent(&node.minimumM, axis,
                   upper ? Component(middle, axis)
                         : Component(parent.minimumM, axis));
      SetComponent(&node.maximumM, axis,
                   upper ? Component(parent.maximumM, axis)
                         : Component(middle, axis));
    }
    children->push_back(node);
  }
}

void RefineRecursive(const LeafBlock& block,
                     const ResolutionConfiguration& configuration,
                     std::vector<LeafBlock>* leaves) {
  const double cellSize = BlockSide(block) / configuration.cellsPerBlockEdge;
  if (block.level < configuration.maximumLevel &&
      cellSize > RequestedInBlock(block, configuration) * (1.0 + 1.0e-12)) {
    std::vector<LeafBlock> children;
    children.reserve(8);
    Split(block, &children);
    for (const LeafBlock& child : children) {
      RefineRecursive(child, configuration, leaves);
    }
  } else {
    leaves->push_back(block);
  }
}

bool IntervalsOverlap(double a0, double a1, double b0, double b1,
                      double tolerance) {
  return std::min(a1, b1) - std::max(a0, b0) > tolerance;
}

bool FaceNeighbours(const LeafBlock& left, const LeafBlock& right,
                    double tolerance) {
  for (int normal = 0; normal < 3; ++normal) {
    const bool touches =
        std::fabs(Component(left.maximumM, normal) -
                  Component(right.minimumM, normal)) <= tolerance ||
        std::fabs(Component(right.maximumM, normal) -
                  Component(left.minimumM, normal)) <= tolerance;
    if (!touches) continue;
    const int tangent0 = (normal + 1) % 3;
    const int tangent1 = (normal + 2) % 3;
    if (IntervalsOverlap(Component(left.minimumM, tangent0),
                         Component(left.maximumM, tangent0),
                         Component(right.minimumM, tangent0),
                         Component(right.maximumM, tangent0), tolerance) &&
        IntervalsOverlap(Component(left.minimumM, tangent1),
                         Component(left.maximumM, tangent1),
                         Component(right.minimumM, tangent1),
                         Component(right.maximumM, tangent1), tolerance)) {
      return true;
    }
  }
  return false;
}

bool SpatialLess(const LeafBlock& left, const LeafBlock& right) {
  if (left.minimumM.z != right.minimumM.z)
    return left.minimumM.z < right.minimumM.z;
  if (left.minimumM.y != right.minimumM.y)
    return left.minimumM.y < right.minimumM.y;
  if (left.minimumM.x != right.minimumM.x)
    return left.minimumM.x < right.minimumM.x;
  return left.level < right.level;
}

Core::Status Solve3(double matrix[3][3], double rhs[3], Core::Vec3* result) {
  if (result == nullptr) return Invalid("gradient output is null");
  double augmented[3][4];
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column)
      augmented[row][column] = matrix[row][column];
    augmented[row][3] = rhs[row];
  }
  for (int pivot = 0; pivot < 3; ++pivot) {
    int best = pivot;
    for (int row = pivot + 1; row < 3; ++row) {
      if (std::fabs(augmented[row][pivot]) >
          std::fabs(augmented[best][pivot])) best = row;
    }
    if (!std::isfinite(augmented[best][pivot]) ||
        std::fabs(augmented[best][pivot]) < 1.0e-30) {
      return Invalid("gradient stencil is rank deficient");
    }
    for (int column = pivot; column < 4; ++column)
      std::swap(augmented[pivot][column], augmented[best][column]);
    const double divisor = augmented[pivot][pivot];
    for (int column = pivot; column < 4; ++column)
      augmented[pivot][column] /= divisor;
    for (int row = 0; row < 3; ++row) {
      if (row == pivot) continue;
      const double factor = augmented[row][pivot];
      for (int column = pivot; column < 4; ++column)
        augmented[row][column] -= factor * augmented[pivot][column];
    }
  }
  *result = Core::Vec3(augmented[0][3], augmented[1][3], augmented[2][3]);
  return Core::Status::OK();
}

}  // namespace

Core::Status Validate(const ResolutionConfiguration& configuration) {
  const double values[] = {
      configuration.innerRadiusM, configuration.outerRadiusM,
      configuration.minimumCellSizeM, configuration.backgroundCellSizeM,
      configuration.tubeCoreRadiusM, configuration.tubeShoulderRadiusM,
      configuration.tubeCellSizeM, configuration.solarWindSpeedMPerS,
      configuration.solarRotationRateRadPerS};
  for (double value : values) {
    if (!std::isfinite(value)) return Invalid("mesh configuration contains a non-finite value");
  }
  if (configuration.innerRadiusM <= 0.0 ||
      configuration.outerRadiusM <= configuration.innerRadiusM) {
    return Invalid("mesh radii must be positive and ordered");
  }
  if (configuration.minimumCellSizeM <= 0.0 ||
      configuration.backgroundCellSizeM < configuration.minimumCellSizeM) {
    return Invalid("mesh cell sizes must be positive and ordered");
  }
  if (configuration.cellsPerBlockEdge == 0 ||
      configuration.maximumLevel > 19) {
    return Invalid("cellsPerBlockEdge must be positive and maximumLevel <= 19");
  }
  if (configuration.memoryBudgetBytes == 0) {
    return Invalid("mesh memory budget must be positive");
  }
  if (configuration.enableTubeRefinement) {
    if ((configuration.tubePolarity != 1 && configuration.tubePolarity != -1) ||
        configuration.tubeCoreRadiusM <= 0.0 ||
        configuration.tubeShoulderRadiusM < configuration.tubeCoreRadiusM ||
        configuration.tubeCellSizeM <= 0.0 ||
        configuration.solarWindSpeedMPerS <= 0.0 ||
        configuration.tubeColatitudeRad < 0.0 ||
        configuration.tubeColatitudeRad > Core::Const::kPi) {
      return Invalid("Parker tube parameters are outside their physical range");
    }
  }
  return Core::Status::OK();
}

DomainBounds MakeDomain(RuntimeModel::DomainPreset preset,
                        double innerRadiusM, double requestedOuterRadiusM) {
  DomainBounds result;
  result.innerRadiusM = innerRadiusM;
  const double presetOuter = preset == RuntimeModel::DomainPreset::Earth
                                 ? Core::Const::AU
                                 : 1.666 * Core::Const::AU;
  result.outerRadiusM = requestedOuterRadiusM > 0.0
                            ? requestedOuterRadiusM : presetOuter;
  result.minimumM = Core::Vec3(-result.outerRadiusM, -result.outerRadiusM,
                               -result.outerRadiusM);
  result.maximumM = Core::Vec3(result.outerRadiusM, result.outerRadiusM,
                               result.outerRadiusM);
  return result;
}

Core::Vec3 ParkerTubeDirection(
    double radiusM, const ResolutionConfiguration& configuration) {
  const double travel = std::max(0.0, radiusM - configuration.innerRadiusM);
  const double winding = configuration.tubePolarity *
      configuration.solarRotationRateRadPerS * travel /
      configuration.solarWindSpeedMPerS;
  const double phi = WrapAngle(configuration.tubeLongitudeRad - winding);
  const double sinTheta = std::sin(configuration.tubeColatitudeRad);
  return Core::Vec3(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
                    std::cos(configuration.tubeColatitudeRad));
}

double TubeDistanceM(const Core::Vec3& positionM,
                     const ResolutionConfiguration& configuration) {
  const double radius = positionM.Norm();
  if (radius <= 0.0 || !std::isfinite(radius)) {
    return std::numeric_limits<double>::infinity();
  }
  const Core::Vec3 direction = positionM.Normalized();
  const Core::Vec3 centreline = ParkerTubeDirection(radius, configuration);
  // atan2(|u x v|,u.v) retains first-order accuracy for nearly coincident
  // directions; acos(u.v) amplifies a one-ulp dot-product error to O(1e-8)
  // radians and would make an exactly generated centreline appear displaced.
  const double sine = direction.Cross(centreline).Norm();
  const double cosine = Clamp(direction.Dot(centreline), -1.0, 1.0);
  return radius * std::atan2(sine, cosine);
}

double RequestedCellSizeM(const Core::Vec3& positionM,
                          const ResolutionConfiguration& configuration) {
  const double radius = positionM.Norm();
  double requested = configuration.backgroundCellSizeM;
  if (configuration.enableRadialRefinement) {
    const double radial = configuration.minimumCellSizeM *
        std::max(radius, configuration.innerRadiusM) /
        configuration.innerRadiusM;
    requested = std::min(requested,
                         Clamp(radial, configuration.minimumCellSizeM,
                               configuration.backgroundCellSizeM));
  }
  if (configuration.enableTubeRefinement) {
    const double distance = TubeDistanceM(positionM, configuration);
    double tube = configuration.backgroundCellSizeM;
    if (distance <= configuration.tubeCoreRadiusM) {
      tube = configuration.tubeCellSizeM;
    } else if (distance < configuration.tubeShoulderRadiusM) {
      const double fraction =
          (distance - configuration.tubeCoreRadiusM) /
          (configuration.tubeShoulderRadiusM -
           configuration.tubeCoreRadiusM);
      tube = configuration.tubeCellSizeM + fraction *
          (configuration.backgroundCellSizeM -
           configuration.tubeCellSizeM);
    }
    requested = std::min(requested, tube);
  }
  return Clamp(requested, configuration.minimumCellSizeM,
               configuration.backgroundCellSizeM);
}

Core::Status StandaloneOctree::Build(
    const DomainBounds& domain,
    const ResolutionConfiguration& resolution,
    const RuntimeModel::StorageLayout& storage, int ownerRanks) {
  const Core::Status valid = Validate(resolution);
  if (!valid.ok()) return valid;
  if (ownerRanks <= 0) return Invalid("ownerRanks must be positive");
  const double sideX = domain.maximumM.x - domain.minimumM.x;
  const double sideY = domain.maximumM.y - domain.minimumM.y;
  const double sideZ = domain.maximumM.z - domain.minimumM.z;
  if (!(sideX > 0.0) || sideX != sideY || sideX != sideZ) {
    return Invalid("standalone octree requires a finite cubic domain");
  }

  std::vector<LeafBlock> candidate;
  LeafBlock root;
  root.minimumM = domain.minimumM;
  root.maximumM = domain.maximumM;
  RefineRecursive(root, resolution, &candidate);

  // Enforce the 2:1 face-neighbour rule.  Refining only the first offending
  // coarse leaf and restarting keeps the result deterministic.
  const double tolerance = sideX * 1.0e-13;
  bool changed = true;
  while (changed) {
    changed = false;
    for (std::size_t i = 0; i < candidate.size() && !changed; ++i) {
      for (std::size_t j = i + 1; j < candidate.size(); ++j) {
        if (!FaceNeighbours(candidate[i], candidate[j], tolerance)) continue;
        const int difference = static_cast<int>(candidate[i].level) -
                               static_cast<int>(candidate[j].level);
        if (std::abs(difference) <= 1) continue;
        const std::size_t coarse = difference < 0 ? i : j;
        if (candidate[coarse].level >= resolution.maximumLevel) {
          return Invalid("maximumLevel prevents required 2:1 balancing");
        }
        const LeafBlock parent = candidate[coarse];
        candidate.erase(candidate.begin() + coarse);
        std::vector<LeafBlock> children;
        Split(parent, &children);
        candidate.insert(candidate.end(), children.begin(), children.end());
        changed = true;
        break;
      }
    }
  }

  std::sort(candidate.begin(), candidate.end(), SpatialLess);
  const std::uint64_t cellsPerBlock =
      static_cast<std::uint64_t>(resolution.cellsPerBlockEdge) *
      resolution.cellsPerBlockEdge * resolution.cellsPerBlockEdge;
  MeshSummary next;
  next.leavesByLevel.assign(resolution.maximumLevel + 1, 0);
  for (std::size_t i = 0; i < candidate.size(); ++i) {
    candidate[i].globalLeaf = i;
    candidate[i].firstCell = static_cast<std::uint64_t>(i) * cellsPerBlock;
    candidate[i].ownerRank = static_cast<int>(i % ownerRanks);
    ++next.leavesByLevel[candidate[i].level];
  }
  next.leafCount = candidate.size();
  next.cellCount = next.leafCount * cellsPerBlock;
  next.estimatedBytes = EstimateMemoryBytes(next, resolution, storage);
  if (next.estimatedBytes > resolution.memoryBudgetBytes) {
    std::ostringstream message;
    message << "estimated mesh memory " << next.estimatedBytes
            << " exceeds budget " << resolution.memoryBudgetBytes;
    return Core::Status(Core::StatusCode::LayoutMismatch, message.str());
  }

  leaves_.swap(candidate);
  summary_ = next;
  return Core::Status::OK();
}

bool StandaloneOctree::IsBalanced() const {
  return AreLeavesBalanced(leaves_);
}

bool AreLeavesBalanced(const std::vector<LeafBlock>& leaves) {
  if (leaves.empty()) return true;
  double scale = 0.0;
  for (const LeafBlock& leaf : leaves) scale = std::max(scale, BlockSide(leaf));
  const double tolerance = scale * 1.0e-13;
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    for (std::size_t j = i + 1; j < leaves.size(); ++j) {
      if (FaceNeighbours(leaves[i], leaves[j], tolerance) &&
          std::abs(static_cast<int>(leaves[i].level) -
                   static_cast<int>(leaves[j].level)) > 1) return false;
    }
  }
  return true;
}

std::size_t EstimateMemoryBytes(
    const MeshSummary& topology,
    const ResolutionConfiguration& resolution,
    const RuntimeModel::StorageLayout& storage) {
  const std::size_t perCell = storage.cellAssociatedBytes +
                              storage.samplingBytesPerCell;
  const long double estimate =
      static_cast<long double>(topology.cellCount) * perCell +
      static_cast<long double>(topology.leafCount) *
          resolution.blockOverheadBytes;
  if (estimate > std::numeric_limits<std::size_t>::max())
    return std::numeric_limits<std::size_t>::max();
  return static_cast<std::size_t>(estimate);
}

Core::Status CellStorage::Allocate(
    const StandaloneOctree& mesh,
    const RuntimeModel::StorageLayout& layout) {
  const std::size_t bytesPerCell = layout.cellAssociatedBytes +
                                   layout.samplingBytesPerCell;
  if (bytesPerCell == 0) return Invalid("cell storage layout is empty");
  if (mesh.summary().cellCount >
      std::numeric_limits<std::size_t>::max() / bytesPerCell) {
    return Invalid("cell storage size overflows size_t");
  }
  std::vector<unsigned char> candidate(
      static_cast<std::size_t>(mesh.summary().cellCount) * bytesPerCell, 0);
  std::vector<int> owners(mesh.summary().cellCount, -1);
  for (const LeafBlock& leaf : mesh.leaves()) {
    const std::uint64_t cells = mesh.summary().cellCount /
                                mesh.summary().leafCount;
    for (std::uint64_t local = 0; local < cells; ++local)
      owners[leaf.firstCell + local] = leaf.ownerRank;
  }
  bytes_.swap(candidate);
  owners_.swap(owners);
  bytesPerCell_ = bytesPerCell;
  return Core::Status::OK();
}

Core::Status CellStorage::Write(std::uint64_t globalCell, int callerRank,
                                std::size_t relativeOffset,
                                const void* source, std::size_t bytes) {
  if (source == nullptr && bytes != 0) return Invalid("cell write source is null");
  if (globalCell >= owners_.size()) return Invalid("cell write ID is out of range");
  if (owners_[globalCell] != callerRank) {
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "only the deterministic owner may fill a cell");
  }
  if (relativeOffset > bytesPerCell_ || bytes > bytesPerCell_ - relativeOffset)
    return Invalid("cell write exceeds the frozen layout");
  std::memcpy(bytes_.data() + globalCell * bytesPerCell_ + relativeOffset,
              source, bytes);
  return Core::Status::OK();
}

Core::Status CellStorage::Read(std::uint64_t globalCell,
                               std::size_t relativeOffset,
                               void* destination, std::size_t bytes) const {
  if (destination == nullptr && bytes != 0) return Invalid("cell read destination is null");
  if (globalCell >= owners_.size()) return Invalid("cell read ID is out of range");
  if (relativeOffset > bytesPerCell_ || bytes > bytesPerCell_ - relativeOffset)
    return Invalid("cell read exceeds the frozen layout");
  std::memcpy(destination,
              bytes_.data() + globalCell * bytesPerCell_ + relativeOffset,
              bytes);
  return Core::Status::OK();
}

int CellStorage::Owner(std::uint64_t globalCell) const {
  return globalCell < owners_.size() ? owners_[globalCell] : -1;
}

Core::Status ReconstructScalarGradient(
    const Core::Vec3& center, double centerValue,
    const std::vector<Core::Vec3>& neighbourPositions,
    const std::vector<double>& neighbourValues, Core::Vec3* gradient) {
  if (gradient == nullptr) return Invalid("gradient output is null");
  if (neighbourPositions.size() != neighbourValues.size() ||
      neighbourPositions.size() < 3) {
    return Invalid("gradient stencil requires at least three paired neighbours");
  }
  double normal[3][3] = {};
  double rhs[3] = {};
  for (std::size_t n = 0; n < neighbourPositions.size(); ++n) {
    const Core::Vec3 displacement = neighbourPositions[n] - center;
    const double d[3] = {displacement.x, displacement.y, displacement.z};
    const double delta = neighbourValues[n] - centerValue;
    const double distance2 = displacement.NormSq();
    if (!(distance2 > 0.0) || !std::isfinite(delta))
      return Invalid("gradient stencil contains a duplicate or non-finite sample");
    const double weight = 1.0 / distance2;
    for (int i = 0; i < 3; ++i) {
      rhs[i] += weight * d[i] * delta;
      for (int j = 0; j < 3; ++j)
        normal[i][j] += weight * d[i] * d[j];
    }
  }
  return Solve3(normal, rhs, gradient);
}

Core::Status ReconstructVectorGradient(
    const Core::Vec3& center, const Core::Vec3& centerValue,
    const std::vector<Core::Vec3>& neighbourPositions,
    const std::vector<Core::Vec3>& neighbourValues,
    Core::Tensor3* gradient) {
  if (gradient == nullptr) return Invalid("vector-gradient output is null");
  if (neighbourPositions.size() != neighbourValues.size())
    return Invalid("vector-gradient position/value counts differ");
  const double centerComponents[3] = {centerValue.x, centerValue.y, centerValue.z};
  for (int component = 0; component < 3; ++component) {
    std::vector<double> values;
    values.reserve(neighbourValues.size());
    for (const Core::Vec3& value : neighbourValues)
      values.push_back(Component(value, component));
    Core::Vec3 row;
    const Core::Status status = ReconstructScalarGradient(
        center, centerComponents[component], neighbourPositions, values, &row);
    if (!status.ok()) return status;
    (*gradient)(component, 0) = row.x;
    (*gradient)(component, 1) = row.y;
    (*gradient)(component, 2) = row.z;
  }
  return Core::Status::OK();
}

}  // namespace Mesh
}  // namespace SEP3D
