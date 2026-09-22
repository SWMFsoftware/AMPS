#include "mesh_model.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
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

double Clamp(double value, double lower, double upper) {
  return std::max(lower, std::min(value, upper));
}

double ProfileFraction(double coordinate,
                       RuntimeModel::RefinementProfile profile,
                       double exponent) {
  const double x = Clamp(coordinate, 0.0, 1.0);
  switch (profile) {
    case RuntimeModel::RefinementProfile::Linear:
      return x;
    case RuntimeModel::RefinementProfile::PowerLaw:
      return std::pow(x, exponent);
    case RuntimeModel::RefinementProfile::Smoothstep: {
      const double smooth = x * x * (3.0 - 2.0 * x);
      return std::pow(smooth, exponent);
    }
  }
  return x;
}

Core::ParkerSpiralGeometry Geometry(
    const ResolutionConfiguration& configuration) {
  Core::ParkerSpiralGeometry geometry;
  geometry.sourceRadiusM = configuration.innerRadiusM;
  geometry.sourceLongitudeRad = configuration.tubeLongitudeRad;
  geometry.sourceColatitudeRad = configuration.tubeColatitudeRad;
  geometry.solarWindSpeedMPerS = configuration.solarWindSpeedMPerS;
  geometry.solarRotationRateRadPerS = configuration.solarRotationRateRadPerS;
  geometry.rotationAxis = configuration.rotationAxis;
  return geometry;
}

std::size_t SaturatingBytes(long double value) {
  if (!(value >= 0.0L) ||
      value > static_cast<long double>(std::numeric_limits<std::size_t>::max()))
    return std::numeric_limits<std::size_t>::max();
  return static_cast<std::size_t>(value);
}

double BlockSide(const LeafBlock& block) {
  return block.maximumM.x - block.minimumM.x;
}

Core::Vec3 BlockCenter(const LeafBlock& block) {
  return 0.5 * (block.minimumM + block.maximumM);
}

double RequestedInBlock(const LeafBlock& block,
                        const ResolutionConfiguration& configuration) {
  // Match the lattice that AMPS uses when it asks localResolution() whether a
  // block needs refinement. Centre-and-corner sampling is not conservative
  // for a thin curved tube: the Parker centreline can cross a coarse block
  // without passing close to any of those nine points.
  double requested = RequestedCellSizeM(BlockCenter(block), configuration);
  const unsigned n = configuration.cellsPerBlockEdge;
  for (unsigned k = 0; k <= n; ++k)
    for (unsigned j = 0; j <= n; ++j)
      for (unsigned i = 0; i <= n; ++i) {
        const double fx = static_cast<double>(i) / n;
        const double fy = static_cast<double>(j) / n;
        const double fz = static_cast<double>(k) / n;
        const Core::Vec3 point(
            block.minimumM.x + fx * (block.maximumM.x - block.minimumM.x),
            block.minimumM.y + fy * (block.maximumM.y - block.minimumM.y),
            block.minimumM.z + fz * (block.maximumM.z - block.minimumM.z));
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
      configuration.originM.x, configuration.originM.y,
      configuration.originM.z, configuration.parkerInitialPointM.x,
      configuration.parkerInitialPointM.y,
      configuration.parkerInitialPointM.z, configuration.parkerLengthM,
      configuration.innerRadiusM, configuration.outerRadiusM,
      configuration.minimumCellSizeM, configuration.backgroundCellSizeM,
      configuration.solarSurfaceCellSizeM,
      configuration.solarRefinementOuterRadiusM,
      configuration.solarRefinementExponent,
      configuration.tubeReferenceRadiusM,
      configuration.tubeRadiusAtReferenceM,
      configuration.tubeCellSizeM,
      configuration.tubeTransverseExponent,
      configuration.solarWindSpeedMPerS,
      configuration.solarRotationRateRadPerS,
      configuration.rotationAxis.x, configuration.rotationAxis.y,
      configuration.rotationAxis.z};
  for (double value : values) {
    if (!std::isfinite(value)) return Invalid("mesh configuration contains a non-finite value");
  }
  if (configuration.innerRadiusM <= 0.0 ||
      configuration.outerRadiusM <= configuration.innerRadiusM) {
    return Invalid("mesh radii must be positive and ordered");
  }
  if (configuration.minimumCellSizeM <= 0.0 ||
      configuration.backgroundCellSizeM < configuration.minimumCellSizeM ||
      configuration.solarSurfaceCellSizeM < configuration.minimumCellSizeM ||
      configuration.solarSurfaceCellSizeM >
          configuration.backgroundCellSizeM ||
      configuration.solarRefinementOuterRadiusM <=
          configuration.innerRadiusM ||
      configuration.solarRefinementExponent <= 0.0 ||
      configuration.tubeTransverseExponent <= 0.0) {
    return Invalid("mesh cell sizes must be positive and ordered");
  }
  if (configuration.cellsPerBlockEdge == 0 ||
      configuration.maximumLevel > 19) {
    return Invalid("cellsPerBlockEdge must be positive and maximumLevel <= 19");
  }
  if (configuration.parkerLengthM <= 0.0 ||
      configuration.parkerPointCount < 2 ||
      configuration.parkerPointCount > 10000000ULL ||
      std::fabs((configuration.parkerInitialPointM - configuration.originM).Norm() -
                configuration.innerRadiusM) >
          1.0e-10 * configuration.innerRadiusM) {
    return Invalid("finite Parker centreline definition is invalid");
  }
  const Core::Vec3 declaredSource = configuration.originM +
      Core::ParkerCurvePoint(configuration.innerRadiusM,
                             Geometry(configuration));
  if ((configuration.parkerInitialPointM - declaredSource).Norm() >
      1.0e-10 * configuration.innerRadiusM) {
    return Invalid("finite Parker initial point disagrees with tube source angles");
  }
  if (configuration.memoryBudgetBytes == 0) {
    return Invalid("mesh memory budget must be positive");
  }
  if (configuration.enableTubeRefinement) {
    if (configuration.tubeReferenceRadiusM <= configuration.innerRadiusM ||
        configuration.tubeRadiusAtReferenceM <= 0.0 ||
        configuration.tubeCellSizeM <= 0.0 ||
        configuration.tubeCellSizeM > configuration.backgroundCellSizeM ||
        configuration.solarWindSpeedMPerS <= 0.0 ||
        configuration.tubeColatitudeRad < 0.0 ||
        configuration.tubeColatitudeRad > Core::Const::kPi ||
        !Core::ValidateParkerGeometry(Geometry(configuration)).ok()) {
      return Invalid("Parker tube parameters are outside their physical range");
    }
  }
  return Core::Status::OK();
}

DomainBounds MakeDomain(
    const RuntimeModel::RunConfiguration3DOptions& configuration) {
  DomainBounds result;
  result.innerRadiusM = configuration.innerRadiusM;
  result.outerRadiusM = configuration.outerRadiusM;
  result.originM = configuration.coordinateOriginM;
  result.innerBoundary = configuration.innerBoundary;
  result.outerBoundary = configuration.outerBoundary;
  result.minimumM = result.originM -
      Core::Vec3(result.outerRadiusM, result.outerRadiusM, result.outerRadiusM);
  result.maximumM = result.originM +
      Core::Vec3(result.outerRadiusM, result.outerRadiusM, result.outerRadiusM);
  return result;
}

double TubeRadiusM(double radiusM,
                   const ResolutionConfiguration& configuration) {
  if (!std::isfinite(radiusM) || radiusM <= 0.0 ||
      configuration.tubeReferenceRadiusM <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (configuration.tubeRadiusMode ==
      RuntimeModel::TubeRadiusMode::PhysicalConstant) {
    return configuration.tubeRadiusAtReferenceM;
  }
  return configuration.tubeRadiusAtReferenceM * radiusM /
         configuration.tubeReferenceRadiusM;
}

Core::Vec3 ParkerTubeDirection(
    double radiusM, const ResolutionConfiguration& configuration) {
  return Core::ParkerCurvePoint(radiusM, Geometry(configuration)).Normalized();
}

Core::Vec3 ParkerTubeTangent(
    double radiusM, const ResolutionConfiguration& configuration) {
  return Core::ParkerCurveTangent(radiusM, Geometry(configuration));
}

double TubeDistanceM(const Core::Vec3& positionM,
                     const ResolutionConfiguration& configuration) {
  const Core::Vec3 relative = positionM - configuration.originM;
  const double radius = relative.Norm();
  if (radius <= 0.0 || !std::isfinite(radius)) {
    return std::numeric_limits<double>::infinity();
  }
  const Core::Vec3 direction = relative.Normalized();
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
  const double radius = (positionM - configuration.originM).Norm();
  double requested = configuration.backgroundCellSizeM;
  if (configuration.enableRadialRefinement) {
    const double fraction = (radius - configuration.innerRadiusM) /
        (configuration.solarRefinementOuterRadiusM -
         configuration.innerRadiusM);
    const double radial = configuration.solarSurfaceCellSizeM +
        ProfileFraction(fraction, configuration.solarRefinementProfile,
                        configuration.solarRefinementExponent) *
        (configuration.backgroundCellSizeM -
         configuration.solarSurfaceCellSizeM);
    requested = std::min(requested, Clamp(
        radial, configuration.solarSurfaceCellSizeM,
        configuration.backgroundCellSizeM));
  }
  if (configuration.enableTubeRefinement) {
    const double distance = TubeDistanceM(positionM, configuration);
    const double boundary = TubeRadiusM(radius, configuration);
    const double fraction = distance / boundary;
    const double tube = configuration.tubeCellSizeM +
        ProfileFraction(fraction, configuration.tubeTransverseProfile,
                        configuration.tubeTransverseExponent) *
        (configuration.backgroundCellSizeM -
         configuration.tubeCellSizeM);

    // AMPS evaluates this point function on a Cartesian lattice. For a curve
    // crossing one lattice cell, the nearest lattice point is at most half a
    // cell diagonal away. Requesting max(h_tube,2*d) at distance d therefore
    // guarantees that a sampled point asks for the next level until the
    // centreline reaches h_tube. This numerical capture envelope fixes
    // sub-cell aliasing without changing the declared physical tube profile.
    constexpr double kStrictRefinement = 2.0 * (1.0 - 1.0e-12);
    const double capture = std::max(
        configuration.tubeCellSizeM, kStrictRefinement * distance);
    requested = std::min(requested, std::min(tube, capture));
  }
  return Clamp(requested, configuration.minimumCellSizeM,
               configuration.backgroundCellSizeM);
}

Core::Status BuildParkerCenterline(
    const ResolutionConfiguration& configuration,
    std::vector<Core::Vec3>* points) {
  if (points == nullptr) return Invalid("Parker centreline output is null");
  const Core::Status valid = Validate(configuration);
  if (!valid.ok()) return valid;

  std::vector<Core::Vec3> candidate;
  candidate.reserve(static_cast<std::size_t>(configuration.parkerPointCount));
  Core::Vec3 point = configuration.parkerInitialPointM;
  candidate.push_back(point);
  const double step = configuration.parkerLengthM /
      static_cast<double>(configuration.parkerPointCount - 1);
  const Core::ParkerSpiralGeometry geometry = Geometry(configuration);
  for (std::uint64_t i = 1; i < configuration.parkerPointCount; ++i) {
    const Core::Vec3 relative = point - configuration.originM;
    const Core::Vec3 first = Core::ParkerLocalTangent(relative, geometry);
    if (first.Norm() == 0.0)
      return Invalid("Parker tangent vanished while sampling centreline");
    const Core::Vec3 midpoint = relative + 0.5 * step * first;
    const Core::Vec3 tangent = Core::ParkerLocalTangent(midpoint, geometry);
    if (tangent.Norm() == 0.0)
      return Invalid("Parker midpoint tangent vanished while sampling centreline");
    point += step * tangent;
    candidate.push_back(point);
  }
  points->swap(candidate);
  return Core::Status::OK();
}

Core::Status WriteParkerCenterlineTecplot(
    const ResolutionConfiguration& configuration, const std::string& path) {
  if (path.empty()) return Invalid("Parker centreline Tecplot path is empty");
  std::vector<Core::Vec3> points;
  const Core::Status built = BuildParkerCenterline(configuration, &points);
  if (!built.ok()) return built;
  std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good())
    return Invalid("cannot open Parker centreline Tecplot file '" + path + "'");
  output << "TITLE=\"srcSEP3D initialized Parker centreline\"\n"
         << "VARIABLES=\"arc_length_m\",\"x_m\",\"y_m\",\"z_m\","
            "\"heliocentric_radius_m\",\"requested_cell_size_m\"\n"
         << "ZONE T=\"parker-centreline\", I=" << points.size()
         << ", F=POINT\n" << std::scientific << std::setprecision(17);
  double arcLengthM = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i != 0) arcLengthM += (points[i] - points[i - 1]).Norm();
    output << arcLengthM << ' ' << points[i].x << ' ' << points[i].y << ' '
           << points[i].z << ' '
           << (points[i] - configuration.originM).Norm() << ' '
           << RequestedCellSizeM(points[i], configuration) << '\n';
  }
  output.flush();
  if (!output.good())
    return Invalid("failed while writing Parker centreline Tecplot file '" +
                   path + "'");
  return Core::Status::OK();
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
  return EstimateWholeRunMemory(topology, resolution, storage).totalBytes;
}

MemoryEstimate EstimateWholeRunMemory(
    const MeshSummary& topology,
    const ResolutionConfiguration& resolution,
    const RuntimeModel::StorageLayout& storage) {
  MemoryEstimate result;
  const long double cells = topology.cellCount;
  const long double leaves = topology.leafCount;
  const RuntimeModel::MemoryModelOptions& model = resolution.memoryModel;
  result.baseCellBytes = SaturatingBytes(cells * model.baseCellBytes);
  result.associatedDataBytes = SaturatingBytes(
      cells * storage.cellAssociatedBytes);
  result.samplingBytes = SaturatingBytes(
      cells * storage.samplingBytesPerCell);
  // A Cartesian cell contributes a bounded number of face/edge/corner nodes;
  // using one full node allocation per cell is deliberately conservative and
  // avoids claiming shared-node savings before the native AMPS calibration.
  result.nodeBytes = SaturatingBytes(cells * model.baseNodeBytes);
  result.blockBytes = SaturatingBytes(leaves *
      (model.blockStructureBytes + resolution.blockOverheadBytes));
  result.particleBytes = SaturatingBytes(
      cells * model.particlesPerCell * model.particleBytes);
  const long double communicationBase = leaves *
      model.communicationBytesPerBlock;
  const long double residentBase =
      static_cast<long double>(result.baseCellBytes) +
      result.associatedDataBytes + result.samplingBytes + result.nodeBytes +
      result.blockBytes + result.particleBytes;
  result.communicationAndHaloBytes = SaturatingBytes(
      communicationBase + model.haloFraction * residentBase);
  result.subtotalBytes = SaturatingBytes(
      residentBase + result.communicationAndHaloBytes);
  result.safetyMarginBytes = SaturatingBytes(
      model.safetyMarginFraction * result.subtotalBytes);
  result.totalBytes = SaturatingBytes(
      static_cast<long double>(result.subtotalBytes) +
      result.safetyMarginBytes);
  return result;
}

Core::Status BuildRefinementPreflight(
    const DomainBounds& domain,
    const ResolutionConfiguration& resolution,
    const RuntimeModel::StorageLayout& storage,
    RefinementPreflight* report) {
  if (report == nullptr) return Invalid("refinement preflight output is null");
  const Core::Status valid = Validate(resolution);
  if (!valid.ok()) return valid;
  RefinementPreflight candidate;
  candidate.minimumRequestedCellM = std::numeric_limits<double>::infinity();
  candidate.maximumRequestedCellM = 0.0;

  // The preflight samples all named limiting manifolds: radial axes, Parker
  // centreline, and the transverse tube boundary.  It does not allocate the
  // AMPS mesh and therefore remains safe for a CLI --dry-run summary.
  constexpr int kRadialSamples = 512;
  for (int i = 0; i <= kRadialSamples; ++i) {
    const double radius = resolution.innerRadiusM +
        (resolution.outerRadiusM - resolution.innerRadiusM) * i /
        static_cast<double>(kRadialSamples);
    const Core::Vec3 probes[] = {
        domain.originM + Core::Vec3(radius, 0.0, 0.0),
        domain.originM + Core::Vec3(-radius, 0.0, 0.0),
        domain.originM + Core::Vec3(0.0, radius, 0.0),
        domain.originM + Core::Vec3(0.0, 0.0, radius),
        domain.originM + radius * ParkerTubeDirection(radius, resolution)};
    for (const Core::Vec3& probe : probes) {
      const double cell = RequestedCellSizeM(probe, resolution);
      if (cell < candidate.minimumRequestedCellM) {
        candidate.minimumRequestedCellM = cell;
        candidate.minimumLocationM = probe;
      }
      if (cell > candidate.maximumRequestedCellM) {
        candidate.maximumRequestedCellM = cell;
        candidate.maximumLocationM = probe;
      }
    }
  }
  candidate.tubeRadiusAtReferenceM =
      TubeRadiusM(resolution.tubeReferenceRadiusM, resolution);

  // Estimate blocks level-by-level from the fraction of sample points that
  // request each level.  Native mesh-count gates compare this planning value
  // with the actual AMPS tree; it is not presented as an exact allocator.
  candidate.estimatedBlocksByLevel.assign(resolution.maximumLevel + 1, 0);
  const double rootCell = 2.0 * domain.outerRadiusM /
      resolution.cellsPerBlockEdge;
  constexpr int kVolumeSamplesPerAxis = 17;
  std::uint64_t sampleCounts[20] = {};
  std::uint64_t totalSamples = 0;
  for (int iz = 0; iz < kVolumeSamplesPerAxis; ++iz) {
    for (int iy = 0; iy < kVolumeSamplesPerAxis; ++iy) {
      for (int ix = 0; ix < kVolumeSamplesPerAxis; ++ix) {
        const auto coordinate = [&](int index) {
          return -domain.outerRadiusM + 2.0 * domain.outerRadiusM * index /
              static_cast<double>(kVolumeSamplesPerAxis - 1);
        };
        const Core::Vec3 point = domain.originM +
            Core::Vec3(coordinate(ix), coordinate(iy), coordinate(iz));
        const double requested = RequestedCellSizeM(point, resolution);
        unsigned level = 0;
        while (level < resolution.maximumLevel &&
               rootCell / std::pow(2.0, level) > requested) ++level;
        ++sampleCounts[level];
        ++totalSamples;
      }
    }
  }
  std::uint64_t leafEstimate = 0;
  for (unsigned level = 0; level <= resolution.maximumLevel; ++level) {
    const long double fullLevelBlocks = std::pow(8.0L, level);
    const std::uint64_t count = static_cast<std::uint64_t>(std::ceil(
        fullLevelBlocks * sampleCounts[level] /
        static_cast<long double>(totalSamples)));
    candidate.estimatedBlocksByLevel[level] = count;
    leafEstimate += count;
  }
  MeshSummary topology;
  topology.leafCount = std::max<std::uint64_t>(1, leafEstimate);
  topology.cellCount = topology.leafCount * resolution.cellsPerBlockEdge *
      resolution.cellsPerBlockEdge * resolution.cellsPerBlockEdge;
  candidate.memory = EstimateWholeRunMemory(topology, resolution, storage);
  if (candidate.memory.totalBytes > resolution.memoryBudgetBytes) {
    std::ostringstream message;
    message << "preflight memory " << candidate.memory.totalBytes
            << " exceeds budget " << resolution.memoryBudgetBytes;
    return Core::Status(Core::StatusCode::LayoutMismatch, message.str());
  }
  *report = candidate;
  return Core::Status::OK();
}

Core::Status ClassifyBoundaryCrossing(const Core::Vec3& previousM,
                                      const Core::Vec3& currentM,
                                      const DomainBounds& domain) {
  const double previousRadius = (previousM - domain.originM).Norm();
  const double currentRadius = (currentM - domain.originM).Norm();
  if (!std::isfinite(previousRadius) || !std::isfinite(currentRadius))
    return Invalid("boundary-crossing position is not finite");
  if (previousRadius >= domain.innerRadiusM &&
      currentRadius < domain.innerRadiusM) {
    return Core::Status(Core::StatusCode::InnerBoundary,
                        "particle crossed inward through the absorbing solar boundary");
  }
  if (previousRadius <= domain.outerRadiusM &&
      currentRadius > domain.outerRadiusM) {
    return Core::Status(Core::StatusCode::DomainExit,
        domain.outerBoundary == RuntimeModel::OuterBoundaryMode::ImportedCoverage
            ? "particle left imported SWMF coverage"
            : "particle escaped through the outer boundary");
  }
  return Core::Status::OK();
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
