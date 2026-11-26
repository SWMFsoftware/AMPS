#include "mesh_model.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <queue>
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

// Exact closest-point sphere/AABB intersection. Unlike center/corner probes,
// this keeps a leaf whose face intersects a small solar neighbourhood.
bool IntersectsSolarNeighbourhood(const Core::Vec3& lower,
                                 const Core::Vec3& upper,
                                 const ResolutionConfiguration& c) {
  if (c.activeSolarSphereRadiusM <= 0.0) return false;
  double distanceSquared = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    const double origin = Component(c.originM, axis);
    const double nearest = Clamp(origin, Component(lower, axis),
                                 Component(upper, axis));
    distanceSquared += (nearest - origin) * (nearest - origin);
  }
  return distanceSquared <=
      c.activeSolarSphereRadiusM * c.activeSolarSphereRadiusM;
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
  if (configuration.enableRadialRefinement &&
      configuration.solarRefinementAnchor ==
          RuntimeModel::SolarRefinementAnchor::Photosphere) {
    // The Sun need not lie on the root's sampling lattice in a corner box.
    // The closest point gives the exact minimum radius in this block; the
    // monotone radial law therefore cannot miss a small photospheric target.
    const Core::Vec3 nearest(
        Clamp(configuration.originM.x, block.minimumM.x, block.maximumM.x),
        Clamp(configuration.originM.y, block.minimumM.y, block.maximumM.y),
        Clamp(configuration.originM.z, block.minimumM.z, block.maximumM.z));
    requested = std::min(requested, RequestedCellSizeM(nearest, configuration));
  }
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
      configuration.activeTubeReferenceRadiusM,
      configuration.activeTubeRadiusAtReferenceM,
      configuration.activeSolarSphereRadiusM,
      configuration.solarWindSpeedMPerS,
      configuration.solarRotationRateRadPerS,
      configuration.rotationAxis.x, configuration.rotationAxis.y,
      configuration.rotationAxis.z};
  for (double value : values) {
    if (!std::isfinite(value)) return Invalid("mesh configuration contains a non-finite value");
  }
  if (configuration.innerRadiusM < Core::Const::R_sun ||
      configuration.outerRadiusM <= configuration.innerRadiusM) {
    return Invalid(
        "mesh radii must be ordered and the Parker/transport inner radius "
        "must not lie below the physical solar surface");
  }
  if (configuration.minimumCellSizeM <= 0.0 ||
      configuration.backgroundCellSizeM < configuration.minimumCellSizeM ||
      configuration.solarSurfaceCellSizeM < configuration.minimumCellSizeM ||
      configuration.solarSurfaceCellSizeM >
          configuration.backgroundCellSizeM ||
      configuration.solarRefinementOuterRadiusM <=
          (configuration.solarRefinementAnchor ==
                  RuntimeModel::SolarRefinementAnchor::Photosphere
              ? Core::Const::R_sun : configuration.innerRadiusM) ||
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
  if (configuration.solarRefinementAnchor !=
          RuntimeModel::SolarRefinementAnchor::SourceShell &&
      configuration.solarRefinementAnchor !=
          RuntimeModel::SolarRefinementAnchor::Photosphere)
    return Invalid("unknown solar refinement anchor");
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
  // An observer flux tube is traced after MPI and the background provider
  // exist, i.e. after this first validation; BuildActiveRegionPlan requires
  // the traced geometry.
  if (configuration.activeRegion ==
      RuntimeModel::ActiveRegionMode::ParkerTube) {
    if (configuration.activeSolarSphereRadiusM < 0.0 ||
        (configuration.activeSolarSphereRadiusM > 0.0 &&
         (configuration.activeSolarSphereRadiusM < configuration.innerRadiusM ||
          configuration.activeSolarSphereRadiusM >= configuration.outerRadiusM)))
      return Invalid("solar neighbourhood must join the finite corridor and lie inside the outer radius");
    if (configuration.activeTubeReferenceRadiusM <= 0.0 ||
        configuration.activeTubeRadiusAtReferenceM <= 0.0 ||
        configuration.solarWindSpeedMPerS <= 0.0 ||
        configuration.tubeColatitudeRad < 0.0 ||
        configuration.tubeColatitudeRad > Core::Const::kPi ||
        !Core::ValidateParkerGeometry(Geometry(configuration)).ok()) {
      return Invalid(
          "active Parker-corridor parameters are outside their physical "
          "range");
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
  // The immutable factory has already validated these bounds. The same
  // resolver is used there for refinement feasibility, here for AMPS init,
  // and by the CLI; a translated corner root cannot use a centered estimate.
  const Core::Status status = RuntimeModel::ResolveDomainBoundsM(
      configuration, &result.minimumM, &result.maximumM);
  if (!status.ok()) {
    const double invalid = std::numeric_limits<double>::quiet_NaN();
    result.minimumM = result.maximumM = Core::Vec3(invalid, invalid, invalid);
  }
  return result;
}

SolarBoundaryGeometry MakeSolarBoundary(
    const RuntimeModel::RunConfiguration3DOptions& configuration) {
  SolarBoundaryGeometry result;
  result.centerM = configuration.coordinateOriginM;
  result.radiusM = Core::Const::R_sun;
  return result;
}

bool AxisAlignedBoxEntirelyInsideSolarBoundary(
    const Core::Vec3& minimumM, const Core::Vec3& maximumM,
    const SolarBoundaryGeometry& boundary) {
  if (!std::isfinite(boundary.centerM.x) ||
      !std::isfinite(boundary.centerM.y) ||
      !std::isfinite(boundary.centerM.z) ||
      !std::isfinite(boundary.radiusM) || boundary.radiusM <= 0.0) {
    return false;
  }

  // The squared distance from the sphere centre is separable by Cartesian
  // axis.  Selecting the farther endpoint on every axis therefore identifies
  // a farthest corner without enumerating all eight corner combinations.
  double farthestSquaredM2 = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    const double lower = Component(minimumM, axis);
    const double upper = Component(maximumM, axis);
    const double center = Component(boundary.centerM, axis);
    if (!std::isfinite(lower) || !std::isfinite(upper) || lower > upper)
      return false;
    const double displacement = std::max(
        std::fabs(lower - center), std::fabs(upper - center));
    farthestSquaredM2 += displacement * displacement;
  }
  return std::isfinite(farthestSquaredM2) &&
      farthestSquaredM2 <= boundary.radiusM * boundary.radiusM;
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

double ActiveTubeRadiusM(double radiusM,
                         const ResolutionConfiguration& configuration) {
  if (!std::isfinite(radiusM) || radiusM <= 0.0 ||
      configuration.activeTubeReferenceRadiusM <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (configuration.activeTubeRadiusMode ==
      RuntimeModel::TubeRadiusMode::PhysicalConstant) {
    return configuration.activeTubeRadiusAtReferenceM;
  }
  return configuration.activeTubeRadiusAtReferenceM * radiusM /
         configuration.activeTubeReferenceRadiusM;
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
    const double anchorRadiusM = configuration.solarRefinementAnchor ==
        RuntimeModel::SolarRefinementAnchor::Photosphere
        ? Core::Const::R_sun : configuration.innerRadiusM;
    const double fraction = (radius - anchorRadiusM) /
        (configuration.solarRefinementOuterRadiusM -
         anchorRadiusM);
    const double radial = configuration.solarSurfaceCellSizeM +
        ProfileFraction(fraction, configuration.solarRefinementProfile,
                        configuration.solarRefinementExponent) *
        (configuration.backgroundCellSizeM -
         configuration.solarSurfaceCellSizeM);
    requested = std::min(requested, Clamp(
        radial, configuration.solarSurfaceCellSizeM,
        configuration.backgroundCellSizeM));
  }
  if (configuration.activeRegion ==
          RuntimeModel::ActiveRegionMode::ObserverFluxTube &&
      configuration.observerTube && configuration.observerTube->index) {
    // Observer flux tube: the radial target inside the tube, graded with
    // distance d outside it (requested <= 2d, the capture envelope below), so
    // the blocks the tube crosses refine to the radial target while the rest
    // of the cube stays coarse.
    //
    // h(x) = clamp(max(h_radial(x), 2 (1 - 1e-12) d(x)), h_min, R_outer)
    //   h_radial  the radial target computed above (at most the global
    //             background cell size);
    //   d         the distance to the nearest tube capsule (0 inside);
    //   R_outer   the domain's outer radius, the coarsest size requested.
    // The factor 2 (1 - 1e-12) is the capture envelope of the analytic
    // corridor below (kStrictRefinement and the comment there).  AMPS
    // samples this function on a lattice, and the nearest sample to a thin
    // tube can sit up to half a cell diagonal away; requesting about 2 d
    // makes that sample still ask for refinement, so a block crossed by the
    // tube is not missed by aliasing.  This branch returns before the
    // analytic-corridor branch (enableTubeRefinement must be off in this
    // mode), and it clamps at R_outer rather than at the background size:
    // far from the tube the blocks may stay coarser than the background
    // cell size.
    constexpr double kStrictCapture = 2.0 * (1.0 - 1.0e-12);
    // Outside the tube the cells grow with the distance d (<= 2d) up to the
    // domain scale, so the cube away from the tube stays at its coarsest
    // levels; inside the tube the radial target (capped at the global size)
    // applies.  The distance query is clipped where 2d reaches that scale.
    const double coarsestM = configuration.outerRadiusM;
    const double distance = DistanceToObserverTubeM(
        *configuration.observerTube, positionM, coarsestM / kStrictCapture);
    const double graded = std::max(requested, kStrictCapture * distance);
    return std::max(configuration.minimumCellSizeM,
                    std::min(coarsestM, graded));
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

namespace {

struct TubeSegment {
  Core::Vec3 firstM;
  Core::Vec3 secondM;
  // The envelope contains the complete curved centreline interval, not only
  // its chord.  See BuildFiniteTubeSegments for the rigorous half-arc bound.
  double envelopeRadiusM = 0.0;
};

bool ValidBox(const Core::Vec3& minimumM, const Core::Vec3& maximumM) {
  for (int axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(Component(minimumM, axis)) ||
        !std::isfinite(Component(maximumM, axis)) ||
        Component(maximumM, axis) <= Component(minimumM, axis)) {
      return false;
    }
  }
  return true;
}

double MinimumSide(const Core::Vec3& minimumM, const Core::Vec3& maximumM) {
  return std::min({maximumM.x - minimumM.x,
                   maximumM.y - minimumM.y,
                   maximumM.z - minimumM.z});
}

double BoxVolume(const LeafBlock& leaf) {
  return (leaf.maximumM.x - leaf.minimumM.x) *
         (leaf.maximumM.y - leaf.minimumM.y) *
         (leaf.maximumM.z - leaf.minimumM.z);
}

double SquaredDistanceToBox(const Core::Vec3& point,
                            const Core::Vec3& minimumM,
                            const Core::Vec3& maximumM) {
  double result = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    const double value = Component(point, axis);
    double difference = 0.0;
    if (value < Component(minimumM, axis))
      difference = value - Component(minimumM, axis);
    else if (value > Component(maximumM, axis))
      difference = value - Component(maximumM, axis);
    result += difference * difference;
  }
  return result;
}

double SquaredDistanceSegmentToBox(
    const Core::Vec3& firstM, const Core::Vec3& secondM,
    const Core::Vec3& minimumM, const Core::Vec3& maximumM) {
  // Distance from p(t)=p0+t*d to a box is a continuous piecewise quadratic.
  // Its active quadratic changes only when one coordinate crosses a box face.
  // Enumerating those at most six breakpoints and the stationary point in
  // each interval gives the exact segment/AABB distance without center/corner
  // sampling or an invalid Lipschitz assumption.
  const Core::Vec3 direction = secondM - firstM;
  std::vector<double> breakpoints;
  breakpoints.reserve(8);
  breakpoints.push_back(0.0);
  breakpoints.push_back(1.0);
  for (int axis = 0; axis < 3; ++axis) {
    const double delta = Component(direction, axis);
    if (delta == 0.0) continue;
    const double first =
        (Component(minimumM, axis) - Component(firstM, axis)) / delta;
    const double second =
        (Component(maximumM, axis) - Component(firstM, axis)) / delta;
    if (first > 0.0 && first < 1.0) breakpoints.push_back(first);
    if (second > 0.0 && second < 1.0) breakpoints.push_back(second);
  }
  std::sort(breakpoints.begin(), breakpoints.end());
  breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()),
                    breakpoints.end());

  auto pointAt = [&](double parameter) {
    return firstM + parameter * direction;
  };
  double best = std::min(SquaredDistanceToBox(firstM, minimumM, maximumM),
                         SquaredDistanceToBox(secondM, minimumM, maximumM));
  for (std::size_t interval = 0; interval + 1 < breakpoints.size();
       ++interval) {
    const double lower = breakpoints[interval];
    const double upper = breakpoints[interval + 1];
    const double middle = 0.5 * (lower + upper);
    const Core::Vec3 middlePoint = pointAt(middle);
    double quadratic = 0.0;
    double linear = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
      const double value = Component(middlePoint, axis);
      double boundary = value;
      if (value < Component(minimumM, axis))
        boundary = Component(minimumM, axis);
      else if (value > Component(maximumM, axis))
        boundary = Component(maximumM, axis);
      else
        continue;
      const double delta = Component(direction, axis);
      const double offset = Component(firstM, axis) - boundary;
      quadratic += delta * delta;
      linear += delta * offset;
    }
    if (quadratic > 0.0) {
      const double stationary = Clamp(-linear / quadratic, lower, upper);
      best = std::min(best, SquaredDistanceToBox(
          pointAt(stationary), minimumM, maximumM));
    }
    best = std::min(best, SquaredDistanceToBox(
        pointAt(lower), minimumM, maximumM));
    best = std::min(best, SquaredDistanceToBox(
        pointAt(upper), minimumM, maximumM));
  }
  return best;
}

Core::Status BuildFiniteTubeSegments(
    const ResolutionConfiguration& configuration, double spatialScaleM,
    std::vector<TubeSegment>* segments, double* effectiveLengthM) {
  if (segments == nullptr || effectiveLengthM == nullptr)
    return Invalid("active-region segment output is null");
  if (!std::isfinite(spatialScaleM) || spatialScaleM <= 0.0)
    return Invalid("active-region spatial scale must be positive");

  const Core::ParkerSpiralGeometry geometry = Geometry(configuration);
  const Core::Status geometryStatus = Core::ValidateParkerGeometry(geometry);
  if (!geometryStatus.ok()) return geometryStatus;
  const double outerLengthM = Core::ParkerCurveArcLengthM(
      configuration.outerRadiusM, geometry);
  if (!std::isfinite(outerLengthM) || outerLengthM <= 0.0)
    return Invalid("physical outer sphere has invalid Parker arc length");
  const double lengthM = std::min(configuration.parkerLengthM, outerLengthM);

  double terminalRadiusM = 0.0;
  const Core::Status terminalStatus = Core::ParkerCurveRadiusAtArcLengthM(
      lengthM, geometry, &terminalRadiusM);
  if (!terminalStatus.ok()) return terminalStatus;
  const double firstRadiusM = ActiveTubeRadiusM(
      configuration.innerRadiusM, configuration);
  const double lastRadiusM = ActiveTubeRadiusM(
      terminalRadiusM, configuration);
  const double narrowestTubeM = std::min(firstRadiusM, lastRadiusM);
  if (!std::isfinite(narrowestTubeM) || narrowestTubeM <= 0.0)
    return Invalid("active Parker tube has an invalid finite-line radius");

  // The mask is independent of the visualization point_count.  A segment is
  // no longer than one quarter of either the narrowest physical tube or the
  // smallest leaf scale.  More importantly, every curved arc of length ds is
  // guaranteed to lie within ds/2 of one of its endpoints (which is on the
  // stored chord). Inflating the chord capsule by ds/2 is therefore a strict
  // no-false-negative bound even without assuming a curvature model.
  const double maximumArcStepM =
      0.25 * std::min(spatialScaleM, narrowestTubeM);
  const double rawSegmentCount = std::ceil(lengthM / maximumArcStepM);
  if (!std::isfinite(rawSegmentCount) || rawSegmentCount < 1.0 ||
      rawSegmentCount > 2000000.0) {
    return Invalid("active Parker tube requires an unreasonable segment count");
  }
  const std::size_t segmentCount =
      static_cast<std::size_t>(rawSegmentCount);
  const double arcStepM = lengthM / static_cast<double>(segmentCount);
  const double roundoffM = 256.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, configuration.outerRadiusM);

  std::vector<TubeSegment> candidate;
  candidate.reserve(segmentCount);
  double firstStationRadiusM = configuration.innerRadiusM;
  Core::Vec3 firstPointM = configuration.originM +
      Core::ParkerCurvePoint(firstStationRadiusM, geometry);
  for (std::size_t index = 0; index < segmentCount; ++index) {
    const double secondArcM =
        (index + 1 == segmentCount) ? lengthM : (index + 1) * arcStepM;
    double secondStationRadiusM = 0.0;
    const Core::Status radiusStatus = Core::ParkerCurveRadiusAtArcLengthM(
        secondArcM, geometry, &secondStationRadiusM);
    if (!radiusStatus.ok()) return radiusStatus;
    const Core::Vec3 secondPointM = configuration.originM +
        Core::ParkerCurvePoint(secondStationRadiusM, geometry);
    TubeSegment segment;
    segment.firstM = firstPointM;
    segment.secondM = secondPointM;
    segment.envelopeRadiusM = std::max(
        ActiveTubeRadiusM(firstStationRadiusM, configuration),
        ActiveTubeRadiusM(secondStationRadiusM, configuration));
    segment.envelopeRadiusM += 0.5 * (secondArcM - index * arcStepM) +
                               roundoffM;
    candidate.push_back(segment);
    firstStationRadiusM = secondStationRadiusM;
    firstPointM = secondPointM;
  }
  segments->swap(candidate);
  *effectiveLengthM = lengthM;
  return Core::Status::OK();
}

bool SegmentEnvelopeIntersectsBox(
    const TubeSegment& segment, const Core::Vec3& minimumM,
    const Core::Vec3& maximumM) {
  for (int axis = 0; axis < 3; ++axis) {
    const double segmentMinimum = std::min(
        Component(segment.firstM, axis), Component(segment.secondM, axis));
    const double segmentMaximum = std::max(
        Component(segment.firstM, axis), Component(segment.secondM, axis));
    if (segmentMaximum + segment.envelopeRadiusM <
            Component(minimumM, axis) ||
        segmentMinimum - segment.envelopeRadiusM >
            Component(maximumM, axis)) {
      return false;
    }
  }
  return SquaredDistanceSegmentToBox(
      segment.firstM, segment.secondM, minimumM, maximumM) <=
      segment.envelopeRadiusM * segment.envelopeRadiusM;
}

bool PointInsideBox(const Core::Vec3& point, const LeafBlock& leaf,
                    double tolerance) {
  for (int axis = 0; axis < 3; ++axis) {
    if (Component(point, axis) < Component(leaf.minimumM, axis) - tolerance ||
        Component(point, axis) > Component(leaf.maximumM, axis) + tolerance)
      return false;
  }
  return true;
}

bool BoxesTouch(const LeafBlock& left, const LeafBlock& right,
                double tolerance) {
  for (int axis = 0; axis < 3; ++axis) {
    if (Component(left.maximumM, axis) <
            Component(right.minimumM, axis) - tolerance ||
        Component(right.maximumM, axis) <
            Component(left.minimumM, axis) - tolerance) {
      return false;
    }
  }
  return true;
}

void AddUndirected(std::size_t left, std::size_t right,
                   std::vector<std::vector<std::size_t>>* rows) {
  if (left == right) return;
  (*rows)[left].push_back(right);
  (*rows)[right].push_back(left);
}

void CanonicalizeRows(std::vector<std::vector<std::size_t>>* rows) {
  for (std::vector<std::size_t>& row : *rows) {
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
  }
}

Core::Status ValidateGraph(const LeafNeighbourGraph& graph,
                           std::size_t leafCount) {
  if (graph.face.size() != leafCount || graph.full.size() != leafCount)
    return Invalid("active-region neighbor graph size does not match leaves");
  for (const auto* rows : {&graph.face, &graph.full}) {
    for (std::size_t from = 0; from < rows->size(); ++from) {
      for (std::size_t to : (*rows)[from]) {
        if (to >= leafCount || to == from)
          return Invalid("active-region neighbor graph contains an invalid index");
      }
    }
  }
  return Core::Status::OK();
}

// Bounding-volume hierarchy over capsule boxes, used for the many traced
// lines of an observer flux tube (the analytic corridor keeps its loop).
//
// An observer tube has many lines, each with many chords.  With a plain loop
// the leaf classification and the refinement distance would cost
// (leaves or lattice points) x (segments).  The hierarchy reduces each query
// to roughly the logarithm of the segment count.
//
// Structure: a binary tree of axis-aligned boxes.  Node k covers the
// segments order[begin, end), and its box [lo, hi] bounds all their capsule
// boxes.  A node with 8 or fewer segments is a leaf (left = right = -1).
// Otherwise its segments are split at the median of the chord midpoints
// along the widest axis of the node's box (std::nth_element), so the tree is
// balanced and the build costs O(n log n).  A capsule lies inside its box,
// so a box test is a conservative filter: a query that misses a node's box
// cannot touch any capsule in it.
struct CapsuleBvh {
  struct Node {
    Core::Vec3 lo, hi;               // box [m] bounding the node's capsules
    int left = -1, right = -1;       // child node indices; -1 for a leaf
    std::size_t begin = 0, end = 0;  // the node's range in `order`
  };
  std::vector<Node> nodes;           // nodes[0] is the root (if any)
  std::vector<std::size_t> order;    // permutation of segment indices, grouped by node
  const std::vector<TubeSegment>* segments = nullptr;   // referenced, not owned

  // Axis-aligned box of one capsule: the chord's bounding box grown by the
  // capsule radius on every side.
  static void Box(const TubeSegment& s, Core::Vec3* lo, Core::Vec3* hi) {
    for (int axis = 0; axis < 3; ++axis) {
      const double a = Component(s.firstM, axis), b = Component(s.secondM, axis);
      SetComponent(lo, axis, std::min(a, b) - s.envelopeRadiusM);
      SetComponent(hi, axis, std::max(a, b) + s.envelopeRadiusM);
    }
  }
  // Build the subtree over order[begin, end) and return its node index.
  // Children are built after the parent is appended, so `nodes` may
  // reallocate during recursion.  The child indices are therefore written
  // through nodes[index] at the end, not through a reference held across
  // the recursive calls.
  int Build(std::size_t begin, std::size_t end) {
    Node node;
    node.begin = begin;
    node.end = end;
    for (int axis = 0; axis < 3; ++axis) {
      SetComponent(&node.lo, axis, std::numeric_limits<double>::infinity());
      SetComponent(&node.hi, axis, -std::numeric_limits<double>::infinity());
    }
    for (std::size_t i = begin; i < end; ++i) {
      Core::Vec3 lo, hi;
      Box((*segments)[order[i]], &lo, &hi);
      for (int axis = 0; axis < 3; ++axis) {
        SetComponent(&node.lo, axis, std::min(Component(node.lo, axis), Component(lo, axis)));
        SetComponent(&node.hi, axis, std::max(Component(node.hi, axis), Component(hi, axis)));
      }
    }
    const int index = static_cast<int>(nodes.size());
    nodes.push_back(node);
    if (end - begin <= 8) return index;
    int axis = 0;
    double widest = -1.0;
    for (int a = 0; a < 3; ++a) {
      const double width = Component(node.hi, a) - Component(node.lo, a);
      if (width > widest) { widest = width; axis = a; }
    }
    // Median split along the widest axis.  The key first + second is twice
    // the chord midpoint (the factor 2 does not change the order).
    const std::size_t middle = begin + (end - begin) / 2;
    std::nth_element(order.begin() + begin, order.begin() + middle,
                     order.begin() + end, [&](std::size_t x, std::size_t y) {
      const TubeSegment& sx = (*segments)[x];
      const TubeSegment& sy = (*segments)[y];
      return Component(sx.firstM, axis) + Component(sx.secondM, axis) <
             Component(sy.firstM, axis) + Component(sy.secondM, axis);
    });
    const int left = Build(begin, middle);
    const int right = Build(middle, end);
    nodes[index].left = left;
    nodes[index].right = right;
    return index;
  }
  explicit CapsuleBvh(const std::vector<TubeSegment>& source) : segments(&source) {
    // `source` must outlive the hierarchy (it is referenced, not copied).
    order.resize(source.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    if (!order.empty()) Build(0, order.size());
  }
  // true when any capsule intersects the box [minimumM, maximumM] (an AMR
  // leaf).  Iterative depth-first descent through the nodes whose boxes
  // overlap the query box.  At a leaf node the exact capsule-box test
  // SegmentEnvelopeIntersectsBox (the one the analytic corridor uses) is
  // applied to each segment.  The result equals that of testing every
  // segment, because the node boxes only prune segments whose capsule boxes
  // miss the query box.
  bool Intersects(const Core::Vec3& minimumM, const Core::Vec3& maximumM) const {
    if (nodes.empty()) return false;
    std::vector<int> stack{0};
    while (!stack.empty()) {
      const Node& node = nodes[stack.back()];
      stack.pop_back();
      bool overlap = true;
      for (int axis = 0; axis < 3 && overlap; ++axis)
        overlap = Component(node.hi, axis) >= Component(minimumM, axis) &&
                  Component(node.lo, axis) <= Component(maximumM, axis);
      if (!overlap) continue;
      if (node.left < 0) {
        for (std::size_t i = node.begin; i < node.end; ++i)
          if (SegmentEnvelopeIntersectsBox((*segments)[order[i]], minimumM, maximumM))
            return true;
        continue;
      }
      stack.push_back(node.left);
      stack.push_back(node.right);
    }
    return false;
  }
};

// Squared Euclidean distance [m^2] from p to the box [lo, hi] (0 inside):
// the sum over axes of the squared gap between p and the box's slab.  A
// lower bound for the distance to any capsule inside the box.
double SquaredDistancePointToBox(const Core::Vec3& p, const Core::Vec3& lo,
                                 const Core::Vec3& hi) {
  double sum = 0.0;
  for (int axis = 0; axis < 3; ++axis) {
    const double x = Component(p, axis);
    const double d = x < Component(lo, axis) ? Component(lo, axis) - x
                   : x > Component(hi, axis) ? x - Component(hi, axis) : 0.0;
    sum += d * d;
  }
  return sum;
}

// Distance [m] from p to the segment [a, b]: the projection parameter
// t = (p - a).(b - a) / |b - a|^2 is clamped to [0, 1].  A degenerate segment
// (a == b) is treated as the point a.
double DistancePointToSegment(const Core::Vec3& p, const Core::Vec3& a,
                              const Core::Vec3& b) {
  const Core::Vec3 ab = b - a;
  const double length2 = ab.Dot(ab);
  const double t = length2 > 0.0
      ? std::max(0.0, std::min(1.0, (p - a).Dot(ab) / length2)) : 0.0;
  return (p - (a + ab * t)).Norm();
}

}  // namespace

// The index owns a copy of the capsules in the internal TubeSegment form
// and the hierarchy built over that copy.  CapsuleBvh holds a pointer to
// `segments`.  The index is only ever held through a shared_ptr and never
// copied or moved, so the pointer stays valid for the index's lifetime.
struct ObserverTubeIndex {
  std::vector<TubeSegment> segments;
  std::unique_ptr<CapsuleBvh> bvh;
};

// Convert the public capsules to TubeSegment (the type that
// SegmentEnvelopeIntersectsBox and the analytic corridor use), build the
// hierarchy over them, and publish the index in geometry->index.
Core::Status IndexObserverTube(ObserverTubeGeometry* geometry) {
  if (!geometry) return Invalid("observer tube geometry is null");
  auto index = std::make_shared<ObserverTubeIndex>();
  index->segments.reserve(geometry->segments.size());
  for (const ActiveTubeSegment& item : geometry->segments) {
    TubeSegment segment;
    segment.firstM = item.firstM;
    segment.secondM = item.secondM;
    segment.envelopeRadiusM = item.envelopeRadiusM;
    index->segments.push_back(segment);
  }
  index->bvh.reset(new CapsuleBvh(index->segments));
  geometry->index = index;
  return Core::Status::OK();
}

double DistanceToObserverTubeM(const ObserverTubeGeometry& geometry,
                               const Core::Vec3& positionM, double capM) {
  if (!geometry.index || geometry.index->bvh->nodes.empty()) return capM;
  const CapsuleBvh& bvh = *geometry.index->bvh;
  double best = capM;
  // Nearest-child-first depth-first search: the first leaves reached give a
  // tight bound, so farther boxes are pruned (cost ~ log of the segments).
  // Stack entries are (squared box distance, node).  A node is skipped when
  // its box is no closer than the best capsule distance found so far.  The
  // box contains its capsules, so the box distance is a lower bound and the
  // pruning never discards the true minimum: the result equals the
  // brute-force minimum over all capsules, clipped at capM.
  std::vector<std::pair<double, int>> stack;
  stack.reserve(64);
  stack.push_back({SquaredDistancePointToBox(positionM, bvh.nodes[0].lo,
                                             bvh.nodes[0].hi), 0});
  while (!stack.empty()) {
    const std::pair<double, int> top = stack.back();
    stack.pop_back();
    if (top.first >= best * best) continue;
    const CapsuleBvh::Node& node = bvh.nodes[top.second];
    if (node.left < 0) {
      for (std::size_t i = node.begin; i < node.end; ++i) {
        const TubeSegment& s = (*bvh.segments)[bvh.order[i]];
        const double d = std::max(0.0,
            DistancePointToSegment(positionM, s.firstM, s.secondM) -
            s.envelopeRadiusM);
        best = std::min(best, d);
      }
      // Inside a capsule: no distance can be smaller.
      if (best == 0.0) return 0.0;
      continue;
    }
    const CapsuleBvh::Node& left = bvh.nodes[node.left];
    const CapsuleBvh::Node& right = bvh.nodes[node.right];
    const double dl = SquaredDistancePointToBox(positionM, left.lo, left.hi);
    const double dr = SquaredDistancePointToBox(positionM, right.lo, right.hi);
    // Push the farther child first so the nearer one is examined next.
    if (dl < dr) {
      stack.push_back({dr, node.right});
      stack.push_back({dl, node.left});
    } else {
      stack.push_back({dl, node.left});
      stack.push_back({dr, node.right});
    }
  }
  return best;
}

const char* ActiveRegionAlgorithmName() {
  return RuntimeModel::kActiveRegionAlgorithmName;
}

const char* ActiveRegionAlgorithmName(const ResolutionConfiguration& configuration) {
  // The observer flux tube has its own algorithm name.  It enters the
  // recorded active-region identity, so a traced-tube plan is never
  // mistaken for an analytic-corridor plan.
  if (configuration.activeRegion ==
      RuntimeModel::ActiveRegionMode::ObserverFluxTube)
    return RuntimeModel::kObserverTubeActiveRegionAlgorithmName;
  return configuration.activeSolarSphereRadiusM > 0.0
      ? RuntimeModel::kSolarSphereActiveRegionAlgorithmName
      : RuntimeModel::kActiveRegionAlgorithmName;
}

bool BlockIntersectsActiveRegion(
    const Core::Vec3& minimumM, const Core::Vec3& maximumM,
    const ResolutionConfiguration& configuration) {
  if (configuration.activeRegion ==
      RuntimeModel::ActiveRegionMode::FullDomain) {
    return true;
  }
  if (!ValidBox(minimumM, maximumM)) {
    // Retain malformed boxes so AMPS' structural validator, rather than the
    // pruning pass, remains the authoritative source of the mesh diagnostic.
    return true;
  }
  if (IntersectsSolarNeighbourhood(minimumM, maximumM, configuration))
    return true;
  std::vector<TubeSegment> segments;
  double effectiveLengthM = 0.0;
  const Core::Status built = BuildFiniteTubeSegments(
      configuration, MinimumSide(minimumM, maximumM), &segments,
      &effectiveLengthM);
  (void)effectiveLengthM;
  if (!built.ok()) return true;
  for (const TubeSegment& segment : segments) {
    if (SegmentEnvelopeIntersectsBox(segment, minimumM, maximumM)) return true;
  }
  return false;
}

Core::Status BuildGeometricLeafNeighbourGraph(
    const std::vector<LeafBlock>& leaves, LeafNeighbourGraph* graph) {
  if (graph == nullptr) return Invalid("leaf-neighbor graph output is null");
  if (leaves.empty()) return Invalid("leaf-neighbor graph requires leaves");
  double scaleM = 0.0;
  for (const LeafBlock& leaf : leaves) {
    if (!ValidBox(leaf.minimumM, leaf.maximumM))
      return Invalid("leaf-neighbor graph contains a malformed box");
    scaleM = std::max(scaleM,
        (leaf.maximumM - leaf.minimumM).Norm());
  }
  const double toleranceM = 128.0 *
      std::numeric_limits<double>::epsilon() * std::max(1.0, scaleM);
  LeafNeighbourGraph candidate;
  candidate.face.resize(leaves.size());
  candidate.full.resize(leaves.size());
  for (std::size_t left = 0; left < leaves.size(); ++left) {
    for (std::size_t right = left + 1; right < leaves.size(); ++right) {
      if (FaceNeighbours(leaves[left], leaves[right], toleranceM))
        AddUndirected(left, right, &candidate.face);
      if (BoxesTouch(leaves[left], leaves[right], toleranceM))
        AddUndirected(left, right, &candidate.full);
    }
  }
  CanonicalizeRows(&candidate.face);
  CanonicalizeRows(&candidate.full);
  *graph = std::move(candidate);
  return Core::Status::OK();
}

Core::Status BuildActiveRegionPlan(
    const std::vector<LeafBlock>& leaves,
    const LeafNeighbourGraph& neighbours,
    const ResolutionConfiguration& configuration,
    ActiveRegionPlan* plan) {
  if (plan == nullptr) return Invalid("active-region plan output is null");
  if (leaves.empty()) return Invalid("active-region plan requires leaves");
  const Core::Status graphStatus = ValidateGraph(neighbours, leaves.size());
  if (!graphStatus.ok()) return graphStatus;

  ActiveRegionPlan candidate;
  candidate.leafClass.assign(leaves.size(), ActiveLeafClass::Inactive);
  double minimumLeafSideM = std::numeric_limits<double>::infinity();
  Core::Vec3 domainMinimum = leaves.front().minimumM;
  Core::Vec3 domainMaximum = leaves.front().maximumM;
  for (const LeafBlock& leaf : leaves) {
    if (!ValidBox(leaf.minimumM, leaf.maximumM))
      return Invalid("active-region plan contains a malformed leaf box");
    minimumLeafSideM = std::min(
        minimumLeafSideM, MinimumSide(leaf.minimumM, leaf.maximumM));
    candidate.totalBlockVolumeM3 += BoxVolume(leaf);
    for (int axis = 0; axis < 3; ++axis) {
      SetComponent(&domainMinimum, axis, std::min(
          Component(domainMinimum, axis), Component(leaf.minimumM, axis)));
      SetComponent(&domainMaximum, axis, std::max(
          Component(domainMaximum, axis), Component(leaf.maximumM, axis)));
    }
  }

  if (configuration.activeRegion ==
      RuntimeModel::ActiveRegionMode::FullDomain) {
    std::fill(candidate.leafClass.begin(), candidate.leafClass.end(),
              ActiveLeafClass::Core);
    candidate.coreLeafCount = leaves.size();
    candidate.activeBlockVolumeM3 = candidate.totalBlockVolumeM3;
    *plan = std::move(candidate);
    return Core::Status::OK();
  }

  std::vector<TubeSegment> segments;
  std::vector<Core::Vec3> sourcePoints, terminalPoints;
  const bool observerTube = configuration.activeRegion ==
      RuntimeModel::ActiveRegionMode::ObserverFluxTube;
  if (observerTube) {
    // Traced lines of the observer flux tube (observer_flux_tube.h).
    // The capsules come ready-made from the application (radius = local
    // spacing + half chord + guard), so the minimum leaf side used by
    // BuildFiniteTubeSegments does not enter.  effectiveLineLengthM is the
    // summed chord length of all lines (the total traced line length), not
    // the length of one line.  Every line contributes its inner end as a
    // possible source point and its outer end as a possible terminal point.
    if (!configuration.observerTube)
      return Invalid("observer flux tube geometry has not been traced");
    for (const ActiveTubeSegment& item : configuration.observerTube->segments) {
      TubeSegment segment;
      segment.firstM = item.firstM;
      segment.secondM = item.secondM;
      segment.envelopeRadiusM = item.envelopeRadiusM;
      candidate.effectiveLineLengthM += (item.secondM - item.firstM).Norm();
      segments.push_back(segment);
    }
    sourcePoints = configuration.observerTube->sourcePointsM;
    terminalPoints = configuration.observerTube->terminalPointsM;
  } else {
    const Core::Status segmentStatus = BuildFiniteTubeSegments(
        configuration, minimumLeafSideM, &segments,
        &candidate.effectiveLineLengthM);
    if (!segmentStatus.ok()) return segmentStatus;
    // Analytic corridor: one line, so one source and one terminal point.
    sourcePoints.push_back(segments.front().firstM);
    terminalPoints.push_back(segments.back().secondM);
  }
  candidate.segmentCount = segments.size();
  // The hierarchy keeps a pointer to its segment list: pass the list itself.
  // Only the observer tube uses it.  The analytic corridor keeps its
  // original segment loop below, so it gets an empty hierarchy (a static
  // empty list, so the pointer is valid) and its classification stays
  // unchanged.
  static const std::vector<TubeSegment> kNoSegments;
  const CapsuleBvh bvh(observerTube ? segments : kNoSegments);
  for (std::size_t leafIndex = 0; leafIndex < leaves.size(); ++leafIndex) {
    if (IntersectsSolarNeighbourhood(leaves[leafIndex].minimumM,
                                    leaves[leafIndex].maximumM, configuration)) {
      candidate.leafClass[leafIndex] = ActiveLeafClass::Core;
      continue;
    }
    // Observer tube: a leaf is Core when any capsule intersects it (the same
    // exact test as the loop below, through the hierarchy).
    if (observerTube) {
      if (bvh.Intersects(leaves[leafIndex].minimumM, leaves[leafIndex].maximumM))
        candidate.leafClass[leafIndex] = ActiveLeafClass::Core;
      continue;
    }
    for (const TubeSegment& segment : segments) {
      if (SegmentEnvelopeIntersectsBox(
              segment, leaves[leafIndex].minimumM,
              leaves[leafIndex].maximumM)) {
        candidate.leafClass[leafIndex] = ActiveLeafClass::Core;
        break;
      }
    }
  }

  std::vector<std::size_t> sourceLeaves;
  std::vector<unsigned char> terminalLeaf(leaves.size(), 0);
  const double pointToleranceM = 256.0 *
      std::numeric_limits<double>::epsilon() *
      std::max(1.0, configuration.outerRadiusM);
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    if (candidate.leafClass[index] != ActiveLeafClass::Core) continue;
    // A Core leaf is a source leaf when it contains any line's inner end,
    // and a terminal leaf when it contains any line's outer end.  For the
    // analytic corridor both lists have one point, as before.  For the
    // observer tube the connectivity search below starts from every source
    // leaf.
    bool source = false, terminal = false;
    for (const Core::Vec3& point : sourcePoints)
      if (PointInsideBox(point, leaves[index], pointToleranceM)) { source = true; break; }
    for (const Core::Vec3& point : terminalPoints)
      if (PointInsideBox(point, leaves[index], pointToleranceM)) { terminal = true; break; }
    if (source) sourceLeaves.push_back(index);
    if (terminal) terminalLeaf[index] = 1;
  }
  if (sourceLeaves.empty() ||
      std::find(terminalLeaf.begin(), terminalLeaf.end(), 1) ==
          terminalLeaf.end()) {
    // The message counts Core, source and terminal leaves and segments, so
    // an under-resolved plan can be diagnosed from the log.  The wording
    // "finite Parker tube" is shared with the analytic corridor and applies
    // to the observer tube as well.
    std::size_t core = 0;
    for (ActiveLeafClass value : candidate.leafClass)
      core += value == ActiveLeafClass::Core ? 1 : 0;
    return Invalid("finite Parker tube does not retain its source/end-point leaf "
                   "(core leaves " + std::to_string(core) + ", source leaves " +
                   std::to_string(sourceLeaves.size()) + ", terminal leaves " +
                   std::to_string(std::count(terminalLeaf.begin(), terminalLeaf.end(), 1)) +
                   ", segments " + std::to_string(segments.size()) + ")");
  }

  // A continuous capsule must form one component in the full (26-neighbor)
  // leaf graph. Rejecting a split here catches an under-resolved segment plan
  // before a partial AMPS domain can be allocated.
  std::vector<unsigned char> reached(leaves.size(), 0);
  std::queue<std::size_t> queue;
  for (std::size_t source : sourceLeaves) {
    if (!reached[source]) {
      reached[source] = 1;
      queue.push(source);
    }
  }
  while (!queue.empty()) {
    const std::size_t from = queue.front();
    queue.pop();
    for (std::size_t to : neighbours.full[from]) {
      if (!reached[to] &&
          candidate.leafClass[to] == ActiveLeafClass::Core) {
        reached[to] = 1;
        queue.push(to);
      }
    }
  }
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    if (candidate.leafClass[index] == ActiveLeafClass::Core &&
        !reached[index]) {
      return Invalid("finite Parker tube core is topologically disconnected");
    }
  }

  // buffer_blocks now means exactly this many complete AMR touching-neighbor
  // layers. It no longer scales with the candidate block diagonal, which was
  // the source of coarse islands separated from fine leaves in the old mask.
  std::vector<std::size_t> frontier;
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    if (candidate.leafClass[index] == ActiveLeafClass::Core)
      frontier.push_back(index);
  }
  for (unsigned layer = 0; layer < configuration.activeTubeBufferBlocks;
       ++layer) {
    std::vector<std::size_t> next;
    for (std::size_t from : frontier) {
      for (std::size_t to : neighbours.full[from]) {
        if (candidate.leafClass[to] == ActiveLeafClass::Inactive) {
          candidate.leafClass[to] = ActiveLeafClass::Halo;
          next.push_back(to);
        }
      }
    }
    std::sort(next.begin(), next.end());
    next.erase(std::unique(next.begin(), next.end()), next.end());
    frontier.swap(next);
    if (frontier.empty()) break;
  }

  // Flood inactive leaves inward from the Cartesian domain boundary. Any
  // inactive leaf not reached is a bounded cavity enclosed by active blocks;
  // retain it as a safety halo. This removes the rectangular internal holes
  // seen in active-mesh plots while leaving the exterior pruned.
  std::fill(reached.begin(), reached.end(), 0);
  const double boundaryToleranceM = 256.0 *
      std::numeric_limits<double>::epsilon() *
      std::max(1.0, (domainMaximum - domainMinimum).Norm());
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    if (candidate.leafClass[index] != ActiveLeafClass::Inactive) continue;
    bool boundary = false;
    for (int axis = 0; axis < 3; ++axis) {
      boundary = boundary || std::fabs(
          Component(leaves[index].minimumM, axis) -
          Component(domainMinimum, axis)) <= boundaryToleranceM;
      boundary = boundary || std::fabs(
          Component(leaves[index].maximumM, axis) -
          Component(domainMaximum, axis)) <= boundaryToleranceM;
    }
    if (boundary) {
      reached[index] = 1;
      queue.push(index);
    }
  }
  while (!queue.empty()) {
    const std::size_t from = queue.front();
    queue.pop();
    for (std::size_t to : neighbours.face[from]) {
      if (!reached[to] &&
          candidate.leafClass[to] == ActiveLeafClass::Inactive) {
        reached[to] = 1;
        queue.push(to);
      }
    }
  }
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    if (candidate.leafClass[index] == ActiveLeafClass::Inactive &&
        !reached[index]) {
      candidate.leafClass[index] = ActiveLeafClass::Halo;
      ++candidate.cavityLeafCount;
    }
  }

  // Face connectivity is the relevant final invariant for a mover crossing
  // block boundaries. Check the complete active set, not only the centreline,
  // so a detached coarse halo island cannot pass initialization silently.
  std::fill(reached.begin(), reached.end(), 0);
  for (std::size_t source : sourceLeaves) {
    reached[source] = 1;
    queue.push(source);
  }
  while (!queue.empty()) {
    const std::size_t from = queue.front();
    queue.pop();
    for (std::size_t to : neighbours.face[from]) {
      if (!reached[to] &&
          candidate.leafClass[to] != ActiveLeafClass::Inactive) {
        reached[to] = 1;
        queue.push(to);
      }
    }
  }
  bool terminalReached = false;
  for (std::size_t index = 0; index < leaves.size(); ++index) {
    if (terminalLeaf[index] && reached[index]) terminalReached = true;
    if (candidate.leafClass[index] != ActiveLeafClass::Inactive &&
        !reached[index]) {
      return Invalid("active Parker tube contains a detached leaf component");
    }
  }
  if (!terminalReached)
    return Invalid("active Parker tube has no face-connected source-to-end path");

  for (std::size_t index = 0; index < leaves.size(); ++index) {
    switch (candidate.leafClass[index]) {
      case ActiveLeafClass::Core: ++candidate.coreLeafCount; break;
      case ActiveLeafClass::Halo: ++candidate.haloLeafCount; break;
      case ActiveLeafClass::Inactive: ++candidate.inactiveLeafCount; break;
    }
    if (candidate.leafClass[index] != ActiveLeafClass::Inactive)
      candidate.activeBlockVolumeM3 += BoxVolume(leaves[index]);
  }
  if (candidate.coreLeafCount == 0)
    return Invalid("active Parker tube retained no core leaves");
  *plan = std::move(candidate);
  return Core::Status::OK();
}

Core::Status BuildParkerCenterline(
    const ResolutionConfiguration& configuration,
    std::vector<Core::Vec3>* points) {
  if (points == nullptr) return Invalid("Parker centreline output is null");
  const Core::Status valid = Validate(configuration);
  if (!valid.ok()) return valid;

  std::vector<Core::Vec3> candidate;
  candidate.reserve(static_cast<std::size_t>(configuration.parkerPointCount));
  const Core::ParkerSpiralGeometry geometry = Geometry(configuration);
  const double arcStepM = configuration.parkerLengthM /
      static_cast<double>(configuration.parkerPointCount - 1);
  for (std::uint64_t i = 0; i < configuration.parkerPointCount; ++i) {
    // Equal-arc stations are inverted against the closed Parker arc-length
    // law and then evaluated on the exact analytic field line. The previous
    // midpoint ODE march accumulated a visible transverse drift and, before
    // the geometry correction, did not even follow the curve used by AMR.
    const double arcLengthM = (i + 1 == configuration.parkerPointCount)
        ? configuration.parkerLengthM
        : static_cast<double>(i) * arcStepM;
    double radiusM = 0.0;
    const Core::Status radiusStatus = Core::ParkerCurveRadiusAtArcLengthM(
        arcLengthM, geometry, &radiusM);
    if (!radiusStatus.ok()) return radiusStatus;
    const Core::Vec3 point = (i == 0)
        ? configuration.parkerInitialPointM
        : configuration.originM +
            Core::ParkerCurvePoint(radiusM, geometry);
    if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
        !std::isfinite(point.z)) {
      return Invalid("Parker centreline produced a non-finite point");
    }
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
  const double cubicTolerance = 128.0 *
      std::numeric_limits<double>::epsilon() * std::max(1.0, sideX);
  if (!std::isfinite(sideX) || !std::isfinite(sideY) || !std::isfinite(sideZ) ||
      !(sideX > 0.0) || std::fabs(sideX - sideY) > cubicTolerance ||
      std::fabs(sideX - sideZ) > cubicTolerance) {
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
    const double anchorRadiusM = resolution.solarRefinementAnchor ==
        RuntimeModel::SolarRefinementAnchor::Photosphere
        ? Core::Const::R_sun : resolution.innerRadiusM;
    const double radius = anchorRadiusM +
        (resolution.outerRadiusM - anchorRadiusM) * i /
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
  const double rootCell = (domain.maximumM.x - domain.minimumM.x) /
      resolution.cellsPerBlockEdge;
  constexpr int kVolumeSamplesPerAxis = 17;
  std::uint64_t sampleCounts[20] = {};
  std::uint64_t totalSamples = 0;
  for (int iz = 0; iz < kVolumeSamplesPerAxis; ++iz) {
    for (int iy = 0; iy < kVolumeSamplesPerAxis; ++iy) {
      for (int ix = 0; ix < kVolumeSamplesPerAxis; ++ix) {
        const auto coordinate = [&](int axis, int index) {
          return Component(domain.minimumM, axis) +
              (Component(domain.maximumM, axis) - Component(domain.minimumM, axis)) * index /
              static_cast<double>(kVolumeSamplesPerAxis - 1);
        };
        const Core::Vec3 point(coordinate(0, ix), coordinate(1, iy), coordinate(2, iz));
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
                        "particle crossed inward through the absorbing "
                        "Parker/CME transport source shell");
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
