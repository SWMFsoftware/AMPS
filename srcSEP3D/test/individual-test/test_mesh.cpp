// Phase M mesh/storage verification.
//
// Every test in this file links only the AMPS-independent mesh model.  The
// same resolution function is called by the production AMPS callback, so the
// deterministic leaf histogram can later be compared without maintaining a
// second implementation in the adapter.

#include "../../core/sep3d_test_registry.h"
#include "../../mesh/mesh_model.h"
#include "../../runtime/run_configuration.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <vector>

namespace {

using SEP3D::Testing::Result;
namespace M = SEP3D::Mesh;
namespace RM = SEP3D::RuntimeModel;

Result Pass(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Pass;
  result.message = message;
  return result;
}
Result Fail(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Fail;
  result.message = message;
  return result;
}

M::ResolutionConfiguration Baseline() {
  M::ResolutionConfiguration configuration;
  configuration.minimumCellSizeM = 0.025 * SEP3D::Core::Const::AU;
  configuration.backgroundCellSizeM = 0.25 * SEP3D::Core::Const::AU;
  configuration.innerRadiusM = 0.1 * SEP3D::Core::Const::AU;
  configuration.outerRadiusM = SEP3D::Core::Const::AU;
  configuration.solarSurfaceCellSizeM = configuration.minimumCellSizeM;
  configuration.solarRefinementOuterRadiusM = 0.5 * SEP3D::Core::Const::AU;
  configuration.solarRefinementProfile = RM::RefinementProfile::Linear;
  configuration.parkerInitialPointM =
      {configuration.innerRadiusM, 0.0, 0.0};
  configuration.parkerLengthM = configuration.outerRadiusM;
  configuration.parkerPointCount = 101;
  configuration.maximumLevel = 3;
  configuration.cellsPerBlockEdge = 4;
  return configuration;
}

RM::StorageLayout Layout() {
  RM::RunConfiguration3DOptions options;
  options.innerRadiusM = 0.1 * SEP3D::Core::Const::AU;
  options.outerRadiusM = SEP3D::Core::Const::AU;
  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  if (!RM::RunConfiguration3D::Create(options, &configuration).ok()) return {};
  return configuration->storage_layout();
}

SEP3D::Core::ParkerSpiralGeometry ParkerGeometry(
    const M::ResolutionConfiguration& configuration) {
  SEP3D::Core::ParkerSpiralGeometry geometry;
  geometry.sourceRadiusM = configuration.innerRadiusM;
  geometry.sourceLongitudeRad = configuration.tubeLongitudeRad;
  geometry.sourceColatitudeRad = configuration.tubeColatitudeRad;
  geometry.solarWindSpeedMPerS = configuration.solarWindSpeedMPerS;
  geometry.solarRotationRateRadPerS =
      configuration.solarRotationRateRadPerS;
  geometry.rotationAxis = configuration.rotationAxis;
  return geometry;
}

bool Active(M::ActiveLeafClass value) {
  return value != M::ActiveLeafClass::Inactive;
}

bool Contains(const M::LeafBlock& leaf, const SEP3D::Core::Vec3& point,
              double tolerance) {
  const double p[] = {point.x, point.y, point.z};
  const double lower[] = {
      leaf.minimumM.x, leaf.minimumM.y, leaf.minimumM.z};
  const double upper[] = {
      leaf.maximumM.x, leaf.maximumM.y, leaf.maximumM.z};
  for (int axis = 0; axis < 3; ++axis) {
    if (p[axis] < lower[axis] - tolerance ||
        p[axis] > upper[axis] + tolerance) return false;
  }
  return true;
}

Result RunMSH3D01() {
  const M::ResolutionConfiguration configuration = Baseline();
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;
  auto uniform = [&state]() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11) *
           (1.0 / 9007199254740992.0);
  };
  double observedMinimum = std::numeric_limits<double>::infinity();
  double observedMaximum = 0.0;
  for (std::size_t i = 0; i < 1000000; ++i) {
    const SEP3D::Core::Vec3 position(
        (2.0 * uniform() - 1.0) * configuration.outerRadiusM,
        (2.0 * uniform() - 1.0) * configuration.outerRadiusM,
        (2.0 * uniform() - 1.0) * configuration.outerRadiusM);
    const double value = M::RequestedCellSizeM(position, configuration);
    if (!std::isfinite(value) || value < configuration.minimumCellSizeM ||
        value > configuration.backgroundCellSizeM) {
      return Fail("resolution left its finite configured bounds");
    }
    observedMinimum = std::min(observedMinimum, value);
    observedMaximum = std::max(observedMaximum, value);
  }
  const SEP3D::Core::Vec3 probes[] = {
      {0.0, 0.0, 0.0}, {configuration.outerRadiusM, 0.0, 0.0},
      {0.0, -configuration.outerRadiusM, 0.0},
      {0.0, 0.0, configuration.outerRadiusM}};
  for (const auto& probe : probes) {
    const double value = M::RequestedCellSizeM(probe, configuration);
    if (!std::isfinite(value)) return Fail("axis/boundary probe is non-finite");
  }
  Result result = Pass("one million deterministic points remain within the finite resolution floor/background bounds");
  result.metrics.push_back({"minimum_m", observedMinimum,
                            configuration.minimumCellSizeM, ">=", "m"});
  result.metrics.push_back({"maximum_m", observedMaximum,
                            configuration.backgroundCellSizeM, "<=", "m"});
  return result;
}

Result RunMSH3D02() {
  const M::ResolutionConfiguration configuration = Baseline();
  const double transition = configuration.solarRefinementOuterRadiusM;
  const double midpoint = 0.5 * (configuration.innerRadiusM + transition);
  const double atSurface = M::RequestedCellSizeM(
      {configuration.innerRadiusM, 0.0, 0.0}, configuration);
  const double atTransition = M::RequestedCellSizeM(
      {transition, 0.0, 0.0}, configuration);
  const double atMidpoint = M::RequestedCellSizeM(
      {midpoint, 0.0, 0.0}, configuration);
  const double expectedMidpoint = 0.5 *
      (configuration.solarSurfaceCellSizeM +
       configuration.backgroundCellSizeM);
  if (atSurface != configuration.solarSurfaceCellSizeM ||
      atTransition != configuration.backgroundCellSizeM ||
      atMidpoint != expectedMidpoint) {
    return Fail("linear near-Sun surface, midpoint, or transition identity is not exact");
  }
  return Pass("linear near-Sun degradation matches its surface, midpoint, and transition values exactly");
}

Result RunMSH3D03() {
  M::ResolutionConfiguration configuration = Baseline();
  configuration.enableTubeRefinement = true;
  double worstRelativeDistance = 0.0;
  for (int i = 1; i <= 200; ++i) {
    const double radius = configuration.innerRadiusM +
        (configuration.outerRadiusM - configuration.innerRadiusM) * i / 200.0;
    const SEP3D::Core::Vec3 point =
        radius * M::ParkerTubeDirection(radius, configuration);
    worstRelativeDistance = std::max(
        worstRelativeDistance, M::TubeDistanceM(point, configuration) / radius);
  }
  if (worstRelativeDistance > 1.0e-9)
    return Fail("analytic Parker centreline is not zero-distance within tolerance");
  return Pass("the polarity-independent analytic Parker centreline has distance below 1e-9 radius");
}

Result RunMSH3D04() {
  const M::ResolutionConfiguration configuration = Baseline();
  const double radius = 0.7 * configuration.outerRadiusM;
  const SEP3D::Core::Vec3 centre =
      M::ParkerTubeDirection(radius, configuration);
  auto error = [&](double angle) {
    // Rotate in the XY plane for this equatorial fixture.  The chord is the
    // fast local approximation; TubeDistanceM is the exact spherical arc.
    const double base = std::atan2(centre.y, centre.x);
    const SEP3D::Core::Vec3 point(
        radius * std::cos(base + angle), radius * std::sin(base + angle), 0.0);
    const double exact = M::TubeDistanceM(point, configuration);
    const double chord = radius * std::sqrt(2.0 * (1.0 - std::cos(angle)));
    return std::fabs(chord - exact);
  };
  const double coarse = error(1.0e-2);
  const double fine = error(5.0e-3);
  const double order = std::log(coarse / fine) / std::log(2.0);
  if (!(order >= 1.9)) return Fail("tube-distance local approximation converges below second order");
  Result result = Pass("fast tube-distance approximation converges to the exact arc above second order");
  result.metrics.push_back({"observed_order", order, 1.9, ">=", ""});
  return result;
}

Result RunMSH3D05() {
  M::ResolutionConfiguration configuration = Baseline();
  configuration.enableTubeRefinement = true;
  configuration.enableRadialRefinement = false;
  configuration.maximumLevel = 4;
  configuration.tubeReferenceRadiusM = SEP3D::Core::Const::AU;
  // Deliberately much narrower than a coarse AMR probe cell. Before the
  // capture-envelope fix this curve could pass through blocks unnoticed and
  // produce the disconnected high-resolution island seen in Tecplot.
  configuration.tubeRadiusAtReferenceM = 0.001 * SEP3D::Core::Const::AU;
  configuration.tubeRadiusMode = RM::TubeRadiusMode::PhysicalConstant;
  configuration.tubeCellSizeM = configuration.minimumCellSizeM;
  M::StandaloneOctree mesh;
  RM::RunConfiguration3DOptions domainOptions;
  domainOptions.innerRadiusM = configuration.innerRadiusM;
  domainOptions.outerRadiusMode = RM::OuterRadiusMode::Explicit;
  domainOptions.outerRadiusM = configuration.outerRadiusM;
  const M::DomainBounds domain = M::MakeDomain(domainOptions);
  if (!mesh.Build(domain, configuration, Layout(), 4).ok() ||
      !mesh.IsBalanced()) return Fail("valid tube profile did not produce a balanced octree");

  const double rootCell = 2.0 * configuration.outerRadiusM /
      configuration.cellsPerBlockEdge;
  const double achievable = std::max(
      configuration.tubeCellSizeM,
      rootCell / std::pow(2.0, configuration.maximumLevel));
  for (int sample = 0; sample <= 80; ++sample) {
    const double radius = configuration.innerRadiusM +
        (configuration.outerRadiusM - configuration.innerRadiusM) *
        sample / 80.0;
    const SEP3D::Core::Vec3 point =
        radius * M::ParkerTubeDirection(radius, configuration);
    bool captured = false;
    for (const M::LeafBlock& leaf : mesh.leaves()) {
      if (point.x < leaf.minimumM.x || point.x > leaf.maximumM.x ||
          point.y < leaf.minimumM.y || point.y > leaf.maximumM.y ||
          point.z < leaf.minimumM.z || point.z > leaf.maximumM.z) continue;
      const double cell = (leaf.maximumM.x - leaf.minimumM.x) /
          configuration.cellsPerBlockEdge;
      if (cell <= achievable * (1.0 + 1.0e-12)) captured = true;
      break;
    }
    if (!captured)
      return Fail("a sub-cell Parker centreline segment escaped AMR capture");
  }

  M::LeafBlock coarse;
  coarse.minimumM = {0.0, 0.0, 0.0};
  coarse.maximumM = {1.0, 1.0, 1.0};
  coarse.level = 0;
  M::LeafBlock fine;
  fine.minimumM = {1.0, 0.0, 0.0};
  fine.maximumM = {1.25, 0.25, 0.25};
  fine.level = 2;
  if (M::AreLeavesBalanced({coarse, fine}))
    return Fail("balance negative control did not detect a two-level jump");
  return Pass("a sub-cell Parker tube remains continuously captured and the 2:1 negative control stays live");
}

Result RunMSH3D06() {
  M::ResolutionConfiguration original = Baseline();
  original.enableTubeRefinement = true;
  original.tubeColatitudeRad = 0.5 * SEP3D::Core::Const::kPi;
  const double radius = 0.5 * original.outerRadiusM;
  const SEP3D::Core::Vec3 point(radius * 0.7, radius * 0.5, 0.0);
  const double angle = 0.731;
  M::ResolutionConfiguration rotated = original;
  rotated.tubeLongitudeRad += angle;
  const SEP3D::Core::Vec3 rotatedPoint(
      std::cos(angle) * point.x - std::sin(angle) * point.y,
      std::sin(angle) * point.x + std::cos(angle) * point.y, point.z);
  const double first = M::RequestedCellSizeM(point, original);
  const double second = M::RequestedCellSizeM(rotatedPoint, rotated);
  if (std::fabs(first - second) > 1.0e-12 * std::max(first, second))
    return Fail("co-rotation changed requested resolution");
  return Pass("co-rotating the point and tube longitude preserves resolution to 1e-12 relative");
}

Result RunMSH3D07() {
  const RM::StorageLayout layout = Layout();
  for (unsigned level = 1; level <= 5; ++level) {
    M::ResolutionConfiguration configuration = Baseline();
    configuration.maximumLevel = level;
    configuration.minimumCellSizeM =
        configuration.backgroundCellSizeM / std::pow(2.0, level);
    configuration.solarSurfaceCellSizeM = configuration.minimumCellSizeM;
    M::StandaloneOctree first;
    M::StandaloneOctree second;
    RM::RunConfiguration3DOptions domainOptions;
    domainOptions.innerRadiusM = configuration.innerRadiusM;
    domainOptions.outerRadiusMode = RM::OuterRadiusMode::Explicit;
    domainOptions.outerRadiusM = configuration.outerRadiusM;
    const M::DomainBounds domain = M::MakeDomain(domainOptions);
    if (!first.Build(domain, configuration, layout, 3).ok() ||
        !second.Build(domain, configuration, layout, 3).ok() ||
        first.summary().leafCount != second.summary().leafCount ||
        first.summary().leavesByLevel != second.summary().leavesByLevel ||
        first.summary().estimatedBytes !=
            M::EstimateMemoryBytes(first.summary(), configuration, layout)) {
      return Fail("octree count, histogram, or memory report is not reproducible");
    }
    if (level == 3) {
      M::CellStorage storage;
      if (!storage.Allocate(first, layout).ok())
        return Fail("standalone cell storage allocation failed");
      const std::uint64_t cell = first.leaves().front().firstCell;
      const int owner = storage.Owner(cell);
      const double value = 17.25;
      if (storage.Write(cell, owner + 1, layout.numberDensityOffset,
                        &value, sizeof(value)).code !=
              SEP3D::Core::StatusCode::ConfigurationConflict ||
          !storage.Write(cell, owner, layout.numberDensityOffset,
                         &value, sizeof(value)).ok()) {
        return Fail("owner-only fill contract was not enforced");
      }
      double recovered = 0.0;
      if (!storage.Read(cell, layout.numberDensityOffset,
                        &recovered, sizeof(recovered)).ok() ||
          recovered != value) return Fail("deterministic cell identity did not round-trip storage");
    }
  }
  return Pass("five octrees reproduce leaf histograms/memory exactly and owner-only cell filling is enforced");
}

Result RunMSH3D08() {
  const double inner = 20.0 * SEP3D::Core::Const::R_sun;
  RM::RunConfiguration3DOptions earthOptions;
  earthOptions.domain = RM::DomainPreset::Earth;
  earthOptions.innerRadiusM = inner;
  RM::RunConfiguration3DOptions marsOptions = earthOptions;
  marsOptions.domain = RM::DomainPreset::Mars;
  marsOptions.maximumMeshLevel = 8;
  std::shared_ptr<const RM::RunConfiguration3D> earthConfiguration;
  std::shared_ptr<const RM::RunConfiguration3D> marsConfiguration;
  if (!RM::RunConfiguration3D::Create(earthOptions, &earthConfiguration).ok() ||
      !RM::RunConfiguration3D::Create(marsOptions, &marsConfiguration).ok())
    return Fail("Earth/Mars preset normalization failed");
  const M::DomainBounds earth = M::MakeDomain(earthConfiguration->options());
  const M::DomainBounds mars = M::MakeDomain(marsConfiguration->options());
  if (earth.innerRadiusM != inner || mars.innerRadiusM != inner ||
      earth.outerRadiusM != SEP3D::Core::Const::AU ||
      mars.outerRadiusM != 1.666 * SEP3D::Core::Const::AU ||
      earth.minimumM.x != -earth.outerRadiusM ||
      mars.maximumM.z != mars.outerRadiusM) {
    return Fail("Earth/Mars preset bounds or shell coverage changed");
  }
  return Pass("Earth and Mars presets exactly enclose their declared outer spheres and retain the inner boundary");
}

Result RunMSH3D09() {
  const SEP3D::Core::Vec3 center(0.0, 0.0, 0.0);
  const SEP3D::Core::Vec3 exact(2.0, -3.0, 0.5);
  const std::vector<SEP3D::Core::Vec3> points = {
      {1.0, 0.0, 0.0}, {-0.5, 0.0, 0.0}, {0.0, 2.0, 0.0},
      {0.0, -1.0, 0.0}, {0.0, 0.0, 0.25}, {0.0, 0.0, -2.0}};
  std::vector<double> values;
  for (const auto& point : points) values.push_back(7.0 + exact.Dot(point));
  SEP3D::Core::Vec3 observed;
  if (!M::ReconstructScalarGradient(center, 7.0, points, values, &observed).ok() ||
      (observed - exact).Norm() > 1.0e-13) {
    return Fail("mixed-spacing refinement-boundary stencil is not linear exact");
  }
  const std::vector<SEP3D::Core::Vec3> deficient = {
      {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
  const std::vector<double> deficientValues = {1.0, 2.0, 3.0};
  if (M::ReconstructScalarGradient(center, 0.0, deficient,
                                   deficientValues, &observed).ok()) {
    return Fail("rank-deficient gradient stencil was silently accepted");
  }
  return Pass("mixed coarse/fine least-squares gradients are linear exact and reject rank-deficient stencils");
}

Result RunMSH3D10() {
  M::ResolutionConfiguration configuration = Baseline();
  configuration.enableTubeRefinement = true;
  configuration.parkerLengthM = 0.7 * SEP3D::Core::Const::AU;
  configuration.parkerPointCount = 401;
  std::vector<SEP3D::Core::Vec3> points;
  if (!M::BuildParkerCenterline(configuration, &points).ok() ||
      points.size() != configuration.parkerPointCount ||
      !(points.front() == configuration.parkerInitialPointM)) {
    return Fail("finite Parker centreline count or initial point is wrong");
  }
  const SEP3D::Core::ParkerSpiralGeometry geometry =
      ParkerGeometry(configuration);
  const double arcStepM = configuration.parkerLengthM /
      static_cast<double>(points.size() - 1);
  for (std::size_t i = 0; i < points.size(); ++i) {
    const double radiusM = (points[i] - configuration.originM).Norm();
    const double observedArcM = SEP3D::Core::ParkerCurveArcLengthM(
        radiusM, geometry);
    const double expectedArcM = (i + 1 == points.size())
        ? configuration.parkerLengthM : i * arcStepM;
    const SEP3D::Core::Vec3 exact = configuration.originM +
        SEP3D::Core::ParkerCurvePoint(radiusM, geometry);
    if (std::fabs(observedArcM - expectedArcM) >
            2.0e-12 * configuration.parkerLengthM ||
        (points[i] - exact).Norm() >
            2.0e-13 * configuration.outerRadiusM) {
      return Fail("finite Parker output is not on its exact equal-arc station");
    }
  }

  const SEP3D::Core::Vec3 shift(3.0e9, -4.0e9, 2.0e9);
  M::ResolutionConfiguration translated = configuration;
  translated.originM += shift;
  translated.parkerInitialPointM += shift;
  const SEP3D::Core::Vec3 probe = points[points.size() / 2];
  const double original = M::RequestedCellSizeM(probe, configuration);
  const double moved = M::RequestedCellSizeM(probe + shift, translated);
  if (std::fabs(original - moved) >
      1.0e-12 * configuration.backgroundCellSizeM) {
    return Fail("translated domain changed the Parker-tube resolution law");
  }
  return Pass("finite Parker sampling uses exact equal-arc field-line stations and mesh refinement is origin-relative");
}

Result RunMSH3D11() {
  const M::ResolutionConfiguration configuration = Baseline();
  const char* path = "test/.initialization-parker-line-test.dat";
  const SEP3D::Core::Status written =
      M::WriteParkerCenterlineTecplot(configuration, path);
  std::ifstream input(path);
  std::ostringstream text;
  text << input.rdbuf();
  const std::string contents = text.str();
  input.close();
  std::remove(path);
  if (!written.ok() || contents.find("TITLE=\"srcSEP3D initialized Parker centreline\"") ==
          std::string::npos ||
      contents.find("ZONE T=\"parker-centreline\", I=101, F=POINT") ==
          std::string::npos ||
      contents.find("x_m") == std::string::npos)
    return Fail("initialization Parker line was not written as unit-labeled Tecplot data");
  return Pass("initialized finite Parker line is written as deterministic unit-labeled Tecplot data");
}

Result RunMSH3D12() {
  M::ResolutionConfiguration configuration = Baseline();
  configuration.activeRegion = RM::ActiveRegionMode::ParkerTube;
  configuration.activeTubeReferenceRadiusM = SEP3D::Core::Const::AU;
  configuration.activeTubeRadiusAtReferenceM =
      0.04 * SEP3D::Core::Const::AU;
  configuration.activeTubeRadiusMode = RM::TubeRadiusMode::PhysicalConstant;
  configuration.activeTubeBufferBlocks = 1;
  const double radiusM = 0.7 * SEP3D::Core::Const::AU;
  const SEP3D::Core::Vec3 centreline =
      radiusM * M::ParkerTubeDirection(radiusM, configuration);
  const double halfSideM = 0.005 * SEP3D::Core::Const::AU;
  const SEP3D::Core::Vec3 half(halfSideM, halfSideM, halfSideM);
  if (!M::BlockIntersectsActiveRegion(
          centreline - half, centreline + half, configuration)) {
    return Fail("a block centred on the Parker line was deactivated");
  }

  const SEP3D::Core::Vec3 opposite = -1.0 * centreline;
  if (M::BlockIntersectsActiveRegion(
          opposite - half, opposite + half, configuration)) {
    return Fail("a remote opposite-longitude block was retained");
  }

  M::ResolutionConfiguration angular = configuration;
  angular.activeTubeRadiusMode = RM::TubeRadiusMode::ConstantAngularWidth;
  const double innerWidth = M::ActiveTubeRadiusM(
      0.5 * SEP3D::Core::Const::AU, angular);
  const double outerWidth = M::ActiveTubeRadiusM(
      SEP3D::Core::Const::AU, angular);
  if (std::fabs(2.0 * innerWidth - outerWidth) >
      1.0e-14 * outerWidth) {
    return Fail("constant-angular active radius did not scale linearly");
  }

  configuration.activeRegion = RM::ActiveRegionMode::FullDomain;
  if (!M::BlockIntersectsActiveRegion(
          opposite - half, opposite + half, configuration)) {
    return Fail("full-domain mode unexpectedly deactivated a block");
  }
  return Pass("finite Parker capsule retains intersected blocks, rejects remote blocks, and preserves full-domain mode");
}

Result RunMSH3D13() {
  M::ResolutionConfiguration configuration = Baseline();
  configuration.tubeLongitudeRad = 0.37;
  configuration.tubeColatitudeRad = 1.13;
  configuration.rotationAxis =
      SEP3D::Core::Vec3(0.21, -0.32, 0.91).Normalized();
  const SEP3D::Core::ParkerSpiralGeometry geometry =
      ParkerGeometry(configuration);
  double worstAngularError = 0.0;
  double worstArcRoundTrip = 0.0;
  for (int i = 1; i <= 200; ++i) {
    const double radiusM = configuration.innerRadiusM +
        (configuration.outerRadiusM - configuration.innerRadiusM) * i / 200.0;
    const double stepM = 1.0e-6 * radiusM;
    const SEP3D::Core::Vec3 derivative =
        (SEP3D::Core::ParkerCurvePoint(radiusM + stepM, geometry) -
         SEP3D::Core::ParkerCurvePoint(radiusM - stepM, geometry)) /
        (2.0 * stepM);
    const SEP3D::Core::Vec3 tangent =
        SEP3D::Core::ParkerCurveTangent(radiusM, geometry);
    const double angularError =
        derivative.Normalized().Cross(tangent).Norm();
    if (derivative.Dot(tangent) <= 0.0)
      return Fail("exact Parker curve derivative points against the IMF tangent");
    worstAngularError = std::max(worstAngularError, angularError);

    const double arcM = SEP3D::Core::ParkerCurveArcLengthM(radiusM, geometry);
    double recoveredRadiusM = 0.0;
    if (!SEP3D::Core::ParkerCurveRadiusAtArcLengthM(
            arcM, geometry, &recoveredRadiusM).ok()) {
      return Fail("Parker arc-length inversion failed");
    }
    worstArcRoundTrip = std::max(
        worstArcRoundTrip, std::fabs(recoveredRadiusM - radiusM));
  }
  if (worstAngularError > 2.0e-10 ||
      worstArcRoundTrip > 2.0e-12 * configuration.outerRadiusM) {
    return Fail("Parker curve, tangent, and arc-length inverse are inconsistent");
  }
  Result result = Pass(
      "exact Parker curve is tangent to the initialized-field law for a rotated axis and arc length round-trips");
  result.metrics.push_back({"maximum_tangent_cross_norm", worstAngularError,
                            2.0e-10, "<=", ""});
  result.metrics.push_back({"maximum_arc_roundtrip_m", worstArcRoundTrip,
                            2.0e-12 * configuration.outerRadiusM, "<=", "m"});
  return result;
}

Result RunMSH3D14() {
  M::ResolutionConfiguration configuration = Baseline();
  configuration.enableTubeRefinement = true;
  configuration.tubeReferenceRadiusM = SEP3D::Core::Const::AU;
  configuration.tubeRadiusAtReferenceM =
      0.03 * SEP3D::Core::Const::AU;
  configuration.tubeRadiusMode = RM::TubeRadiusMode::ConstantAngularWidth;
  configuration.tubeCellSizeM = configuration.minimumCellSizeM;
  configuration.parkerLengthM = 2.0 * SEP3D::Core::Const::AU;
  configuration.parkerPointCount = 401;
  configuration.activeRegion = RM::ActiveRegionMode::ParkerTube;
  configuration.activeTubeReferenceRadiusM = SEP3D::Core::Const::AU;
  configuration.activeTubeRadiusAtReferenceM =
      0.05 * SEP3D::Core::Const::AU;
  configuration.activeTubeRadiusMode =
      RM::TubeRadiusMode::ConstantAngularWidth;
  configuration.activeTubeBufferBlocks = 0;

  RM::RunConfiguration3DOptions domainOptions;
  domainOptions.innerRadiusM = configuration.innerRadiusM;
  domainOptions.outerRadiusMode = RM::OuterRadiusMode::Explicit;
  domainOptions.outerRadiusM = configuration.outerRadiusM;
  const M::DomainBounds domain = M::MakeDomain(domainOptions);
  M::StandaloneOctree mesh;
  if (!mesh.Build(domain, configuration, Layout(), 4).ok())
    return Fail("could not construct active-mask regression octree");
  M::LeafNeighbourGraph graph;
  if (!M::BuildGeometricLeafNeighbourGraph(mesh.leaves(), &graph).ok())
    return Fail("could not construct active-mask regression neighbor graph");

  M::ActiveRegionPlan plans[3];
  for (unsigned layers = 0; layers < 3; ++layers) {
    configuration.activeTubeBufferBlocks = layers;
    const SEP3D::Core::Status status = M::BuildActiveRegionPlan(
        mesh.leaves(), graph, configuration, &plans[layers]);
    if (!status.ok())
      return Fail("hole-free active-region planner rejected the AMR fixture: " +
                  status.message);
    if (plans[layers].cavityLeafCount != 0)
      return Fail("well-resolved Parker fixture unexpectedly required cavity filling");
  }
  const std::size_t total = mesh.leaves().size();
  if (plans[0].coreLeafCount == 0 || plans[0].inactiveLeafCount == 0 ||
      plans[0].coreLeafCount >= total)
    return Fail("physical Parker tube did not prune a strict subset of leaves");

  // Exact neighbor-layer semantics: every newly active leaf at N=1 is a
  // touching neighbor of the N=0 set, and similarly for N=2. This catches the
  // former candidate-local block-diagonal inflation at coarse/fine interfaces.
  for (unsigned layers = 1; layers < 3; ++layers) {
    for (std::size_t index = 0; index < total; ++index) {
      if (!Active(plans[layers].leafClass[index]) ||
          Active(plans[layers - 1].leafClass[index])) continue;
      bool hasPreviousNeighbour = false;
      for (std::size_t neighbour : graph.full[index]) {
        if (Active(plans[layers - 1].leafClass[neighbour])) {
          hasPreviousNeighbour = true;
          break;
        }
      }
      if (!hasPreviousNeighbour)
        return Fail("buffer_blocks retained a leaf beyond its topological layer");
    }
    for (std::size_t index = 0; index < total; ++index) {
      if (Active(plans[layers - 1].leafClass[index]) &&
          !Active(plans[layers].leafClass[index]))
        return Fail("increasing buffer_blocks removed an active leaf");
    }
  }

  // A dense, independent sampling of the exact field line must always land
  // in an active physical-core leaf. This is the direct no-hole invariant.
  const SEP3D::Core::ParkerSpiralGeometry geometry =
      ParkerGeometry(configuration);
  const double effectiveLengthM = std::min(
      configuration.parkerLengthM,
      SEP3D::Core::ParkerCurveArcLengthM(
          configuration.outerRadiusM, geometry));
  const double toleranceM = 1.0e-12 * configuration.outerRadiusM;
  for (int station = 0; station <= 2000; ++station) {
    double radiusM = 0.0;
    if (!SEP3D::Core::ParkerCurveRadiusAtArcLengthM(
            effectiveLengthM * station / 2000.0,
            geometry, &radiusM).ok())
      return Fail("dense Parker coverage oracle could not invert arc length");
    const SEP3D::Core::Vec3 point = configuration.originM +
        SEP3D::Core::ParkerCurvePoint(radiusM, geometry);
    bool covered = false;
    for (std::size_t index = 0; index < total; ++index) {
      if (plans[0].leafClass[index] == M::ActiveLeafClass::Core &&
          Contains(mesh.leaves()[index], point, toleranceM)) {
        covered = true;
        break;
      }
    }
    if (!covered)
      return Fail("finite Parker centreline crosses an inactive AMR leaf");
  }

  // The finite line ends at the physical outer sphere in this fixture. A box
  // on the analytic continuation well beyond that endpoint must not be kept;
  // this proves the Cartesian cube corners cannot grow an unintended branch.
  const double remoteRadiusM = 1.35 * configuration.outerRadiusM;
  const SEP3D::Core::Vec3 remote = configuration.originM +
      SEP3D::Core::ParkerCurvePoint(remoteRadiusM, geometry);
  const SEP3D::Core::Vec3 remoteHalf(
      0.005 * SEP3D::Core::Const::AU,
      0.005 * SEP3D::Core::Const::AU,
      0.005 * SEP3D::Core::Const::AU);
  if (M::BlockIntersectsActiveRegion(
          remote - remoteHalf, remote + remoteHalf, configuration))
    return Fail("finite active tube retained its analytic continuation");

  Result result = Pass(
      "finite Parker capsule covers every dense line station, prunes the cube, and grows by exact AMR-neighbor layers without cavities");
  result.metrics.push_back({"total_leaves", static_cast<double>(total),
                            1.0, ">=", "leaves"});
  result.metrics.push_back({"core_leaves",
                            static_cast<double>(plans[0].coreLeafCount),
                            1.0, ">=", "leaves"});
  result.metrics.push_back({"inactive_leaves",
                            static_cast<double>(plans[0].inactiveLeafCount),
                            1.0, ">=", "leaves"});
  return result;
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterMeshTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using R = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* name, const char* description,
                 SEP3D::Testing::TestCallback callback) {
    D d;
    d.id = id; d.name = name; d.group = "MSH3D"; d.description = description;
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = R::Routine; d.seedPolicy = "deterministic";
    d.stateIsolation = "fresh mesh configuration and octree per test";
    d.callback = std::move(callback); return d;
  };
  return {
      make("MSH3D01", "Resolution bounds", "One million deterministic resolution probes.", RunMSH3D01),
      make("MSH3D02", "Radial closed forms", "Surface, midpoint, and transition identities.", RunMSH3D02),
      make("MSH3D03", "Tube centreline", "Analytic Parker centreline distance.", RunMSH3D03),
      make("MSH3D04", "Tube distance convergence", "Fast chord converges to exact arc.", RunMSH3D04),
      make("MSH3D05", "Tube profile and balance", "Composite tube refinement with a live 2:1 negative control.", RunMSH3D05),
      make("MSH3D06", "Rotation invariance", "Co-rotation leaves resolution invariant.", RunMSH3D06),
      make("MSH3D07", "Octree budget and ownership", "Reproducible memory and owner-only fills.", RunMSH3D07),
      make("MSH3D08", "Earth and Mars presets", "Exact preset bounds and shell coverage.", RunMSH3D08),
      make("MSH3D09", "Refinement gradients", "Mixed-spacing gradient reconstruction.", RunMSH3D09),
      make("MSH3D10", "Finite Parker initialization", "Point-count, arc-length, and translated-origin identities.", RunMSH3D10),
      make("MSH3D11", "Initialization Tecplot", "Finite Parker-line visualization output.", RunMSH3D11),
      make("MSH3D12", "Active Parker corridor", "Finite capsule single-block classifier contract.", RunMSH3D12),
      make("MSH3D13", "Parker geometry authority", "Exact curve/tangent/arc-length identity for a rotated axis.", RunMSH3D13),
      make("MSH3D14", "Hole-free active mask", "Finite-tube coverage, pruning, cavity, and topological halo contract.", RunMSH3D14),
  };
}
