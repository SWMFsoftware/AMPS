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
  configuration.tubeReferenceRadiusM = SEP3D::Core::Const::AU;
  configuration.tubeRadiusAtReferenceM = 0.12 * SEP3D::Core::Const::AU;
  configuration.tubeCellSizeM = configuration.minimumCellSizeM;
  M::StandaloneOctree mesh;
  RM::RunConfiguration3DOptions domainOptions;
  domainOptions.innerRadiusM = configuration.innerRadiusM;
  domainOptions.outerRadiusMode = RM::OuterRadiusMode::Explicit;
  domainOptions.outerRadiusM = configuration.outerRadiusM;
  const M::DomainBounds domain = M::MakeDomain(domainOptions);
  if (!mesh.Build(domain, configuration, Layout(), 4).ok() ||
      !mesh.IsBalanced()) return Fail("valid tube profile did not produce a balanced octree");

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
  return Pass("tube profile produces a 2:1 mesh and the negative control detects an illegal level jump");
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
  double polylineLength = 0.0;
  for (std::size_t i = 1; i < points.size(); ++i)
    polylineLength += (points[i] - points[i - 1]).Norm();
  if (std::fabs(polylineLength - configuration.parkerLengthM) >
      1.0e-12 * configuration.parkerLengthM) {
    return Fail("sampled Parker polyline did not preserve the configured arc length");
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
  return Pass("finite Parker sampling preserves point count/length and mesh refinement is origin-relative");
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
  };
}
