// Portable acceptance probe for the actual production geometry and mesh
// kernels. It needs no replacement AMPS/SWCME headers or numerical mocks.
#include "../../mesh/mesh_model.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace C = SEP3D::Core;
namespace M = SEP3D::Mesh;
namespace R = SEP3D::RuntimeModel;

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}
bool Inside(const C::Vec3& p, const C::Vec3& a, const C::Vec3& b,
            double pad = 0.0) {
  return p.x >= a.x + pad && p.x <= b.x - pad &&
         p.y >= a.y + pad && p.y <= b.y - pad &&
         p.z >= a.z + pad && p.z <= b.z - pad;
}
R::RunConfiguration3DOptions Options() {
  R::RunConfiguration3DOptions c;
  c.innerRadiusM = 0.1 * C::Const::AU;
  c.outerRadiusM = C::Const::AU;
  c.outerRadiusMode = R::OuterRadiusMode::FieldLineEndpoint;
  c.parkerSpiralEndRadiusM = c.outerRadiusM;
  c.domainBoxGeometry = R::DomainBoxGeometry::FieldLineCornerCube;
  c.domainCornerMarginM = 0.01 * C::Const::AU;
  c.activeRegion = R::ActiveRegionMode::ParkerTube;
  c.activeSolarSphereRadiusM = 0.2 * C::Const::AU;
  c.activeTubeRadiusAtReferenceM = 0.05 * C::Const::AU;
  return c;
}
C::ParkerSpiralGeometry Geometry(const R::RunConfiguration3DOptions& c) {
  C::ParkerSpiralGeometry g;
  g.sourceRadiusM = c.innerRadiusM;
  g.sourceLongitudeRad = c.tubeLongitudeRad;
  g.sourceColatitudeRad = c.tubeColatitudeRad;
  g.solarWindSpeedMPerS = c.parker.solarWindSpeedMPerS;
  g.solarRotationRateRadPerS = c.parker.solarRotationRateRadPerS;
  g.rotationAxis = c.parker.rotationAxis;
  return g;
}
M::ResolutionConfiguration Resolution(const R::RunConfiguration3DOptions& o) {
  M::ResolutionConfiguration c;
  c.originM = o.coordinateOriginM;
  c.rotationAxis = o.parker.rotationAxis;
  c.solarWindSpeedMPerS = o.parker.solarWindSpeedMPerS;
  c.solarRotationRateRadPerS = o.parker.solarRotationRateRadPerS;
  c.tubeLongitudeRad = o.tubeLongitudeRad;
  c.tubeColatitudeRad = o.tubeColatitudeRad;
  c.innerRadiusM = o.innerRadiusM;
  c.outerRadiusM = o.outerRadiusM;
  c.parkerInitialPointM = c.originM + C::ParkerCurvePoint(c.innerRadiusM, Geometry(o));
  c.parkerLengthM = C::ParkerCurveArcLengthM(c.outerRadiusM, Geometry(o));
  c.parkerPointCount = 21;
  c.activeRegion = o.activeRegion;
  c.activeSolarSphereRadiusM = o.activeSolarSphereRadiusM;
  c.activeTubeRadiusAtReferenceM = o.activeTubeRadiusAtReferenceM;
  c.activeTubeBufferBlocks = 1;
  c.solarRefinementAnchor = R::SolarRefinementAnchor::Photosphere;
  c.solarRefinementOuterRadiusM = 0.3 * C::Const::AU;
  c.minimumCellSizeM = c.solarSurfaceCellSizeM = 0.06 * C::Const::AU;
  c.backgroundCellSizeM = 0.3 * C::Const::AU;
  c.maximumLevel = 3;
  c.memoryBudgetBytes = std::size_t{64} * 1024 * 1024 * 1024;
  return c;
}

void DomainBounds() {
  auto o = Options();
  // Arbitrary axes, translations, weak/strong winding and all eight corners
  // exercise the actual complete-curve enclosure, not an endpoint AABB.
  for (double spin : {0.0, 2.86533e-6, 1.4e-5}) {
    o.parker.solarRotationRateRadPerS = spin;
    o.parker.rotationAxis = C::Vec3(1.0, 2.0, 3.0);
    o.coordinateOriginM = C::Vec3(0.03, -0.04, 0.02) * C::Const::AU;
    for (int corner = 0; corner < 9; ++corner) {
      o.domainCornerDirection = corner == 8 ? C::Vec3() :
          C::Vec3(corner & 1 ? 1 : -1, corner & 2 ? 1 : -1, corner & 4 ? 1 : -1);
      const auto d = M::MakeDomain(o);
      const double side = d.maximumM.x - d.minimumM.x;
      Require(std::fabs(side - (d.maximumM.y - d.minimumM.y)) < 1e-12 * side,
              "root is not cubic");
      Require(Inside(o.coordinateOriginM, d.minimumM, d.maximumM,
                     o.activeSolarSphereRadiusM), "full near-Sun sphere clipped");
      const auto geometry = Geometry(o);
      for (int i = 0; i <= 4000; ++i) {
        const double radius = o.innerRadiusM +
            (o.outerRadiusM - o.innerRadiusM) * i / 4000.0;
        const auto p = o.coordinateOriginM + C::ParkerCurvePoint(radius, geometry);
        Require(Inside(p, d.minimumM, d.maximumM,
                       o.activeTubeRadiusAtReferenceM * radius / C::Const::AU),
                "a bent corridor or cross-section was clipped");
      }
      o.parkerSpiralPointCount = 2;
      const auto sparse = M::MakeDomain(o);
      o.parkerSpiralPointCount = 10001;
      const auto dense = M::MakeDomain(o);
      Require((sparse.minimumM - dense.minimumM).Norm() == 0.0 &&
              (sparse.maximumM - dense.maximumM).Norm() == 0.0,
              "plot point count changes domain size");
    }
  }
  // A broad corridor on the opposite side of the chosen corner must not
  // push the near-Sun sphere through the far face. Check both faces even
  // when the padded line's farthest coordinate is smaller than the sphere.
  o = Options();
  o.parker.solarRotationRateRadPerS = 0.0;
  o.tubeLongitudeRad = C::Const::kPi;
  o.activeTubeRadiusMode = R::TubeRadiusMode::PhysicalConstant;
  o.activeTubeRadiusAtReferenceM = 0.25 * C::Const::AU;
  for (int corner = 0; corner < 8; ++corner) {
    o.domainCornerDirection = C::Vec3(corner & 1 ? 1 : -1,
        corner & 2 ? 1 : -1, corner & 4 ? 1 : -1);
    const auto broad = M::MakeDomain(o);
    Require(Inside(o.coordinateOriginM, broad.minimumM, broad.maximumM,
                   o.activeSolarSphereRadiusM),
            "broad opposite corridor clips the solar sphere at a far face");
  }
  o = Options();
  const double small = (M::MakeDomain(o).maximumM - M::MakeDomain(o).minimumM).x;
  o.outerRadiusM *= 1.3;
  const double large = (M::MakeDomain(o).maximumM - M::MakeDomain(o).minimumM).x;
  Require(large > small, "endpoint distance does not control domain size");
  o.domainCornerDirection.x = 2.0;
  C::Vec3 a, b;
  Require(!R::ResolveDomainBoundsM(o, &a, &b).ok(), "invalid corner sign accepted");
  o = Options();
  o.activeSolarSphereRadiusM = o.outerRadiusM;
  Require(!R::ResolveDomainBoundsM(o, &a, &b).ok(), "invalid sphere accepted");
  o = Options();
  o.activeSolarSphereRadiusM = 0.5 * o.innerRadiusM;
  Require(!R::ResolveDomainBoundsM(o, &a, &b).ok(), "disconnected sphere bounds accepted");
  o = Options();
  o.domainBoxGeometry = R::DomainBoxGeometry::SunCenteredCube;
  const auto legacy = M::MakeDomain(o);
  Require(legacy.minimumM.x == -C::Const::AU &&
          legacy.maximumM.z == C::Const::AU, "legacy centered bounds changed");
}

void SphereAndResolution() {
  auto c = Resolution(Options());
  const C::Vec3 point(-0.15 * C::Const::AU, 0.0, 0.0);
  const C::Vec3 half(0.005 * C::Const::AU, 0.005 * C::Const::AU, 0.005 * C::Const::AU);
  Require(M::BlockIntersectsActiveRegion(point - half, point + half, c),
          "off-corridor solar-neighbourhood leaf was removed");
  c.activeSolarSphereRadiusM = 0.0;
  Require(!M::BlockIntersectsActiveRegion(point - half, point + half, c),
          "disabled sphere still retains off-corridor leaves");
  c.activeSolarSphereRadiusM = 0.2 * C::Const::AU;
  const C::Vec3 tangent(-0.2 * C::Const::AU, 0.0, 0.0);
  Require(M::BlockIntersectsActiveRegion(tangent - half, tangent, c),
          "sphere/box face tangency was missed");
  const double r0 = C::Const::R_sun;
  const double r1 = c.solarRefinementOuterRadiusM;
  for (auto profile : {R::RefinementProfile::Linear,
                       R::RefinementProfile::PowerLaw,
                       R::RefinementProfile::Smoothstep}) {
    c.solarRefinementProfile = profile;
    c.solarRefinementExponent = 2.0;
    const double middle = M::RequestedCellSizeM(C::Vec3(0.5 * (r0 + r1), 0.0, 0.0), c);
    const double fraction = profile == R::RefinementProfile::Linear ? 0.5 : 0.25;
    const double expected = c.solarSurfaceCellSizeM + fraction *
        (c.backgroundCellSizeM - c.solarSurfaceCellSizeM);
    Require(std::fabs(middle - expected) < 1e-12 * expected, "radial coarsening law differs");
    Require(M::RequestedCellSizeM(C::Vec3(r0, 0.0, 0.0), c) == c.solarSurfaceCellSizeM,
            "photosphere does not have the surface target");
    Require(M::RequestedCellSizeM(C::Vec3(r1, 0.0, 0.0), c) == c.backgroundCellSizeM,
            "transition does not reach global resolution");
    double previous = 0.0;
    for (int i = 0; i <= 500; ++i) {
      const double cell = M::RequestedCellSizeM(C::Vec3(r0 + (r1-r0)*i/500.0, 0.0, 0.0), c);
      Require(cell >= previous, "radial coarsening is not monotone");
      previous = cell;
    }
  }
  c.solarRefinementAnchor = R::SolarRefinementAnchor::SourceShell;
  Require(M::RequestedCellSizeM(C::Vec3(c.innerRadiusM, 0.0, 0.0), c) == c.solarSurfaceCellSizeM,
          "legacy source-shell anchor changed");
  c.activeSolarSphereRadiusM = 0.05 * C::Const::AU;
  Require(!M::Validate(c).ok(), "disconnected solar sphere accepted");
}

void MeshAllocation(bool centerZ = false) {
  auto o = Options();
  if (centerZ) o.domainBoxGeometry = R::DomainBoxGeometry::FieldLineXYCornerCube;
  auto c = Resolution(o);
  // Resolve several coarse leaves across the root so one topological halo
  // cannot legitimately consume the entire tiny planning fixture.
  c.backgroundCellSizeM = 0.12 * C::Const::AU;
  c.minimumCellSizeM = c.solarSurfaceCellSizeM = 0.03 * C::Const::AU;
  c.maximumLevel = 4;
  const auto d = M::MakeDomain(o);
  M::StandaloneOctree mesh;
  const auto built = mesh.Build(d, c, R::StorageLayout(), 4);
  Require(built.ok(), "corner octree build failed: " + built.message);
  M::LeafNeighbourGraph graph;
  Require(M::BuildGeometricLeafNeighbourGraph(mesh.leaves(), &graph).ok(), "neighbor graph failed");
  M::ActiveRegionPlan plan;
  const auto planned = M::BuildActiveRegionPlan(mesh.leaves(), graph, c, &plan);
  Require(planned.ok(), "sphere/corridor union is disconnected: " + planned.message);
  Require(plan.coreLeafCount > 0 && plan.inactiveLeafCount > 0, "active union does not prune the root");
  if (centerZ) {
    // In the equatorial fixture both the sphere and corridor are symmetric
    // about the Sun's z plane. Verify the real balanced mesh and active plan,
    // not just that the enclosing box has symmetric scalar bounds.
    const double tolerance = 1e-12 * (d.maximumM.z - d.minimumM.z);
    for (std::size_t i = 0; i < mesh.leaves().size(); ++i) {
      const auto& leaf = mesh.leaves()[i];
      const C::Vec3 mirrorMin(leaf.minimumM.x, leaf.minimumM.y, -leaf.maximumM.z);
      const C::Vec3 mirrorMax(leaf.maximumM.x, leaf.maximumM.y, -leaf.minimumM.z);
      bool found = false;
      for (std::size_t j = 0; j < mesh.leaves().size(); ++j) {
        const auto& mirror = mesh.leaves()[j];
        if ((mirror.minimumM - mirrorMin).Norm() <= tolerance &&
            (mirror.maximumM - mirrorMax).Norm() <= tolerance) {
          Require(plan.leafClass[i] == plan.leafClass[j],
                  "active solar neighbourhood/corridor is asymmetric in z");
          found = true;
          break;
        }
      }
      Require(found, "balanced mesh lacks a reflected z leaf");
    }
  }
  const auto coreCovers = [&](const C::Vec3& p) {
    for (std::size_t i = 0; i < mesh.leaves().size(); ++i)
      if (plan.leafClass[i] == M::ActiveLeafClass::Core &&
          Inside(p, mesh.leaves()[i].minimumM, mesh.leaves()[i].maximumM)) return true;
    return false;
  };
  for (int i = 0; i <= 1000; ++i) {
    const double radius = c.innerRadiusM + (c.outerRadiusM-c.innerRadiusM)*i/1000.0;
    Require(coreCovers(C::ParkerCurvePoint(radius, Geometry(o))), "finite line crosses an inactive leaf");
  }
  for (int i = 0; i < 512; ++i) {
    const double z = 1.0 - 2.0*(i+0.5)/512.0;
    const double phi = i * 2.399963229728653;
    const double transverse = std::sqrt(1.0-z*z);
    const C::Vec3 direction(transverse*std::cos(phi), transverse*std::sin(phi), z);
    Require(coreCovers(c.activeSolarSphereRadiusM*direction), "near-Sun sphere has an allocation hole");
  }
  M::RefinementPreflight preflight;
  Require(M::BuildRefinementPreflight(d, c, R::StorageLayout(), &preflight).ok(), "corner preflight failed");
  Require(preflight.minimumRequestedCellM == c.solarSurfaceCellSizeM,
          "preflight misses photospheric resolution");
  Require(std::fabs(preflight.minimumLocationM.Norm()-C::Const::R_sun) < 1e-10*C::Const::R_sun,
          "preflight still anchors at the source shell");
  std::cout << "leaves=" << mesh.leaves().size() << " core=" << plan.coreLeafCount
            << " inactive=" << plan.inactiveLeafCount << '\n';
}

void XYCornerCenteredZ() {
  auto o = Options();
  o.domainBoxGeometry = R::DomainBoxGeometry::FieldLineXYCornerCube;
  o.coordinateOriginM = C::Vec3(0.03, -0.04, 0.02) * C::Const::AU;
  // Check automatic and all four x/y corners, including a tilted line and
  // both polar directions that force enlargement of a centered z extent.
  for (int shape = 0; shape < 4; ++shape) {
    o.tubeColatitudeRad = shape == 2 ? 0.0 :
        (shape == 3 ? C::Const::kPi : 0.5 * C::Const::kPi);
    o.parker.rotationAxis = shape == 1 ? C::Vec3(1.0, 2.0, 3.0) : C::Vec3(0.0, 0.0, 1.0);
    for (int corner = 0; corner < 5; ++corner) {
      o.domainCornerDirection = corner == 4 ? C::Vec3() :
          C::Vec3(corner & 1 ? 1.0 : -1.0, corner & 2 ? 1.0 : -1.0, 0.0);
      const auto d = M::MakeDomain(o);
      const auto sun = M::MakeSolarBoundary(o);
      const double side = d.maximumM.x - d.minimumM.x;
      Require((sun.centerM - o.coordinateOriginM).Norm() == 0.0,
              "solar boundary moved away from the heliocentric origin");
      Require(std::fabs(d.minimumM.z + d.maximumM.z - 2.0 * sun.centerM.z) < 1e-12 * side,
              "Sun is not on the root z midplane");
      Require(std::fabs(side - (d.maximumM.z - d.minimumM.z)) < 1e-12 * side &&
              std::fabs(side - (d.maximumM.y - d.minimumM.y)) < 1e-12 * side,
              "x-y corner root is not cubic");
      Require(Inside(sun.centerM, d.minimumM, d.maximumM,
                     o.activeSolarSphereRadiusM + o.domainCornerMarginM),
              "centered-z layout clips the complete solar neighbourhood");
      const auto geometry = Geometry(o);
      for (int i = 0; i <= 2000; ++i) {
        const double radius = o.innerRadiusM + (o.outerRadiusM - o.innerRadiusM) * i / 2000.0;
        const auto p = o.coordinateOriginM + C::ParkerCurvePoint(radius, geometry);
        Require(Inside(p, d.minimumM, d.maximumM,
                       o.activeTubeRadiusAtReferenceM * radius / C::Const::AU),
                "centered-z layout clips a tilted/polar corridor cross-section");
      }
      if (shape >= 2)
        Require(side >= 2.0 * o.outerRadiusM,
                "polar corridor failed to enlarge the symmetric z span");
    }
  }
  o.domainCornerDirection.z = 1.0;
  C::Vec3 a, b;
  Require(!R::ResolveDomainBoundsM(o, &a, &b).ok(),
          "ambiguous z corner selection was accepted in centered-z mode");
  MeshAllocation(true);
}

int main(int argc, char** argv) {
  try {
    Require(argc == 2, "expected one test ID");
    const std::string id = argv[1];
    if (id == "DOM3D01") DomainBounds();
    else if (id == "DOM3D02") SphereAndResolution();
    else if (id == "DOM3D03") MeshAllocation();
    else if (id == "DOM3D04") XYCornerCenteredZ();
    else throw std::runtime_error("unknown geometry test");
    std::cout << '[' << id << "] PASS\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
