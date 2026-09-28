#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/ellipsoid_geometry.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a, double b, double relative = 1.0e-8,
           double absolute = 1.0e-10) {
  return std::abs(a - b) <= absolute +
      relative * std::max(std::abs(a), std::abs(b));
}

RadialPrincipalBasis Basis() {
  auto basis = BuildRadialPrincipalBasis(0.2, 0.4, 0.3);
  Require(basis.ok(), basis.status.message);
  return basis.value;
}

FixedOrientationEllipsoid Ellipsoid(
    double center = 5.0, double ar = 2.0, double a1 = 3.0,
    double a2 = 4.0, double centerRate = 0.7, double arRate = 0.2,
    double a1Rate = 0.3, double a2Rate = 0.4,
    double solarRadius = 0.0) {
  EllipsoidKinematics k{{center, centerRate, 0.0},
                        {ar, arRate, 0.0},
                        {a1, a1Rate, 0.0},
                        {a2, a2Rate, 0.0}};
  auto ellipsoid = FixedOrientationEllipsoid::FromCenter(
      Basis(), k, solarRadius);
  Require(ellipsoid.ok(), ellipsoid.status.message);
  return ellipsoid.value;
}

void ELL3D01() {
  auto sphere = Ellipsoid(5.0, 2.0, 2.0, 2.0);
  const Vec3 point = sphere.Point(0.8, 1.2);
  auto evaluated = sphere.Evaluate(point);
  const Vec3 center = 5.0 * sphere.Basis().radial;
  Require(evaluated.ok() && Close(evaluated.value.implicitValue, 0.0) &&
          Close(Norm(evaluated.value.outwardNormal - Unit(point - center)),
                0.0),
          "sphere implicit surface/normal is incorrect");
  auto triaxial = Ellipsoid();
  auto triEval = triaxial.Evaluate(triaxial.Point(1.1, 0.7));
  Require(triEval.ok() && Close(triEval.value.implicitValue, 0.0),
          "triaxial point does not satisfy its implicit surface");
}

void ELL3D02() {
  auto ellipsoid = Ellipsoid();
  const double polar = 1.0, azimuth = 0.6;
  const Vec3 point = ellipsoid.Point(polar, azimuth);
  auto evaluation = ellipsoid.Evaluate(point);
  const double projected = Dot(ellipsoid.SurfaceVelocity(polar, azimuth),
                               evaluation.value.outwardNormal);
  Require(evaluation.ok() &&
          Close(evaluation.value.normalSpeedMPerS, projected),
          "implicit Q-dot normal speed disagrees with surface motion");
}

void ELL3D03() {
  auto sphere = Ellipsoid(5.0, 2.0, 2.0, 2.0);
  double previousError = 1.0e100;
  for (int resolution : {8, 16, 32, 64}) {
    auto patches = sphere.Tessellate(resolution, 2 * resolution);
    Require(patches.ok(), patches.status.message);
    const double area = std::accumulate(patches.value.begin(),
        patches.value.end(), 0.0,
        [](double total, const SurfacePatch& patch) {
          return total + patch.areaM2;
        });
    const double error = std::abs(area - 16.0 * Constants::kPi);
    Require(error < previousError, "sphere area quadrature did not converge");
    previousError = error;
  }
  Require(previousError / (16.0 * Constants::kPi) < 2.0e-4,
          "sphere area quadrature missed its analytic value");
}

void ELL3D04() {
  // The back of this ellipsoid lies below the unit solar sphere, forcing a
  // real dome clipping while retaining the outward nose.
  auto dome = Ellipsoid(1.2, 0.8, 0.7, 0.7, 0.0, 0.0, 0.0, 0.0, 1.0);
  auto patches = dome.Tessellate(60, 120);
  Require(patches.ok(), patches.status.message);
  bool hasIdGap = false;
  for (std::size_t index = 0; index < patches.value.size(); ++index) {
    Require(Norm(patches.value[index].centerM) >= 1.0,
            "solar-clipped surface contains a sub-surface patch center");
    if (index > 0 && patches.value[index].physicalId >
        patches.value[index - 1].physicalId + 1) hasIdGap = true;
  }
  Require(hasIdGap || patches.value.front().physicalId != 0 ||
              patches.value.back().physicalId != 60U * 120U - 1U,
          "solar clipping did not preserve skipped physical IDs");
}

void ELL3D05() {
  auto surface = Ellipsoid();
  auto first = surface.Tessellate(12, 24);
  auto second = surface.Tessellate(12, 24);
  Require(first.ok() && second.ok() && first.value.size() == second.value.size(),
          "deterministic tessellation size changed");
  for (std::size_t index = 0; index < first.value.size(); ++index) {
    Require(first.value[index].physicalId == second.value[index].physicalId,
            "physical patch ID changed with traversal/decomposition");
  }
  PatchSpatialIndex spatial(first.value);
  auto nearby = spatial.Query(first.value.front().centerM, 0.0);
  Require(nearby.size() == 1 && nearby.front() == first.value.front().physicalId,
          "surface spatial index lost stable patch identity");
}

void ELL3D06() {
  const auto basis = Basis();
  EllipsoidKinematics k{{5.0, 0.7, 0.1}, {2.0, 0.2, 0.03},
                        {3.0, 0.3, 0.04}, {4.0, 0.4, 0.05}};
  auto center = FixedOrientationEllipsoid::FromCenter(basis, k);
  KinematicValue apex{7.0, 0.9, 0.13};
  auto fromApex = FixedOrientationEllipsoid::FromApex(
      basis, apex, k.radialSemiAxisM, k.firstLateralSemiAxisM,
      k.secondLateralSemiAxisM);
  Require(center.ok() && fromApex.ok() &&
          Close(Norm(center.value.Point(1.0, 0.8) -
                     fromApex.value.Point(1.0, 0.8)), 0.0) &&
          Close(center.value.Kinematics().centerDistanceM.firstDerivative,
                fromApex.value.Kinematics().centerDistanceM.firstDerivative),
          "center/apex parameterizations did not construct one surface");
}

void ELL3D07() {
  auto before = EvaluateSmoothKinematics(1.0, 0.0, 10.0, 2.0, 4.0,
                                         1.0, 3.0);
  auto middle = EvaluateSmoothKinematics(4.0, 0.0, 10.0, 2.0, 4.0,
                                         1.0, 3.0);
  auto after = EvaluateSmoothKinematics(8.0, 0.0, 10.0, 2.0, 4.0,
                                        1.0, 3.0);
  Require(before.ok() && middle.ok() && after.ok() &&
          Close(before.value.firstDerivative, 1.0) &&
          Close(middle.value.firstDerivative, 2.0) &&
          Close(after.value.firstDerivative, 3.0) &&
          middle.value.secondDerivative > 0.0,
          "independent quintic rate history is incorrect");
}

void ELL3D08() {
  const std::vector<HermiteKnot> knots{{0.0, 1.0, 0.5},
                                        {2.0, 3.0, 1.5},
                                        {5.0, 4.0, -0.2}};
  auto knot = EvaluateCubicHermiteHistory(knots, 2.0);
  auto left = EvaluateCubicHermiteHistory(knots, 2.0 - 1.0e-7);
  auto right = EvaluateCubicHermiteHistory(knots, 2.0 + 1.0e-7);
  Require(knot.ok() && left.ok() && right.ok() &&
          Close(knot.value.value, 3.0) &&
          Close(knot.value.firstDerivative, 1.5) &&
          Close(left.value.value, right.value.value, 1.0e-6) &&
          Close(left.value.firstDerivative, right.value.firstDerivative,
                1.0e-6),
          "tabulated Hermite history is not C1 at a knot");
  Require(!EvaluateCubicHermiteHistory(knots, 6.0).ok(),
          "tabulated history extrapolated beyond coverage");
}

void ELL3D09() {
  const auto basis = Basis();
  Require(Close(Dot(Cross(basis.radial, basis.firstLateral),
                    basis.secondLateral), 1.0),
          "radial principal basis is not right handed");
  auto ellipsoid = Ellipsoid();
  const double polar = 1.2, azimuth = 2.1;
  auto evaluation = ellipsoid.Evaluate(ellipsoid.Point(polar, azimuth));
  Require(evaluation.ok() && evaluation.value.meanCurvaturePerM > 0.0 &&
          evaluation.value.gaussianCurvaturePerM2 > 0.0 &&
          Close(evaluation.value.normalSpeedMPerS,
                Dot(ellipsoid.SurfaceVelocity(polar, azimuth),
                    evaluation.value.outwardNormal)),
          "basis/Q-dot/curvature regression failed");
}

void ELL3D10() {
  auto front = Ellipsoid(5.0, 3.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0);
  auto piston = Ellipsoid(5.0, 2.0, 3.0, 3.0, 0.0, 0.0, 0.0, 0.0);
  auto nested = CheckPistonNesting(front, piston, 0.5, 30, 60);
  Require(nested.ok() && nested.value.minimumRadialSeparationM >= 0.5,
          "independently prescribed piston is not nested");
  auto invalid = Ellipsoid(5.0, 3.2, 4.2, 4.2, 0.0, 0.0, 0.0, 0.0);
  Require(!CheckPistonNesting(front, invalid, 0.0, 20, 40).ok(),
          "piston outside the candidate front was accepted");
}

}  // namespace

void RegisterStage5(Registry* tests) {
  (*tests)["ELL3D01"] = ELL3D01; (*tests)["ELL3D02"] = ELL3D02;
  (*tests)["ELL3D03"] = ELL3D03; (*tests)["ELL3D04"] = ELL3D04;
  (*tests)["ELL3D05"] = ELL3D05; (*tests)["ELL3D06"] = ELL3D06;
  (*tests)["ELL3D07"] = ELL3D07; (*tests)["ELL3D08"] = ELL3D08;
  (*tests)["ELL3D09"] = ELL3D09; (*tests)["ELL3D10"] = ELL3D10;
}

}  // namespace SCCMTest
