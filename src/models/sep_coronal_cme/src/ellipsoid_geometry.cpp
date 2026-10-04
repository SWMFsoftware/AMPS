#include "sep_coronal_cme/ellipsoid_geometry.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronalCME { namespace {

bool Finite(double value) { return std::isfinite(value); }

bool FiniteKinematics(const KinematicValue& value) {
  return Finite(value.value) && Finite(value.firstDerivative) &&
      Finite(value.secondDerivative);
}

double SmoothRate(double tau) {
  return 10.0 * std::pow(tau, 3) - 15.0 * std::pow(tau, 4) +
      6.0 * std::pow(tau, 5);
}

double SmoothIntegral(double tau) {
  return 2.5 * std::pow(tau, 4) - 3.0 * std::pow(tau, 5) +
      std::pow(tau, 6);
}

Vec3 LocalToGlobal(const RadialPrincipalBasis& basis,
                   double radial, double first, double second) {
  return radial * basis.radial + first * basis.firstLateral +
      second * basis.secondLateral;
}

bool SameBasis(const RadialPrincipalBasis& a,
               const RadialPrincipalBasis& b) {
  return Norm(a.radial - b.radial) < 1.0e-12 &&
      Norm(a.firstLateral - b.firstLateral) < 1.0e-12 &&
      Norm(a.secondLateral - b.secondLateral) < 1.0e-12;
}

bool ValidBasis(const RadialPrincipalBasis& basis) {
  // Factories are public and may receive a manually assembled record rather
  // than BuildRadialPrincipalBasis's output.  Recheck SO(3) membership here so
  // Q remains a physical symmetric-positive-definite ellipsoid tensor.
  return std::abs(Norm(basis.radial) - 1.0) < 1.0e-12 &&
      std::abs(Norm(basis.firstLateral) - 1.0) < 1.0e-12 &&
      std::abs(Norm(basis.secondLateral) - 1.0) < 1.0e-12 &&
      std::abs(Dot(basis.radial, basis.firstLateral)) < 1.0e-12 &&
      std::abs(Dot(basis.radial, basis.secondLateral)) < 1.0e-12 &&
      std::abs(Dot(basis.firstLateral, basis.secondLateral)) < 1.0e-12 &&
      Dot(Cross(basis.radial, basis.firstLateral), basis.secondLateral) >
          1.0 - 1.0e-12;
}

}  // namespace

Core::Result<KinematicValue> EvaluateSmoothKinematics(
    double timeS, double referenceTimeS, double qReference,
    double transitionStartS, double transitionDurationS,
    double initialRatePerS, double finalRatePerS) {
  if (!(Finite(timeS) && Finite(referenceTimeS) && Finite(qReference) &&
        Finite(transitionStartS) && transitionDurationS > 0.0 &&
        Finite(initialRatePerS) && Finite(finalRatePerS))) {
    return Core::Result<KinematicValue>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "smooth kinematics requires finite values and transition duration>0");
  }
  const double startValue = qReference + initialRatePerS *
      (transitionStartS - referenceTimeS);
  if (timeS < transitionStartS) {
    return Core::Result<KinematicValue>::Success(
        {qReference + initialRatePerS * (timeS - referenceTimeS),
         initialRatePerS, 0.0});
  }
  const double endTime = transitionStartS + transitionDurationS;
  const double deltaRate = finalRatePerS - initialRatePerS;
  if (timeS <= endTime) {
    const double elapsed = timeS - transitionStartS;
    const double tau = elapsed / transitionDurationS;
    const double acceleration = deltaRate / transitionDurationS *
        30.0 * tau * tau * std::pow(1.0 - tau, 2);
    return Core::Result<KinematicValue>::Success({
        startValue + initialRatePerS * elapsed +
            deltaRate * transitionDurationS * SmoothIntegral(tau),
        initialRatePerS + deltaRate * SmoothRate(tau), acceleration});
  }
  const double endValue = startValue +
      initialRatePerS * transitionDurationS +
      deltaRate * transitionDurationS * SmoothIntegral(1.0);
  return Core::Result<KinematicValue>::Success(
      {endValue + finalRatePerS * (timeS - endTime), finalRatePerS, 0.0});
}

Core::Result<KinematicValue> EvaluateCubicHermiteHistory(
    const std::vector<HermiteKnot>& knots, double timeS) {
  if (knots.size() < 2 || !Finite(timeS)) {
    return Core::Result<KinematicValue>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "tabulated C1 history requires at least two knots and a finite time");
  }
  for (std::size_t index = 0; index < knots.size(); ++index) {
    if (!(Finite(knots[index].timeS) && Finite(knots[index].value) &&
          Finite(knots[index].ratePerS)) ||
        (index > 0 && knots[index].timeS <= knots[index - 1].timeS)) {
      return Core::Result<KinematicValue>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "tabulated C1 knots must be finite and strictly time ordered");
    }
  }
  if (timeS < knots.front().timeS || timeS > knots.back().timeS) {
    return Core::Result<KinematicValue>::Failure(
        Core::StatusCode::OutOfDomain,
        "tabulated C1 history does not extrapolate");
  }
  std::size_t right = 1;
  while (right < knots.size() && timeS > knots[right].timeS) ++right;
  if (right == knots.size()) right = knots.size() - 1;
  const HermiteKnot& leftKnot = knots[right - 1];
  const HermiteKnot& rightKnot = knots[right];
  const double duration = rightKnot.timeS - leftKnot.timeS;
  const double t = (timeS - leftKnot.timeS) / duration;
  const double h00 = 2.0 * t * t * t - 3.0 * t * t + 1.0;
  const double h10 = t * t * t - 2.0 * t * t + t;
  const double h01 = -2.0 * t * t * t + 3.0 * t * t;
  const double h11 = t * t * t - t * t;
  const double value = h00 * leftKnot.value +
      h10 * duration * leftKnot.ratePerS + h01 * rightKnot.value +
      h11 * duration * rightKnot.ratePerS;
  const double derivative =
      ((6.0 * t * t - 6.0 * t) * leftKnot.value +
       (-6.0 * t * t + 6.0 * t) * rightKnot.value) / duration +
      (3.0 * t * t - 4.0 * t + 1.0) * leftKnot.ratePerS +
      (3.0 * t * t - 2.0 * t) * rightKnot.ratePerS;
  const double second =
      ((12.0 * t - 6.0) * leftKnot.value +
       (-12.0 * t + 6.0) * rightKnot.value) /
          (duration * duration) +
      (6.0 * t - 4.0) * leftKnot.ratePerS / duration +
      (6.0 * t - 2.0) * rightKnot.ratePerS / duration;
  return Core::Result<KinematicValue>::Success({value, derivative, second});
}

Core::Result<RadialPrincipalBasis> BuildRadialPrincipalBasis(
    double latitudeRad, double longitudeRad, double lateralTiltRad) {
  if (!(Finite(latitudeRad) && Finite(longitudeRad) &&
        Finite(lateralTiltRad) &&
        latitudeRad >= -0.5 * Constants::kPi &&
        latitudeRad <= 0.5 * Constants::kPi)) {
    return Core::Result<RadialPrincipalBasis>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "radial basis requires finite longitude/tilt and physical latitude");
  }
  const double cl = std::cos(latitudeRad);
  const double sl = std::sin(latitudeRad);
  const double cp = std::cos(longitudeRad);
  const double sp = std::sin(longitudeRad);
  const Vec3 radial{cl * cp, cl * sp, sl};
  const Vec3 east{-sp, cp, 0.0};
  const Vec3 north{-sl * cp, -sl * sp, cl};
  const Vec3 first = std::cos(lateralTiltRad) * east +
      std::sin(lateralTiltRad) * north;
  const Vec3 second = -std::sin(lateralTiltRad) * east +
      std::cos(lateralTiltRad) * north;
  if (Dot(Cross(radial, first), second) < 1.0 - 1.0e-12) {
    return Core::Result<RadialPrincipalBasis>::Failure(
        Core::StatusCode::InvalidState,
        "radial principal basis is not right handed");
  }
  return Core::Result<RadialPrincipalBasis>::Success(
      {radial, first, second});
}

Core::Result<FixedOrientationEllipsoid> FixedOrientationEllipsoid::FromCenter(
    const RadialPrincipalBasis& basis,
    const EllipsoidKinematics& kinematics, double solarRadiusM) {
  if (!(FiniteKinematics(kinematics.centerDistanceM) &&
        FiniteKinematics(kinematics.radialSemiAxisM) &&
        FiniteKinematics(kinematics.firstLateralSemiAxisM) &&
        FiniteKinematics(kinematics.secondLateralSemiAxisM) &&
        kinematics.centerDistanceM.value > 0.0 &&
        kinematics.radialSemiAxisM.value > 0.0 &&
        kinematics.firstLateralSemiAxisM.value > 0.0 &&
        kinematics.secondLateralSemiAxisM.value > 0.0 &&
        solarRadiusM >= 0.0 && ValidBasis(basis))) {
    return Core::Result<FixedOrientationEllipsoid>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "ellipsoid requires positive center/semiaxes and finite kinematics");
  }
  FixedOrientationEllipsoid result;
  result.basis_ = basis;
  result.kinematics_ = kinematics;
  result.solarRadiusM_ = solarRadiusM;
  return Core::Result<FixedOrientationEllipsoid>::Success(result);
}

Core::Result<FixedOrientationEllipsoid> FixedOrientationEllipsoid::FromApex(
    const RadialPrincipalBasis& basis, KinematicValue apexRadiusM,
    KinematicValue radialSemiAxisM,
    KinematicValue firstLateralSemiAxisM,
    KinematicValue secondLateralSemiAxisM, double solarRadiusM) {
  KinematicValue center;
  center.value = apexRadiusM.value - radialSemiAxisM.value;
  center.firstDerivative = apexRadiusM.firstDerivative -
      radialSemiAxisM.firstDerivative;
  center.secondDerivative = apexRadiusM.secondDerivative -
      radialSemiAxisM.secondDerivative;
  return FromCenter(basis, {center, radialSemiAxisM,
      firstLateralSemiAxisM, secondLateralSemiAxisM}, solarRadiusM);
}

Vec3 FixedOrientationEllipsoid::Point(
    double polarParameterRad, double azimuthParameterRad) const {
  const double sinPolar = std::sin(polarParameterRad);
  const double localRadial = kinematics_.radialSemiAxisM.value *
      std::cos(polarParameterRad);
  const double localFirst = kinematics_.firstLateralSemiAxisM.value *
      sinPolar * std::cos(azimuthParameterRad);
  const double localSecond = kinematics_.secondLateralSemiAxisM.value *
      sinPolar * std::sin(azimuthParameterRad);
  return kinematics_.centerDistanceM.value * basis_.radial +
      LocalToGlobal(basis_, localRadial, localFirst, localSecond);
}

Vec3 FixedOrientationEllipsoid::SurfaceVelocity(
    double polarParameterRad, double azimuthParameterRad) const {
  const double sinPolar = std::sin(polarParameterRad);
  return kinematics_.centerDistanceM.firstDerivative * basis_.radial +
      LocalToGlobal(basis_,
          kinematics_.radialSemiAxisM.firstDerivative *
              std::cos(polarParameterRad),
          kinematics_.firstLateralSemiAxisM.firstDerivative *
              sinPolar * std::cos(azimuthParameterRad),
          kinematics_.secondLateralSemiAxisM.firstDerivative *
              sinPolar * std::sin(azimuthParameterRad));
}

Core::Result<SurfaceEvaluation> FixedOrientationEllipsoid::Evaluate(
    Vec3 positionM) const {
  const Vec3 center = kinematics_.centerDistanceM.value * basis_.radial;
  const Vec3 relative = positionM - center;
  const double coordinates[3] = {
      Dot(relative, basis_.radial),
      Dot(relative, basis_.firstLateral),
      Dot(relative, basis_.secondLateral)};
  const double axes[3] = {
      kinematics_.radialSemiAxisM.value,
      kinematics_.firstLateralSemiAxisM.value,
      kinematics_.secondLateralSemiAxisM.value};
  const double rates[3] = {
      kinematics_.radialSemiAxisM.firstDerivative,
      kinematics_.firstLateralSemiAxisM.firstDerivative,
      kinematics_.secondLateralSemiAxisM.firstDerivative};
  const Vec3 basis[3] = {basis_.radial, basis_.firstLateral,
                         basis_.secondLateral};
  double implicit = -1.0;
  Vec3 qRelative;
  double q2 = 0.0;
  double q3 = 0.0;
  double traceQ = 0.0;
  double axisTimeTerm = 0.0;
  for (int index = 0; index < 3; ++index) {
    const double inverseAxisSquared = 1.0 / (axes[index] * axes[index]);
    implicit += coordinates[index] * coordinates[index] *
        inverseAxisSquared;
    qRelative = qRelative +
        coordinates[index] * inverseAxisSquared * basis[index];
    q2 += coordinates[index] * coordinates[index] *
        inverseAxisSquared * inverseAxisSquared;
    q3 += coordinates[index] * coordinates[index] *
        inverseAxisSquared * inverseAxisSquared * inverseAxisSquared;
    traceQ += inverseAxisSquared;
    axisTimeTerm += -2.0 * rates[index] *
        coordinates[index] * coordinates[index] /
        (axes[index] * axes[index] * axes[index]);
  }
  const double gradientMagnitude = 2.0 * std::sqrt(q2);
  if (!(gradientMagnitude > 0.0)) {
    return Core::Result<SurfaceEvaluation>::Failure(
        Core::StatusCode::InvalidState,
        "ellipsoid normal is undefined at its center");
  }
  const Vec3 centerVelocity =
      kinematics_.centerDistanceM.firstDerivative * basis_.radial;
  const double timeDerivative = -2.0 * Dot(centerVelocity, qRelative) +
      axisTimeTerm;
  SurfaceEvaluation result;
  result.implicitValue = implicit;
  result.outwardNormal = Unit(qRelative);
  result.normalSpeedMPerS = -timeDerivative / gradientMagnitude;
  result.meanCurvaturePerM = (q2 * traceQ - q3) /
      (2.0 * std::pow(q2, 1.5));
  result.gaussianCurvaturePerM2 = 1.0 /
      (axes[0] * axes[0] * axes[1] * axes[1] *
       axes[2] * axes[2] * q2 * q2);
  return Core::Result<SurfaceEvaluation>::Success(result);
}

Core::Result<std::vector<SurfacePatch>>
FixedOrientationEllipsoid::Tessellate(
    int polarCells, int azimuthCells) const {
  if (polarCells < 2 || azimuthCells < 3) {
    return Core::Result<std::vector<SurfacePatch>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "ellipsoid tessellation requires >=2 polar and >=3 azimuth cells");
  }
  const double dPolar = Constants::kPi / polarCells;
  const double dAzimuth = 2.0 * Constants::kPi / azimuthCells;
  std::vector<SurfacePatch> patches;
  for (int polarIndex = 0; polarIndex < polarCells; ++polarIndex) {
    const double polar = (polarIndex + 0.5) * dPolar;
    const double sinPolar = std::sin(polar);
    const double cosPolar = std::cos(polar);
    for (int azimuthIndex = 0; azimuthIndex < azimuthCells;
         ++azimuthIndex) {
      const double azimuth = (azimuthIndex + 0.5) * dAzimuth;
      const Vec3 center = Point(polar, azimuth);
      if (solarRadiusM_ > 0.0 && Norm(center) < solarRadiusM_) continue;
      const Vec3 dPolarVector = LocalToGlobal(basis_,
          -kinematics_.radialSemiAxisM.value * sinPolar,
          kinematics_.firstLateralSemiAxisM.value * cosPolar *
              std::cos(azimuth),
          kinematics_.secondLateralSemiAxisM.value * cosPolar *
              std::sin(azimuth));
      const Vec3 dAzimuthVector = LocalToGlobal(basis_, 0.0,
          -kinematics_.firstLateralSemiAxisM.value * sinPolar *
              std::sin(azimuth),
          kinematics_.secondLateralSemiAxisM.value * sinPolar *
              std::cos(azimuth));
      auto evaluation = Evaluate(center);
      if (!evaluation.ok()) {
        return Core::Result<std::vector<SurfacePatch>>::Failure(
            evaluation.status.code, evaluation.status.message);
      }
      patches.push_back({
          1U + static_cast<std::uint64_t>(polarIndex) *
              static_cast<std::uint64_t>(azimuthCells) +
              static_cast<std::uint64_t>(azimuthIndex),
          center, evaluation.value.outwardNormal,
          Norm(Cross(dPolarVector, dAzimuthVector)) * dPolar * dAzimuth,
          polar, azimuth});
    }
  }
  if (patches.empty()) {
    return Core::Result<std::vector<SurfacePatch>>::Failure(
        Core::StatusCode::InvalidState,
        "solar clipping removed the complete candidate surface");
  }
  return Core::Result<std::vector<SurfacePatch>>::Success(std::move(patches));
}

std::vector<std::uint64_t> PatchSpatialIndex::Query(
    Vec3 pointM, double radiusM) const {
  std::vector<std::uint64_t> result;
  if (!(radiusM >= 0.0 && Finite(radiusM))) return result;
  for (const auto& patch : patches_) {
    if (Norm(patch.centerM - pointM) <= radiusM) {
      result.push_back(patch.physicalId);
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

Core::Result<PistonNestingDiagnostics> CheckPistonNesting(
    const FixedOrientationEllipsoid& front,
    const FixedOrientationEllipsoid& piston,
    double requiredMinimumSeparationM, int polarSamples,
    int azimuthSamples) {
  if (!(requiredMinimumSeparationM >= 0.0 && polarSamples >= 2 &&
        azimuthSamples >= 3 && SameBasis(front.Basis(), piston.Basis()))) {
    return Core::Result<PistonNestingDiagnostics>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "piston nesting requires a common fixed basis and valid sampling");
  }
  PistonNestingDiagnostics result;
  result.minimumRadialSeparationM = std::numeric_limits<double>::infinity();
  for (int i = 0; i < polarSamples; ++i) {
    const double polar = (i + 0.5) * Constants::kPi / polarSamples;
    for (int j = 0; j < azimuthSamples; ++j) {
      const double azimuth = (j + 0.5) * 2.0 * Constants::kPi /
          azimuthSamples;
      const Vec3 pistonPoint = piston.Point(polar, azimuth);
      auto inside = front.Evaluate(pistonPoint);
      if (!inside.ok() || inside.value.implicitValue > 1.0e-10) {
        return Core::Result<PistonNestingDiagnostics>::Failure(
            Core::StatusCode::InvalidState,
            "enabled piston is not inside the candidate front");
      }
      // Compare both surfaces on the same heliocentric ray.  Equal ellipsoid
      // parameters are not generally collinear with the Sun when the center
      // is displaced, so subtracting two parameterized radii would implement
      // the wrong standoff definition.
      const Vec3 ray = Unit(pistonPoint);
      const auto& fk = front.Kinematics();
      const auto& fb = front.Basis();
      const double direction[3] = {Dot(ray, fb.radial),
                                   Dot(ray, fb.firstLateral),
                                   Dot(ray, fb.secondLateral)};
      const double center[3] = {fk.centerDistanceM.value, 0.0, 0.0};
      const double axes[3] = {fk.radialSemiAxisM.value,
                              fk.firstLateralSemiAxisM.value,
                              fk.secondLateralSemiAxisM.value};
      double quadratic = 0.0, linear = 0.0, constant = -1.0;
      for (int component = 0; component < 3; ++component) {
        const double inverseAxisSquared = 1.0 /
            (axes[component] * axes[component]);
        quadratic += direction[component] * direction[component] *
            inverseAxisSquared;
        linear += -2.0 * direction[component] * center[component] *
            inverseAxisSquared;
        constant += center[component] * center[component] *
            inverseAxisSquared;
      }
      const double discriminant = linear * linear -
          4.0 * quadratic * constant;
      if (!(quadratic > 0.0 && discriminant >= 0.0)) {
        return Core::Result<PistonNestingDiagnostics>::Failure(
            Core::StatusCode::InvalidState,
            "piston ray does not intersect the candidate front");
      }
      const double frontRadius = (-linear + std::sqrt(discriminant)) /
          (2.0 * quadratic);
      const double separation = frontRadius - Norm(pistonPoint);
      if (separation < result.minimumRadialSeparationM) {
        result.minimumRadialSeparationM = separation;
        result.limitingPatchId = static_cast<std::uint64_t>(i) *
            static_cast<std::uint64_t>(azimuthSamples) +
            static_cast<std::uint64_t>(j);
      }
    }
  }
  if (result.minimumRadialSeparationM < requiredMinimumSeparationM) {
    return Core::Result<PistonNestingDiagnostics>::Failure(
        Core::StatusCode::InvalidState,
        "piston/front separation is below the configured minimum");
  }
  return Core::Result<PistonNestingDiagnostics>::Success(result);
}

} }  // namespace SEP::CoronalCME
