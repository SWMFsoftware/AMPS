#include "parker_geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Core {
namespace {

Vec3 Rotate(const Vec3& vector, const Vec3& unitAxis, double angle) {
  // Rodrigues' formula works for arbitrary rotation-axis orientation and
  // avoids a hidden assumption that the solar axis is exactly the mesh Z axis.
  return std::cos(angle) * vector +
         std::sin(angle) * unitAxis.Cross(vector) +
         (1.0 - std::cos(angle)) * unitAxis.Dot(vector) * unitAxis;
}

Vec3 SourceDirection(const ParkerSpiralGeometry& geometry) {
  const double sine = std::sin(geometry.sourceColatitudeRad);
  return {sine * std::cos(geometry.sourceLongitudeRad),
          sine * std::sin(geometry.sourceLongitudeRad),
          std::cos(geometry.sourceColatitudeRad)};
}

double TransverseWindingPerM(const ParkerSpiralGeometry& geometry) {
  // The rotation axis is configurable, so sin(theta) must be measured from
  // that actual axis rather than inferred from the +Z-based input angle.
  // |a_hat x r_hat| is invariant under the subsequent Parker rotation.
  const Vec3 axis = geometry.rotationAxis.Normalized();
  const double sine = axis.Cross(SourceDirection(geometry).Normalized()).Norm();
  return std::fabs(geometry.solarRotationRateRadPerS) * sine /
         geometry.solarWindSpeedMPerS;
}

double ArcLengthFromRadialOffset(
    double radialOffsetM, const ParkerSpiralGeometry& geometry) {
  const double windingPerM = TransverseWindingPerM(geometry);
  if (windingPerM == 0.0) return radialOffsetM;

  // ds/dr = sqrt(1 + [Omega sin(theta) (r-r0) / V]^2).  Its
  // antiderivative is evaluated with asinh, which remains stable for both
  // weakly and strongly wound lines and has no subtraction cancellation.
  const double scaled = windingPerM * radialOffsetM;
  return 0.5 * (radialOffsetM * std::sqrt(1.0 + scaled * scaled) +
                std::asinh(scaled) / windingPerM);
}

}  // namespace

Status ValidateParkerGeometry(const ParkerSpiralGeometry& geometry) {
  const double values[] = {
      geometry.sourceRadiusM, geometry.sourceLongitudeRad,
      geometry.sourceColatitudeRad, geometry.solarWindSpeedMPerS,
      geometry.solarRotationRateRadPerS, geometry.rotationAxis.x,
      geometry.rotationAxis.y, geometry.rotationAxis.z};
  for (double value : values) {
    if (!std::isfinite(value))
      return Status(StatusCode::InvalidInput,
                    "Parker geometry contains a non-finite value");
  }
  if (geometry.sourceRadiusM <= 0.0 || geometry.solarWindSpeedMPerS <= 0.0 ||
      geometry.sourceColatitudeRad < 0.0 ||
      geometry.sourceColatitudeRad > Const::kPi ||
      geometry.rotationAxis.Norm() <= 0.0) {
    return Status(StatusCode::InvalidInput,
                  "Parker geometry is outside its physical range");
  }
  return Status::OK();
}

Status ParkerCurveBoundsM(double endRadiusM,
                          const ParkerSpiralGeometry& geometry,
                          double paddingM, Vec3* minimumM, Vec3* maximumM) {
  const Status valid = ValidateParkerGeometry(geometry);
  if (!valid.ok()) return valid;
  if (!minimumM || !maximumM || !std::isfinite(endRadiusM) ||
      endRadiusM < geometry.sourceRadiusM || !std::isfinite(paddingM) ||
      paddingM < 0.0)
    return Status(StatusCode::InvalidInput, "invalid Parker curve enclosure");
  Vec3 lower(std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity());
  Vec3 upper = -1.0 * lower;
  const double winding = TransverseWindingPerM(geometry);
  constexpr int intervals = 512;
  for (int i = 0; i < intervals; ++i) {
    const double a = geometry.sourceRadiusM +
        (endRadiusM - geometry.sourceRadiusM) * i / intervals;
    const double b = geometry.sourceRadiusM +
        (endRadiusM - geometry.sourceRadiusM) * (i + 1) / intervals;
    const Vec3 middle = ParkerCurvePoint(0.5 * (a + b), geometry);
    // The speed in radial coordinates is monotone on the outward branch.
    // The mean-value inequality bounds every point by this midpoint ball.
    const double speed = std::hypot(1.0,
        winding * (b - geometry.sourceRadiusM));
    const double roundoff = 128.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, endRadiusM);
    const double width = paddingM + 0.5 * (b - a) * speed + roundoff;
    if (!std::isfinite(width))
      return Status(StatusCode::InvalidInput, "Parker enclosure overflow");
    lower.x = std::min(lower.x, middle.x - width);
    lower.y = std::min(lower.y, middle.y - width);
    lower.z = std::min(lower.z, middle.z - width);
    upper.x = std::max(upper.x, middle.x + width);
    upper.y = std::max(upper.y, middle.y + width);
    upper.z = std::max(upper.z, middle.z + width);
  }
  *minimumM = lower;
  *maximumM = upper;
  return Status::OK();
}

Vec3 ParkerCurvePoint(double radiusM,
                      const ParkerSpiralGeometry& geometry) {
  if (!ValidateParkerGeometry(geometry).ok() || !std::isfinite(radiusM) ||
      radiusM <= 0.0) return {};
  const Vec3 axis = geometry.rotationAxis.Normalized();
  const double travel = std::max(0.0, radiusM - geometry.sourceRadiusM);

  // SWCME initializes B_phi/B_r = -Omega (r-r0) sin(theta) / V.  A field
  // line therefore satisfies dphi/dr = -Omega (r-r0)/(V r), whose exact
  // integral is the expression below.  The former -Omega(r-r0)/V formula
  // described a different Archimedean curve and made the refinement/mask
  // centreline diverge from both the exported line and the initialized IMF.
  const double radialIntegral = travel - geometry.sourceRadiusM *
      std::log1p(travel / geometry.sourceRadiusM);
  const double angle = -geometry.solarRotationRateRadPerS * radialIntegral /
                       geometry.solarWindSpeedMPerS;
  return radiusM * Rotate(SourceDirection(geometry), axis, angle).Normalized();
}

double ParkerCurveArcLengthM(
    double radiusM, const ParkerSpiralGeometry& geometry) {
  if (!ValidateParkerGeometry(geometry).ok() || !std::isfinite(radiusM) ||
      radiusM < geometry.sourceRadiusM) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return ArcLengthFromRadialOffset(radiusM - geometry.sourceRadiusM,
                                   geometry);
}

Status ParkerCurveRadiusAtArcLengthM(
    double arcLengthM, const ParkerSpiralGeometry& geometry,
    double* radiusM) {
  if (radiusM == nullptr) {
    return Status(StatusCode::InvalidInput,
                  "Parker radius output pointer is null");
  }
  const Status valid = ValidateParkerGeometry(geometry);
  if (!valid.ok()) return valid;
  if (!std::isfinite(arcLengthM) || arcLengthM < 0.0) {
    return Status(StatusCode::InvalidInput,
                  "Parker arc length must be finite and non-negative");
  }
  if (arcLengthM == 0.0) {
    *radiusM = geometry.sourceRadiusM;
    return Status::OK();
  }
  if (!std::isfinite(geometry.sourceRadiusM + arcLengthM)) {
    return Status(StatusCode::InvalidInput,
                  "Parker arc length overflows the radial coordinate");
  }

  // ds/dr >= 1, hence the radial offset is always in [0, arcLengthM].
  // Newton gives rapid convergence in normal heliospheric configurations;
  // the bracket is updated first and rejects any Newton step that could leave
  // the monotone interval, making the solve deterministic even at the poles.
  const double windingPerM = TransverseWindingPerM(geometry);
  double lower = 0.0;
  double upper = arcLengthM;
  double offset = 0.5 * (lower + upper);
  for (int iteration = 0; iteration < 96; ++iteration) {
    const double residual =
        ArcLengthFromRadialOffset(offset, geometry) - arcLengthM;
    if (residual > 0.0) upper = offset;
    else lower = offset;

    const double scaled = windingPerM * offset;
    const double derivative = std::sqrt(1.0 + scaled * scaled);
    const double newton = offset - residual / derivative;
    offset = (std::isfinite(newton) && newton > lower && newton < upper)
        ? newton : 0.5 * (lower + upper);

    const double scale = std::max(
        geometry.sourceRadiusM, geometry.sourceRadiusM + upper);
    if (upper - lower <=
        16.0 * std::numeric_limits<double>::epsilon() * scale) {
      offset = 0.5 * (lower + upper);
      break;
    }
  }

  const double candidate = geometry.sourceRadiusM + offset;
  if (!std::isfinite(candidate)) {
    return Status(StatusCode::InvalidInput,
                  "Parker arc-length inversion produced a non-finite radius");
  }
  *radiusM = candidate;
  return Status::OK();
}

Vec3 ParkerLocalTangent(const Vec3& positionM,
                        const ParkerSpiralGeometry& geometry) {
  const double radius = positionM.Norm();
  if (!ValidateParkerGeometry(geometry).ok() || !std::isfinite(radius) ||
      radius <= 0.0) return {};
  const Vec3 radial = positionM / radius;
  const Vec3 axis = geometry.rotationAxis.Normalized();
  const double spiral = geometry.solarRotationRateRadPerS *
      std::max(0.0, radius - geometry.sourceRadiusM) /
      geometry.solarWindSpeedMPerS;
  return (radial - spiral * axis.Cross(radial)).Normalized();
}

Vec3 ParkerCurveTangent(double radiusM,
                        const ParkerSpiralGeometry& geometry) {
  return ParkerLocalTangent(ParkerCurvePoint(radiusM, geometry), geometry);
}

}  // namespace Core
}  // namespace SEP3D
