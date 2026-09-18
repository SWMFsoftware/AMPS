#include "parker_geometry.h"

#include <algorithm>
#include <cmath>

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

Vec3 ParkerCurvePoint(double radiusM,
                      const ParkerSpiralGeometry& geometry) {
  if (!ValidateParkerGeometry(geometry).ok() || !std::isfinite(radiusM) ||
      radiusM <= 0.0) return {};
  const Vec3 axis = geometry.rotationAxis.Normalized();
  const double travel = std::max(0.0, radiusM - geometry.sourceRadiusM);
  const double angle = -geometry.solarRotationRateRadPerS * travel /
                       geometry.solarWindSpeedMPerS;
  return radiusM * Rotate(SourceDirection(geometry), axis, angle).Normalized();
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
