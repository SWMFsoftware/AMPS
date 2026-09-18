// ============================================================================
// Polarity-independent Parker-spiral geometry
//
// This AMPS-independent helper is the single geometric authority shared by
// mesh refinement and the analytic Parker background.  Magnetic polarity is
// intentionally absent: polarity reverses a vector field, not the spatial
// curve on which the field lies.  Keeping the sign out of this type prevents
// a polarity change from mirroring the high-resolution mesh corridor.
// ============================================================================

#ifndef SEP3D_CORE_PARKER_GEOMETRY_H
#define SEP3D_CORE_PARKER_GEOMETRY_H

#include "sep3d_types.h"

namespace SEP3D {
namespace Core {

struct ParkerSpiralGeometry {
  double sourceRadiusM = 20.0 * Const::R_sun;
  double sourceLongitudeRad = 0.0;
  double sourceColatitudeRad = 0.5 * Const::kPi;
  double solarWindSpeedMPerS = Const::V_sw_default;
  double solarRotationRateRadPerS = Const::Omega_sun;
  Vec3 rotationAxis = {0.0, 0.0, 1.0};
};

Status ValidateParkerGeometry(const ParkerSpiralGeometry& geometry);

// Point on the field line that leaves the source sphere at the configured
// longitude/colatitude.  The returned vector has norm radiusM.
Vec3 ParkerCurvePoint(double radiusM, const ParkerSpiralGeometry& geometry);

// Unit tangent in the increasing-radius direction.  This direction is the
// positive-polarity Parker field direction; callers apply magnetic polarity
// only after this geometric quantity has been computed.
Vec3 ParkerCurveTangent(double radiusM,
                        const ParkerSpiralGeometry& geometry);

// Unit Parker tangent through an arbitrary heliocentric position.  On the
// configured centreline it is identical to ParkerCurveTangent().
Vec3 ParkerLocalTangent(const Vec3& positionM,
                        const ParkerSpiralGeometry& geometry);

}  // namespace Core
}  // namespace SEP3D

#endif  // SEP3D_CORE_PARKER_GEOMETRY_H
