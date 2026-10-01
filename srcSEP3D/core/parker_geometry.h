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
// longitude/colatitude.  The returned vector has norm radiusM.  The azimuth
// is the analytic integral of the same B_phi/B_r law used by SWCME and by
// ParkerLocalTangent(); this is important because a merely Archimedean curve
// is not an integral curve when the field contains the usual (r-r0) source-
// surface correction.
Vec3 ParkerCurvePoint(double radiusM, const ParkerSpiralGeometry& geometry);

// Arc length measured outwards from the source sphere to radiusM.  Invalid
// geometry or a radius inside the source sphere returns NaN.  Keeping this
// operation in the geometry authority lets mesh masks and visualization use
// exact equal-arc stations without integrating a second, drifting curve.
double ParkerCurveArcLengthM(
    double radiusM, const ParkerSpiralGeometry& geometry);

// Invert ParkerCurveArcLengthM on the monotonically increasing outward
// branch.  A safeguarded Newton/bisection solve is used so zero rotation,
// polar field lines, and tightly wound equatorial lines share one contract.
// The result is committed only on success.
Status ParkerCurveRadiusAtArcLengthM(
    double arcLengthM, const ParkerSpiralGeometry& geometry,
    double* radiusM);

// Unit tangent in the increasing-radius direction.  This direction is the
// positive-polarity Parker field direction; callers apply magnetic polarity
// only after this geometric quantity has been computed.
Vec3 ParkerCurveTangent(double radiusM,
                        const ParkerSpiralGeometry& geometry);

// Certified enclosure of the complete outward curve and a surrounding ball
// of radius paddingM. Bounds use the analytic |dx/dr| on each radial interval,
// not the diagnostic point count, so a bend between plotted points is retained.
Status ParkerCurveBoundsM(double endRadiusM,
                          const ParkerSpiralGeometry& geometry,
                          double paddingM, Vec3* minimumM, Vec3* maximumM);

// Unit Parker tangent through an arbitrary heliocentric position.  On the
// configured centreline it is identical to ParkerCurveTangent().
Vec3 ParkerLocalTangent(const Vec3& positionM,
                        const ParkerSpiralGeometry& geometry);

}  // namespace Core
}  // namespace SEP3D

#endif  // SEP3D_CORE_PARKER_GEOMETRY_H
