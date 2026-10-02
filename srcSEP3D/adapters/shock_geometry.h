#ifndef SEP3D_ADAPTERS_SHOCK_GEOMETRY_H
#define SEP3D_ADAPTERS_SHOCK_GEOMETRY_H

// Geometry carried across the provider/AMPS boundary. This record contains no
// plasma or acceleration model: SWCME remains the authority for local Mach
// number, Rankine-Hugoniot states and source weights. All lengths are SI.
#include "../core/sep3d_types.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace SEP3D { namespace Adapters {
enum class ShockGeometryKind { Sphere = 0, FiniteSSE = 1 };

struct ExpandingShock {
  // centerM is the heliocentric origin, NOT the translated SSE sphere center.
  // radiusAtStepStartM and radialSpeedMPerS mean apex distance/speed for SSE.
  Core::Vec3 centerM;
  double radiusAtStepStartM = 0.0;
  double radialSpeedMPerS = 0.0;
  std::uint64_t generation = 0;
  bool active = false;
  ShockGeometryKind geometry = ShockGeometryKind::Sphere;
  Core::Vec3 cmeDirection = {1.0, 0.0, 0.0};
  double halfWidthRad = Core::Const::kPi / 2.0;
};

// Preserve source compatibility for existing spherical hosts and tests. New
// code uses ExpandingShock; the old spelling no longer discards finite shape.
using ExpandingSphericalShock = ExpandingShock;

inline Core::Status ValidateShockGeometry(const ExpandingShock& shock) {
  const auto finite=[](const Core::Vec3& v) {
    return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);
  };
  if (!finite(shock.centerM)||!std::isfinite(shock.radiusAtStepStartM)||
      shock.radiusAtStepStartM<=0||!std::isfinite(shock.radialSpeedMPerS))
    return Core::Status(Core::StatusCode::InvalidInput,"invalid shock origin/apex kinematics");
  if (shock.geometry==ShockGeometryKind::Sphere) return Core::Status::OK();
  if (shock.geometry!=ShockGeometryKind::FiniteSSE||!finite(shock.cmeDirection)||
      std::fabs(shock.cmeDirection.Norm()-1.0)>1e-12||
      !std::isfinite(shock.halfWidthRad)||shock.halfWidthRad<=0||
      shock.halfWidthRad>Core::Const::kPi/2.0)
    return Core::Status(Core::StatusCode::InvalidInput,"finite SSE requires a unit axis and half width in (0,pi/2]");
  return Core::Status::OK();
}

struct GeneratingSphere {
  Core::Vec3 centerM, centerVelocityMPerS;
  double radiusM = 0.0, radiusSpeedMPerS = 0.0;
};

// The SSE generating sphere translates AND grows. Its distance from the Sun
// is c=R_apex/(1+sin(lambda)), and its radius is a=c*sin(lambda). Treating it
// as a fixed-center sphere would give wrong particle crossings at the flanks.
inline GeneratingSphere ShockGeneratingSphere(const ExpandingShock& shock) {
  GeneratingSphere out;
  out.centerM=shock.centerM;
  out.radiusM=shock.radiusAtStepStartM;
  out.radiusSpeedMPerS=shock.radialSpeedMPerS;
  if (shock.geometry==ShockGeometryKind::FiniteSSE) {
    const double s=std::sin(shock.halfWidthRad), factor=1.0/(1.0+s);
    out.centerM+=shock.cmeDirection*(shock.radiusAtStepStartM*factor);
    out.centerVelocityMPerS=shock.cmeDirection*(shock.radialSpeedMPerS*factor);
    out.radiusM=shock.radiusAtStepStartM*s*factor;
    out.radiusSpeedMPerS=shock.radialSpeedMPerS*s*factor;
  }
  return out;
}

// A full generating sphere includes a rear surface that is NOT the SSE front.
// At an outward ray intersection, n dot e_r >= 0 selects the outer root.
// The angular test separately enforces the finite cap. The tolerance accepts
// a true tangent boundary without extending the cap by a macroscopic angle.
inline bool OnOutwardShockCap(const ExpandingShock& shock,
                             const Core::Vec3& positionM,
                             const GeneratingSphere& sphere) {
  if (shock.geometry==ShockGeometryKind::Sphere) return true;
  const auto radial=positionM-shock.centerM;
  const double r=radial.Norm();
  if (!(r>0)) return false;
  const double tolerance=256*std::numeric_limits<double>::epsilon();
  if (radial.Dot(shock.cmeDirection)/r<std::cos(shock.halfWidthRad)-tolerance)
    return false;
  return radial.Dot(positionM-sphere.centerM)>=-tolerance*r*sphere.radiusM;
}

struct ShockSurfacePoint {
  bool exists = false;
  double radiusM = 0.0, normalSpeedMPerS = 0.0;
  Core::Vec3 positionM, outwardNormal;
};

// Directional geometry diagnostics and an oracle for the mover handoff.
// A surface is a GEOMETRIC front; the existence of a fast shock must still be
// queried from SWCME. Directions outside the cap return exists=false.
inline Core::Status EvaluateShockSurface(const ExpandingShock& shock,
    const Core::Vec3& direction, ShockSurfacePoint* output) {
  auto status=ValidateShockGeometry(shock);
  const double norm=direction.Norm();
  if (!status.ok()) return status;
  if (!output||!std::isfinite(norm)||norm<=0)
    return Core::Status(Core::StatusCode::InvalidInput,"invalid shock surface direction/output");
  ShockSurfacePoint result;
  const auto u=direction/norm;
  const auto sphere=ShockGeneratingSphere(shock);
  double radius=shock.radiusAtStepStartM;
  if (shock.geometry==ShockGeometryKind::FiniteSSE) {
    const double cosine=u.Dot(shock.cmeDirection);
    const double tolerance=256*std::numeric_limits<double>::epsilon();
    if (cosine<std::cos(shock.halfWidthRad)-tolerance) { *output=result; return Core::Status::OK(); }
    const double c=shock.radiusAtStepStartM/(1+std::sin(shock.halfWidthRad));
    double d=sphere.radiusM*sphere.radiusM-c*c*std::max(0.0,1-cosine*cosine);
    if (d<-tolerance*c*c) { *output=result; return Core::Status::OK(); }
    radius=c*cosine+std::sqrt(std::max(0.0,d));
  }
  result.exists=radius>0;
  result.radiusM=radius;
  result.positionM=shock.centerM+u*radius;
  result.outwardNormal=(result.positionM-sphere.centerM).Normalized();
  // Self-similar motion x(t)=R_apex(t)*x(0)/R_apex(0) gives the local normal
  // speed. It vanishes at a tangent flank; no empirical cosine-speed factor.
  result.normalSpeedMPerS=shock.radialSpeedMPerS*radius/shock.radiusAtStepStartM*
      u.Dot(result.outwardNormal);
  *output=result;
  return Core::Status::OK();
}
} }
#endif
