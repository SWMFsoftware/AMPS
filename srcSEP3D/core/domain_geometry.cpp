#include "domain_geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D { namespace Core {
Status BuildDomainBoundsM(const DomainGeometryParameters& p,
                          Vec3* minimumM, Vec3* maximumM) {
  if (!minimumM || !maximumM || !std::isfinite(p.endpointRadiusM) ||
      p.endpointRadiusM <= 0.0 || !std::isfinite(p.originM.Norm()))
    return Status(StatusCode::InvalidInput, "invalid domain bounds inputs");
  if (!p.cornerCube) {
    const Vec3 half(p.endpointRadiusM, p.endpointRadiusM, p.endpointRadiusM);
    *minimumM = p.originM - half;
    *maximumM = p.originM + half;
    return Status::OK();
  }
  if (!std::isfinite(p.solarSphereRadiusM) ||
      p.solarSphereRadiusM < Const::R_sun ||
      p.solarSphereRadiusM < p.parker.sourceRadiusM ||
      p.solarSphereRadiusM >= p.endpointRadiusM ||
      !std::isfinite(p.cornerMarginM) || p.cornerMarginM < 0.0)
    return Status(StatusCode::InvalidInput, "corner domain requires a full solar sphere inside its endpoint radius");
  Vec3 curveMinimum, curveMaximum;
  const Status enclosed = ParkerCurveBoundsM(p.endpointRadiusM, p.parker,
      p.corridorPaddingM, &curveMinimum, &curveMaximum);
  if (!enclosed.ok()) return enclosed;
  const Vec3 endpoint = ParkerCurvePoint(p.endpointRadiusM, p.parker);
  const double low[] = {curveMinimum.x, curveMinimum.y, curveMinimum.z};
  const double high[] = {curveMaximum.x, curveMaximum.y, curveMaximum.z};
  const double end[] = {endpoint.x, endpoint.y, endpoint.z};
  const double chosen[] = {p.cornerDirection.x, p.cornerDirection.y,
                            p.cornerDirection.z};
  const double origin[] = {p.originM.x, p.originM.y, p.originM.z};
  double direction[3], inset[3];
  // The endpoint distance sets a common cubic root scale. Increase it only
  // where the full winding, corridor width or solar sphere requires room.
  // A tightly wound line may need a larger inset on an axis; clipping it to
  // enforce a cosmetic corner placement would lose the requested field line.
  double side = p.endpointRadiusM +
      2.0 * (p.solarSphereRadiusM + p.cornerMarginM);
  for (int axis = 0; axis < 3; ++axis) {
    if (chosen[axis] != -1.0 && chosen[axis] != 0.0 && chosen[axis] != 1.0)
      return Status(StatusCode::InvalidInput, "corner directions must be -1, 0 or +1");
    if (p.centerZ && axis == 2) {
      if (chosen[axis] != 0.0)
        return Status(StatusCode::InvalidInput,
                      "x-y corner geometry requires corner_direction_z=0");
      // The z midplane passes through the Sun. Certified whole-curve bounds
      // include the corridor cross-section, so a polar/tilted line enlarges
      // BOTH z half-extents instead of shifting the Sun off that midplane.
      // Keeping the root cubic also enlarges x/y when this sets the scale.
      const double extent = std::max(p.solarSphereRadiusM,
          std::max(std::fabs(low[axis]), std::fabs(high[axis])));
      side = std::max(side, 2.0 * (extent + p.cornerMarginM));
      direction[axis] = 0.0;
      inset[axis] = 0.0;
      continue;
    }
    const double signProbe = std::fabs(end[axis]) >
        128.0 * std::numeric_limits<double>::epsilon() * p.endpointRadiusM
        ? end[axis] : high[axis] + low[axis];
    direction[axis] = chosen[axis] == 0.0
        ? (signProbe < 0.0 ? -1.0 : 1.0) : chosen[axis];
    const double signedLow = direction[axis] > 0.0 ? low[axis] : -high[axis];
    const double signedHigh = direction[axis] > 0.0 ? high[axis] : -low[axis];
    inset[axis] = std::max(p.solarSphereRadiusM, -signedLow) + p.cornerMarginM;
    // The opposite face must contain both members of the active union. A
    // wide corridor aimed away from the selected corner can require a large
    // inset while ending short of the sphere on the far side of that axis.
    side = std::max(side, inset[axis] +
        std::max(p.solarSphereRadiusM, signedHigh) + p.cornerMarginM);
  }
  if (!std::isfinite(side))
    return Status(StatusCode::InvalidInput, "corner domain extent overflow");
  double lower[3], upper[3];
  for (int axis = 0; axis < 3; ++axis) {
    if (p.centerZ && axis == 2) {
      lower[axis] = origin[axis] - 0.5 * side;
      upper[axis] = origin[axis] + 0.5 * side;
    } else {
      lower[axis] = origin[axis] + (direction[axis] > 0.0
          ? -inset[axis] : inset[axis] - side);
      upper[axis] = lower[axis] + side;
    }
    if (!std::isfinite(lower[axis]) || !std::isfinite(upper[axis]))
      return Status(StatusCode::InvalidInput, "corner domain coordinate overflow");
  }
  *minimumM = Vec3(lower[0], lower[1], lower[2]);
  *maximumM = Vec3(upper[0], upper[1], upper[2]);
  return Status::OK();
}
} }  // namespace SEP3D::Core
