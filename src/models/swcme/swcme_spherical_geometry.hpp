#pragma once

#include <array>
#include <cmath>
#include <limits>

#include "swcme_status.hpp"

namespace swcme {
namespace geometry {

// Canonical Sun-centred spherical limit shared by the legacy SWCME3D model
// and the separately typed coronal-ambient SHOCK_ONLY event. The caller owns
// the actual apex history and supplies an already normalized HCI ray. This
// geometry has no Leblanc density, wind, source, sheath or ejecta inputs.
// Preserving the supplied unit components avoids changing the legacy model's
// normalization arithmetic when it delegates its existing sphere branch.
struct SphereSurface {
  double radius_m = 0;
  std::array<double, 3> position_m{}, outward_normal{};
  double normal_speed_m_s = 0;
};
inline ModelStatus spherical_surface(double radius_m,
                                     const std::array<double, 3>& unit_ray,
                                     double normal_speed_m_s, SphereSurface& out) {
  double squared = 0;
  for (double x : unit_ray) {
    if (!std::isfinite(x))
      return ModelStatus::make(StatusCode::DegenerateVector,
                               "spherical front unit ray");
    squared += x * x;
  }
  if (!std::isfinite(radius_m) || radius_m <= 0 || !std::isfinite(normal_speed_m_s) ||
      normal_speed_m_s < 0 || !std::isfinite(squared) ||
      std::abs(squared - 1) > 64 * std::numeric_limits<double>::epsilon())
    return ModelStatus::make(StatusCode::GeometryFailure,
                             "spherical front requires positive SI radius, nonnegative "
                             "speed and unit HCI ray");
  SphereSurface candidate;
  candidate.radius_m = radius_m;
  candidate.outward_normal = unit_ray;
  candidate.normal_speed_m_s = normal_speed_m_s;
  for (unsigned i = 0; i < 3; ++i)
    candidate.position_m[i] = radius_m * unit_ray[i];
  out = candidate;
  return ModelStatus::success();
}

}  // namespace geometry
}  // namespace swcme
