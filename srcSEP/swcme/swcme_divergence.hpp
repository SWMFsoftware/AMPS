#pragma once

// ============================================================================
// swcme_divergence.hpp
// ----------------------------------------------------------------------------
// Shared velocity-divergence mathematics used by the 1-D and 3-D SWCME
// interfaces.
//
// The central design rule is that the numerical operator must match the actual
// dimensionality of the velocity field:
//   * purely radial profiles use the exact spherical identity
//         div(V) = 2 V_r/r + dV_r/dr ;
//   * the 3-D FULL_ICME field may contain non-radial Rankine-Hugoniot velocity
//     components and angular variation of the shock/sheath surfaces, so its
//     divergence is evaluated from the full Cartesian Jacobian
//         dVx/dx + dVy/dy + dVz/dz.
//
// The previous 3-D implementation always differentiated only r^2 V_r along a
// ray.  That expression is exact only for a radial vector field and therefore
// silently omitted angular/tangential contributions in the resolved ICME.
// ============================================================================

#include "swcme_status.hpp"

#include <array>
#include <cmath>

namespace swcme {
namespace divergence {

// Decomposed radial result used by diagnostics/tests.  Keeping the geometric
// 2V/r term separate from dV/dr makes it straightforward to identify whether a
// discrepancy originates in the spherical geometry or in the profile
// derivative itself.
struct RadialTerms {
  double geometric_s_inv = 0.0;
  double derivative_s_inv = 0.0;
  double total_s_inv = 0.0;
};

inline RadialTerms radial_terms(double radius_m,
                                double radial_velocity_m_s,
                                double d_radial_velocity_dr_s_inv) {
  RadialTerms out;
  out.geometric_s_inv = 2.0 * radial_velocity_m_s / radius_m;
  out.derivative_s_inv = d_radial_velocity_dr_s_inv;
  out.total_s_inv = out.geometric_s_inv + out.derivative_s_inv;
  return out;
}

// Generic second-order Cartesian divergence operator.  `evaluate_velocity`
// must have signature
//
//   ModelStatus(const std::array<double,3>& x_m,
//               std::array<double,3>& velocity_m_s)
//
// and return OUTSIDE_MODEL_DOMAIN when a stencil point lies outside its valid
// domain.  Centered differences are used whenever possible.  If exactly one
// side leaves the domain, the operator switches explicitly to the matching
// second-order one-sided stencil (-3,4,-1)/(2h) or its backward analogue.
// This is a documented boundary rule, not a hidden radius clip.
template <class VelocityEvaluator>
ModelStatus cartesian_second_order(
    const std::array<double,3>& point_m,
    double step_m,
    VelocityEvaluator&& evaluate_velocity,
    double& divergence_s_inv) {
  if (!std::isfinite(step_m) || !(step_m > 0.0)) {
    return ModelStatus::make_value(
        StatusCode::InvalidNumericalStep,
        "swcme::divergence::cartesian_second_order step",step_m);
  }
  for (double q : point_m) {
    if (!std::isfinite(q)) {
      return ModelStatus::make(
          StatusCode::NonFiniteInput,
          "swcme::divergence::cartesian_second_order coordinate");
    }
  }

  std::array<double,3> center_v{{0.0,0.0,0.0}};
  bool have_center=false;
  double div=0.0;

  // Array subscripts use size_type directly.  Besides satisfying the strict
  // sign-conversion gate, this prevents a future signed axis value from being
  // silently converted to a very large index before bounds checking.
  for (std::size_t axis=0;axis<point_m.size();++axis) {
    std::array<double,3> plus=point_m;
    std::array<double,3> minus=point_m;
    plus[axis]+=step_m;
    minus[axis]-=step_m;

    std::array<double,3> vp{{0.0,0.0,0.0}};
    std::array<double,3> vm{{0.0,0.0,0.0}};
    ModelStatus sp=evaluate_velocity(plus,vp);
    ModelStatus sm=evaluate_velocity(minus,vm);

    double deriv=0.0;
    if (sp.ok() && sm.ok()) {
      deriv=(vp[axis]-vm[axis])/(2.0*step_m);
    } else if (sm.code==StatusCode::OutsideModelDomain && sp.ok()) {
      // Forward one-sided second-order derivative.  Evaluate the center and
      // x+2h only when the centered stencil actually crosses the domain edge.
      if (!have_center) {
        ModelStatus s0=evaluate_velocity(point_m,center_v);
        if (!s0.ok()) return s0;
        have_center=true;
      }
      std::array<double,3> plus2=point_m;
      plus2[axis]+=2.0*step_m;
      std::array<double,3> vp2{{0.0,0.0,0.0}};
      ModelStatus sp2=evaluate_velocity(plus2,vp2);
      if (!sp2.ok()) return sp2;
      deriv=(-3.0*center_v[axis]+4.0*vp[axis]-vp2[axis])/(2.0*step_m);
    } else if (sp.code==StatusCode::OutsideModelDomain && sm.ok()) {
      // Backward one-sided second-order derivative (included for completeness
      // even though SWCME's inner spherical boundary normally triggers the
      // forward form on the locally radial coordinate).
      if (!have_center) {
        ModelStatus s0=evaluate_velocity(point_m,center_v);
        if (!s0.ok()) return s0;
        have_center=true;
      }
      std::array<double,3> minus2=point_m;
      minus2[axis]-=2.0*step_m;
      std::array<double,3> vm2{{0.0,0.0,0.0}};
      ModelStatus sm2=evaluate_velocity(minus2,vm2);
      if (!sm2.ok()) return sm2;
      deriv=(3.0*center_v[axis]-4.0*vm[axis]+vm2[axis])/(2.0*step_m);
    } else {
      // Preserve the first real evaluator failure.  If both sides are merely
      // outside the domain there is no valid second-order stencil at this
      // point, which is reported explicitly as an invalid numerical step.
      if (sp.failure() && sp.code!=StatusCode::OutsideModelDomain) return sp;
      if (sm.failure() && sm.code!=StatusCode::OutsideModelDomain) return sm;
      return ModelStatus::make(
          StatusCode::InvalidNumericalStep,
          "swcme::divergence::cartesian_second_order stencil");
    }

    if (!std::isfinite(deriv)) {
      return ModelStatus::make(
          StatusCode::NonFiniteResult,
          "swcme::divergence::cartesian_second_order derivative");
    }
    div+=deriv;
  }

  if (!std::isfinite(div)) {
    return ModelStatus::make(
        StatusCode::NonFiniteResult,
        "swcme::divergence::cartesian_second_order divergence");
  }
  divergence_s_inv=div;
  return ModelStatus::success();
}

} // namespace divergence
} // namespace swcme
