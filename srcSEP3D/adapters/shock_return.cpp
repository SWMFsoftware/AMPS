// ============================================================================
// srcSEP3D/adapters/shock_return.cpp
//
// Absorbing-shock test of one transport substep; see shock_return.h for the
// contract and the caller.  The four cases below are tried in order.
// ============================================================================

#include "shock_return.h"

#include <cmath>

namespace SEP3D {
namespace Adapters {

ShockReturnTest TestShockReturn(const Core::Vec3& x0, const Core::Vec3& x1,
                                double dtS, const ExpandingShock& shock,
                                double kappaNormalDtM2,
                                Transport::KeyedRandomStream* bridge) {
  ShockReturnTest test;
  // Generating sphere at the substep start: centre c0 [m] moving at c' [m/s]
  // (zero for the heliocentric sphere), radius a0 [m] growing at a' [m/s].
  const GeneratingSphere sphere = ShockGeneratingSphere(shock);
  // Particle position relative to the moving centre, linear in tau:
  //   q(tau) = (x0 - c0) + tau [(x1 - x0) - c' dt] = q0 + tau dq,
  // and the radius a(tau) = a0 + tau da with da = a' dt.
  const Core::Vec3 q0 = x0 - sphere.centerM;
  const Core::Vec3 dq = (x1 - x0) - sphere.centerVelocityMPerS * dtS;
  const double a0 = sphere.radiusM, da = sphere.radiusSpeedMPerS * dtS;
  // g(tau) [m]: signed distance of the particle from the moving sphere,
  // > 0 upstream (outside), <= 0 downstream (inside).
  const auto g = [&](double tau) { return (q0 + dq * tau).Norm() - (a0 + da * tau); };
  const double g0 = g(0.0), g1 = g(1.0);
  // (1) Downstream start (g0 <= 0), e.g. overtaken by the front between
  // substeps: hit at tau = 0.
  if (g0 <= 0.0) { test.hit = true; test.fraction = 0.0; return test; }
  // (2) Downstream end: the zero of g interpolated linearly between the ends.
  if (g1 <= 0.0) { test.hit = true; test.fraction = g0 / (g0 - g1); return test; }
  // (3) Diffusive (Parker) substep with both ends upstream.  Near the front
  // the normal coordinate g performs, to first order, a one-dimensional
  // Brownian motion with variance sigma^2 = 2 kappa_nn dt over the substep,
  // plus a drift (advection and the front's own motion).  Conditioned on its
  // end values g0 > 0 and g1 > 0 it is a Brownian bridge, whose law does not
  // depend on a constant drift, and the bridge touches g = 0 with probability
  //   P = exp(-2 g0 g1 / sigma^2) = exp(-g0 g1 / (kappa_nn dt)).
  // The test draws u in (0,1) and hits when log(u) < -g0 g1/(kappa_nn dt),
  // i.e. u < P; the exponent is formed as (g0/sqrt(kappa_nn dt)) times
  // (g1/sqrt(kappa_nn dt)).  This treats the front as planar and kappa_nn
  // as constant over the diffusion length sqrt(2 kappa_nn dt).  The hit
  // fraction g0/(g0 + g1) is where the straight line from g0 to the mirror end -g1 of
  // the reflection principle crosses zero; it locates the hit for the caller
  // and is not a sample of the first-passage time.  A diffusive substep is
  // never also given the chord test below.
  if (kappaNormalDtM2 > 0.0 && bridge != nullptr) {
    const double u = bridge->UniformOpen01();
    if (std::log(u) < -(g0 / std::sqrt(kappaNormalDtM2)) *
                       (g1 / std::sqrt(kappaNormalDtM2))) {
      test.hit = true;
      test.fraction = g0 / (g0 + g1);
    }
    return test;
  }
  // (4) Deterministic (focused) substep with both ends upstream.  The
  // straight chord between the ends can still dip into the sphere.  g is
  // convex in tau (the norm of an affine function minus an affine function),
  // so it dips below zero only if its minimum does:
  //   g'(tau) = q(tau).dq / |q(tau)| - da.
  // g'(0) >= 0: g increases on [0,1], its minimum is g0 > 0, no hit.
  // g'(1) <= 0: g decreases on [0,1], its minimum is g1 > 0, no hit.
  // Convex g: locate the minimum from the sign of g'.
  const auto slope = [&](double tau) {
    const Core::Vec3 q = q0 + dq * tau;
    const double n = q.Norm();
    return (n > 0.0 ? q.Dot(dq) / n : 0.0) - da;
  };
  if (slope(0.0) >= 0.0 || slope(1.0) <= 0.0) return test;
  // Bisection on the sign of g'; 60 halvings resolve tau to 2^-60 ~ 1e-18,
  // below double resolution of [0,1].
  double lo = 0.0, hi = 1.0;
  for (int i = 0; i < 60; ++i) {
    const double m = 0.5 * (lo + hi);
    if (slope(m) < 0.0) lo = m; else hi = m;
  }
  const double tauMin = 0.5 * (lo + hi);
  if (g(tauMin) > 0.0) return test;   // the chord stays upstream
  // g(0) > 0 >= g(tauMin) and g decreases on [0, tauMin]: bisect for the
  // first root, which is where the particle reaches the front.  hi keeps
  // g(hi) <= 0, so the reported fraction is on or just past the front.
  lo = 0.0; hi = tauMin;   // first root on [0, tauMin]
  for (int i = 0; i < 60; ++i) {
    const double m = 0.5 * (lo + hi);
    if (g(m) > 0.0) lo = m; else hi = m;
  }
  test.hit = true;
  test.fraction = hi;
  return test;
}

}  // namespace Adapters
}  // namespace SEP3D
