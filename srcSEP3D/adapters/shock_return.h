// ============================================================================
// srcSEP3D/adapters/shock_return.h
//
// Absorbing-shock test of one transport substep.
//
// WHY THIS EXISTS
// source.shock.returning_particles = absorb removes a particle that returns
// from upstream into the downstream side of an accepted shock.  The 1-D net
// source already counts only the particles that escape upstream (srcSEP 1-D
// README section 2.2), so a returning particle would otherwise be counted
// twice (srcSEP3D/test/event_2013_04_11/README.md, section 2.5).  The 1-D
// Parker mover asks the same question with a Brownian bridge
// (srcSEP/parker_mover.cpp, ReducedParticles::Native::CheckBoundary).
//
// WHERE IT FITS
// adapters/transport_adapter.cpp calls TestShockReturn after every mover
// substep.  On a hit it places the crossing at x0 + fraction (x1 - x0),
// projects it radially from the generating-sphere centre onto the front,
// requires that point to lie on the outward SSE cap (OnOutwardShockCap,
// adapters/shock_geometry.h), and asks the run's ShockSurface whether the
// front there is an accepted shock (main_lib.cpp ShockAcceptedAt).  Only then
// is the particle absorbed.  This function is geometry and probability only.
// ============================================================================
#ifndef SEP3D_ADAPTERS_SHOCK_RETURN_H
#define SEP3D_ADAPTERS_SHOCK_RETURN_H

#include "../transport/keyed_random.h"
#include "shock_geometry.h"

namespace SEP3D {
namespace Adapters {

// Absorbing-shock test of one substep.  The front is the shock's generating
// sphere (adapters/shock_geometry.h, ShockGeneratingSphere) with centre and
// radius moving linearly over the substep from their substep-start values:
// the caller supplies the exact radius and speed at the step start and
// advances the radius linearly to each substep start.  The sphere is the
// whole front for the heliocentric shock; for the finite SSE cap the caller
// additionally requires the hit to lie on the outward cap, and the cap axis
// is held at its step-start direction.  With g(tau) = |q(tau)| - a(tau), q
// the particle position relative to the moving centre and a the radius,
// tau in [0,1]:
//  * g(0) <= 0: the particle starts downstream (overtaken by the front);
//  * g(1) <= 0: the particle ends downstream;
//  * diffusive (Parker) motion with both ends upstream: the Brownian bridge
//    of the 1-D boundary, hit probability exp(-g0 g1 / (kappa_r dt)) with
//    kappa_r = kappa_par (b.e_r)^2 the diffusion normal to the front;
//  * deterministic (focused) motion: the straight chord can dip below the
//    front between positive ends; g is convex, so its minimum decides.
//
// Arguments (all in the computation frame):
//   x0, x1           particle position at the substep start and end [m];
//   dtS              substep length [s];
//   shock            mover-side front record at the substep start;
//   kappaNormalDtM2  kappa_nn dt [m^2], the diffusion coefficient normal to
//                    the front times dtS.  > 0 selects the Brownian bridge
//                    (the caller passes kappa_par (b.n)^2 dt for the Parker
//                    mover, with b the field direction and n the front normal
//                    at the substep start), 0 the straight-chord test;
//   bridge           keyed random stream for the bridge decision (purpose
//                    ShockReturnBridge), so a rerun takes the same decision;
//                    nullptr also selects the straight-chord test.
// Result: hit, and the substep fraction of the crossing (see below).
struct ShockReturnTest {
  bool hit = false;       // the path reached the downstream side
  // Substep fraction tau of the crossing: 0 for a downstream start, the zero
  // of g interpolated between the ends for a downstream end, g0/(g0 + g1)
  // for a bridge hit, and the first root of g for a chord dip.
  double fraction = 0.0;
};

ShockReturnTest TestShockReturn(const Core::Vec3& x0, const Core::Vec3& x1,
                                double dtS, const ExpandingShock& shock,
                                double kappaNormalDtM2,
                                Transport::KeyedRandomStream* bridge);

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_SHOCK_RETURN_H
