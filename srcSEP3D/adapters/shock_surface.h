// ============================================================================
// srcSEP3D/adapters/shock_surface.h
//
// Provider-neutral shock surface for srcSEP3D particle sources and boundaries.
//
// The 3-D transport needs four things from a shock model, whatever produces
// it (the reduced shock-front core with either of its geometries today; other
// core models, or a front extracted from an MHD solution, later):
//
//   1. a level function L(x,t): L > 0 upstream, L < 0 in the geometric
//      interior (downstream), L = 0 on the front;
//   2. the crossings of a straight segment with the front at a fixed epoch;
//   3. the local front state at a point of the front: support, admission as a
//      fast shock, outward normal, normal speed, shock-frame inflow, upstream
//      proton density and field, compression ratio and theta_Bn;
//   4. the release ("free-escape") surface at a distance L_ref = f * r_front
//      upstream of the front, f = release_distance_shock_radius_fraction,
//      r_front the heliocentric distance of the front point (for the
//      heliocentric sphere, L_ref = f R_sh as in the 1-D srcSEP source).
//
// All positions, vectors and velocities are expressed in the run's
// computation frame (core/computation_frame.h).  An implementation converts
// to and from its own frame internally, so a consumer never sees a provider
// frame.  Epochs are the model's elapsed SI seconds.
//
// A new shock model is added by implementing this interface; the particle
// source (release on the tube field lines), the absorbing shock boundary and
// the time-step limit consume nothing else.
//
// IMPLEMENTATIONS
//   ReducedFrontShockSurface (adapters/reduced_front_shock_surface.*): the
//   shared reduced shock-front core, v1.3 heliocentric sphere or finite SSE
//   cap (srcSEP3D/test/event_2013_04_11/README.md, section 2.2).
//
// CONSUMERS IN THIS TREE
//   adapters/tube_line_source.cpp  Level and ReleaseLevel (where a tube line
//                                  crosses the front and the release surface),
//                                  FrontAt (upstream state of the net source);
//   main_lib.cpp ShockAcceptedAt   FrontAt: is the crossing point of the
//                                  absorbing-shock test an accepted front?
//   main_lib.cpp start-up log      GeometryName.
// The absorbing-shock crossing itself (adapters/shock_return.*) is computed on
// the mover-side ExpandingShock record (adapters/shock_geometry.h), which
// carries the exact front radius and speed at each step start.  Intersect,
// SupportBeginS/SupportEndS, Identity and MaximumNormalSpeedMPerS have no
// caller in this tree yet.
//
// STATUS CONVENTION
// Every query returns Core::Status.  An error status means the query itself
// failed (bad arguments, provider failure).  A point where the model has no
// front is not an error: FrontAt returns OK with state->supported = false,
// and Intersect omits crossings on unsupported parts of the surface.
// ============================================================================

#ifndef SEP3D_ADAPTERS_SHOCK_SURFACE_H
#define SEP3D_ADAPTERS_SHOCK_SURFACE_H

#include "../core/computation_frame.h"
#include "../core/sep3d_types.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace SEP3D {
namespace Adapters {

// Local state of the front at one front point, in the computation frame.
// Only `supported` is meaningful when it is false; every other field then
// keeps its default.  The source and the absorbing boundary act only where
// `accepted` is true.
struct ShockFrontState {
  bool supported = false;       // the point lies on the physical front
  bool accepted = false;        // admitted fast-mode shock with a valid jump
  Core::Vec3 positionM;         // on the front (frame) [m]
  Core::Vec3 outwardNormal;     // unit, from downstream toward upstream
  double normalSpeedMPerS = 0;  // front speed along outwardNormal (frame) [m/s]
  // w1 = V_n - U1.n [m/s], shock-frame normal inflow.  The same in every
  // frame: a frame change shifts V_n and U1.n by the same frame velocity.
  double inflowMPerS = 0;
  double compressionRatio = 1;  // rho2/rho1 of the accepted jump; 1 otherwise
  // Angle [rad] between upstream B and the normal, folded into [0, pi/2]
  // (the field polarity does not enter).
  double thetaBnRad = 0;
  double upstreamProtonDensityM3 = 0;  // n1 [m^-3]
  Core::Vec3 upstreamMagneticFieldT;   // B1, frame vector [T]
  Core::Vec3 upstreamVelocityMPerS;    // U1, frame velocity [m/s]
  std::uint64_t generation = 0;        // provider epoch generation (0 = none)
  std::string classification;          // provider status name (diagnostics)
};

// One crossing of a straight segment with the front (Intersect).
struct ShockSegmentCrossing {
  double fraction = 0;     // x = x0 + fraction (x1 - x0), fraction in [0,1]
  int direction = 0;       // +1 enters the downstream side, -1 leaves it
  ShockFrontState front;   // FrontAt at the crossing point
};

// Abstract shock surface.  It is shared read-only (ShockSurfacePointer) and
// every query is a const method that takes its epoch explicitly.  Arguments
// and results are in the computation frame, lengths in m, epochs in model
// elapsed s.
class ShockSurface {
 public:
  virtual ~ShockSurface() = default;

  // Short geometry name for receipts, e.g. "heliocentric-sphere".
  virtual std::string GeometryName() const = 0;
  // Identity of the physical event behind the surface (fingerprint text).
  virtual std::string Identity() const = 0;
  // Epoch support [begin, end] of the surface history.
  virtual double SupportBeginS() const = 0;
  virtual double SupportEndS() const = 0;

  // L(x,t) > 0 upstream, < 0 downstream.  For the heliocentric sphere this
  // is the exact signed distance |x| - R(t); for other geometries any
  // continuous level that is a distance to first order near the front.
  // *level [m] is written only on success.
  virtual Core::Status Level(const Core::Vec3& positionM, double epochS,
                             double* level) const = 0;

  // Crossings of the segment [x0,x1] with the front frozen at epochS,
  // sorted by fraction.  A tangent touch is not a crossing.  *crossings is
  // cleared first; only supported front points are reported.
  virtual Core::Status Intersect(const Core::Vec3& x0, const Core::Vec3& x1,
                                 double epochS,
                                 std::vector<ShockSegmentCrossing>* crossings)
      const = 0;

  // Front state at the front point nearest to positionM along the
  // provider's natural projection (radial for a heliocentric sphere).
  // positionM need not lie on the front.  *state is reset first; OK with
  // state->supported = false means the projection meets no physical front.
  virtual Core::Status FrontAt(const Core::Vec3& positionM, double epochS,
                               ShockFrontState* state) const = 0;

  // Level of the release surface: > 0 beyond it (farther upstream),
  // < 0 between it and the front.  fraction = f of L_ref = f r_front.
  // A negative or non-finite fraction is an error.  The level is also < 0
  // downstream of the front.
  virtual Core::Status ReleaseLevel(const Core::Vec3& positionM,
                                    double epochS, double fraction,
                                    double* level) const = 0;

  // Upper bound of |dL/dt| over [t0,t1]: the largest front normal speed.
  // Used only to limit substeps and to prune candidates.
  virtual double MaximumNormalSpeedMPerS(double t0S, double t1S) const = 0;
};

// Shared, read-only handle; main_lib.cpp holds the run's surface.
using ShockSurfacePointer = std::shared_ptr<const ShockSurface>;

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_SHOCK_SURFACE_H
