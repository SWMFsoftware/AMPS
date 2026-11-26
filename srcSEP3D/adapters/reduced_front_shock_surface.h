// ============================================================================
// srcSEP3D/adapters/reduced_front_shock_surface.h
//
// ShockSurface over the shared reduced shock-front core
// (src/models/sep_corona_swcme/shock_front): the v1.3 closed heliocentric
// sphere and the finite fixed-width SSE cap.  The core works in HCI; this
// adapter converts every query from and every result to the run's
// computation frame.  It adds no shock physics: support, admission, the
// Rankine-Hugoniot jump and the upstream state are the core's.
//
// GEOMETRY
// Both core geometries are described by a generating sphere at each epoch,
// Provider::Sphere(t), with centre c(t) and radius a(t) in HCI:
//   heliocentric sphere  c = 0, a = R_sh(t) (the front IS the sphere);
//   finite SSE cap       c = R_apex/(1 + sin(lambda)) e_cme,
//                        a = |c| sin(lambda), lambda the half width
//                        (the front is the outward cap of that sphere; its
//                        rear and flanks are geometry only).
// Level, Intersect and ReleaseLevel are computed on the generating sphere;
// FrontAt and Intersect ask the core whether a point is on the physical
// front (EvaluateRay / EvaluateFrontPoint).
//
// FRAMES
// A frame point x is converted to HCI as Rz(theta) (x - origin), where origin
// is the Sun's position in mesh coordinates (RunOptions::coordinateOriginM)
// and Rz(theta) the computation-frame rotation (core/computation_frame.h).
// Results go back the same way.  Beyond the rotation of every vector, only
// the front's normal speed and the upstream velocity lose a frame-velocity
// term in a rotating frame; see FrontAt.
//
// Used by the tube-line source and the absorbing shock
// (srcSEP3D/test/event_2013_04_11/README.md, sections 2.2, 2.4 and 2.5).
// ============================================================================

#ifndef SEP3D_ADAPTERS_REDUCED_FRONT_SHOCK_SURFACE_H
#define SEP3D_ADAPTERS_REDUCED_FRONT_SHOCK_SURFACE_H

#include "shock_surface.h"

#include <memory>

namespace SEP { namespace CoronaSwcme { namespace ShockFront {
class Provider;
} } }

namespace SEP3D {
namespace Adapters {

class ReducedFrontShockSurface final : public ShockSurface {
 public:
  // Build the surface over `provider` (the run's single shock-front
  // authority, ShockFrontBackgroundAdapter::SharedProvider()) for the run's
  // computation `frame`; originM [m] is the Sun's position in mesh
  // coordinates.  Errors: null provider or output, an invalid frame, or an
  // event whose coordinate frame is not HCI.  The equality of the frame rate
  // and the event's ambient rotation rate is checked by
  // ShockFrontBackgroundAdapter::Create, not here.
  static Core::Status Create(
      std::shared_ptr<const SEP::CoronaSwcme::ShockFront::Provider> provider,
      const Core::ComputationFrame& frame, const Core::Vec3& originM,
      ShockSurfacePointer* output);

  // "heliocentric-sphere" or "finite-sse-fixed-width".
  std::string GeometryName() const override;
  // The event's physics fingerprint.
  std::string Identity() const override;
  // The event's ambient support interval (ambient.support) [model s].
  double SupportBeginS() const override;
  double SupportEndS() const override;
  Core::Status Level(const Core::Vec3& positionM, double epochS,
                     double* level) const override;
  Core::Status Intersect(const Core::Vec3& x0, const Core::Vec3& x1,
                         double epochS,
                         std::vector<ShockSegmentCrossing>* crossings)
      const override;
  Core::Status FrontAt(const Core::Vec3& positionM, double epochS,
                       ShockFrontState* state) const override;
  Core::Status ReleaseLevel(const Core::Vec3& positionM, double epochS,
                            double fraction, double* level) const override;
  double MaximumNormalSpeedMPerS(double t0S, double t1S) const override;

 private:
  ReducedFrontShockSurface() = default;   // construct through Create()
  // True for the v1.3 closed heliocentric sphere, false for the SSE cap.
  bool Spherical() const;
  // Generating sphere in HCI at epochS (centre [m], radius [m]); fails when
  // the core's trajectory has no state at epochS.
  Core::Status GeneratingSphere(double epochS, Core::Vec3* centreHci,
                                double* radiusM) const;
  std::shared_ptr<const SEP::CoronaSwcme::ShockFront::Provider> provider_;
  Core::ComputationFrame frame_;   // run's computation frame
  Core::Vec3 originM_;             // Sun's position in mesh coordinates [m]
};

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_REDUCED_FRONT_SHOCK_SURFACE_H
