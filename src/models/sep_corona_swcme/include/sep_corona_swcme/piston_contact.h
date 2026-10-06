#ifndef SEP_CORONA_SWCME_PISTON_CONTACT_H
#define SEP_CORONA_SWCME_PISTON_CONTACT_H

#include "sep_corona_swcme/cme_event.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace SEP { namespace CoronaSwcme {

class AmbientModel;

// This phase belongs to the independently prescribed material contact, not to
// EventConfiguration::At(), whose phase still describes the Level-A shock
// history.  Keeping the types distinct is a compile-time guard against using
// the observed/prescribed front as the Level-B piston.
enum class PistonContactPhase {
  StartupRamp,
  CoronalAnalytic,
  HandoffTransition,
  SwcmeOuter
};

struct PistonContactState {
  double timeS = 0.0;
  double handoffWeight = 0.0;
  PistonContactPhase phase = PistonContactPhase::StartupRamp;
  CoronalCME::EllipsoidKinematics ellipsoid;
  CoronalCME::KinematicValue apexRadiusM;
};

enum class PistonRayDisposition {
  Supported,
  NoIntersection,
  NonpositiveIntersection,
  UnsupportedFlank,
  BelowAmbientSupport
};

// Ray/contact data are evaluated from one implicit ellipsoid.  The radial
// speed and acceleration are analytic implicit derivatives; they are not
// finite differences of heliocentric positions.  This avoids catastrophic
// cancellation in the contact-leakage check which affected the diagnostic
// relaxation map.
struct PistonRayState {
  PistonRayDisposition disposition = PistonRayDisposition::NoIntersection;
  double timeS = 0.0;
  double radiusM = 0.0;
  double radialSpeedMPerS = 0.0;
  double radialAccelerationMPerS2 = 0.0;
  double incidence = 0.0; // q_hat dot n_c; positive on the outer intersection
  double normalSpeedMPerS = 0.0;
  CoronalCME::Vec3 positionM;
  CoronalCME::Vec3 outwardNormal;

  bool Supported() const noexcept {
    return disposition==PistonRayDisposition::Supported;
  }
};

struct PistonStartupDiagnostics {
  std::size_t supportedRays = 0;
  std::size_t unsupportedRays = 0;
  std::size_t incompatibleRays = 0;
  std::size_t limitingRay = 0;
  double maximumAmbientRadialMach = 0.0;
  bool compatible = false;
};

// The generic qualification contact is a prescribed ejecta boundary, not an
// MHD solution.  Its job is to provide one differentiable contact authority
// to geometry, the Level-B piston boundary, future BG3D-5 matching and native
// publication.  The downstream plasma and shock will be solved by the
// separate per-ray Lagrangian system.
class PistonContactModel final {
 public:
  static Core::Result<std::shared_ptr<const PistonContactModel>> Create(
      std::shared_ptr<const EventConfiguration> event);

  Core::Result<PistonContactState> At(double timeS) const;
  Core::Result<PistonRayState> EvaluateRay(
      CoronalCME::Vec3 unitDirection,double timeS) const;

  // The contact starts at rest because all global parameter rates use the
  // quintic startup ramp.  Compatibility additionally requires the ambient
  // radial flow to be small compared with its physical fast speed on every
  // supported ray.  This check is deliberately ray-set dependent and must be
  // rerun for each checksummed production quadrature.
  Core::Result<PistonStartupDiagnostics> CheckStartupCompatibility(
      const AmbientModel& ambient,
      const std::vector<CoronalCME::Vec3>& unitDirections) const;

  double HandoffTimeS() const noexcept { return handoffTimeS_; }
  const CoronalCME::RadialPrincipalBasis& Basis() const noexcept {
    return basis_;
  }
  const EventConfiguration& Event() const noexcept { return *event_; }

 private:
  std::shared_ptr<const EventConfiguration> event_;
  CoronalCME::RadialPrincipalBasis basis_;
  double handoffTimeS_ = 0.0;
};

const char* Name(PistonContactPhase phase) noexcept;
const char* Name(PistonRayDisposition disposition) noexcept;

} } // namespace SEP::CoronaSwcme

#endif
