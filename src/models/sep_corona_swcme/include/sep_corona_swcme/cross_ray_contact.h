#ifndef SEP_CORONA_SWCME_CROSS_RAY_CONTACT_H
#define SEP_CORONA_SWCME_CROSS_RAY_CONTACT_H

#include "sep_corona_swcme/cross_ray_transport.h"
#include "sep_corona_swcme/piston_contact.h"

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

// New, explicitly selected geometry increment; it does not change the fixed-
// HCI legacy contact. The full base center AND axes rotate about the Sun's HCI
// z axis. At startup omega equals the retained ambient rotation rate; a
// quintic RATE transition reduces omega to zero over an EXPLICIT duration.
// Its integral leaves theta_final=Omega*T/2, preserving continuous geometry.
// The user selected T=1200 s for the future expanded production fixture.
// No duration is inferred from a legacy asset or silently defaulted here.
class CorotatingPistonContact final {
 public:
  static Core::Result<std::shared_ptr<const CorotatingPistonContact>> Create(
      std::shared_ptr<const PistonContactModel> reference,double transitionDurationS);
  Core::Result<CoronalCME::KinematicValue> RotationAt(double epochS) const;
  Core::Result<CoronalCME::Vec3> SurfacePoint(
      double polarRad,double azimuthRad,double epochS) const;
  Core::Result<CoronalCME::Vec3> MaterialVelocity(
      double polarRad,double azimuthRad,double epochS) const;
  Core::Result<CoronalCME::SurfaceEvaluation> Evaluate(
      CoronalCME::Vec3 positionHciM,double epochS) const;
  Core::Result<PistonRayState> EvaluateRay(
      CoronalCME::Vec3 unitHciDirection,double epochS) const;
  const std::string& Fingerprint() const noexcept {return fingerprint_;}
  const std::shared_ptr<const PistonContactModel>& Reference() const noexcept {return reference_;}

 private:
  // Only a validated factory may create an authority; an empty public default
  // object would bypass mesh/epoch/positivity checks before a later query.
  CorotatingPistonContact()=default;
  std::shared_ptr<const PistonContactModel> reference_;
  double durationS_=0;
  std::string fingerprint_;
};

struct IdealContactResiduals {
  double leftRelativeMassFluxKgM2S=0,rightRelativeMassFluxKgM2S=0;
  double normalMagneticJumpT=0,pressureJumpPa=0;
  CoronalCME::Vec3 velocityJumpMPerS,magneticJumpT;
  CoronalCME::Vec3 tractionJumpPa;
  CoronalCME::Vec3 tangentialElectricJumpVPerM;
  double relativeTotalEnergyFluxJumpWPerM2=0;
};

// Dimensional diagnostics for the explicitly selected continuous-field ideal
// material contact. With nonzero Bn, zero relative normal mass flux requires
// continuous v, B and p (rho/entropy may jump). We retain that same continuous
// branch through Bn=0; there is no raw Bn!=0 switch or division by Bn, and no
// claim that a discontinuous tangential-contact branch is implemented.
//
// n points left/ejecta to right/sheath; wn is the physical normal surface
// speed in m/s. Maxwell traction is (p+B^2/2mu)n-(Bn/mu)B in Pa. Tangential
// E is evaluated in a frame moving at wn*n, while total energy remains HCI
// and its interface flux is F_E-wn*E. No tolerances/admissibility verdicts are
// invented here: production qualification applies its unchanged gates to the
// raw signed SI residuals. This is not an ejecta-state generator.
Core::Result<IdealContactResiduals> EvaluateIdealContactResiduals(
    const Primitive& ejecta,const Primitive& sheath,
    const FaceGeometry& contact,const Thermodynamics& gas);

} } } // namespace SEP::CoronaSwcme::CrossRay

#endif
