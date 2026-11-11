#include "sep_corona_swcme/cross_ray_contact.h"

#include "sep_coronal_cme/configuration_parser.h"

#include <cmath>
#include <limits>
#include <locale>
#include <sstream>

namespace SEP { namespace CoronaSwcme { namespace CrossRay { namespace {

using CoronalCME::Cross;
using CoronalCME::Dot;
using CoronalCME::Vec3;

Vec3 Rotate(Vec3 x,double angle) {
  const double c=std::cos(angle),s=std::sin(angle);
  return {c*x.x-s*x.y,s*x.x+c*x.y,x.z};
}

Core::Result<CoronalCME::FixedOrientationEllipsoid> ReferenceShape(
    const PistonContactModel& model,double t) {
  const auto state=model.At(t);
  if(!state.ok())return Core::Result<CoronalCME::FixedOrientationEllipsoid>::Failure(
      state.status.code,state.status.message);
  return CoronalCME::FixedOrientationEllipsoid::FromCenter(model.Basis(),
      state.value.ellipsoid,model.Event().support.solarRadiusM);
}

bool Finite(Vec3 v) {
  return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);
}

} // namespace

Core::Result<std::shared_ptr<const CorotatingPistonContact>> CorotatingPistonContact::Create(
    std::shared_ptr<const PistonContactModel> reference,double duration) {
  using Return=Core::Result<std::shared_ptr<const CorotatingPistonContact>>;
  if(!reference||!std::isfinite(duration)||!(duration>0))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"corotating contact requires an explicit positive transition duration");
  const auto& e=reference->Event();
  if(e.pistonContact.startS!=e.support.startS||
      !std::isfinite(e.ambient.rotationRateRadPerS)||
      !(e.support.startS+duration>e.support.startS)||
      !std::isfinite(e.support.startS+duration)||
      e.support.startS+duration>reference->HandoffTimeS())return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
          "corotation must start with the event and finish before the declared handoff");
  std::shared_ptr<CorotatingPistonContact> out(new CorotatingPistonContact);
  out->reference_=std::move(reference);out->durationS_=duration;
  std::ostringstream bytes;bytes.imbue(std::locale::classic());
  bytes<<"corotating-startup-continuous-field-contact-v1\n"<<e.physicsFingerprint<<'\n'<<
      std::hexfloat<<e.support.startS<<','<<duration<<','<<e.ambient.rotationRateRadPerS<<'\n';
  out->fingerprint_=CoronalCME::ComputeContentChecksum(bytes.str());
  return Return::Success(std::move(out));
}

Core::Result<CoronalCME::KinematicValue> CorotatingPistonContact::RotationAt(double t) const {
  const auto& e=reference_->Event();
  if(!std::isfinite(t)||t<e.support.startS||t>e.support.endS)
    return Core::Result<CoronalCME::KinematicValue>::Failure(
        Core::StatusCode::OutOfDomain,"corotation epoch outside owning event support");
  // theta'=Omega*(1-S5(tau)), S5=10tau^3-15tau^4+6tau^5.
  // Reuse the qualified analytic integral with full first/second derivatives;
  // rotating only axes or blending already computed velocities omits center
  // motion and produces false normal leakage at startup or during transition.
  return CoronalCME::EvaluateSmoothKinematics(t,e.support.startS,0,e.support.startS,
      durationS_,e.ambient.rotationRateRadPerS,0);
}

Core::Result<Vec3> CorotatingPistonContact::SurfacePoint(double polar,double azimuth,double t) const {
  const auto rotation=RotationAt(t);const auto shape=ReferenceShape(*reference_,t);
  if(!rotation.ok())return Core::Result<Vec3>::Failure(rotation.status.code,rotation.status.message);
  if(!shape.ok())return Core::Result<Vec3>::Failure(shape.status.code,shape.status.message);
  if(!std::isfinite(polar)||!std::isfinite(azimuth))return Core::Result<Vec3>::Failure(
      Core::StatusCode::InvalidConfiguration,"contact material parameters are nonfinite");
  return Core::Result<Vec3>::Success(Rotate(shape.value.Point(polar,azimuth),rotation.value.value));
}

Core::Result<Vec3> CorotatingPistonContact::MaterialVelocity(double polar,double azimuth,double t) const {
  const auto point=SurfacePoint(polar,azimuth,t);const auto rotation=RotationAt(t);
  const auto shape=ReferenceShape(*reference_,t);
  if(!point.ok())return Core::Result<Vec3>::Failure(point.status.code,point.status.message);
  if(!rotation.ok())return Core::Result<Vec3>::Failure(rotation.status.code,rotation.status.message);
  if(!shape.ok())return Core::Result<Vec3>::Failure(shape.status.code,shape.status.message);
  const Vec3 omega{0,0,rotation.value.firstDerivative};
  return Core::Result<Vec3>::Success(Rotate(shape.value.SurfaceVelocity(polar,azimuth),
      rotation.value.value)+Cross(omega,point.value));
}

Core::Result<CoronalCME::SurfaceEvaluation> CorotatingPistonContact::Evaluate(Vec3 x,double t) const {
  using Return=Core::Result<CoronalCME::SurfaceEvaluation>;
  const auto rotation=RotationAt(t);const auto shape=ReferenceShape(*reference_,t);
  if(!rotation.ok())return Return::Failure(rotation.status.code,rotation.status.message);
  if(!shape.ok())return Return::Failure(shape.status.code,shape.status.message);
  if(!Finite(x))return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "contact query position is nonfinite");
  const auto value=shape.value.Evaluate(Rotate(x,-rotation.value.value));
  if(!value.ok())return value;
  auto out=value.value;out.outwardNormal=Rotate(out.outwardNormal,rotation.value.value);
  out.normalSpeedMPerS+=Dot(Cross({0,0,rotation.value.firstDerivative},x),out.outwardNormal);
  if(!std::isfinite(out.normalSpeedMPerS))return Return::Failure(
      Core::StatusCode::NumericalFailure,"rotating contact normal speed overflows SI output");
  return Return::Success(out);
}

Core::Result<PistonRayState> CorotatingPistonContact::EvaluateRay(Vec3 q,double t) const {
  using Return=Core::Result<PistonRayState>;
  const auto angle=RotationAt(t);
  if(!angle.ok())return Return::Failure(angle.status.code,angle.status.message);
  const auto original=reference_->EvaluateRay(Rotate(q,-angle.value.value),t);
  if(!original.ok())return original;
  auto out=original.value;
  // Exclusion receipts may retain a positive intersection and its normal
  // (below plasma support or at a grazing flank). Those geometry values still
  // belong in HCI even when no physical state/derivative is eligible. Do not
  // leak reference-body vectors into unsupported-area/coverage accounting.
  out.positionM=Rotate(out.positionM,angle.value.value);
  out.outwardNormal=Rotate(out.outwardNormal,angle.value.value);
  if(!out.Supported())return Return::Success(out);
  const auto state=reference_->At(t);
  if(!state.ok())return Return::Failure(state.status.code,state.status.message);
  const auto& k=state.value.ellipsoid;const auto& b=reference_->Basis();
  const Vec3 axes[]={Rotate(b.radial,angle.value.value),
      Rotate(b.firstLateral,angle.value.value),Rotate(b.secondLateral,angle.value.value)};
  const CoronalCME::KinematicValue lengths[]={k.radialSemiAxisM,
      k.firstLateralSemiAxisM,k.secondLateralSemiAxisM};
  const Vec3 omega{0,0,angle.value.firstDerivative},alpha{0,0,angle.value.secondDerivative};
  long double fr=0,ft=0,frr=0,frt=0,ftt=0;
  for(int i=0;i<3;++i) {
    // F=sum_i z_i^2-1, z_i=(r*q.e_i-c*delta_i0)/a_i. Derivatives below
    // are partials at FIXED HCI ray q and fixed r. Rotating the query ray
    // without differentiating it would return the legacy speed incorrectly.
    const long double s=Dot(q,axes[i]),sd=Dot(q,Cross(omega,axes[i]));
    const long double sdd=Dot(q,Cross(alpha,axes[i])+Cross(omega,Cross(omega,axes[i])));
    const long double a=lengths[i].value,ad=lengths[i].firstDerivative,
        add=lengths[i].secondDerivative;
    const long double y=out.radiusM*s-(i==0?k.centerDistanceM.value:0);
    const long double yd=out.radiusM*sd-(i==0?k.centerDistanceM.firstDerivative:0);
    const long double ydd=out.radiusM*sdd-(i==0?k.centerDistanceM.secondDerivative:0);
    const long double z=y/a,zd=yd/a-y*ad/(a*a);
    const long double zdd=ydd/a-2*yd*ad/(a*a)-y*add/(a*a)+2*y*ad*ad/(a*a*a);
    fr+=2*z*s/a;ft+=2*z*zd;frr+=2*s*s/(a*a);
    frt+=2*(zd*s/a+z*(sd/a-s*ad/(a*a)));ftt+=2*(zd*zd+z*zdd);
  }
  if(!std::isfinite(fr)||!(fr>0))return Return::Failure(Core::StatusCode::NumericalFailure,
      "rotating outer contact ray has unresolved implicit radial derivative");
  const long double speed=-ft/fr;
  out.radialSpeedMPerS=static_cast<double>(speed);
  out.radialAccelerationMPerS2=static_cast<double>(-(ftt+2*frt*speed+frr*speed*speed)/fr);
  out.normalSpeedMPerS=out.radialSpeedMPerS*out.incidence;
  if(!std::isfinite(out.radialSpeedMPerS)||!std::isfinite(out.radialAccelerationMPerS2)||
      !std::isfinite(out.normalSpeedMPerS))return Return::Failure(
          Core::StatusCode::NumericalFailure,"rotating contact ray derivatives overflow SI output");
  return Return::Success(out);
}

Core::Result<IdealContactResiduals> EvaluateIdealContactResiduals(
    const Primitive& left,const Primitive& right,const FaceGeometry& face,const Thermodynamics& gas) {
  using Return=Core::Result<IdealContactResiduals>;
  const auto ul=ToConserved(left,gas),ur=ToConserved(right,gas);
  if(!ul.ok())return Return::Failure(ul.status.code,ul.status.message);
  if(!ur.ok())return Return::Failure(ur.status.code,ur.status.message);
  if(!Finite(face.normalLeftToRight)||!std::isfinite(face.normalSpeedMPerS)||
      std::abs(Dot(face.normalLeftToRight,face.normalLeftToRight)-1)>
          128*std::numeric_limits<double>::epsilon())return Return::Failure(
              Core::StatusCode::InvalidConfiguration,"contact diagnostics require a finite physical unit normal/speed");
  const auto n=face.normalLeftToRight;const double wn=face.normalSpeedMPerS;
  const double mu=gas.magneticPermeabilityNPerA2;
  const auto traction=[&](const Primitive& p) {
    return (p.pressurePa+Dot(p.magneticFieldT,p.magneticFieldT)/(2*mu))*n-
        (Dot(p.magneticFieldT,n)/mu)*p.magneticFieldT;
  };
  const auto energy=[&](const Primitive& p,const ConservedDensity& u) {
    const double pt=p.pressurePa+Dot(p.magneticFieldT,p.magneticFieldT)/(2*mu);
    return (u.totalEnergyDensityJPerM3+pt)*Dot(p.velocityMPerS,n)-
        Dot(p.magneticFieldT,p.velocityMPerS)*Dot(p.magneticFieldT,n)/mu-
        wn*u.totalEnergyDensityJPerM3;
  };
  IdealContactResiduals out;
  out.leftRelativeMassFluxKgM2S=left.densityKgM3*(Dot(left.velocityMPerS,n)-wn);
  out.rightRelativeMassFluxKgM2S=right.densityKgM3*(Dot(right.velocityMPerS,n)-wn);
  out.velocityJumpMPerS=right.velocityMPerS-left.velocityMPerS;
  out.magneticJumpT=right.magneticFieldT-left.magneticFieldT;
  out.normalMagneticJumpT=Dot(out.magneticJumpT,n);
  out.pressureJumpPa=right.pressurePa-left.pressurePa;
  out.tractionJumpPa=traction(right)-traction(left);
  const auto el=-1.0*Cross(left.velocityMPerS-wn*n,left.magneticFieldT);
  const auto er=-1.0*Cross(right.velocityMPerS-wn*n,right.magneticFieldT);
  out.tangentialElectricJumpVPerM=(er-el)-Dot(er-el,n)*n;
  out.relativeTotalEnergyFluxJumpWPerM2=energy(right,ur.value)-energy(left,ul.value);
  if(!std::isfinite(out.leftRelativeMassFluxKgM2S)||!std::isfinite(out.rightRelativeMassFluxKgM2S)||
      !Finite(out.velocityJumpMPerS)||!Finite(out.magneticJumpT)||
      !std::isfinite(out.normalMagneticJumpT)||!std::isfinite(out.pressureJumpPa)||
      !Finite(out.tractionJumpPa)||!Finite(out.tangentialElectricJumpVPerM)||
      !std::isfinite(out.relativeTotalEnergyFluxJumpWPerM2))return Return::Failure(
          Core::StatusCode::NumericalFailure,"ideal-contact dimensional residuals overflow SI output");
  return Return::Success(out);
}

} } } // namespace SEP::CoronaSwcme::CrossRay
