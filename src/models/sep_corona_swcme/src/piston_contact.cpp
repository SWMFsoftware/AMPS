#include "sep_corona_swcme/piston_contact.h"

#include "sep_corona_swcme/ambient_state.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace {

using CoronalCME::Dot;
using CoronalCME::EllipsoidKinematics;
using CoronalCME::KinematicValue;
using CoronalCME::Norm;
using CoronalCME::Vec3;

constexpr double kAstronomicalUnitM=149597870700.0;

template<class T> Core::Result<T> Invalid(const std::string& message) {
  return Core::Result<T>::Failure(Core::StatusCode::InvalidConfiguration,message);
}

KinematicValue Apex(const EllipsoidKinematics& k) {
  return {k.centerDistanceM.value+k.radialSemiAxisM.value,
      k.centerDistanceM.firstDerivative+
          k.radialSemiAxisM.firstDerivative,
      k.centerDistanceM.secondDerivative+
          k.radialSemiAxisM.secondDerivative};
}

// Integral and first two time derivatives of the requested startup law.
// q_dot=S(tau) q_dot_presc and
// integral_0^tau S(z) dz = 2.5 tau^4-3 tau^5+tau^6.
// The integral is one half at tau=1, so the post-ramp linear continuation is
// q=q0+qdot_presc*(t-t0-T/2), with no position or velocity reset.
KinematicValue Ramped(double initial,double rate,double elapsed,double duration) {
  if(elapsed<=0)return {initial,0,0};
  if(elapsed>=duration)return {
      initial+rate*(elapsed-0.5*duration),rate,0};
  const double u=elapsed/duration;
  const double u2=u*u,u3=u2*u,u4=u3*u,u5=u4*u,u6=u5*u;
  const double integral=2.5*u4-3*u5+u6;
  const double s=10*u3-15*u4+6*u5;
  const double ds=30*u2-60*u3+30*u4;
  return {initial+rate*duration*integral,rate*s,rate*ds/duration};
}

EllipsoidKinematics Early(const PistonContactInput& c,double timeS) {
  const double elapsed=timeS-c.startS;
  return {
      Ramped(c.initialCenterDistanceM,c.centerRateMPerS,elapsed,
          c.startupRampDurationS),
      Ramped(c.initialRadialSemiAxisM,c.radialRateMPerS,elapsed,
          c.startupRampDurationS),
      Ramped(c.initialFirstLateralSemiAxisM,c.firstLateralRateMPerS,elapsed,
          c.startupRampDurationS),
      Ramped(c.initialSecondLateralSemiAxisM,c.secondLateralRateMPerS,elapsed,
          c.startupRampDurationS)};
}

Core::Result<KinematicValue> Dbm(const PistonContactInput& law,
    const KinematicValue& crossing,double elapsedS) {
  using Return=Core::Result<KinematicValue>;
  const double difference=crossing.firstDerivative-
      law.outerAmbientSpeedMPerS;
  const double magnitude=std::abs(difference);
  if(!(elapsedS>=0)||!std::isfinite(elapsedS))return Return::Failure(
      Core::StatusCode::OutOfDomain,"piston DBM elapsed time is invalid");
  if(law.outerDragCoefficientPerM==0||magnitude==0)return Return::Success({
      std::fma(crossing.firstDerivative,elapsedS,crossing.value),
      crossing.firstDerivative,0});
  const double x=law.outerDragCoefficientPerM*magnitude*elapsedS;
  if(!std::isfinite(x)||!(1+x>0))return Return::Failure(
      Core::StatusCode::OutOfDomain,"piston DBM left its finite branch");
  const double velocityDifference=difference/(1+x);
  const double speed=law.outerAmbientSpeedMPerS+velocityDifference;
  const double distance=std::copysign(
      std::log1p(x)/law.outerDragCoefficientPerM,difference);
  return Return::Success({crossing.value+law.outerAmbientSpeedMPerS*elapsedS+
      distance,speed,-law.outerDragCoefficientPerM*velocityDifference*
      std::abs(velocityDifference)});
}

KinematicValue Scale(const KinematicValue& scale,double reference) {
  return {reference*scale.value,reference*scale.firstDerivative,
      reference*scale.secondDerivative};
}

KinematicValue Blend(const KinematicValue& a,const KinematicValue& b,
    double w,double dw,double ddw) {
  // This is one differentiated position curve.  The terms involving w_dot
  // and w_ddot are essential: blending precomputed velocities would not be
  // the derivative of the published contact and would reintroduce interface
  // leakage at the corona/SWCME handoff.
  return {(1-w)*a.value+w*b.value,
      (1-w)*a.firstDerivative+w*b.firstDerivative+dw*(b.value-a.value),
      (1-w)*a.secondDerivative+w*b.secondDerivative+
          2*dw*(b.firstDerivative-a.firstDerivative)+ddw*(b.value-a.value)};
}

std::array<double,3> Local(const CoronalCME::RadialPrincipalBasis& b,Vec3 x) {
  return {Dot(x,b.radial),Dot(x,b.firstLateral),Dot(x,b.secondLateral)};
}

} // namespace

const char* Name(PistonContactPhase phase) noexcept {
  switch(phase) {
    case PistonContactPhase::StartupRamp:return "startup-ramp";
    case PistonContactPhase::CoronalAnalytic:return "coronal-analytic";
    case PistonContactPhase::HandoffTransition:return "handoff-transition";
    case PistonContactPhase::SwcmeOuter:return "swcme-outer";
  }
  return "unknown";
}

const char* Name(PistonRayDisposition disposition) noexcept {
  switch(disposition) {
    case PistonRayDisposition::Supported:return "supported";
    case PistonRayDisposition::NoIntersection:return "no-intersection";
    case PistonRayDisposition::NonpositiveIntersection:
      return "nonpositive-intersection";
    case PistonRayDisposition::UnsupportedFlank:return "unsupported-flank";
    case PistonRayDisposition::BelowAmbientSupport:
      return "below-ambient-support";
  }
  return "unknown";
}

Core::Result<std::shared_ptr<const PistonContactModel>>
PistonContactModel::Create(std::shared_ptr<const EventConfiguration> event) {
  using Return=Core::Result<std::shared_ptr<const PistonContactModel>>;
  if(!event||!event->pistonContact.enabled||
      event->sheathModel!="per-ray-lagrangian-piston-v1")
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "piston contact requires the Level-B event selector and contact asset");
  const auto basis=CoronalCME::BuildRadialPrincipalBasis(
      event->pistonContact.latitudeRad,event->pistonContact.longitudeRad,
      event->pistonContact.lateralTiltRad);
  if(!basis.ok())return Return::Failure(basis.status.code,basis.status.message);
  std::shared_ptr<PistonContactModel> model(new PistonContactModel);
  model->event_=std::move(event);
  model->basis_=basis.value;

  const auto& c=model->event_->pistonContact;
  const double startApex=Apex(Early(c,c.startS)).value;
  const double endApex=Apex(Early(c,model->event_->support.endS)).value;
  if(!(startApex<c.handoffApexRadiusM&&
      c.handoffApexRadiusM<endApex))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
          "contact handoff radius is not crossed inside event support");
  double low=c.startS,high=model->event_->support.endS;
  // Monotone analytic motion makes bisection an unambiguous setup operation;
  // it is not used for contact speed or acceleration.  The bracket is reduced
  // well below the configured time root tolerance.
  for(int i=0;i<100;++i) {
    const double mid=0.5*(low+high);
    if(Apex(Early(c,mid)).value<c.handoffApexRadiusM)low=mid;else high=mid;
  }
  model->handoffTimeS_=0.5*(low+high);
  if(model->handoffTimeS_+c.handoffTransitionDurationS>=
      model->event_->support.endS)return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
          "contact handoff transition exceeds event support");

  const auto initial=model->At(c.startS);
  const auto final=model->At(model->event_->support.endS);
  if(!initial.ok()||!final.ok())return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "contact endpoints cannot be evaluated");
  const auto initialExtent=EvaluateRadialExtent(initial.value.ellipsoid,
      model->event_->support.solarRadiusM);
  const auto finalExtent=EvaluateRadialExtent(final.value.ellipsoid,
      model->event_->support.solarRadiusM);
  if(!initialExtent.ok()||!finalExtent.ok()||
      !initialExtent.value.intersectsSolarSurface||
      finalExtent.value.minimumRadiusM<=model->event_->support.solarRadiusM)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "generic contact does not realize declared attachment then detachment");
  if(final.value.apexRadiusM.value<kAstronomicalUnitM||
      final.value.apexRadiusM.value>model->event_->support.coverageRadiusM)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "generic contact does not propagate to 1 AU inside ambient coverage");

  const double exactRampMaximum=1.875*
      (c.centerRateMPerS+c.radialRateMPerS)/c.startupRampDurationS;
  if(exactRampMaximum>c.maximumApexAccelerationMPerS2*(1+1e-13))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "startup ramp exceeds the preregistered apex acceleration");
  return Return::Success(std::move(model));
}

Core::Result<PistonContactState> PistonContactModel::At(double timeS) const {
  using Return=Core::Result<PistonContactState>;
  const auto& support=event_->support;
  const auto& c=event_->pistonContact;
  if(!std::isfinite(timeS)||timeS<support.startS||timeS>support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "piston contact query is outside event support");
  const EllipsoidKinematics early=Early(c,timeS);
  if(timeS<=handoffTimeS_) {
    const auto phase=timeS<c.startS+c.startupRampDurationS?
        PistonContactPhase::StartupRamp:PistonContactPhase::CoronalAnalytic;
    return Return::Success({timeS,0,phase,early,Apex(early)});
  }

  const EllipsoidKinematics reference=Early(c,handoffTimeS_);
  const KinematicValue crossing=Apex(reference);
  const auto apex=Dbm(c,crossing,timeS-handoffTimeS_);
  if(!apex.ok())return Return::Failure(apex.status.code,apex.status.message);
  const KinematicValue scale={apex.value.value/crossing.value,
      apex.value.firstDerivative/crossing.value,
      apex.value.secondDerivative/crossing.value};
  const EllipsoidKinematics outer={
      Scale(scale,reference.centerDistanceM.value),
      Scale(scale,reference.radialSemiAxisM.value),
      Scale(scale,reference.firstLateralSemiAxisM.value),
      Scale(scale,reference.secondLateralSemiAxisM.value)};
  const double transitionEnd=handoffTimeS_+c.handoffTransitionDurationS;
  if(timeS>=transitionEnd)return Return::Success({timeS,1,
      PistonContactPhase::SwcmeOuter,outer,apex.value});

  const double duration=c.handoffTransitionDurationS;
  const double u=(timeS-handoffTimeS_)/duration;
  const double u2=u*u,u3=u2*u,u4=u3*u,u5=u4*u;
  const double w=10*u3-15*u4+6*u5;
  const double dw=(30*u2-60*u3+30*u4)/duration;
  const double ddw=(60*u-180*u2+120*u3)/(duration*duration);
  const EllipsoidKinematics state={
      Blend(early.centerDistanceM,outer.centerDistanceM,w,dw,ddw),
      Blend(early.radialSemiAxisM,outer.radialSemiAxisM,w,dw,ddw),
      Blend(early.firstLateralSemiAxisM,outer.firstLateralSemiAxisM,w,dw,ddw),
      Blend(early.secondLateralSemiAxisM,outer.secondLateralSemiAxisM,w,dw,ddw)};
  return Return::Success({timeS,w,PistonContactPhase::HandoffTransition,
      state,Apex(state)});
}

Core::Result<PistonRayState> PistonContactModel::EvaluateRay(
    Vec3 q,double timeS) const {
  using Return=Core::Result<PistonRayState>;
  const double magnitude=Norm(q);
  if(!std::isfinite(magnitude)||std::abs(magnitude-1)>1e-12)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "piston ray direction must be a finite HCI unit vector");
  const auto state=At(timeS);
  if(!state.ok())return Return::Failure(state.status.code,state.status.message);
  const auto local=Local(basis_,q);
  const auto& k=state.value.ellipsoid;
  const double a[3]={k.radialSemiAxisM.value,
      k.firstLateralSemiAxisM.value,k.secondLateralSemiAxisM.value};
  const double ad[3]={k.radialSemiAxisM.firstDerivative,
      k.firstLateralSemiAxisM.firstDerivative,
      k.secondLateralSemiAxisM.firstDerivative};
  const double add[3]={k.radialSemiAxisM.secondDerivative,
      k.firstLateralSemiAxisM.secondDerivative,
      k.secondLateralSemiAxisM.secondDerivative};
  const double center[3]={k.centerDistanceM.value,0,0};
  const double centerRate[3]={k.centerDistanceM.firstDerivative,0,0};
  const double centerAcceleration[3]={k.centerDistanceM.secondDerivative,0,0};

  double alpha=0,beta=0,gamma=-1;
  for(int i=0;i<3;++i) {
    const double inverse=1/(a[i]*a[i]);
    alpha+=local[i]*local[i]*inverse;
    beta-=2*local[i]*center[i]*inverse;
    gamma+=center[i]*center[i]*inverse;
  }
  PistonRayState out;
  out.timeS=timeS;
  double discriminant=beta*beta-4*alpha*gamma;
  const double roundoff=128*std::numeric_limits<double>::epsilon()*
      (beta*beta+std::abs(4*alpha*gamma));
  if(discriminant<0&&discriminant>=-roundoff)discriminant=0;
  if(discriminant<0) {
    out.disposition=PistonRayDisposition::NoIntersection;
    return Return::Success(out);
  }
  const double root=std::sqrt(discriminant);
  // The larger quadratic root is the unique outermost intersection requested
  // by the closure.  The smaller root may be another positive crossing of the
  // same convex ellipsoid and is intentionally not used as the piston face.
  const double radius=(-beta+root)/(2*alpha);
  if(!(radius>0)) {
    out.disposition=PistonRayDisposition::NonpositiveIntersection;
    return Return::Success(out);
  }

  double y[3]={radius*local[0]-center[0],radius*local[1],radius*local[2]};
  Vec3 gradient;
  for(int i=0;i<3;++i) {
    const double component=2*y[i]/(a[i]*a[i]);
    gradient=gradient+component*(i==0?basis_.radial:
        (i==1?basis_.firstLateral:basis_.secondLateral));
  }
  const double gradientMagnitude=Norm(gradient);
  if(!(gradientMagnitude>0))return Return::Failure(
      Core::StatusCode::NumericalFailure,
      "contact outer intersection has an unresolved normal");
  out.outwardNormal=gradient/gradientMagnitude;
  out.incidence=Dot(q,out.outwardNormal);
  out.radiusM=radius;
  out.positionM=radius*q;
  if(out.incidence<event_->pistonContact.minimumContactIncidence) {
    out.disposition=PistonRayDisposition::UnsupportedFlank;
    return Return::Success(out);
  }
  if(radius<event_->support.firstValidPlasmaRadiusM) {
    out.disposition=PistonRayDisposition::BelowAmbientSupport;
    return Return::Success(out);
  }

  // Differentiate F(R(t),t)=0.  F_t and F_tt below hold r fixed; this
  // distinction is what makes R_dot=-F_t/F_r an independent analytic contact
  // velocity rather than a velocity copied from a sheath cell.
  double fr=0,frr=0,ft=0,frt=0,ftt=0;
  for(int i=0;i<3;++i) {
    const double ai=a[i],yi=y[i],qi=local[i],ci1=centerRate[i],
        ci2=centerAcceleration[i],ai1=ad[i],ai2=add[i];
    const double ai2p=ai*ai,ai3=ai2p*ai,ai4=ai3*ai;
    fr+=2*yi*qi/ai2p;
    frr+=2*qi*qi/ai2p;
    ft+=-2*yi*ci1/ai2p-2*yi*yi*ai1/ai3;
    frt+=-2*ci1*qi/ai2p-4*yi*qi*ai1/ai3;
    ftt+=2*ci1*ci1/ai2p-2*yi*ci2/ai2p+
        8*yi*ci1*ai1/ai3-2*yi*yi*ai2/ai3+
        6*yi*yi*ai1*ai1/ai4;
  }
  if(!(fr>0))return Return::Failure(Core::StatusCode::NumericalFailure,
      "contact outer intersection violates positive radial incidence");
  out.radialSpeedMPerS=-ft/fr;
  out.radialAccelerationMPerS2=-(ftt+2*frt*out.radialSpeedMPerS+
      frr*out.radialSpeedMPerS*out.radialSpeedMPerS)/fr;
  out.normalSpeedMPerS=out.radialSpeedMPerS*out.incidence;
  out.disposition=PistonRayDisposition::Supported;
  return Return::Success(out);
}

Core::Result<PistonStartupDiagnostics>
PistonContactModel::CheckStartupCompatibility(const AmbientModel& ambient,
    const std::vector<Vec3>& directions) const {
  using Return=Core::Result<PistonStartupDiagnostics>;
  if(directions.empty())return Invalid<PistonStartupDiagnostics>(
      "startup compatibility requires a nonempty ray set");
  PistonStartupDiagnostics out;
  const double time=event_->pistonContact.startS;
  for(std::size_t i=0;i<directions.size();++i) {
    const auto contact=EvaluateRay(directions[i],time);
    if(!contact.ok())return Return::Failure(
        contact.status.code,contact.status.message);
    if(!contact.value.Supported()) {
      ++out.unsupportedRays;
      continue;
    }
    ++out.supportedRays;
    const auto state=ambient.Evaluate(contact.value.positionM,time);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    const auto& plasma=state.value.plasma;
    const double magnetic=Norm(state.value.magneticFieldT);
    const double cosine=std::abs(Dot(state.value.magneticFieldT,directions[i]))/
        magnetic;
    const double a2=plasma.soundSpeedMPerS*plasma.soundSpeedMPerS;
    const double va2=plasma.alfvenSpeedMPerS*plasma.alfvenSpeedMPerS;
    const double sum=a2+va2;
    const double fastSquared=0.5*(sum+
        std::sqrt(std::max(0.0,sum*sum-4*a2*va2*cosine*cosine)));
    const double mach=std::abs(Dot(state.value.velocityMPerS,directions[i]))/
        std::sqrt(fastSquared);
    if(mach>out.maximumAmbientRadialMach) {
      out.maximumAmbientRadialMach=mach;
      out.limitingRay=i;
    }
    if(mach>event_->pistonContact.startupMachTolerance)
      ++out.incompatibleRays;
  }
  if(out.supportedRays==0)return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "startup ray set contains no supported piston direction");
  out.compatible=out.incompatibleRays==0;
  return Return::Success(out);
}

} } // namespace SEP::CoronaSwcme
