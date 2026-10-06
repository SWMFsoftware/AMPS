#include "sep_corona_swcme/piston_ambient.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/vector_math.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace {

using CoronalCME::Cross;
using CoronalCME::Dot;
using CoronalCME::Norm;
using CoronalCME::Unit;
using CoronalCME::Vec3;

struct ProjectedPrimitive {
  double rho=0,p=0,e=0,u=0,br=0;
  std::array<double,2> ut{},bt{};
  AmbientRegion region=AmbientRegion::PfssClosed;
  int sector=0;
};

bool SameBranch(const ProjectedPrimitive& a,const ProjectedPrimitive& b) {
  return a.region==b.region&&a.sector==b.sector;
}

ProjectedPrimitive Combine(const ProjectedPrimitive& a,double wa,
    const ProjectedPrimitive& b,double wb,const ProjectedPrimitive& c,
    double wc,const ProjectedPrimitive& d,double wd,double denominator) {
  ProjectedPrimitive out;
  const auto sum=[&](double ProjectedPrimitive::*member) {
    return (wa*a.*member+wb*b.*member+wc*c.*member+wd*d.*member)/denominator;
  };
  out.rho=sum(&ProjectedPrimitive::rho);out.p=sum(&ProjectedPrimitive::p);
  out.e=sum(&ProjectedPrimitive::e);out.u=sum(&ProjectedPrimitive::u);
  out.br=sum(&ProjectedPrimitive::br);
  for(int i=0;i<2;++i) {
    out.ut[i]=(wa*a.ut[i]+wb*b.ut[i]+wc*c.ut[i]+wd*d.ut[i])/denominator;
    out.bt[i]=(wa*a.bt[i]+wb*b.bt[i]+wc*c.bt[i]+wd*d.bt[i])/denominator;
  }
  return out;
}

} // namespace

Core::Result<std::shared_ptr<const PistonAmbientProjection>>
PistonAmbientProjection::Create(std::shared_ptr<const AmbientModel> ambient,
    Vec3 rayDirection,double frozenEpochS,double derivativeRelativeStep) {
  using Return=Core::Result<std::shared_ptr<const PistonAmbientProjection>>;
  if(!ambient||!std::isfinite(rayDirection.x)||
      !std::isfinite(rayDirection.y)||!std::isfinite(rayDirection.z)||
      !(Norm(rayDirection)>0)||!std::isfinite(frozenEpochS)||
      !(derivativeRelativeStep>0&&derivativeRelativeStep<=1e-2))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "piston ambient projection requires an ambient, finite ray/epoch and "
        "a derivative step in (0,1e-2]");
  if(frozenEpochS<ambient->Event().support.startS||
      frozenEpochS>ambient->Event().support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "piston ambient projection epoch is outside event support");

  std::shared_ptr<PistonAmbientProjection> out(new PistonAmbientProjection);
  out->ambient_=std::move(ambient);out->ray_=Unit(rayDirection);
  out->epochS_=frozenEpochS;out->relativeStep_=derivativeRelativeStep;
  // The reference axis changes only near the HCI pole and is selected once,
  // never as a function of radius.  Thus signed transverse components cannot
  // flip because a sampling stencil crossed an arbitrary basis branch.
  const Vec3 reference=std::abs(out->ray_.z)<0.9?Vec3{0,0,1}:Vec3{1,0,0};
  out->tangent1_=Unit(Cross(reference,out->ray_));
  out->tangent2_=Cross(out->ray_,out->tangent1_);
  return Return::Success(std::move(out));
}

Core::Result<PistonAmbientRayState> PistonAmbientProjection::Evaluate(
    double radius) const {
  using Return=Core::Result<PistonAmbientRayState>;
  const auto sample=[&](double r)->Core::Result<ProjectedPrimitive> {
    const auto primitive=ambient_->Evaluate(r*ray_,epochS_);
    if(!primitive.ok())return Core::Result<ProjectedPrimitive>::Failure(
        primitive.status.code,primitive.status.message);
    ProjectedPrimitive out;
    out.rho=primitive.value.plasma.massDensityKgM3;
    out.p=primitive.value.plasma.pressurePa;
    out.e=out.p/((ambient_->Event().composition.gammaAdiabatic-1)*out.rho);
    out.u=Dot(primitive.value.velocityMPerS,ray_);
    out.ut={Dot(primitive.value.velocityMPerS,tangent1_),
        Dot(primitive.value.velocityMPerS,tangent2_)};
    out.br=Dot(primitive.value.magneticFieldT,ray_);
    out.bt={Dot(primitive.value.magneticFieldT,tangent1_),
        Dot(primitive.value.magneticFieldT,tangent2_)};
    out.region=primitive.value.region;out.sector=primitive.value.magneticSector;
    return Core::Result<ProjectedPrimitive>::Success(out);
  };
  if(!std::isfinite(radius)||
      radius<ambient_->Event().support.firstValidPlasmaRadiusM||
      radius>ambient_->Event().support.coverageRadiusM)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "piston ambient radius is outside declared plasma coverage");
  const auto center=sample(radius);
  if(!center.ok())return Return::Failure(center.status.code,center.status.message);

  // A four-point one-sided formula is used when a central stencil would cross
  // a PFSS/Parker or topology/sector interface.  The step is reduced before
  // choosing a side.  No derivative ever averages two physical branches.
  ProjectedPrimitive derivative;
  bool resolved=false;
  double step=std::max(1.0,relativeStep_*radius);
  for(int attempt=0;attempt<=12&&!resolved;++attempt,step*=0.5) {
    const auto m2=sample(radius-2*step),m1=sample(radius-step),
        p1=sample(radius+step),p2=sample(radius+2*step);
    if(m2.ok()&&m1.ok()&&p1.ok()&&p2.ok()&&
        SameBranch(center.value,m2.value)&&SameBranch(center.value,m1.value)&&
        SameBranch(center.value,p1.value)&&SameBranch(center.value,p2.value)) {
      derivative=Combine(m2.value,1,m1.value,-8,p1.value,8,p2.value,-1,
          12*step);
      resolved=true;
      break;
    }
    const auto p3=sample(radius+3*step);
    if(p1.ok()&&p2.ok()&&p3.ok()&&SameBranch(center.value,p1.value)&&
        SameBranch(center.value,p2.value)&&SameBranch(center.value,p3.value)) {
      derivative=Combine(center.value,-11,p1.value,18,p2.value,-9,p3.value,2,
          6*step);
      resolved=true;
      break;
    }
    const auto m3=sample(radius-3*step);
    if(m1.ok()&&m2.ok()&&m3.ok()&&SameBranch(center.value,m1.value)&&
        SameBranch(center.value,m2.value)&&SameBranch(center.value,m3.value)) {
      derivative=Combine(center.value,11,m1.value,-18,m2.value,9,m3.value,-2,
          6*step);
      resolved=true;
    }
  }
  if(!resolved)return Return::Failure(Core::StatusCode::NumericalFailure,
      "piston ambient radial derivative cannot remain on one physical branch");

  const auto& a=center.value;
  const double mu0=CoronalCME::Constants::kVacuumPermeabilityHPerM;
  const double gamma=ambient_->Event().composition.gammaAdiabatic;
  const double bt2=a.bt[0]*a.bt[0]+a.bt[1]*a.bt[1];
  const double b2=bt2+a.br*a.br;
  const double dMagneticPressure=(a.bt[0]*derivative.bt[0]+
      a.bt[1]*derivative.bt[1])/mu0;
  PistonAmbientRayState out;
  out.radiusM=radius;out.densityKgM3=a.rho;out.pressurePa=a.p;
  out.specificInternalEnergyJPerKg=a.e;out.radialVelocityMPerS=a.u;
  out.radialVelocityGradientPerS=derivative.u;
  out.materialRadialAccelerationMPerS2=a.u*derivative.u;
  out.transverseVelocityMPerS=a.ut;out.radialMagneticFieldT=a.br;
  out.transverseMagneticFieldT=a.bt;out.region=a.region;
  out.magneticSector=a.sector;
  out.reducedFastSpeedMPerS=std::sqrt((gamma*a.p+bt2/mu0)/a.rho);
  const double sound2=gamma*a.p/a.rho,alfven2=b2/(mu0*a.rho),
      radialAlfven2=a.br*a.br/(mu0*a.rho);
  const double discriminant=std::max(0.0,
      (sound2+alfven2)*(sound2+alfven2)-4*sound2*radialAlfven2);
  out.canonicalFastSpeedMPerS=std::sqrt(
      0.5*(sound2+alfven2+std::sqrt(discriminant)));

  const double pressureGradient=derivative.p+dMagneticPressure;
  const double netAcceleration=a.u*derivative.u+pressureGradient/a.rho+
      bt2/(mu0*a.rho*radius);
  out.gravityAccelerationMPerS2=-
      CoronalCME::Constants::kSolarGravitationalParameterM3PerS2/(radius*radius);
  out.ambientMaintainingAccelerationMPerS2=
      netAcceleration-out.gravityAccelerationMPerS2;
  out.source.radialAccelerationMPerS2=netAcceleration;
  out.source.heatingWPerM3=a.rho*a.u*(derivative.e-
      a.p*derivative.rho/(a.rho*a.rho));
  for(int i=0;i<2;++i)out.source.transverseInvariantRate[i]=a.u*(
      derivative.bt[i]/(a.rho*radius)-
      a.bt[i]*derivative.rho/(a.rho*a.rho*radius)-
      a.bt[i]/(a.rho*radius*radius));

  const double values[]={out.densityKgM3,out.pressurePa,
      out.specificInternalEnergyJPerKg,out.reducedFastSpeedMPerS,
      out.canonicalFastSpeedMPerS,out.source.radialAccelerationMPerS2,
      out.source.heatingWPerM3,out.source.transverseInvariantRate[0],
      out.source.transverseInvariantRate[1]};
  for(double value:values)if(!std::isfinite(value))return Return::Failure(
      Core::StatusCode::NumericalFailure,
      "piston ambient projection produced a nonfinite primitive or source");
  return Return::Success(out);
}

Core::Result<PistonAmbientRayState> PistonAmbientProjection::EvaluatePrimitive(
    double radius) const {
  using Return=Core::Result<PistonAmbientRayState>;
  if(!std::isfinite(radius)||
      radius<ambient_->Event().support.firstValidPlasmaRadiusM||
      radius>ambient_->Event().support.coverageRadiusM)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "piston ambient primitive radius is outside plasma coverage");
  const auto primitive=ambient_->Evaluate(radius*ray_,epochS_);
  if(!primitive.ok())return Return::Failure(
      primitive.status.code,primitive.status.message);
  PistonAmbientRayState out;
  out.radiusM=radius;
  out.densityKgM3=primitive.value.plasma.massDensityKgM3;
  out.pressurePa=primitive.value.plasma.pressurePa;
  out.specificInternalEnergyJPerKg=out.pressurePa/
      ((ambient_->Event().composition.gammaAdiabatic-1)*out.densityKgM3);
  out.radialVelocityMPerS=Dot(primitive.value.velocityMPerS,ray_);
  out.transverseVelocityMPerS={Dot(primitive.value.velocityMPerS,tangent1_),
      Dot(primitive.value.velocityMPerS,tangent2_)};
  out.radialMagneticFieldT=Dot(primitive.value.magneticFieldT,ray_);
  out.transverseMagneticFieldT={Dot(primitive.value.magneticFieldT,tangent1_),
      Dot(primitive.value.magneticFieldT,tangent2_)};
  out.region=primitive.value.region;out.magneticSector=primitive.value.magneticSector;
  const double mu0=CoronalCME::Constants::kVacuumPermeabilityHPerM;
  const double bt2=out.transverseMagneticFieldT[0]*
      out.transverseMagneticFieldT[0]+out.transverseMagneticFieldT[1]*
      out.transverseMagneticFieldT[1];
  const double b2=bt2+out.radialMagneticFieldT*out.radialMagneticFieldT;
  const double sound2=ambient_->Event().composition.gammaAdiabatic*
      out.pressurePa/out.densityKgM3;
  out.reducedFastSpeedMPerS=std::sqrt((sound2*out.densityKgM3+bt2/mu0)/
      out.densityKgM3);
  const double alfven2=b2/(mu0*out.densityKgM3);
  const double radialAlfven2=out.radialMagneticFieldT*
      out.radialMagneticFieldT/(mu0*out.densityKgM3);
  const double discriminant=std::max(0.0,
      (sound2+alfven2)*(sound2+alfven2)-4*sound2*radialAlfven2);
  out.canonicalFastSpeedMPerS=std::sqrt(
      0.5*(sound2+alfven2+std::sqrt(discriminant)));
  return Return::Success(out);
}

Core::Result<std::shared_ptr<const PistonAmbientSourceTable>>
PistonAmbientSourceTable::Create(
    std::shared_ptr<const PistonAmbientProjection> projection,
    double minimumRadius,double maximumRadius,int points,
    bool logarithmicSpacing) {
  using Return=Core::Result<std::shared_ptr<const PistonAmbientSourceTable>>;
  if(!projection||!std::isfinite(minimumRadius)||
      !std::isfinite(maximumRadius)||!(maximumRadius>minimumRadius)||points<17)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "piston ambient source table requires a projection, ordered radii "
        "and at least 17 samples");
  std::shared_ptr<PistonAmbientSourceTable> out(new PistonAmbientSourceTable);
  out->projection_=projection;
  out->radius_.reserve(points);out->source_.reserve(points);
  out->region_.reserve(points);out->sector_.reserve(points);
  for(int i=0;i<points;++i) {
    const double fraction=static_cast<double>(i)/(points-1);
    const double radius=logarithmicSpacing?
        std::exp((1-fraction)*std::log(minimumRadius)+
            fraction*std::log(maximumRadius)):
        minimumRadius+(maximumRadius-minimumRadius)*fraction;
    const auto sample=projection->Evaluate(radius);
    if(!sample.ok())return Return::Failure(
        sample.status.code,sample.status.message);
    out->radius_.push_back(radius);out->source_.push_back(sample.value.source);
    out->region_.push_back(sample.value.region);
    out->sector_.push_back(sample.value.magneticSector);
    const double gravity=std::abs(sample.value.gravityAccelerationMPerS2);
    out->maximumAmbientForceOverGravity_=std::max(
        out->maximumAmbientForceOverGravity_,
        std::abs(sample.value.ambientMaintainingAccelerationMPerS2)/gravity);
    out->maximumAbsoluteHeatingWPerM3_=std::max(
        out->maximumAbsoluteHeatingWPerM3_,
        std::abs(sample.value.source.heatingWPerM3));
  }
  return Return::Success(std::move(out));
}

Core::Result<PistonVolumeSource> PistonAmbientSourceTable::Evaluate(
    double radius) const {
  using Return=Core::Result<PistonVolumeSource>;
  if(!std::isfinite(radius)||radius<radius_.front()||radius>radius_.back())
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "piston ambient source request is outside the frozen radial table");
  auto upper=std::upper_bound(radius_.begin(),radius_.end(),radius);
  if(upper==radius_.begin())return Return::Success(source_.front());
  if(upper==radius_.end())return Return::Success(source_.back());
  const std::size_t right=static_cast<std::size_t>(upper-radius_.begin());
  const std::size_t left=right-1;
  if(region_[left]!=region_[right]||sector_[left]!=sector_[right]) {
    // Do not smear a source-surface/topology/current-sheet jump merely to use
    // a fast table.  Transition intervals are rare; resolve the authority on
    // demand with its one-sided physical-branch derivative instead.
    const auto resolved=projection_->Evaluate(radius);
    if(!resolved.ok())return Return::Failure(
        resolved.status.code,resolved.status.message);
    return Return::Success(resolved.value.source);
  }
  const double fraction=(radius-radius_[left])/(radius_[right]-radius_[left]);
  PistonVolumeSource out;
  const auto blend=[&](double a,double b) { return (1-fraction)*a+fraction*b; };
  out.radialAccelerationMPerS2=blend(
      source_[left].radialAccelerationMPerS2,
      source_[right].radialAccelerationMPerS2);
  out.heatingWPerM3=blend(source_[left].heatingWPerM3,
      source_[right].heatingWPerM3);
  for(int i=0;i<2;++i)out.transverseInvariantRate[i]=blend(
      source_[left].transverseInvariantRate[i],
      source_[right].transverseInvariantRate[i]);
  return Return::Success(out);
}

Core::Result<std::shared_ptr<const PistonAmbientTrajectory>>
PistonAmbientTrajectory::Create(
    std::shared_ptr<const PistonAmbientProjection> projection,
    double initialRadius,double start,double end,double maximumStep) {
  using Return=Core::Result<std::shared_ptr<const PistonAmbientTrajectory>>;
  if(!projection||!std::isfinite(initialRadius)||!std::isfinite(start)||
      !std::isfinite(end)||!std::isfinite(maximumStep)||!(initialRadius>0)||
      !(end>start)||!(maximumStep>0))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
          "ambient material trajectory requires finite ordered support");
  std::shared_ptr<PistonAmbientTrajectory> out(new PistonAmbientTrajectory);
  double time=start,radius=initialRadius;
  while(true) {
    const auto state=projection->Evaluate(radius);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    out->time_.push_back(time);out->radius_.push_back(radius);
    out->velocity_.push_back(state.value.radialVelocityMPerS);
    out->acceleration_.push_back(state.value.materialRadialAccelerationMPerS2);
    if(time>=end)break;
    const double dt=std::min(maximumStep,end-time);
    const auto velocity=[&](double r)->Core::Result<double> {
      const auto sample=projection->Evaluate(r);
      if(!sample.ok())return Core::Result<double>::Failure(
          sample.status.code,sample.status.message);
      return Core::Result<double>::Success(sample.value.radialVelocityMPerS);
    };
    const auto k1=velocity(radius);
    if(!k1.ok())return Return::Failure(k1.status.code,k1.status.message);
    const auto k2=velocity(radius+0.5*dt*k1.value);
    if(!k2.ok())return Return::Failure(k2.status.code,k2.status.message);
    const auto k3=velocity(radius+0.5*dt*k2.value);
    if(!k3.ok())return Return::Failure(k3.status.code,k3.status.message);
    const auto k4=velocity(radius+dt*k3.value);
    if(!k4.ok())return Return::Failure(k4.status.code,k4.status.message);
    radius+=dt*(k1.value+2*k2.value+2*k3.value+k4.value)/6;
    time+=dt;
  }
  return Return::Success(std::move(out));
}

Core::Result<CoronalCME::KinematicValue> PistonAmbientTrajectory::Evaluate(
    double time) const {
  using Return=Core::Result<CoronalCME::KinematicValue>;
  if(!std::isfinite(time)||time<time_.front()||time>time_.back())
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "ambient material trajectory query is outside time support");
  auto upper=std::upper_bound(time_.begin(),time_.end(),time);
  if(upper==time_.begin())return Return::Success(
      {radius_.front(),velocity_.front(),acceleration_.front()});
  if(upper==time_.end())return Return::Success(
      {radius_.back(),velocity_.back(),acceleration_.back()});
  const std::size_t right=static_cast<std::size_t>(upper-time_.begin());
  const std::size_t left=right-1;
  const double h=time_[right]-time_[left];
  const double s=(time-time_[left])/h;
  const double c0=radius_[left],c1=h*velocity_[left],
      c2=0.5*h*h*acceleration_[left];
  const double d=radius_[right]-c0-c1-c2;
  const double v=h*velocity_[right]-c1-2*c2;
  const double a=h*h*acceleration_[right]-2*c2;
  const double c3=10*d-4*v+0.5*a,c4=-15*d+7*v-a,
      c5=6*d-3*v+0.5*a;
  const double value=c0+s*(c1+s*(c2+s*(c3+s*(c4+s*c5))));
  const double first=(c1+s*(2*c2+s*(3*c3+s*(4*c4+s*5*c5))))/h;
  const double second=(2*c2+s*(6*c3+s*(12*c4+s*20*c5)))/(h*h);
  return Return::Success({value,first,second});
}

} } // namespace SEP::CoronaSwcme
