// BG3D-2 composition of maintained coronal kernels. This code is adapted from
// the independently audited Mars1 background prototype, but is owned here by
// the particle-free composite event and revalidated in Mars2.
#include "sep_corona_swcme/ambient_state.h"

#include "sep_coronal_cme/closed_field_plasma.h"
#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/flux_tube_wind.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace SEP { namespace CoronaSwcme { namespace {

using CoronalCME::FieldLineTopology;
using CoronalCME::HarmonicCoefficient;
using CoronalCME::Norm;
using CoronalCME::Unit;
using CoronalCME::Vec3;

Vec3 Rotate(Vec3 vector,double angle) {
  const double cosine=std::cos(angle),sine=std::sin(angle);
  return {cosine*vector.x-sine*vector.y,sine*vector.x+cosine*vector.y,vector.z};
}

Core::Result<std::vector<HarmonicCoefficient>> Harmonics(const std::string& text) {
  using Return=Core::Result<std::vector<HarmonicCoefficient>>;
  std::istringstream rows(text);
  std::string row;
  std::vector<HarmonicCoefficient> result;
  while(std::getline(rows,row,';')) {
    if(row.empty())continue;
    std::replace(row.begin(),row.end(),',',' ');
    HarmonicCoefficient coefficient;
    std::istringstream input(row);
    if(!(input>>coefficient.degree>>coefficient.order>>coefficient.cosineT>>
        coefficient.sineT))return Return::Failure(
            Core::StatusCode::InvalidConfiguration,
            "ambient harmonics require l,m,cosine_T,sine_T rows separated by semicolons");
    std::string trailing;
    if(input>>trailing||coefficient.degree>32)return Return::Failure(
        Core::StatusCode::InvalidConfiguration,"invalid ambient harmonic row");
    result.push_back(coefficient);
  }
  if(result.empty())return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "ambient harmonic table is empty");
  return Return::Success(std::move(result));
}

bool SameBranch(const AmbientPrimitive& a,const AmbientPrimitive& b) {
  return a.region==b.region&&a.magneticSector==b.magneticSector;
}

} // namespace

const char* Name(AmbientRegion region) noexcept {
  switch(region) {
    case AmbientRegion::PfssOpen:return "pfss-open";
    case AmbientRegion::PfssClosed:return "pfss-closed";
    case AmbientRegion::ParkerExterior:return "parker-exterior";
  }
  return "unknown";
}

bool ValidAmbientRegionSector(AmbientRegion region,int sector) noexcept {
  if(region==AmbientRegion::PfssClosed)return sector==0;
  return (region==AmbientRegion::PfssOpen||region==AmbientRegion::ParkerExterior)&&
      (sector==-1||sector==1);
}

Core::Result<std::shared_ptr<const AmbientModel>> AmbientModel::Create(
    std::shared_ptr<const EventConfiguration> event) {
  using Return=Core::Result<std::shared_ptr<const AmbientModel>>;
  if(!event)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "ambient model requires an immutable event");
  const auto valid=ValidateEventConfiguration(*event);
  if(!valid.ok())return Return::Failure(valid.code,valid.message);
  const auto harmonics=Harmonics(event->ambient.harmonics);
  if(!harmonics.ok())return Return::Failure(harmonics.status.code,harmonics.status.message);
  const auto pfss=CoronalCME::PfssHarmonics::Create(event->support.solarRadiusM,
      event->ambient.sourceSurfaceRadiusM,harmonics.value);
  if(!pfss.ok())return Return::Failure(pfss.status.code,pfss.status.message);

  std::shared_ptr<AmbientModel> model(new AmbientModel);
  model->event_=std::move(event);
  model->pfss_=pfss.value;
  const auto& c=model->event_->composition;
  model->ions_={{"H+",1,1,CoronalCME::Constants::kProtonMassKg,c.protonTemperatureK}};
  if(c.alphaToProtonNumberRatio>0)model->ions_.push_back({"He++",
      c.alphaToProtonNumberRatio,2,CoronalCME::Constants::kAlphaMassKg,
      c.alphaTemperatureK});
  const auto reference=CoronalCME::EvaluatePlasmaFromElectronDensity(
      model->event_->ambient.electronDensityAtReferenceM3,c.electronTemperatureK,
      model->ions_,c.includeElectronMass,c.gammaAdiabatic,0);
  if(!reference.ok())return Return::Failure(reference.status.code,reference.status.message);
  model->soundSquared_=reference.value.pressurePa/reference.value.massDensityKgM3;
  model->criticalRadiusM_=CoronalCME::Constants::kSolarGravitationalParameterM3PerS2/
      (2*model->soundSquared_);
  const auto referenceWind=CoronalCME::SolveRadialIsothermalParker(
      {model->event_->ambient.referenceRadiusM},std::sqrt(model->soundSquared_),
      model->criticalRadiusM_,1);
  if(!referenceWind.ok())return Return::Failure(
      referenceWind.status.code,referenceWind.status.message);
  const double referenceRadius=model->event_->ambient.referenceRadiusM;
  model->massFluxPerSr_=reference.value.massDensityKgM3*
      referenceWind.value[0].speedMPerS*referenceRadius*referenceRadius;

  const int count=model->event_->ambient.windTablePoints;
  const double begin=std::log(model->event_->support.solarRadiusM);
  const double end=std::log(model->event_->support.coverageRadiusM);
  std::vector<double> radii;
  radii.reserve(static_cast<std::size_t>(count));
  for(int i=0;i<count;++i)radii.push_back(std::exp(begin+(end-begin)*i/(count-1)));
  radii.front()=model->event_->support.solarRadiusM;
  radii.back()=model->event_->support.coverageRadiusM;
  const auto wind=CoronalCME::SolveRadialIsothermalParker(radii,
      std::sqrt(model->soundSquared_),model->criticalRadiusM_,model->massFluxPerSr_);
  if(!wind.ok())return Return::Failure(wind.status.code,wind.status.message);
  for(const auto& point:wind.value) {
    const double y=point.speedMPerS*point.speedMPerS/model->soundSquared_;
    model->logRadius_.push_back(std::log(point.radiusM));
    model->logSpeed_.push_back(std::log(point.speedMPerS));
    model->slope_.push_back(
        std::abs(point.radiusM/model->criticalRadiusM_-1)<1e-7?1:
        2*(1-model->criticalRadiusM_/point.radiusM)/(y-1));
  }
  model->winding_.resize(radii.size());
  for(std::size_t i=1;i<radii.size();++i)model->winding_[i]=
      model->winding_[i-1]+model->IntegrateWinding(radii[i-1],radii[i]);
  return Return::Success(std::move(model));
}

double AmbientModel::WindSpeed(double radiusM) const {
  const double x=std::log(radiusM);
  auto upper=std::upper_bound(logRadius_.begin(),logRadius_.end(),x);
  const std::size_t i=std::min(logRadius_.size()-2,
      upper==logRadius_.begin()?std::size_t(0):
      static_cast<std::size_t>(upper-logRadius_.begin()-1));
  const double interval=logRadius_[i+1]-logRadius_[i];
  const double t=(x-logRadius_[i])/interval;
  return std::exp((2*t*t*t-3*t*t+1)*logSpeed_[i]+
      (t*t*t-2*t*t+t)*interval*slope_[i]+
      (-2*t*t*t+3*t*t)*logSpeed_[i+1]+
      (t*t*t-t*t)*interval*slope_[i+1]);
}

double AmbientModel::IntegrateWinding(double beginM,double endM) const {
  const double source=event_->ambient.sourceSurfaceRadiusM;
  if(endM<=source)return 0;
  beginM=std::max(beginM,source);
  static const double nodes[]={-0.8611363115940526,-0.3399810435848563,
      0.3399810435848563,0.8611363115940526};
  static const double weights[]={0.3478548451374538,0.6521451548625461,
      0.6521451548625461,0.3478548451374538};
  double integral=0;
  for(int i=0;i<4;++i) {
    const double radius=0.5*(beginM+endM)+0.5*(endM-beginM)*nodes[i];
    integral+=weights[i]*event_->ambient.rotationRateRadPerS*
        (1-std::pow(source/radius,2))/WindSpeed(radius);
  }
  return 0.5*(endM-beginM)*integral;
}

double AmbientModel::Winding(double radiusM) const {
  const auto upper=std::upper_bound(logRadius_.begin(),logRadius_.end(),
      std::log(radiusM));
  const std::size_t i=std::min(logRadius_.size()-2,
      static_cast<std::size_t>(upper-logRadius_.begin()-1));
  return winding_[i]+IntegrateWinding(std::exp(logRadius_[i]),radiusM);
}

Core::Result<std::pair<FieldLineTopology,Vec3>> AmbientModel::Trace(Vec3 start) const {
  using Return=Core::Result<std::pair<FieldLineTopology,Vec3>>;
  const auto& a=event_->ambient;
  const double solar=event_->support.solarRadiusM;
  const double epsilon=1e-6*solar;
  for(double sign:{1.0,-1.0}) {
    const double seedRadius=Norm(start);
    Vec3 point=seedRadius<=solar+epsilon?((solar+2*epsilon)/seedRadius)*start:start;
    bool ended=false;
    for(int i=0;i<a.traceMaximumSteps;++i) {
      const double radius=Norm(point);
      if(radius>=a.sourceSurfaceRadiusM-epsilon) {
        // PfssHarmonics owns the closed spherical shell. Roundoff in scaling a
        // Cartesian direction can otherwise put a nominal boundary point a
        // few ulps outside it. Evaluate the explicit one-sided coronal limit;
        // this is a coordinate guard, not a magnetic-field floor.
        const double inside=a.sourceSurfaceRadiusM*
            (1-64*std::numeric_limits<double>::epsilon());
        return Return::Success({FieldLineTopology::OpenToOuterBoundary,
            (inside/radius)*point});
      }
      if(radius<=solar+epsilon){ended=true;break;}
      const double step=std::min({a.traceStepM,0.4*(radius-solar),
          0.4*(a.sourceSurfaceRadiusM-radius)});
      const auto direction=[&](Vec3 query)->Core::Result<Vec3> {
        const auto field=pfss_.EvaluateCartesian(query);
        if(!field.ok())return Core::Result<Vec3>::Failure(
            field.status.code,field.status.message);
        if(Norm(field.value)<=a.minimumMagneticFieldT)return Core::Result<Vec3>::Failure(
            Core::StatusCode::NumericalFailure,
            "ambient topology trace encountered an unresolved magnetic null");
        return Core::Result<Vec3>::Success(sign*Unit(field.value));
      };
      const auto k1=direction(point);
      if(!k1.ok())return Return::Failure(k1.status.code,k1.status.message);
      const auto k2=direction(point+0.5*step*k1.value);
      if(!k2.ok())return Return::Failure(k2.status.code,k2.status.message);
      const auto k3=direction(point+0.5*step*k2.value);
      if(!k3.ok())return Return::Failure(k3.status.code,k3.status.message);
      const auto k4=direction(point+step*k3.value);
      if(!k4.ok())return Return::Failure(k4.status.code,k4.status.message);
      point=point+(step/6)*(k1.value+2*k2.value+2*k3.value+k4.value);
    }
    if(!ended)return Return::Failure(Core::StatusCode::NumericalFailure,
        "ambient topology trace exceeded trace_maximum_steps");
  }
  return Return::Success({FieldLineTopology::ClosedBelowOuterBoundary,{}});
}

Core::Result<AmbientPrimitive> AmbientModel::Evaluate(Vec3 position,double epochS) const {
  using Return=Core::Result<AmbientPrimitive>;
  const auto& support=event_->support;
  const auto& ambient=event_->ambient;
  const auto& composition=event_->composition;
  const double radius=Norm(position);
  if(!std::isfinite(position.x)||!std::isfinite(position.y)||
      !std::isfinite(position.z)||!std::isfinite(epochS)||
      epochS<support.startS||epochS>support.endS||
      radius<support.firstValidPlasmaRadiusM||radius>support.coverageRadiusM)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "ambient position or epoch is outside declared coverage");
  const double rotation=ambient.rotationRateRadPerS*epochS;
  const Vec3 body=Rotate(position,-rotation);
  AmbientPrimitive out;
  double density=0;
  const double speed=WindSpeed(radius);
  if(radius>=ambient.sourceSurfaceRadiusM) {
    const Vec3 direction=Rotate(Unit(body),Winding(radius));
    const double inside=ambient.sourceSurfaceRadiusM*
        (1-64*std::numeric_limits<double>::epsilon());
    const Vec3 foot=(inside/Norm(direction))*direction;
    const auto source=pfss_.EvaluateCartesian(foot);
    if(!source.ok())return Return::Failure(source.status.code,source.status.message);
    const double radialField=CoronalCME::Dot(source.value,Unit(foot))*
        std::pow(ambient.sourceSurfaceRadiusM/radius,2);
    const Vec3 radial=position/radius;
    const Vec3 azimuth={-radial.y,radial.x,0};
    out.magneticFieldT=radialField*(radial-
        (ambient.rotationRateRadPerS*(radius-ambient.sourceSurfaceRadiusM*
        ambient.sourceSurfaceRadiusM/radius)/speed)*azimuth);
    out.velocityMPerS=speed*radial+(ambient.rotationRateRadPerS*
        ambient.sourceSurfaceRadiusM*ambient.sourceSurfaceRadiusM/radius)*azimuth;
    density=massFluxPerSr_/(speed*radius*radius);
    out.region=AmbientRegion::ParkerExterior;
    out.magneticSector=radialField>0?1:-1;
  } else {
    const auto field=pfss_.EvaluateCartesian(body);
    if(!field.ok())return Return::Failure(field.status.code,field.status.message);
    out.magneticFieldT=Rotate(field.value,rotation);
    const auto topology=Trace(body);
    if(!topology.ok())return Return::Failure(topology.status.code,topology.status.message);
    const Vec3 corotation={-ambient.rotationRateRadPerS*position.y,
        ambient.rotationRateRadPerS*position.x,0};
    if(topology.value.first==FieldLineTopology::OpenToOuterBoundary) {
      const auto source=pfss_.EvaluateCartesian(topology.value.second);
      if(!source.ok())return Return::Failure(source.status.code,source.status.message);
      const double sourceRadial=CoronalCME::Dot(source.value,Unit(topology.value.second));
      if(std::abs(sourceRadial)<=ambient.minimumMagneticFieldT)
        return Return::Failure(Core::StatusCode::NumericalFailure,
            "open tube has unresolved source-surface flux");
      const double massPerFlux=massFluxPerSr_/(ambient.sourceSurfaceRadiusM*
          ambient.sourceSurfaceRadiusM*std::abs(sourceRadial));
      density=massPerFlux*Norm(field.value)/speed;
      out.magneticSector=sourceRadial>0?1:-1;
      out.region=AmbientRegion::PfssOpen;
      out.velocityMPerS=corotation+(out.magneticSector*speed)*Unit(out.magneticFieldT);
    } else {
      const auto base=CoronalCME::EvaluatePlasmaFromElectronDensity(
          ambient.closedBaseElectronDensityM3,composition.electronTemperatureK,
          ions_,composition.includeElectronMass,composition.gammaAdiabatic,0);
      if(!base.ok())return Return::Failure(base.status.code,base.status.message);
      const double potential=CoronalCME::RotatingEffectivePotential(position,
          {0,0,ambient.rotationRateRadPerS})+
          CoronalCME::Constants::kSolarGravitationalParameterM3PerS2/
          support.solarRadiusM;
      const auto hydro=CoronalCME::IsothermalHydrostatic(
          base.value.massDensityKgM3,base.value.pressurePa,potential);
      if(!hydro.ok())return Return::Failure(hydro.status.code,hydro.status.message);
      density=hydro.value.densityKgM3;
      out.velocityMPerS=corotation;
      out.region=AmbientRegion::PfssClosed;
      out.magneticSector=0;
    }
  }
  const double magneticMagnitude=Norm(out.magneticFieldT);
  if(!(magneticMagnitude>ambient.minimumMagneticFieldT))return Return::Failure(
      Core::StatusCode::NumericalFailure,
      "ambient sample has an unresolved magnetic null; no floor is applied");
  const auto unit=CoronalCME::EvaluatePlasmaFromElectronDensity(1,
      composition.electronTemperatureK,ions_,composition.includeElectronMass,
      composition.gammaAdiabatic,magneticMagnitude);
  if(!unit.ok())return Return::Failure(unit.status.code,unit.status.message);
  const auto plasma=CoronalCME::EvaluatePlasmaFromElectronDensity(
      density/unit.value.massDensityKgM3,composition.electronTemperatureK,ions_,
      composition.includeElectronMass,composition.gammaAdiabatic,magneticMagnitude);
  if(!plasma.ok())return Return::Failure(plasma.status.code,plasma.status.message);
  out.plasma=plasma.value;
  out.protonTemperatureK=composition.protonTemperatureK;
  out.electronTemperatureK=composition.electronTemperatureK;
  out.electronPressurePa=out.plasma.electronNumberDensityM3*
      CoronalCME::Constants::kBoltzmannJPerK*composition.electronTemperatureK;
  if(!ValidAmbientRegionSector(out.region,out.magneticSector))return Return::Failure(
      Core::StatusCode::DataIntegrityFailure,"ambient region/sector category is invalid");
  return Return::Success(out);
}

Core::Result<AmbientState> AmbientModel::EvaluateWithDerivatives(
    Vec3 point,double epochS,std::uint64_t generation) const {
  using Return=Core::Result<AmbientState>;
  if(generation==0)return Return::Failure(Core::StatusCode::InvalidState,
      "ambient generation zero is reserved");
  const auto center=Evaluate(point,epochS);
  if(!center.ok())return Return::Failure(center.status.code,center.status.message);
  AmbientState out;
  out.primitive=center.value;
  out.epochS=epochS;
  out.generation=generation;
  out.eventIdentity=event_->physicsFingerprint;
  const double radius=Norm(point);
  for(int axis=0;axis<3;++axis) {
    double step=std::max(1.0,event_->ambient.gradientRelativeStep*radius);
    bool resolved=false;
    for(int attempt=0;attempt<=12&&!resolved;++attempt,step*=0.5) {
      Vec3 delta;
      if(axis==0)delta.x=step;else if(axis==1)delta.y=step;else delta.z=step;
      const auto plus=Evaluate(point+delta,epochS);
      const auto minus=Evaluate(point-delta,epochS);
      const bool goodPlus=plus.ok()&&SameBranch(out.primitive,plus.value);
      const bool goodMinus=minus.ok()&&SameBranch(out.primitive,minus.value);
      Vec3 magneticDerivative,velocityDerivative;
      if(goodPlus&&goodMinus) {
        magneticDerivative=(plus.value.magneticFieldT-minus.value.magneticFieldT)/(2*step);
        velocityDerivative=(plus.value.velocityMPerS-minus.value.velocityMPerS)/(2*step);
        out.stencils[axis]=DifferenceStencil::CentralSecondOrder;
        resolved=true;
      } else for(int sign:{1,-1}) {
        if((sign==1&&!goodPlus)||(sign==-1&&!goodMinus))continue;
        const auto twice=Evaluate(point+2*sign*delta,epochS);
        if(!twice.ok()||!SameBranch(out.primitive,twice.value))continue;
        const auto& first=sign==1?plus.value:minus.value;
        magneticDerivative=(-3*out.primitive.magneticFieldT+
            4*first.magneticFieldT-twice.value.magneticFieldT)/(2*sign*step);
        velocityDerivative=(-3*out.primitive.velocityMPerS+
            4*first.velocityMPerS-twice.value.velocityMPerS)/(2*sign*step);
        out.stencils[axis]=sign==1?DifferenceStencil::ForwardSecondOrder:
            DifferenceStencil::BackwardSecondOrder;
        resolved=true;
        break;
      }
      if(resolved) {
        const double magnetic[]={magneticDerivative.x,magneticDerivative.y,
            magneticDerivative.z};
        const double velocity[]={velocityDerivative.x,velocityDerivative.y,
            velocityDerivative.z};
        for(int component=0;component<3;++component) {
          if(!std::isfinite(magnetic[component])||!std::isfinite(velocity[component]))
            return Return::Failure(Core::StatusCode::NumericalFailure,
                "ambient derivative is nonfinite");
          out.gradientB[3*component+axis]=magnetic[component];
          out.gradientU[3*component+axis]=velocity[component];
        }
      }
    }
    if(!resolved)return Return::Failure(Core::StatusCode::NumericalFailure,
        "ambient derivative cannot resolve a one-sided physical branch");
  }
  const Vec3 magnetic=out.primitive.magneticFieldT;
  const double squared=CoronalCME::Dot(magnetic,magnetic);
  if(!(squared>0))return Return::Failure(Core::StatusCode::NumericalFailure,
      "ambient logarithmic magnetic gradient lacks support");
  out.gradientLogBPerM={
      (magnetic.x*out.gradientB[0]+magnetic.y*out.gradientB[3]+
       magnetic.z*out.gradientB[6])/squared,
      (magnetic.x*out.gradientB[1]+magnetic.y*out.gradientB[4]+
       magnetic.z*out.gradientB[7])/squared,
      (magnetic.x*out.gradientB[2]+magnetic.y*out.gradientB[5]+
       magnetic.z*out.gradientB[8])/squared};
  out.derivativesValid=true;
  return Return::Success(std::move(out));
}

} } // namespace SEP::CoronaSwcme
