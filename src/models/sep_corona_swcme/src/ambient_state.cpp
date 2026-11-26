// BG3D-2 composition of maintained coronal kernels. This code is adapted from
// the independently audited Mars1 background prototype, but is owned here by
// the particle-free composite event and revalidated in Mars2.
#include "sep_corona_swcme/ambient_state.h"

#include "ambient_derivatives.h"
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

// Subsonic root 0<y<1 of the isothermal Parker/Bernoulli relation
// y - ln y = rhs. g(z)=e^z-z-rhs, z=ln y, is strictly decreasing for z<0 with
// g(-rhs-1)=e^(-rhs-1)+1>0 and g(0)=1-rhs<0, so bisection in z brackets the
// unique root for every rhs>1 without the radial solver's y floor.
Core::Result<double> SubsonicParkerY(double rhs) {
  using Return=Core::Result<double>;
  if(!std::isfinite(rhs)||!(rhs>1))return Return::Failure(
      Core::StatusCode::NumericalFailure,
      "open-tube isothermal wind has no subsonic solution here: the PFSS tube "
      "would need a sonic point inside the source surface");
  double low=-rhs-1,high=0;
  for(int i=0;i<300&&high-low>4*std::numeric_limits<double>::epsilon()*(1-low);++i) {
    const double middle=0.5*(low+high);
    if(std::exp(middle)-middle-rhs>0)low=middle;else high=middle;
  }
  const double y=std::exp(0.5*(low+high));
  if(!(y>0)||!(y<1))return Return::Failure(Core::StatusCode::NumericalFailure,
      "open-tube isothermal wind root left the subsonic branch");
  return Return::Success(y);
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

namespace {

Core::Status ValidateDefinition(const AmbientDefinition& definition) {
  const auto& c=definition.composition;
  const auto& a=definition.ambient;
  const auto& s=definition.support;
  const bool composition=std::isfinite(c.gammaAdiabatic)&&c.gammaAdiabatic>1&&
      std::isfinite(c.electronTemperatureK)&&c.electronTemperatureK>0&&
      std::isfinite(c.protonTemperatureK)&&c.protonTemperatureK>0&&
      std::isfinite(c.alphaTemperatureK)&&c.alphaTemperatureK>0&&
      std::isfinite(c.alphaToProtonNumberRatio)&&
      c.alphaToProtonNumberRatio>=0;
  const bool ambient=std::isfinite(a.sourceSurfaceRadiusM)&&
      std::isfinite(a.referenceRadiusM)&&
      std::isfinite(a.electronDensityAtReferenceM3)&&
      std::isfinite(a.closedBaseElectronDensityM3)&&
      std::isfinite(a.rotationRateRadPerS)&&
      std::isfinite(a.minimumMagneticFieldT)&&
      std::isfinite(a.traceStepM)&&std::isfinite(a.gradientRelativeStep)&&
      a.sourceSurfaceRadiusM>s.solarRadiusM&&
      a.sourceSurfaceRadiusM<s.coverageRadiusM&&
      a.referenceRadiusM>=a.sourceSurfaceRadiusM&&
      a.referenceRadiusM<=s.coverageRadiusM&&
      a.electronDensityAtReferenceM3>0&&a.closedBaseElectronDensityM3>0&&
      a.minimumMagneticFieldT>=0&&a.traceStepM>0&&
      a.gradientRelativeStep>0&&a.traceMaximumSteps>0&&
      a.windTablePoints>=17&&!a.harmonics.empty();
  const bool support=std::isfinite(s.startS)&&std::isfinite(s.endS)&&
      std::isfinite(s.solarRadiusM)&&
      std::isfinite(s.firstValidPlasmaRadiusM)&&
      std::isfinite(s.coverageRadiusM)&&s.endS>s.startS&&
      s.solarRadiusM>0&&s.firstValidPlasmaRadiusM>=s.solarRadiusM&&
      s.firstValidPlasmaRadiusM<a.sourceSurfaceRadiusM&&
      s.coverageRadiusM>a.sourceSurfaceRadiusM;
  // Optional v1.3 closures. Earlier profiles leave both at their inactive
  // defaults; an active closure must be finite and lie inside its support.
  const double calibrationRadius=CoronalCME::Norm(a.magneticCalibrationPositionM);
  const bool wind=a.windSpeedAtReferenceMPerS==0||
      (std::isfinite(a.windSpeedAtReferenceMPerS)&&a.windSpeedAtReferenceMPerS>0);
  const bool calibration=!a.magneticCalibration||
      (std::isfinite(calibrationRadius)&&
       calibrationRadius>=s.firstValidPlasmaRadiusM&&
       calibrationRadius<=s.coverageRadiusM&&
       std::isfinite(a.magneticCalibrationEpochS)&&
       a.magneticCalibrationEpochS>=s.startS&&a.magneticCalibrationEpochS<=s.endS&&
       std::isfinite(a.magneticCalibrationMagnitudeT)&&
       a.magneticCalibrationMagnitudeT>0);
  // The tube closure replaces the open-shell density law; it cannot be
  // combined with the radial rho U r^2 alternative.
  const bool openClosure=!a.isothermalFluxTubeOpenWind||!a.sphericalOpenMassFlux;
  if(!composition||!ambient||!support||!wind||!calibration||!openClosure||
      definition.coordinateFrame.empty()||
      definition.physicsFingerprint.empty())
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
        "canonical ambient definition is incomplete or outside its physical domain");
  return Core::Status::Success();
}

} // namespace

Core::Result<std::shared_ptr<const AmbientModel>> AmbientModel::CreateResolved(
    const AmbientDefinition& definition,
    std::shared_ptr<const EventConfiguration> event) {
  using Return=Core::Result<std::shared_ptr<const AmbientModel>>;
  const auto valid=ValidateDefinition(definition);
  if(!valid.ok())return Return::Failure(valid.code,valid.message);
  const auto harmonics=Harmonics(definition.ambient.harmonics);
  if(!harmonics.ok())return Return::Failure(harmonics.status.code,harmonics.status.message);
  const auto pfss=CoronalCME::PfssHarmonics::Create(definition.support.solarRadiusM,
      definition.ambient.sourceSurfaceRadiusM,harmonics.value);
  if(!pfss.ok())return Return::Failure(pfss.status.code,pfss.status.message);

  std::shared_ptr<AmbientModel> model(new AmbientModel);
  model->event_=std::move(event);
  model->definition_=definition;
  model->pfss_=pfss.value;
  // A Parker-exterior query on the exact rotation axis repeatedly asks PFSS
  // for one of these identical source points, at every trajectory RK stage.
  // Compute its original canonical result once per immutable model. This is
  // memoization, not an axisymmetric field approximation. Keep the failure
  // result and signed-zero azimuthal limits; other source points still execute
  // EvaluateCartesian normally. No query outside its actual domain is repaired.
  model->CacheAxialSourceFields();
  const auto& c=model->definition_.composition;
  model->ions_={{"H+",1,1,CoronalCME::Constants::kProtonMassKg,c.protonTemperatureK}};
  if(c.alphaToProtonNumberRatio>0)model->ions_.push_back({"He++",
      c.alphaToProtonNumberRatio,2,CoronalCME::Constants::kAlphaMassKg,
      c.alphaTemperatureK});
  const auto reference=CoronalCME::EvaluatePlasmaFromElectronDensity(
      model->definition_.ambient.electronDensityAtReferenceM3,c.electronTemperatureK,
      model->ions_,c.includeElectronMass,c.gammaAdiabatic,0);
  if(!reference.ok())return Return::Failure(reference.status.code,reference.status.message);
  const double referenceRadius=model->definition_.ambient.referenceRadiusM;
  const double gm=CoronalCME::Constants::kSolarGravitationalParameterM3PerS2;
  const auto referenceSpeed=[&](double sound)->Core::Result<double> {
    const auto point=CoronalCME::SolveRadialIsothermalParker({referenceRadius},
        sound,gm/(2*sound*sound),1);
    if(!point.ok())return Core::Result<double>::Failure(point.status.code,point.status.message);
    return Core::Result<double>::Success(point.value[0].speedMPerS);
  };
  const double targetWind=model->definition_.ambient.windSpeedAtReferenceMPerS;
  if(targetWind>0) {
    // Parker (1958) transonic isothermal wind: U/c_w=F(r/r_c), r_c=GM/(2c_w^2),
    // with F increasing. U(r_ref)=c_w F(2 r_ref c_w^2/GM) is therefore strictly
    // increasing in c_w and the measured U(r_ref) has exactly one root. The
    // logarithmic bracket runs from r_c=2 r_ref (deeply subsonic at r_ref, yet
    // inside the Parker solver's resolved branch) to r_c=1e-4 r_ref; a target
    // outside it is a configuration failure, never a clamped wind.
    const double lowSound=std::sqrt(gm/(4*referenceRadius));
    double low=std::log(lowSound),high=std::log(100*lowSound);
    const auto lower=referenceSpeed(std::exp(low)),upper=referenceSpeed(std::exp(high));
    if(!lower.ok()||!upper.ok()||!(lower.value<targetWind)||!(upper.value>targetWind))
      return Return::Failure(Core::StatusCode::InvalidConfiguration,
          "reference wind speed is outside the isothermal Parker bracket");
    for(int i=0;i<200&&high-low>4*std::numeric_limits<double>::epsilon()*std::abs(high);++i) {
      const double middle=0.5*(low+high);
      const auto speed=referenceSpeed(std::exp(middle));
      if(!speed.ok())return Return::Failure(speed.status.code,speed.status.message);
      if(speed.value<targetWind)low=middle;else high=middle;
    }
    const double sound=std::exp(0.5*(low+high));
    const auto fitted=referenceSpeed(sound);
    if(!fitted.ok()||std::abs(fitted.value-targetWind)>1e-9*targetWind)
      return Return::Failure(Core::StatusCode::NumericalFailure,
          "isothermal Parker wind fit did not reproduce the reference speed");
    model->soundSquared_=sound*sound;
  } else {
    model->soundSquared_=reference.value.pressurePa/reference.value.massDensityKgM3;
  }
  model->criticalRadiusM_=gm/(2*model->soundSquared_);
  // The open-tube closure continues the exterior's subsonic branch down each
  // PFSS tube. If the exterior is already supersonic at the source surface,
  // the tube would need its own sonic point inside the shell, which this
  // closure does not solve; reject rather than mix branches.
  if(model->definition_.ambient.isothermalFluxTubeOpenWind&&
      !(model->criticalRadiusM_>model->definition_.ambient.sourceSurfaceRadiusM))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "isothermal flux-tube open wind requires the Parker critical radius above the source surface");
  const auto referenceWind=referenceSpeed(std::sqrt(model->soundSquared_));
  if(!referenceWind.ok())return Return::Failure(
      referenceWind.status.code,referenceWind.status.message);
  model->massFluxPerSr_=reference.value.massDensityKgM3*
      referenceWind.value*referenceRadius*referenceRadius;

  const int count=model->definition_.ambient.windTablePoints;
  const double begin=std::log(model->definition_.support.solarRadiusM);
  const double end=std::log(model->definition_.support.coverageRadiusM);
  std::vector<double> radii;
  radii.reserve(static_cast<std::size_t>(count));
  for(int i=0;i<count;++i)radii.push_back(std::exp(begin+(end-begin)*i/(count-1)));
  radii.front()=model->definition_.support.solarRadiusM;
  radii.back()=model->definition_.support.coverageRadiusM;
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
  const auto& a=model->definition_.ambient;
  if(a.magneticCalibration) {
    // B is linear in the harmonic coefficients in both the PFSS shell and the
    // Parker continuation, while U, rho and topology do not depend on their
    // amplitude. One factor therefore reproduces the declared |B| exactly.
    const auto unit=model->Evaluate(a.magneticCalibrationPositionM,
        a.magneticCalibrationEpochS);
    if(!unit.ok())return Return::Failure(unit.status.code,
        "magnetic calibration point: "+unit.status.message);
    const double magnitude=Norm(unit.value.magneticFieldT);
    model->magneticFactor_=a.magneticCalibrationMagnitudeT/magnitude;
    if(!std::isfinite(model->magneticFactor_)||!(model->magneticFactor_>0))
      return Return::Failure(Core::StatusCode::NumericalFailure,
          "magnetic calibration factor is not finite and positive");
    model->pfss_=model->pfss_.Scaled(model->magneticFactor_);
    model->CacheAxialSourceFields();
    const auto check=model->Evaluate(a.magneticCalibrationPositionM,
        a.magneticCalibrationEpochS);
    if(!check.ok()||std::abs(Norm(check.value.magneticFieldT)-
        a.magneticCalibrationMagnitudeT)>1e-12*a.magneticCalibrationMagnitudeT)
      return Return::Failure(Core::StatusCode::NumericalFailure,
          "calibrated magnetogram does not reproduce the declared |B|");
  }
  return Return::Success(std::move(model));
}

void AmbientModel::CacheAxialSourceFields() {
  const double axialInside=definition_.ambient.sourceSurfaceRadiusM*
      (1-64*std::numeric_limits<double>::epsilon());
  axialSourceRadiiM_[2]=axialInside;
  for(unsigned i=1;i<=2;++i) {
    axialSourceRadiiM_[2-i]=std::nextafter(axialSourceRadiiM_[3-i],0.0);
    axialSourceRadiiM_[2+i]=std::nextafter(axialSourceRadiiM_[1+i],
        std::numeric_limits<double>::infinity());
  }
  for(unsigned j=0;j<axialSourceRadiiM_.size();++j)
    for(unsigned i=0;i<8;++i)
      axialSourceFields_[8*j+i]=Field({
          i&1?-0.0:0.0,i&2?-0.0:0.0,i&4?-axialSourceRadiiM_[j]:axialSourceRadiiM_[j]});
}

Core::Result<std::shared_ptr<const AmbientModel>> AmbientModel::Create(
    std::shared_ptr<const EventConfiguration> event) {
  using Return=Core::Result<std::shared_ptr<const AmbientModel>>;
  if(!event)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "ambient model requires an immutable event");
  const auto valid=ValidateEventConfiguration(*event);
  if(!valid.ok())return Return::Failure(valid.code,valid.message);
  AmbientDefinition definition;
  definition.composition=event->composition;
  definition.ambient=event->ambient;
  definition.support=event->support;
  definition.coordinateFrame=event->coordinateFrame;
  definition.physicsFingerprint=event->physicsFingerprint;
  return CreateResolved(definition,std::move(event));
}

Core::Result<std::shared_ptr<const AmbientModel>> AmbientModel::Create(
    const AmbientDefinition& definition) {
  return CreateResolved(definition,nullptr);
}

double AmbientModel::SpeedBoundMPerS() const {
  double logMaximum=logSpeed_.empty()?0:logSpeed_.front();
  for(double value:logSpeed_)logMaximum=std::max(logMaximum,value);
  // Hermite interpolation between monotone nodes can exceed a node value by
  // a small fraction of one interval's increase; the 1% margin covers it.
  return 1.01*std::max(std::exp(logMaximum),std::sqrt(soundSquared_))+
      std::abs(definition_.ambient.rotationRateRadPerS)*definition_.ambient.sourceSurfaceRadiusM;
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
  const double source=definition_.ambient.sourceSurfaceRadiusM;
  if(endM<=source)return 0;
  beginM=std::max(beginM,source);
  static const double nodes[]={-0.8611363115940526,-0.3399810435848563,
      0.3399810435848563,0.8611363115940526};
  static const double weights[]={0.3478548451374538,0.6521451548625461,
      0.6521451548625461,0.3478548451374538};
  double integral=0;
  for(int i=0;i<4;++i) {
    const double radius=0.5*(beginM+endM)+0.5*(endM-beginM)*nodes[i];
    integral+=weights[i]*definition_.ambient.rotationRateRadPerS*
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

Core::Result<Vec3> AmbientModel::Field(Vec3 body) const {
  return definition_.ambient.isothermalFluxTubeOpenWind?pfss_.EvaluateFieldCartesian(body):
      pfss_.EvaluateCartesian(body);
}

Core::Result<std::pair<FieldLineTopology,Vec3>> AmbientModel::TraceOpenByRadius(Vec3 start) const {
  using Return=Core::Result<std::pair<FieldLineTopology,Vec3>>;
  const auto& a=definition_.ambient;
  const double inside=a.sourceSurfaceRadiusM*(1-64*std::numeric_limits<double>::epsilon());
  const double r0=Norm(start);
  if(!(r0>definition_.support.solarRadiusM)||!(r0<inside))return Trace(start);
  const auto initial=Field(start);
  if(!initial.ok())return Return::Failure(initial.status.code,initial.status.message);
  const double sign=CoronalCME::Dot(initial.value,start)>=0?1.0:-1.0;
  constexpr double kMinimumRadialCosine=0.2;
  bool valid=true;
  const auto slope=[&](Vec3 x)->Vec3 {
    if(!valid)return {};
    const double r=Norm(x);
    // RK stages of the last step can exceed the shell by roundoff only.
    if(r>inside)x=(inside/r)*x;
    const auto b=Field(x);
    if(!b.ok()){valid=false;return {};}
    const double radial=CoronalCME::Dot(b.value,x)/Norm(x);
    if(!(sign*radial>=kMinimumRadialCosine*Norm(b.value))){valid=false;return {};}
    return b.value/radial;
  };
  const int steps=std::max(1,static_cast<int>(std::ceil((inside-r0)/a.traceStepM)));
  const double h=(inside-r0)/steps;
  Vec3 x=start;
  for(int i=0;i<steps&&valid;++i) {
    const Vec3 k1=slope(x),k2=slope(x+0.5*h*k1),k3=slope(x+0.5*h*k2),k4=slope(x+h*k3);
    if(valid)x=x+(h/6)*(k1+2*k2+2*k3+k4);
  }
  if(!valid||!std::isfinite(Norm(x)))return Trace(start);
  return Return::Success({FieldLineTopology::OpenToOuterBoundary,(inside/Norm(x))*x});
}

Core::Result<std::pair<FieldLineTopology,Vec3>> AmbientModel::Trace(Vec3 start) const {
  using Return=Core::Result<std::pair<FieldLineTopology,Vec3>>;
  const auto& a=definition_.ambient;
  const double solar=definition_.support.solarRadiusM;
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
        const auto field=Field(query);
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
  const auto& support=definition_.support;
  const auto& ambient=definition_.ambient;
  const auto& composition=definition_.composition;
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
    // On the exact Cartesian rotation axis, every azimuthal rotation leaves
    // the position direction unchanged. Avoid integrating the Parker winding
    // angle for that identity operation: moving polar material lines query
    // this primitive at every RK stage and particle-boundary refinement.
    // This is exact for any PFSS magnetogram, not an axisymmetric-field or
    // constant-wind approximation. Off-axis points retain the full winding
    // integration. In both branches the same PFSS, magnetic-null, coverage,
    // plasma and finite-value checks below remain authoritative; no validity
    // test or physical component of the background is bypassed.
    const Vec3 direction=(body.x==0&&body.y==0)?Unit(body):
        Rotate(Unit(body),Winding(radius));
    const double inside=ambient.sourceSurfaceRadiusM*
        (1-64*std::numeric_limits<double>::epsilon());
    const Vec3 foot=(inside/Norm(direction))*direction;
    // Floating-point normalization can put an axial foot one ULP away from
    // inside. Reuse a result ONLY for an exactly equal source point, including
    // the signs of its zero coordinates. Every cached radius has its own
    // independently evaluated canonical result; an unmatched foot uses PFSS.
    // no coordinate quantization, radius tolerance or field interpolation is
    // introduced. All subsequent null, plasma and finite-value checks remain.
    const unsigned axialKey=(std::signbit(foot.x)?1U:0U)|
        (std::signbit(foot.y)?2U:0U)|(std::signbit(foot.z)?4U:0U);
    std::size_t axialIndex=axialSourceFields_.size();
    if(foot.x==0&&foot.y==0)for(std::size_t j=0;j<axialSourceRadiiM_.size();++j)
      if(std::abs(foot.z)==axialSourceRadiiM_[j]){axialIndex=8*j+axialKey;break;}
    const auto source=axialIndex<axialSourceFields_.size()?
        axialSourceFields_[axialIndex]:Field(foot);
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
    const auto field=Field(body);
    if(!field.ok())return Return::Failure(field.status.code,field.status.message);
    out.magneticFieldT=Rotate(field.value,rotation);
    const auto topology=ambient.isothermalFluxTubeOpenWind?TraceOpenByRadius(body):Trace(body);
    if(!topology.ok())return Return::Failure(topology.status.code,topology.status.message);
    const Vec3 corotation={-ambient.rotationRateRadPerS*position.y,
        ambient.rotationRateRadPerS*position.x,0};
    if(topology.value.first==FieldLineTopology::OpenToOuterBoundary) {
      const auto source=Field(topology.value.second);
      if(!source.ok())return Return::Failure(source.status.code,source.status.message);
      const double sourceRadial=CoronalCME::Dot(source.value,Unit(topology.value.second));
      if(std::abs(sourceRadial)<=ambient.minimumMagneticFieldT)
        return Return::Failure(Core::StatusCode::NumericalFailure,
            "open tube has unresolved source-surface flux");
      double alongSpeed=speed;
      if(ambient.isothermalFluxTubeOpenWind) {
        // Steady corotating-frame isothermal wind along this open tube (see
        // AmbientInput). The tube area is A=Phi/|B|; the Bernoulli constant
        // and the mass flux per unit flux K are those of the exterior Parker
        // wind continued from the tube's source-surface foot, so u, rho and
        // B are continuous at Rss. The rotating frame's centrifugal term is
        // omitted exactly as in the exterior Parker solution: below Rss it is
        // at most Omega^2 Rss^3/(2 GM) ~ 2e-4 of the gravitational potential.
        const double magnitude=Norm(field.value);
        const double rss=ambient.sourceSurfaceRadiusM;
        const double rc=criticalRadiusM_;
        const auto y=SubsonicParkerY(4*std::log(rss/rc)-3+
            2*std::log(std::abs(sourceRadial)/magnitude)+4*rc/radius);
        if(!y.ok())return Return::Failure(y.status.code,y.status.message);
        alongSpeed=std::sqrt(soundSquared_*y.value);
        density=massFluxPerSr_/(rss*rss*std::abs(sourceRadial))*magnitude/alongSpeed;
      } else if(ambient.sphericalOpenMassFlux) {
        // Radial Parker reference: rho U r^2 is the prescribed invariant.
        // This is a kinematic plasma closure and does not claim that a radial
        // wind is force-balanced across a nonradial PFSS field.
        density=massFluxPerSr_/(speed*radius*radius);
      } else {
        const double massPerFlux=massFluxPerSr_/(ambient.sourceSurfaceRadiusM*
            ambient.sourceSurfaceRadiusM*std::abs(sourceRadial));
        density=massPerFlux*Norm(field.value)/speed;
      }
      out.magneticSector=sourceRadial>0?1:-1;
      out.region=AmbientRegion::PfssOpen;
      out.velocityMPerS=corotation+(out.magneticSector*alongSpeed)*Unit(out.magneticFieldT);
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
  if(composition.alphaToProtonNumberRatio>0)out.alphaTemperatureK=composition.alphaTemperatureK;
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
  out.eventIdentity=definition_.physicsFingerprint;
  return Detail::DifferentiateAmbient(point,std::move(out),definition_.ambient.gradientRelativeStep,
      [&](Vec3 p){return Evaluate(p,epochS);});
}

Core::Result<AmbientState> Detail::DifferentiateAmbient(Vec3 point,AmbientState out,
    double relativeStep,const std::function<Core::Result<AmbientPrimitive>(Vec3)>& evaluate,
    bool requireConvergedNeighbors) {
  using Return=Core::Result<AmbientState>;
  // For the resolved wind, only a declared geometric boundary permits a
  // one-sided stencil. An unresolved sonic solution inside the supported
  // population is a failed physics/numerical query, not a boundary.
  const auto unresolved=[&](const Core::Result<AmbientPrimitive>& value) {
    return requireConvergedNeighbors&&!value.ok()&&value.status.code!=Core::StatusCode::OutOfDomain;
  };
  const double radius=Norm(point);
  for(int axis=0;axis<3;++axis) {
    double step=std::max(1.0,relativeStep*radius);
    bool resolved=false;
    for(int attempt=0;attempt<=12&&!resolved;++attempt,step*=0.5) {
      Vec3 delta;
      if(axis==0)delta.x=step;else if(axis==1)delta.y=step;else delta.z=step;
      const auto plus=evaluate(point+delta);
      const auto minus=evaluate(point-delta);
      for(const auto* neighbor:{&plus,&minus})if(unresolved(*neighbor))
        return Return::Failure(neighbor->status.code,neighbor->status.message);
      const bool goodPlus=plus.ok()&&SameBranch(out.primitive,plus.value);
      const bool goodMinus=minus.ok()&&SameBranch(out.primitive,minus.value);
      Vec3 magneticDerivative,velocityDerivative;
      double densityDerivative=0,pressureDerivative=0;
      if(goodPlus&&goodMinus) {
        magneticDerivative=(plus.value.magneticFieldT-minus.value.magneticFieldT)/(2*step);
        velocityDerivative=(plus.value.velocityMPerS-minus.value.velocityMPerS)/(2*step);
        densityDerivative=(plus.value.plasma.massDensityKgM3-
            minus.value.plasma.massDensityKgM3)/(2*step);
        pressureDerivative=(plus.value.plasma.pressurePa-
            minus.value.plasma.pressurePa)/(2*step);
        out.stencils[axis]=DifferenceStencil::CentralSecondOrder;
        resolved=true;
      } else for(int sign:{1,-1}) {
        if((sign==1&&!goodPlus)||(sign==-1&&!goodMinus))continue;
        const auto twice=evaluate(point+2*sign*delta);
        if(unresolved(twice))return Return::Failure(twice.status.code,twice.status.message);
        if(!twice.ok()||!SameBranch(out.primitive,twice.value))continue;
        const auto& first=sign==1?plus.value:minus.value;
        magneticDerivative=(-3*out.primitive.magneticFieldT+
            4*first.magneticFieldT-twice.value.magneticFieldT)/(2*sign*step);
        velocityDerivative=(-3*out.primitive.velocityMPerS+
            4*first.velocityMPerS-twice.value.velocityMPerS)/(2*sign*step);
        densityDerivative=(-3*out.primitive.plasma.massDensityKgM3+
            4*first.plasma.massDensityKgM3-twice.value.plasma.massDensityKgM3)/(2*sign*step);
        pressureDerivative=(-3*out.primitive.plasma.pressurePa+
            4*first.plasma.pressurePa-twice.value.plasma.pressurePa)/(2*sign*step);
        out.stencils[axis]=sign==1?DifferenceStencil::ForwardSecondOrder:
            DifferenceStencil::BackwardSecondOrder;
        resolved=true;
        break;
      }
      if(resolved) {
        if(!std::isfinite(densityDerivative)||!std::isfinite(pressureDerivative))
          return Return::Failure(Core::StatusCode::NumericalFailure,
              "ambient scalar derivative is nonfinite");
        // Populate the same Cartesian coordinate column as grad(B), grad(U).
        // These scalar derivatives retain kg/m^4 and Pa/m independently;
        // a divergence or inferred density-history scalar cannot replace them.
        if(axis==0) {
          out.gradientDensityKgM4.x=densityDerivative;
          out.gradientPressurePaPerM.x=pressureDerivative;
        } else if(axis==1) {
          out.gradientDensityKgM4.y=densityDerivative;
          out.gradientPressurePaPerM.y=pressureDerivative;
        } else {
          out.gradientDensityKgM4.z=densityDerivative;
          out.gradientPressurePaPerM.z=pressureDerivative;
        }
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
