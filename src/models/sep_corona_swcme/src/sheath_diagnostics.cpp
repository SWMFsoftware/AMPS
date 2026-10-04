#include "sep_corona_swcme/sheath_diagnostics.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace SEP { namespace CoronaSwcme {
namespace {

using CoronalCME::Cross;
using CoronalCME::Dot;
using CoronalCME::Norm;
using CoronalCME::Vec3;

constexpr double kPi=3.141592653589793238462643383279502884;
constexpr double kMu0=1.25663706212e-6;

double Determinant(const std::array<Vec3,3>& columns) {
  return Dot(columns[0],Cross(columns[1],columns[2]));
}

template<class Value>
Value Fourth(const Value& minus2,const Value& minus,const Value& plus,
    const Value& plus2,double step) {
  return (minus2-8*minus+8*plus-plus2)/(12*step);
}

std::array<Vec3,3> SpatialGradient(const std::array<Vec3,3>& coordinateDerivatives,
    const std::array<Vec3,3>& valueDerivatives) {
  const double determinant=Determinant(coordinateDerivatives);
  const auto derivative=[&](Vec3 direction) {
    const double c0=Dot(direction,Cross(coordinateDerivatives[1],
        coordinateDerivatives[2]))/determinant;
    const double c1=Dot(direction,Cross(coordinateDerivatives[2],
        coordinateDerivatives[0]))/determinant;
    const double c2=Dot(direction,Cross(coordinateDerivatives[0],
        coordinateDerivatives[1]))/determinant;
    return c0*valueDerivatives[0]+c1*valueDerivatives[1]+c2*valueDerivatives[2];
  };
  return {derivative({1,0,0}),derivative({0,1,0}),derivative({0,0,1})};
}

Vec3 SpatialScalarGradient(const std::array<Vec3,3>& coordinateDerivatives,
    const std::array<double,3>& valueDerivatives) {
  const double determinant=Determinant(coordinateDerivatives);
  const double c0x=Dot(Vec3{1,0,0},Cross(coordinateDerivatives[1],
      coordinateDerivatives[2]))/determinant;
  const double c1x=Dot(Vec3{1,0,0},Cross(coordinateDerivatives[2],
      coordinateDerivatives[0]))/determinant;
  const double c2x=Dot(Vec3{1,0,0},Cross(coordinateDerivatives[0],
      coordinateDerivatives[1]))/determinant;
  const double c0y=Dot(Vec3{0,1,0},Cross(coordinateDerivatives[1],
      coordinateDerivatives[2]))/determinant;
  const double c1y=Dot(Vec3{0,1,0},Cross(coordinateDerivatives[2],
      coordinateDerivatives[0]))/determinant;
  const double c2y=Dot(Vec3{0,1,0},Cross(coordinateDerivatives[0],
      coordinateDerivatives[1]))/determinant;
  const double c0z=Dot(Vec3{0,0,1},Cross(coordinateDerivatives[1],
      coordinateDerivatives[2]))/determinant;
  const double c1z=Dot(Vec3{0,0,1},Cross(coordinateDerivatives[2],
      coordinateDerivatives[0]))/determinant;
  const double c2z=Dot(Vec3{0,0,1},Cross(coordinateDerivatives[0],
      coordinateDerivatives[1]))/determinant;
  return {c0x*valueDerivatives[0]+c1x*valueDerivatives[1]+c2x*valueDerivatives[2],
      c0y*valueDerivatives[0]+c1y*valueDerivatives[1]+c2y*valueDerivatives[2],
      c0z*valueDerivatives[0]+c1z*valueDerivatives[1]+c2z*valueDerivatives[2]};
}

double EnergyDensity(const SheathMappedState& state,double gamma) {
  const auto& primitive=state.primitive;
  return 0.5*primitive.massDensityKgM3*Dot(primitive.velocityMPerS,
      primitive.velocityMPerS)+primitive.pressurePa/(gamma-1)+
      Dot(primitive.magneticFieldT,primitive.magneticFieldT)/(2*kMu0);
}

Vec3 MaterialBoundaryEnergyFlux(const SheathMappedState& state) {
  const auto& primitive=state.primitive;
  const double totalPressure=primitive.pressurePa+
      Dot(primitive.magneticFieldT,primitive.magneticFieldT)/(2*kMu0);
  return totalPressure*primitive.velocityMPerS-
      (Dot(primitive.velocityMPerS,primitive.magneticFieldT)/kMu0)*
          primitive.magneticFieldT;
}

} // namespace

Core::Result<SheathResidualSample> EvaluateSheathResidual(
    const ShockFedSheathModel& model,const SheathMaterialLabel& label,
    double epochS,double gamma,const SheathDiagnosticOptions& options) {
  using Return=Core::Result<SheathResidualSample>;
  if(!(std::isfinite(epochS)&&std::isfinite(gamma)&&gamma>1&&
      options.angularStepRad>0&&options.admissionTimeStepS>0&&
      options.evolutionTimeStepS>0&&options.gravitationalParameterM3S2>=0&&
      label.polarRad-2*options.angularStepRad>0&&
      label.polarRad+2*options.angularStepRad<kPi&&
      label.azimuthRad-2*options.angularStepRad>=0&&
      label.azimuthRad+2*options.angularStepRad<2*kPi&&
      label.crossingTimeS-2*options.admissionTimeStepS>=0&&
      epochS-2*options.evolutionTimeStepS>=label.crossingTimeS))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "residual stencil is outside its smooth material-label support");

  const auto evaluate=[&](const SheathMaterialLabel& q,double time) {
    return model.Evaluate(q,time);
  };
  const auto central=evaluate(label,epochS);
  if(!central.ok())return Return::Failure(central.status.code,central.status.message);

  std::array<std::array<SheathMappedState,4>,3> neighbor;
  for(int dimension=0;dimension<3;++dimension)for(int k=0;k<4;++k) {
    const double multiplier=k<2?static_cast<double>(k-2):static_cast<double>(k-1);
    SheathMaterialLabel q=label;
    if(dimension==0)q.polarRad+=multiplier*options.angularStepRad;
    if(dimension==1)q.azimuthRad+=multiplier*options.angularStepRad;
    if(dimension==2)q.crossingTimeS+=multiplier*options.admissionTimeStepS;
    const auto value=evaluate(q,epochS);
    if(!value.ok())return Return::Failure(Core::StatusCode::UnsupportedCapability,
        "distributed residual stencil crosses an unsupported material branch: "+
        value.status.message);
    neighbor[dimension][k]=value.value;
  }
  const double steps[3]={options.angularStepRad,options.angularStepRad,
      options.admissionTimeStepS};
  std::array<Vec3,3> dx,du,db,dg;
  std::array<double,3> dp;
  for(int dimension=0;dimension<3;++dimension) {
    const auto& s=neighbor[dimension];
    dx[dimension]=Fourth(s[0].positionM,s[1].positionM,s[2].positionM,
        s[3].positionM,steps[dimension]);
    du[dimension]=Fourth(s[0].primitive.velocityMPerS,s[1].primitive.velocityMPerS,
        s[2].primitive.velocityMPerS,s[3].primitive.velocityMPerS,steps[dimension]);
    db[dimension]=Fourth(s[0].primitive.magneticFieldT,s[1].primitive.magneticFieldT,
        s[2].primitive.magneticFieldT,s[3].primitive.magneticFieldT,steps[dimension]);
    dp[dimension]=Fourth(s[0].primitive.pressurePa,s[1].primitive.pressurePa,
        s[2].primitive.pressurePa,s[3].primitive.pressurePa,steps[dimension]);
    dg[dimension]=Fourth(MaterialBoundaryEnergyFlux(s[0]),
        MaterialBoundaryEnergyFlux(s[1]),MaterialBoundaryEnergyFlux(s[2]),
        MaterialBoundaryEnergyFlux(s[3]),steps[dimension]);
  }
  const double coordinateDeterminant=Determinant(dx);
  const double coordinateScale=std::max(1.0,Norm(dx[0])*Norm(dx[1])*Norm(dx[2]));
  if(!(std::isfinite(coordinateDeterminant)&&
      std::abs(coordinateDeterminant)>1e-14*coordinateScale))
    return Return::Failure(Core::StatusCode::NumericalFailure,
        "distributed residual stencil has a singular Cartesian map");
  const auto gradientU=SpatialGradient(dx,du);
  const auto gradientB=SpatialGradient(dx,db);
  const auto gradientG=SpatialGradient(dx,dg);
  const Vec3 gradientP=SpatialScalarGradient(dx,dp);
  const double divergenceU=gradientU[0].x+gradientU[1].y+gradientU[2].z;
  const double divergenceG=gradientG[0].x+gradientG[1].y+gradientG[2].z;
  const Vec3 curlB{gradientB[1].z-gradientB[2].y,
      gradientB[2].x-gradientB[0].z,gradientB[0].y-gradientB[1].x};

  std::array<SheathMappedState,4> temporal;
  for(int k=0;k<4;++k) {
    const double multiplier=k<2?static_cast<double>(k-2):static_cast<double>(k-1);
    const auto value=evaluate(label,epochS+multiplier*options.evolutionTimeStepS);
    if(!value.ok())return Return::Failure(Core::StatusCode::UnsupportedCapability,
        "distributed residual time stencil leaves material support: "+
        value.status.message);
    temporal[k]=value.value;
  }
  const Vec3 acceleration=Fourth(temporal[0].primitive.velocityMPerS,
      temporal[1].primitive.velocityMPerS,temporal[2].primitive.velocityMPerS,
      temporal[3].primitive.velocityMPerS,options.evolutionTimeStepS);
  const double materialEnergyRate=Fourth(EnergyDensity(temporal[0],gamma),
      EnergyDensity(temporal[1],gamma),EnergyDensity(temporal[2],gamma),
      EnergyDensity(temporal[3],gamma),options.evolutionTimeStepS);

  const auto& primitive=central.value.primitive;
  const double radius=Norm(central.value.positionM);
  if(!(std::isfinite(radius)&&radius>0))return Return::Failure(
      Core::StatusCode::InvalidState,"distributed residual has zero radius");
  const Vec3 gravitationalAcceleration=(-options.gravitationalParameterM3S2/
      (radius*radius*radius))*central.value.positionM;
  SheathResidualSample result;
  result.label=label;result.epochS=epochS;result.positionM=central.value.positionM;
  result.inertiaPaPerM=primitive.massDensityKgM3*acceleration;
  result.pressureGradientPaPerM=gradientP;
  result.lorentzPaPerM=(1/kMu0)*Cross(curlB,primitive.magneticFieldT);
  result.gravityPaPerM=primitive.massDensityKgM3*gravitationalAcceleration;
  result.momentumResidualPaPerM=result.inertiaPaPerM+
      result.pressureGradientPaPerM-result.lorentzPaPerM-result.gravityPaPerM;
  result.momentumScalePaPerM=Norm(result.inertiaPaPerM)+
      Norm(result.pressureGradientPaPerM)+Norm(result.lorentzPaPerM)+
      Norm(result.gravityPaPerM);
  result.localMomentumRatio=Norm(result.momentumResidualPaPerM)/
      std::max(result.momentumScalePaPerM,std::numeric_limits<double>::min());
  const double energy=EnergyDensity(central.value,gamma);
  result.gravityWorkWPerM3=Dot(result.gravityPaPerM,primitive.velocityMPerS);
  result.energyResidualWPerM3=materialEnergyRate+energy*divergenceU+
      divergenceG-result.gravityWorkWPerM3;
  result.energyScaleWPerM3=std::abs(materialEnergyRate)+
      std::abs(energy*divergenceU)+std::abs(divergenceG)+
      std::abs(result.gravityWorkWPerM3);
  result.inertiaWorkWPerM3=Dot(result.inertiaPaPerM,primitive.velocityMPerS);
  result.pressureGradientWorkWPerM3=Dot(result.pressureGradientPaPerM,
      primitive.velocityMPerS);
  result.lorentzWorkWPerM3=Dot(result.lorentzPaPerM,primitive.velocityMPerS);
  result.residualForceWorkWPerM3=Dot(result.momentumResidualPaPerM,
      primitive.velocityMPerS);
  return Return::Success(std::move(result));
}

Core::Result<SheathResidualReport> IntegrateSheathResiduals(
    const ShockFedSheathModel& model,
    const std::vector<WeightedSheathDiagnosticPoint>& points,double epochS,
    double gamma,const SheathDiagnosticOptions& options) {
  using Return=Core::Result<SheathResidualReport>;
  if(points.empty())return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "distributed residual integration requires material points");
  SheathResidualReport report;
  std::vector<std::pair<double,double>> ratios;
  for(const auto& point:points) {
    if(!(std::isfinite(point.volumeM3)&&point.volumeM3>0))return Return::Failure(
        Core::StatusCode::InvalidConfiguration,
        "distributed residual point has nonpositive physical volume");
    const auto sample=EvaluateSheathResidual(model,point.label,epochS,gamma,options);
    if(!sample.ok()) {
      if(sample.status.code==Core::StatusCode::UnsupportedCapability) {
        ++report.unsupportedSamples;continue;
      }
      return Return::Failure(sample.status.code,sample.status.message);
    }
    const double volume=point.volumeM3;
    report.sampledVolumeM3+=volume;
    report.integratedAbsoluteMomentumResidualN+=
        Norm(sample.value.momentumResidualPaPerM)*volume;
    report.integratedMomentumScaleN+=sample.value.momentumScalePaPerM*volume;
    report.signedEnergyResidualW+=sample.value.energyResidualWPerM3*volume;
    report.absoluteEnergyResidualW+=std::abs(sample.value.energyResidualWPerM3)*volume;
    report.integratedEnergyScaleW+=sample.value.energyScaleWPerM3*volume;
    report.signedResidualForceWorkW+=sample.value.residualForceWorkWPerM3*volume;
    report.absoluteResidualForceWorkW+=
        std::abs(sample.value.residualForceWorkWPerM3)*volume;
    report.signedInertiaWorkW+=sample.value.inertiaWorkWPerM3*volume;
    report.signedPressureGradientWorkW+=
        sample.value.pressureGradientWorkWPerM3*volume;
    report.signedLorentzWorkW+=sample.value.lorentzWorkWPerM3*volume;
    report.signedGravityWorkW+=sample.value.gravityWorkWPerM3*volume;
    ratios.push_back({sample.value.localMomentumRatio,volume});
    report.samples.push_back(sample.value);
  }
  if(report.samples.empty())return Return::Failure(Core::StatusCode::UnsupportedCapability,
      "no distributed residual sample has a complete smooth stencil");
  report.integratedMomentumRatio=report.integratedAbsoluteMomentumResidualN/
      std::max(report.integratedMomentumScaleN,std::numeric_limits<double>::min());
  report.forceWorkRatio=report.absoluteResidualForceWorkW/
      std::max(report.integratedEnergyScaleW,std::numeric_limits<double>::min());
  std::sort(ratios.begin(),ratios.end(),[](const auto& left,const auto& right) {
    return left.first<right.first;
  });
  const double target=0.99*report.sampledVolumeM3;
  double cumulative=0;
  for(const auto& item:ratios) {
    cumulative+=item.second;
    report.volumeWeightedLocalMomentumP99=item.first;
    if(cumulative>=target)break;
  }
  return Return::Success(std::move(report));
}

} } // namespace SEP::CoronaSwcme
