#include "sep_corona_swcme/open_tube_isopoly.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP::CoronaSwcme {
namespace {

using Core::Result;
using Core::StatusCode;
using CoronalCME::Vec3;
namespace C=CoronalCME;
constexpr double gm=C::Constants::kSolarGravitationalParameterM3PerS2;

bool Positive(double x) { return std::isfinite(x)&&x>0; }

struct MeridionalNode { double radius, theta, thetaPrime; };

double Colatitude(const MeridionalNode& lo,const MeridionalNode& hi,double r) {
  const double h=hi.radius-lo.radius;
  const double z=(r-lo.radius)/h;
  const double d=(hi.theta-lo.theta)/h;
  return lo.theta+h*z*(lo.thetaPrime+z*(3*d-2*lo.thetaPrime-hi.thetaPrime+
      z*(lo.thetaPrime+hi.thetaPrime-2*d)));
}

}  // namespace

struct OpenTubeIsoPolyWind::Geometry {
  C::PfssHarmonics pfss;
  AmbientDefinition ambient;
  std::vector<MeridionalNode> nodes;
  double sourceTheta=0, sourceBr=0, angularMomentum=0;

  Result<double> Direction(double r,double theta) const {
    const auto field=pfss.Evaluate(r,theta,0);
    if(!field.ok())return Result<double>::Failure(field.status.code,field.status.message);
    const auto& b=field.value;
    // In an axisymmetric open tube the meridional tangent obeys dtheta/dr
    // = B_theta/(r B_r), irrespective of polarity. The radial chart stops
    // at a null, polarity change or unresolved radial turning point.
    if(std::abs(b.brT)<=ambient.ambient.minimumMagneticFieldT||
        std::signbit(b.brT)!=std::signbit(sourceBr))return Result<double>::Failure(
        StatusCode::OutOfDomain,"open-tube wind trace leaves its monotone radial magnetic sector");
    if(theta==0||theta==C::Constants::kPi)return Result<double>::Success(0);
    return Result<double>::Success(b.bThetaT/(r*b.brT));
  }

  Result<double> TraceStep(double r,double theta,double end) const {
    const double h=end-r, mid=r+h/2;
    const auto a=Direction(r,theta);
    if(!a.ok())return a;
    const auto b=Direction(mid,theta+h*a.value/2);
    if(!b.ok())return b;
    const auto c=Direction(mid,theta+h*b.value/2);
    if(!c.ok())return c;
    const auto d=Direction(end,theta+h*c.value);
    if(!d.ok())return d;
    return Result<double>::Success(theta+h*(a.value+2*b.value+2*c.value+d.value)/6);
  }

  Core::Status Trace(double tolerance) {
    const auto& a=ambient.ambient;
    const double target=ambient.support.firstValidPlasmaRadiusM;
    double r=a.sourceSurfaceRadiusM, theta=sourceTheta, step=a.traceStepM;
    const auto first=Direction(r,theta);
    if(!first.ok())return first.status;
    nodes.push_back({r,theta,first.value});
    for(int i=0;i<a.traceMaximumSteps;++i) {
      const double end=r-target<=step?target:r-step;
      const double mid=r+(end-r)/2;
      if(mid==r||mid==end)return Core::Status::Failure(StatusCode::NumericalFailure,
          "open-tube wind could not refine its angular trace in the radial chart");
      const auto coarse=TraceStep(r,theta,end);
      const auto half=TraceStep(r,theta,mid);
      const auto fine=half.ok()?TraceStep(mid,half.value,end):half;
      const auto slope=fine.ok()?Direction(end,fine.value):fine;
      if(!coarse.ok()||!fine.ok()||!slope.ok()) { step*=0.5; continue; }
      const MeridionalNode next{end,fine.value,slope.value};
      // Resolve both the endpoint and dense trajectory. The geometry solver
      // will query the tube between trace nodes while locating sonic points.
      const double error=std::max(std::abs(coarse.value-fine.value)/15,
          std::abs(half.value-Colatitude(next,nodes.back(),mid)));
      if(error>tolerance) { step*=0.5; continue; }
      nodes.push_back(next);
      if(end==target) {
        std::reverse(nodes.begin(),nodes.end());
        return Core::Status::Success();
      }
      r=end;
      theta=fine.value;
      step=std::min(a.traceStepM,1.75*step);
    }
    return Core::Status::Failure(StatusCode::NumericalFailure,
        "open-tube wind angular trace did not cover the declared inner guard within the work limit");
  }

  double Theta(double r) const {
    if(r>=pfss.SourceSurfaceRadiusM())return sourceTheta;
    const auto hi=std::upper_bound(nodes.begin(),nodes.end(),r,
        [](double key,const MeridionalNode& node){return key<node.radius;});
    return Colatitude(*(hi-1),*hi,r);
  }

  Result<C::TubeGeometryPoint> At(double r) const {
    using Return=Result<C::TubeGeometryPoint>;
    const double source=pfss.SourceSurfaceRadiusM();
    const double omega=ambient.ambient.rotationRateRadPerS;
    if(r>=source) {
      // Outside the source surface use radial speed and AREA NORMAL TO r.
      // The axisymmetric azimuthal flow has r U_phi constant. Substitution of
      // U-relative along the spiral into the rotating projected momentum
      // equation cancels the apparent Omega^2 r acceleration and gives
      // dPhi/dr = GM/r^2 - U_phi^2/r. Choose the additive constant to match
      // the interior potential at the source surface. Here 1 Wb merely
      // normalizes tube area and mass flux by the same constant.
      const double uPhi=angularMomentum/r;
      const double area=std::pow(r/source,2)/std::abs(sourceBr);
      const double offset=omega*source*std::sin(sourceTheta);
      return Return::Success({r,r,area,-gm/r+0.5*uPhi*uPhi-offset*offset,
          2/r,gm/(r*r)-uPhi*uPhi/r});
    }
    const double theta=Theta(r);
    const auto field=pfss.Evaluate(r,theta,0);
    if(!field.ok())return Return::Failure(field.status.code,field.status.message);
    const auto& b=field.value;
    if(std::abs(b.brT)<=ambient.ambient.minimumMagneticFieldT||
        std::signbit(b.brT)!=std::signbit(sourceBr))return Return::Failure(
        StatusCode::OutOfDomain,"open-tube wind chart loses its monotone radial magnetic sector");
    // Interior speed follows the field tangent while r parametrizes the
    // monotone curve. All terms in the projected equation are multiplied
    // by ds/dr consistently. Area is perpendicular to B on this side.
    const double thetaPrime=b.bThetaT/(r*b.brT);
    const double b2=b.brT*b.brT+b.bThetaT*b.bThetaT;
    const double dLogArea=-(b.brT*(b.derivative[0][0]+thetaPrime*b.derivative[0][1])+
        b.bThetaT*(b.derivative[1][0]+thetaPrime*b.derivative[1][1]))/b2;
    const double sin=std::sin(theta), cos=std::cos(theta);
    return Return::Success({r,r,1/std::sqrt(b2),-gm/r-0.5*std::pow(omega*r*sin,2),
        dLogArea,gm/(r*r)-omega*omega*r*sin*(sin+r*cos*thetaPrime)});
  }
};

Result<OpenTubeIsoPolyWind> OpenTubeIsoPolyWind::Prepare(std::shared_ptr<const AmbientModel> ambient,
    const OpenTubeIsoPolyRequest& request) {
  using Return=Result<OpenTubeIsoPolyWind>;
  if(!ambient||!(request.sourceColatitudeRad>=0&&request.sourceColatitudeRad<=C::Constants::kPi)||
      !Positive(request.referenceDensityKgM3)||!Positive(request.massLoadingKgPerSWb)||
      !Positive(request.angularTraceToleranceRad)||request.angularTraceToleranceRad>=1)
    return Return::Failure(StatusCode::InvalidConfiguration,
        "open-tube iso-poly preparation requires a source angle, thermal/mass anchors and trace accuracy");
  for(const auto& coefficient:ambient->pfss_.Coefficients())
    if(coefficient.order!=0&&(coefficient.cosineT!=0||coefficient.sineT!=0))
      return Return::Failure(StatusCode::UnsupportedCapability,
          "tube-dependent nonaxisymmetric Parker winding requires a conservative longitude map");

  auto g=std::make_shared<Geometry>();
  g->pfss=ambient->pfss_;
  g->ambient=ambient->Definition();
  g->sourceTheta=request.sourceColatitudeRad;
  const double source=g->pfss.SourceSurfaceRadiusM();
  const auto field=g->pfss.Evaluate(source,g->sourceTheta,0);
  if(!field.ok())return Return::Failure(field.status.code,field.status.message);
  g->sourceBr=field.value.brT;
  if(std::abs(g->sourceBr)<=g->ambient.ambient.minimumMagneticFieldT)
    return Return::Failure(StatusCode::OutOfDomain,"open-tube source lies on an unresolved magnetic neutral line");
  g->angularMomentum=g->ambient.ambient.rotationRateRadPerS*source*source*
      std::sin(g->sourceTheta);
  const auto trace=g->Trace(request.angularTraceToleranceRad);
  if(!trace.ok())return Return::Failure(trace.code,trace.message);

  const auto& c=g->ambient.composition;
  const auto unit=C::EvaluatePlasmaFromElectronDensity(1,c.electronTemperatureK,
      ambient->ions_,c.includeElectronMass,c.gammaAdiabatic,0);
  if(!unit.ok())return Return::Failure(unit.status.code,unit.status.message);
  // Use the same charge-neutral EOS conversion as ambient/jump export.
  // A unit of electron number density fixes every species per kg exactly.
  const double particles=1/unit.value.massDensityKgM3;
  const double protons=particles*unit.value.protonNumberDensityM3;
  C::IsoPolyTubeInput input;
  input.geometry=[geometry=std::shared_ptr<const Geometry>(g)](double r){return geometry->At(r);};
  input.beginM=g->ambient.support.firstValidPlasmaRadiusM;
  input.endM=g->ambient.support.coverageRadiusM;
  input.massFluxKgPerS=request.massLoadingKgPerSWb; // normalized to 1 Wb
  input.criticalScanIntervals=request.criticalScanIntervals;
  input.relativeTolerance=request.relativeWindTolerance;
  const auto addSpecies=[&](double concentration,double temperature,const OpenTubeThermalBand& band) {
    input.species.push_back({concentration,temperature,request.referenceDensityKgM3,
        band.polytropicIndex,band.blendBeginRadiusM,band.blendEndRadiusM,c.gammaAdiabatic});
  };
  addSpecies(particles,c.electronTemperatureK,request.electron);
  addSpecies(protons,c.protonTemperatureK,request.proton);
  if(c.alphaToProtonNumberRatio>0)
    addSpecies(protons*c.alphaToProtonNumberRatio,c.alphaTemperatureK,request.alpha);
  auto wind=C::IsoPolyTubeWind::Solve(input);
  if(!wind.ok())return Return::Failure(wind.status.code,wind.status.message);
  OpenTubeIsoPolyWind answer;
  answer.ambient_=std::move(ambient);
  answer.geometry_=std::move(g);
  answer.wind_=std::move(wind.value);
  return Return::Success(std::move(answer));
}

Result<OpenTubeIsoPolySample> OpenTubeIsoPolyWind::Evaluate(double r,double phi,double epoch) const {
  using Return=Result<OpenTubeIsoPolySample>;
  if(!geometry_||!std::isfinite(phi)||!std::isfinite(epoch)||
      epoch<geometry_->ambient.support.startS||epoch>geometry_->ambient.support.endS)
    return Return::Failure(StatusCode::OutOfDomain,"open-tube wind query has no declared epoch or azimuth");
  auto balance=wind_.Evaluate(r);
  if(!balance.ok())return Return::Failure(balance.status.code,balance.status.message);
  const auto& g=*geometry_;
  const auto& c=g.ambient.composition;
  const double source=g.pfss.SourceSurfaceRadiusM();
  const double theta=g.Theta(r), sin=std::sin(theta), cos=std::cos(theta);
  const Vec3 eR={sin*std::cos(phi),sin*std::sin(phi),cos};
  const Vec3 eTheta={cos*std::cos(phi),cos*std::sin(phi),-sin};
  const Vec3 ePhi={-std::sin(phi),std::cos(phi),0};
  const double speed=balance.value.wind.speedMPerS;
  OpenTubeIsoPolySample out;
  out.positionM=r*eR;
  auto& p=out.primitive;
  p.magneticSector=g.sourceBr>0?1:-1;
  const double corotation=g.ambient.ambient.rotationRateRadPerS*r*sin;
  if(r<source) {
    const auto field=g.pfss.Evaluate(r,theta,0);
    if(!field.ok())return Return::Failure(field.status.code,field.status.message);
    p.magneticFieldT=field.value.brT*eR+field.value.bThetaT*eTheta;
    p.velocityMPerS=p.magneticSector*speed*C::Unit(p.magneticFieldT)+corotation*ePhi;
    p.region=AmbientRegion::PfssOpen;
    // The momentum equation is reparametrized by r while its speed is the
    // magnitude along B. The generic heat expression multiplied its radial
    // gradients by that magnitude. Material advection across a radial shell
    // instead uses U_r=u*dr/ds (corotation has zero radial component).
    const double radialProjection=p.magneticSector*field.value.brT/C::Norm(p.magneticFieldT);
    for(double& rate:balance.value.netHeatWPerM3)rate*=radialProjection;
  } else {
    const double br=g.sourceBr*std::pow(source/r,2);
    const double uPhi=g.angularMomentum/r;
    p.magneticFieldT=br*(eR+((uPhi-corotation)/speed)*ePhi);
    p.velocityMPerS=speed*eR+uPhi*ePhi;
    p.region=AmbientRegion::ParkerExterior;
  }
  p.electronTemperatureK=balance.value.temperatureK[0];
  p.protonTemperatureK=balance.value.temperatureK[1];
  p.electronPressurePa=balance.value.partialPressurePa[0];
  std::vector<C::IonSpecies> ions={{"H+",1,1,C::Constants::kProtonMassKg,p.protonTemperatureK}};
  if(c.alphaToProtonNumberRatio>0)ions.push_back({"He++",c.alphaToProtonNumberRatio,2,
      C::Constants::kAlphaMassKg,balance.value.temperatureK[2]});
  const double fieldMagnitude=C::Norm(p.magneticFieldT);
  if(!std::isfinite(fieldMagnitude)||fieldMagnitude<=g.ambient.ambient.minimumMagneticFieldT)
    return Return::Failure(StatusCode::OutOfDomain,
        "resolved open-tube wind does not cover this magnetic field magnitude");
  const auto unit=C::EvaluatePlasmaFromElectronDensity(1,p.electronTemperatureK,ions,
      c.includeElectronMass,c.gammaAdiabatic,0);
  if(!unit.ok())return Return::Failure(unit.status.code,unit.status.message);
  const auto plasma=C::EvaluatePlasmaFromElectronDensity(balance.value.wind.densityKgM3/
      unit.value.massDensityKgM3,p.electronTemperatureK,ions,c.includeElectronMass,
      c.gammaAdiabatic,fieldMagnitude);
  if(!plasma.ok())return Return::Failure(plasma.status.code,plasma.status.message);
  p.plasma=plasma.value;
  out.balance=std::move(balance.value);
  return Return::Success(std::move(out));
}

}  // namespace SEP::CoronaSwcme
