#include "sep_coronal_cme/flux_tube_wind.h"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronalCME {
// These routines operate on already traced geometry. They neither infer a
// magnetic tube nor introduce an event normalization: geometry and one mass-
// flux authority must be supplied by their owning providers.
double FluxTubeArea(double flux,double field) {
  return flux>0.0&&field>0.0 ? flux/field : std::numeric_limits<double>::quiet_NaN();
}
double EffectivePotential(double r,double cylindrical,double omega) {
  return -Constants::kSolarGravitationalParameterM3PerS2/r-
      0.5*omega*omega*cylindrical*cylindrical;
}

namespace {
Core::Result<double> ParkerY(double radius,double critical) {
  const double rhs=4.0*std::log(radius/critical)+4.0*critical/radius-3.0;
  if(std::abs(radius-critical)<=1e-13*critical)
    return Core::Result<double>::Success(1.0);
  auto f=[&](double y){return y-std::log(y)-rhs;};
  double low,high;
  if(radius<critical) { low=1e-15; high=1.0-1e-13; }
  else { low=1.0+1e-13; high=std::max(4.0,rhs+std::log(std::max(2.0,rhs))+2.0); while(f(high)<0.0)high*=2.0; }
  if(f(low)*f(high)>0.0) return Core::Result<double>::Failure(
      Core::StatusCode::NumericalFailure,"Parker branch could not be bracketed");
  for(int i=0;i<200;++i) { const double mid=0.5*(low+high); if(f(low)*f(mid)<=0)high=mid;else low=mid; }
  return Core::Result<double>::Success(0.5*(low+high));
}
double Relative(double value,double reference) {
  return std::abs(value-reference)/std::max({1.0,std::abs(value),std::abs(reference)});
}
}

Core::Result<std::vector<WindPoint>> SolveRadialIsothermalParker(
    const std::vector<double>& radii,double sound,double critical,double massFlux,
    double areaAtOne) {
  if(radii.empty()||!(sound>0.0&&critical>0.0&&massFlux>0.0&&areaAtOne>0.0))
    return Core::Result<std::vector<WindPoint>>::Failure(Core::StatusCode::InvalidConfiguration,
        "isothermal Parker solve requires positive inputs");
  std::vector<WindPoint> result; result.reserve(radii.size());
  for(double r:radii) {
    if(!(r>0.0)) return Core::Result<std::vector<WindPoint>>::Failure(
        Core::StatusCode::InvalidConfiguration,"Parker radius must be positive");
    const auto y=ParkerY(r,critical); if(!y.ok()) return Core::Result<std::vector<WindPoint>>::Failure(y.status.code,y.status.message);
    WindPoint p; p.radiusM=r; p.areaM2=areaAtOne*r*r;
    p.speedMPerS=sound*std::sqrt(y.value); p.densityKgM3=massFlux/(p.speedMPerS*p.areaM2);
    p.pressurePa=sound*sound*p.densityKgM3;p.soundSpeedSquaredM2S2=sound*sound;
    p.effectivePotentialM2S2=-2.0*sound*sound*critical/r;
    result.push_back(p);
  }
  return Core::Result<std::vector<WindPoint>>::Success(std::move(result));
}

Core::Status CheckWindInvariants(const std::vector<WindPoint>& points,double gamma,
    double tolerance,double* massResidual,double* bernoulliResidual,double* entropyResidual) {
  if(points.size()<2||!(gamma>=1.0&&tolerance>=0.0)) return Core::Status::Failure(
      Core::StatusCode::InvalidConfiguration,"invariant check requires points and gamma>=1");
  const auto& p0=points.front(); const double mass0=p0.densityKgM3*p0.speedMPerS*p0.areaM2;
  const bool isothermal=std::abs(gamma-1.0)<1e-14;
  const double entropy0=isothermal?p0.pressurePa/p0.densityKgM3:
      p0.pressurePa/std::pow(p0.densityKgM3,gamma);
  const double bernoulli0=isothermal?
      0.5*p0.speedMPerS*p0.speedMPerS+p0.soundSpeedSquaredM2S2*std::log(p0.densityKgM3)+p0.effectivePotentialM2S2:
      0.5*p0.speedMPerS*p0.speedMPerS+p0.soundSpeedSquaredM2S2/(gamma-1.0)+p0.effectivePotentialM2S2;
  double mr=0,br=0,er=0;
  for(const auto& p:points) {
    const double mass=p.densityKgM3*p.speedMPerS*p.areaM2;
    const double entropy=isothermal?p.pressurePa/p.densityKgM3:
        p.pressurePa/std::pow(p.densityKgM3,gamma);
    const double bernoulli=isothermal?
        0.5*p.speedMPerS*p.speedMPerS+p.soundSpeedSquaredM2S2*std::log(p.densityKgM3)+p.effectivePotentialM2S2:
        0.5*p.speedMPerS*p.speedMPerS+p.soundSpeedSquaredM2S2/(gamma-1.0)+p.effectivePotentialM2S2;
    mr=std::max(mr,Relative(mass,mass0));br=std::max(br,Relative(bernoulli,bernoulli0));er=std::max(er,Relative(entropy,entropy0));
  }
  if(massResidual) *massResidual=mr;
  if(bernoulliResidual) *bernoulliResidual=br;
  if(entropyResidual) *entropyResidual=er;
  if(mr>tolerance||br>tolerance||er>tolerance) return Core::Status::Failure(
      Core::StatusCode::NumericalFailure,"wind invariant exceeds tolerance");
  return Core::Status::Success();
}

Core::Result<std::vector<CriticalCandidate>> FindCriticalCandidates(
    const std::vector<double>& s,const std::vector<double>& area,
    const std::vector<double>& potential) {
  if(s.size()<2||s.size()!=area.size()||s.size()!=potential.size())
    return Core::Result<std::vector<CriticalCandidate>>::Failure(
        Core::StatusCode::InvalidConfiguration,"critical-point arrays are inconsistent");
  std::vector<CriticalCandidate> result;
  // At a regular nozzle critical point a_c^2=Phi'/[ln(A)]'. Positive ratios
  // are retained as candidates; the global branch selector needs to see all
  // of them and therefore no nearest-root shortcut is used here.
  for(std::size_t i=0;i+1<s.size();++i) {
    if(!(s[i+1]>s[i])) return Core::Result<std::vector<CriticalCandidate>>::Failure(
        Core::StatusCode::InvalidConfiguration,"tube coordinate must increase");
    const double q0=area[i]!=0.0?potential[i]/area[i]:-1.0;
    const double q1=area[i+1]!=0.0?potential[i+1]/area[i+1]:-1.0;
    if(q0>0.0)result.push_back({s[i],q0,false});
    if(i+1==s.size()-1&&q1>0.0)result.push_back({s[i+1],q1,false});
  }
  if(result.empty()) return Core::Result<std::vector<CriticalCandidate>>::Failure(
      Core::StatusCode::NumericalFailure,"no physically positive critical candidate");
  // The outermost positive candidate is the global accelerating-branch choice
  // for this sampled topology; every candidate remains reported.
  result.back().globallyAdmissible=true;
  return Core::Result<std::vector<CriticalCandidate>>::Success(std::move(result));
}

Core::Result<PolytropicWindSolution> SolvePolytropicTube(
    const std::vector<TubeGeometryPoint>& geometry,double gamma,double entropy,
    double rhoCritical,std::size_t selected,double criticalTolerance) {
  if(geometry.size()<3||selected==0||selected+1>=geometry.size()||
      !(gamma>1.0&&gamma<1.5&&entropy>0.0&&rhoCritical>0.0&&criticalTolerance>0.0))
    return Core::Result<PolytropicWindSolution>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "polytropic tube requires 1<gamma_w<3/2, positive K/rho_c, and an interior critical index");
  std::vector<double> coordinate,areaDerivative,potentialDerivative;
  coordinate.reserve(geometry.size());areaDerivative.reserve(geometry.size());
  potentialDerivative.reserve(geometry.size());
  for(std::size_t i=0;i<geometry.size();++i) {
    const auto& g=geometry[i];
    if(!(g.areaM2>0.0&&g.radiusM>0.0&&std::isfinite(g.effectivePotentialM2S2)&&
         std::isfinite(g.dLogAreaDsPerM)&&std::isfinite(g.dPotentialDsMPerS2)&&
         (i==0||g.coordinateM>geometry[i-1].coordinateM)))
      return Core::Result<PolytropicWindSolution>::Failure(
          Core::StatusCode::InvalidConfiguration,"polytropic tube geometry is invalid");
    coordinate.push_back(g.coordinateM);areaDerivative.push_back(g.dLogAreaDsPerM);
    potentialDerivative.push_back(g.dPotentialDsMPerS2);
  }
  auto candidates=FindCriticalCandidates(coordinate,areaDerivative,potentialDerivative);
  if(!candidates.ok())return Core::Result<PolytropicWindSolution>::Failure(
      candidates.status.code,candidates.status.message);
  const auto& critical=geometry[selected];
  if(!(critical.dLogAreaDsPerM!=0.0))return Core::Result<PolytropicWindSolution>::Failure(
      Core::StatusCode::NumericalFailure,"selected critical point has zero area derivative");
  const double geometricSound2=critical.dPotentialDsMPerS2/critical.dLogAreaDsPerM;
  const double thermodynamicSound2=gamma*entropy*std::pow(rhoCritical,gamma-1.0);
  if(!(geometricSound2>0.0)||Relative(geometricSound2,thermodynamicSound2)>criticalTolerance)
    return Core::Result<PolytropicWindSolution>::Failure(
        Core::StatusCode::NumericalFailure,"selected critical point violates nozzle regularity");
  const double criticalSpeed=std::sqrt(thermodynamicSound2);
  const double massFlux=rhoCritical*criticalSpeed*critical.areaM2;
  const double bernoulli=0.5*thermodynamicSound2+
      thermodynamicSound2/(gamma-1.0)+critical.effectivePotentialM2S2;

  PolytropicWindSolution solution;solution.criticalCandidates=std::move(candidates.value);
  solution.selectedCriticalIndex=selected;solution.massFluxKgPerS=massFlux;
  solution.bernoulliM2S2=bernoulli;solution.polytropicConstantSI=entropy;
  solution.points.reserve(geometry.size());
  for(std::size_t index=0;index<geometry.size();++index) {
    const auto& g=geometry[index];
    double density=rhoCritical;
    if(index!=selected) {
      // At fixed A and Phi, Bernoulli as a function of log(rho) tends to
      // infinity at both ends. A valid transonic state has two roots: choose
      // the subsonic root inward and the supersonic root outward.
      auto residual=[&](double logDensity) {
        const double rho=std::exp(logDensity);
        const double u=massFlux/(rho*g.areaM2);
        const double a2=gamma*entropy*std::pow(rho,gamma-1.0);
        return 0.5*u*u+a2/(gamma-1.0)+g.effectivePotentialM2S2-bernoulli;
      };
      std::vector<std::pair<double,double>> brackets;
      double left=std::log(rhoCritical)-40.0,fl=residual(left);
      const int intervals=1600;
      for(int n=1;n<=intervals;++n) {
        const double right=std::log(rhoCritical)-40.0+80.0*n/intervals;
        const double fr=residual(right);
        if(fl==0.0||fl*fr<0.0)brackets.emplace_back(left,right);
        left=right;fl=fr;
      }
      bool found=false;
      for(const auto& bracket:brackets) {
        double low=bracket.first,high=bracket.second,fLow=residual(low);
        for(int n=0;n<120;++n){const double mid=0.5*(low+high),fm=residual(mid);
          if(fLow*fm<=0.0)high=mid;else{low=mid;fLow=fm;}}
        const double candidate=std::exp(0.5*(low+high));
        const double speed=massFlux/(candidate*g.areaM2);
        const double sound=std::sqrt(gamma*entropy*std::pow(candidate,gamma-1.0));
        const bool desired=index<selected ? speed<sound : speed>sound;
        if(desired){density=candidate;found=true;break;}
      }
      if(!found)return Core::Result<PolytropicWindSolution>::Failure(
          Core::StatusCode::NumericalFailure,
          "global polytropic transonic branch does not exist at every tube point");
    }
    const double speed=massFlux/(density*g.areaM2);
    const double pressure=entropy*std::pow(density,gamma);
    const double sound2=gamma*pressure/density;
    solution.points.push_back({g.radiusM,g.areaM2,speed,density,pressure,sound2,
                               g.effectivePotentialM2S2});
  }
  return Core::Result<PolytropicWindSolution>::Success(std::move(solution));
}

Core::Result<double> InvertTargetSpeed(const std::function<Core::Result<double>(double)>& model,
    double target,double low,double high,double tolerance) {
  if(!(target>0&&low>0&&high>low&&tolerance>0)) return Core::Result<double>::Failure(
      Core::StatusCode::InvalidConfiguration,"invalid target-speed inversion bracket");
  auto a=model(low),b=model(high);if(!a.ok()||!b.ok())return Core::Result<double>::Failure(
      Core::StatusCode::NumericalFailure,"target model failed at bracket endpoint");
  double fa=a.value-target,fb=b.value-target;
  if(fa*fb>0)return Core::Result<double>::Failure(Core::StatusCode::NumericalFailure,
      "target speed is not bracketed");
  for(int n=0;n<200;++n) { double mid=0.5*(low+high);auto m=model(mid);if(!m.ok())return Core::Result<double>::Failure(m.status.code,m.status.message);
    const double fm=m.value-target;if(std::abs(fm)<=tolerance*target)return Core::Result<double>::Success(mid);
    if(fa*fm<=0){high=mid;fb=fm;}else{low=mid;fa=fm;}}
  return Core::Result<double>::Failure(Core::StatusCode::NumericalFailure,"target inversion did not converge");
}

Core::Result<double> ResolveMassLoading(double rho,double speed,double field) {
  if(!(rho>0&&speed>0&&field>0))return Core::Result<double>::Failure(
      Core::StatusCode::InvalidConfiguration,"mass loading requires positive co-located rho,u,|B|");
  return Core::Result<double>::Success(rho*speed/field);
}
Core::Result<double> ResolveMassLoadingFromRadialFlux(double flux,double br) {
  if(!(flux>0&&std::abs(br)>0))return Core::Result<double>::Failure(
      Core::StatusCode::InvalidConfiguration,"radial mass loading requires positive flux and nonzero mapped Br");
  return Core::Result<double>::Success(flux/std::abs(br));
}
} }
