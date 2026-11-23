#include "sep_coronal_cme/tube_isopoly_wind.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace SEP::CoronalCME {
namespace {

using Core::Result;
using Core::StatusCode;
constexpr double roundoff=64*std::numeric_limits<double>::epsilon();

bool Positive(double a) { return std::isfinite(a)&&a>0; }

struct SpeciesState {
  double temperature, pressure, densityExponent, explicitLogGradient;
};

SpeciesState Thermal(const IsoPolySpecies& species,double s,double rho) {
  double weight=0, weightDerivative=0, logTemperatureRatio=0;
  if(species.polytropicIndex!=1&&s>species.blendBeginM) {
    const double width=species.blendEndM-species.blendBeginM;
    const double z=(s-species.blendBeginM)/width;
    if(z>=1)weight=1;
    else {
      // The quintic has vanishing first and second derivatives at both ends.
      // Blending log(T) preserves positive temperature. The explicit gradient
      // at fixed density also contributes to the numerator of the wind ODE.
      weight=z*z*z*(10+z*(-15+6*z));
      weightDerivative=30*z*z*(1-z)*(1-z)/width;
    }
    logTemperatureRatio=(species.polytropicIndex-1)*
        (std::log(rho)-std::log(species.referenceDensityKgM3));
  }
  const double temperature=species.referenceTemperatureK*
      std::exp(weight*logTemperatureRatio);
  return {temperature,species.particlesPerKg*rho*Constants::kBoltzmannJPerK*
      temperature,weight*(species.polytropicIndex-1),
      weightDerivative*logTemperatureRatio};
}

struct Equation {
  double density=0, pressure=0, soundSquared=0, explicitPressureGradient=0;
  double denominator=0, numerator=0;
};

struct Problem {
  const IsoPolyTubeInput& input;
  double scale, extent;

  // Logarithmic distance resolves the inner corona without requiring the
  // scan or integrator to use the same metre spacing at several AU. x is
  // dimensionless. Thermal/geometry callbacks always receive metres.
  double Coordinate(double x) const {
    if(x==0)return input.beginM;
    if(x==1)return input.endM;
    return input.beginM+scale*std::expm1(x*extent);
  }
  double X(double s) const { return std::log1p((s-input.beginM)/scale)/extent; }
  double Jacobian(double x) const { return extent*scale*std::exp(x*extent); }

  Result<TubeGeometryPoint> GeometryAtCoordinate(double s) const {
    const auto g=input.geometry(s);
    if(!g.ok())return g;
    if(g.value.coordinateM!=s||!Positive(g.value.radiusM)||!Positive(g.value.areaM2)||
        !std::isfinite(g.value.effectivePotentialM2S2)||
        !std::isfinite(g.value.dLogAreaDsPerM)||!std::isfinite(g.value.dPotentialDsMPerS2))
      return Result<TubeGeometryPoint>::Failure(StatusCode::InvalidConfiguration,
          "tube-isopoly geometry is nonfinite, nonpositive, or evaluated at a different coordinate");
    return g;
  }

  Result<TubeGeometryPoint> Geometry(double x) const {
    return GeometryAtCoordinate(Coordinate(x));
  }

  Result<Equation> At(double x,const TubeGeometryPoint& g,double logSpeed) const {
    Equation eq;
    const double u=std::exp(logSpeed);
    eq.density=input.massFluxKgPerS/(u*g.areaM2);
    if(!Positive(eq.density)||!Positive(u))return Result<Equation>::Failure(
        StatusCode::NumericalFailure,"tube-isopoly density or speed became unrepresentable");
    for(const auto& species:input.species) {
      const auto q=Thermal(species,g.coordinateM,eq.density);
      if(!Positive(q.temperature)||!Positive(q.pressure))return Result<Equation>::Failure(
          StatusCode::NumericalFailure,"tube-isopoly species temperature or pressure is unrepresentable");
      eq.pressure+=q.pressure;
      eq.soundSquared+=(1+q.densityExponent)*q.pressure/eq.density;
      eq.explicitPressureGradient+=q.pressure*q.explicitLogGradient;
    }
    // Steady momentum after substituting rho = massFlux/(A*u). In log(u)
    // the sonic denominator is u^2-c_*^2. Both numerator and denominator
    // vanish at a regular sonic point. No sustaining force is added.
    eq.denominator=u*u-eq.soundSquared;
    eq.numerator=Jacobian(x)*(eq.soundSquared*g.dLogAreaDsPerM-
        eq.explicitPressureGradient/eq.density-g.dPotentialDsMPerS2);
    if(!Positive(eq.pressure)||!Positive(eq.soundSquared)||
        !std::isfinite(eq.denominator)||!std::isfinite(eq.numerator))
      return Result<Equation>::Failure(StatusCode::NumericalFailure,
          "tube-isopoly thermal or momentum state is unrepresentable");
    return Result<Equation>::Success(eq);
  }

  Result<Equation> At(double x,double v) const {
    const auto g=Geometry(x);
    if(!g.ok())return Result<Equation>::Failure(g.status.code,g.status.message);
    return At(x,g.value,v);
  }

  Result<double> SonicSpeed(double x,const TubeGeometryPoint& g) const {
    // For the supported gamma >= 1 the derivative of u^2-c_*^2 with
    // respect to log(u) is positive. Bracket its only zero at fixed x.
    double isothermal=0;
    for(const auto& q:input.species)
      isothermal+=q.particlesPerKg*Constants::kBoltzmannJPerK*q.referenceTemperatureK;
    const double guess=0.5*std::log(isothermal);
    for(double width=1;width<=64;width*=2) {
      double lo=guess-width, hi=guess+width;
      const auto left=At(x,g,lo), right=At(x,g,hi);
      if(!left.ok()||!right.ok())return Result<double>::Failure(StatusCode::NumericalFailure,
          "tube-isopoly cannot evaluate the sonic-speed bracket");
      if(left.value.denominator>0||right.value.denominator<0)continue;
      for(int i=0;i<64&&hi-lo>roundoff;++i) {
        const double mid=lo+(hi-lo)/2;
        const auto m=At(x,g,mid);
        if(!m.ok())return Result<double>::Failure(m.status.code,m.status.message);
        if(m.value.denominator<0)lo=mid;
        else hi=mid;
      }
      return Result<double>::Success(lo+(hi-lo)/2);
    }
    return Result<double>::Failure(StatusCode::NumericalFailure,
        "tube-isopoly could not bracket the sonic speed");
  }

  Result<std::pair<double,double>> Sonic(double x) const {
    using Return=Result<std::pair<double,double>>;
    const auto g=Geometry(x);
    if(!g.ok())return Return::Failure(g.status.code,g.status.message);
    const auto v=SonicSpeed(x,g.value);
    if(!v.ok())return Return::Failure(v.status.code,v.status.message);
    const auto eq=At(x,g.value,v.value);
    if(!eq.ok())return Return::Failure(eq.status.code,eq.status.message);
    return Return::Success({v.value,eq.value.numerator});
  }

  Result<double> Slope(double x,double v,int side) const {
    const auto eq=At(x,v);
    if(!eq.ok())return Result<double>::Failure(eq.status.code,eq.status.message);
    if(side*eq.value.denominator<=roundoff*eq.value.soundSquared)
      return Result<double>::Failure(StatusCode::NoPhysicalSolution,
          "tube-isopoly integration reached a second sonic point or the wrong branch");
    return Result<double>::Success(eq.value.numerator/eq.value.denominator);
  }
};

struct ShotNode { double x, logSpeed, slope; };

// Differentiating D*v'=N at the sonic point yields a quadratic for v'.
// A wind from subsonic to supersonic must have dD/dx>0 there. A degenerate
// root or two admissible slopes do not specify a unique regular crossing.
Result<double> CriticalSlope(const Problem& p,double x,double v) {
  using Return=Result<double>;
  const double h=std::min({1e-6,x/8,(1-x)/8});
  if(h<roundoff)return Return::Failure(StatusCode::NoPhysicalSolution,
      "tube-isopoly sonic point is unresolved at a declared endpoint");
  const double hv=1e-6;
  const auto left=p.At(x-h,v),right=p.At(x+h,v),slow=p.At(x,v-hv),fast=p.At(x,v+hv);
  for(const auto* value:{&left,&right,&slow,&fast})
    if(!value->ok())return Return::Failure(value->status.code,value->status.message);
  const double dx=(right.value.denominator-left.value.denominator)/(2*h);
  const double dv=(fast.value.denominator-slow.value.denominator)/(2*hv);
  const double nx=(right.value.numerator-left.value.numerator)/(2*h);
  const double nv=(fast.value.numerator-slow.value.numerator)/(2*hv);
  const double b=dx-nv, discriminant=b*b+4*dv*nx;
  if(!Positive(dv)||!Positive(discriminant))return Return::Failure(
      StatusCode::NoPhysicalSolution,"tube-isopoly critical point has no regular real slope");
  const double q=-0.5*(b+std::copysign(std::sqrt(discriminant),b));
  const std::array<double,2> slopes={q/dv,-nx/q};
  double chosen=0;
  unsigned accelerating=0;
  for(const double slope:slopes)if(std::isfinite(slope)&&
      dx+dv*slope>1e-8*(std::abs(dx)+std::abs(dv*slope))) {
    chosen=slope;
    ++accelerating;
  }
  if(accelerating!=1)return Return::Failure(StatusCode::NoPhysicalSolution,
      "tube-isopoly critical point has no unique subsonic-to-supersonic slope");
  return Return::Success(chosen);
}

Result<double> RK4(const Problem& p,double x,double v,double h,double end,int side) {
  const auto k1=p.Slope(x,v,side);
  if(!k1.ok())return k1;
  const auto k2=p.Slope(x+h/2,v+h*k1.value/2,side);
  if(!k2.ok())return k2;
  const auto k3=p.Slope(x+h/2,v+h*k2.value/2,side);
  if(!k3.ok())return k3;
  const auto k4=p.Slope(end,v+h*k3.value,side);
  if(!k4.ok())return k4;
  return Result<double>::Success(v+h*(k1.value+2*k2.value+2*k3.value+k4.value)/6);
}

Result<std::vector<ShotNode>> Shoot(const Problem& p,ShotNode sonic,int side) {
  using Return=Result<std::vector<ShotNode>>;
  const double target=side<0?0:1;
  const double distance=std::abs(target-sonic.x);
  const double launch=std::min({distance/4,1e-6,
      std::sqrt(p.input.relativeTolerance)*0.05/(1+std::abs(sonic.slope))});
  if(launch<roundoff)return Return::Failure(StatusCode::NoPhysicalSolution,
      "tube-isopoly cannot start a resolvable critical-point shot");
  // Start either side from the regular slope. The launch truncation is
  // second order in distance. Step doubling controls both the endpoint
  // integration error and the cubic-Hermite dense interpolant's midpoint.
  double x=sonic.x+side*launch, v=sonic.logSpeed+side*launch*sonic.slope;
  auto current=p.Slope(x,v,side);
  if(!current.ok())return Return::Failure(current.status.code,current.status.message);
  std::vector<ShotNode> nodes={{x,v,current.value}};
  const double stepLimit=1.0/p.input.criticalScanIntervals;
  double step=launch;
  for(unsigned i=0;i<p.input.maximumIntegrationAttempts;++i) {
    const double remaining=std::abs(target-x);
    const double end=remaining<=step?target:x+side*step;
    const double h=end-x;
    if(std::abs(h)<roundoff)return Return::Failure(StatusCode::NoPhysicalSolution,
        "tube-isopoly shot cannot resolve the requested branch across the support");
    const auto full=RK4(p,x,v,h,end,side);
    const auto mid=RK4(p,x,v,h/2,x+h/2,side);
    const auto fine=mid.ok()?RK4(p,x+h/2,mid.value,h/2,end,side):mid;
    const auto next=fine.ok()?p.Slope(end,fine.value,side):fine;
    bool complete=true;
    for(const auto* state:{&full,&mid,&fine,&next})if(!state->ok()) {
      if(state->status.code!=StatusCode::NoPhysicalSolution)
        return Return::Failure(state->status.code,state->status.message);
      complete=false;
    }
    const double error=complete?std::max(std::abs(fine.value-full.value)/15,
        std::abs(mid.value-0.5*(v+fine.value)-h*(current.value-next.value)/8)):
        std::numeric_limits<double>::infinity();
    if(error>p.input.relativeTolerance) { step*=0.5; continue; }
    x=end;
    v=fine.value;
    current=next;
    nodes.push_back({x,v,current.value});
    if(x==target)return Return::Success(std::move(nodes));
    // Avoid leaving a terminal interval shorter than a representable step.
    step=std::min({stepLimit,1.75*step,std::abs(target-x)});
  }
  return Return::Failure(StatusCode::NumericalFailure,
      "tube-isopoly shot exhausted the configured integration attempts");
}

}  // namespace

Result<IsoPolyTubeWind> IsoPolyTubeWind::Solve(const IsoPolyTubeInput& input) {
  using Return=Result<IsoPolyTubeWind>;
  if(!input.geometry||input.species.empty()||!std::isfinite(input.beginM)||
      !std::isfinite(input.endM)||!(input.endM>input.beginM)||
      !Positive(input.massFluxKgPerS)||input.criticalScanIntervals<16||
      input.maximumIntegrationAttempts==0||!std::isfinite(input.relativeTolerance)||
      input.relativeTolerance<roundoff||input.relativeTolerance>1e-3)
    return Return::Failure(StatusCode::InvalidConfiguration,
        "tube-isopoly requires finite covered geometry, a mass anchor and valid numerical controls");
  for(const auto& q:input.species)if(!Positive(q.particlesPerKg)||
      !Positive(q.referenceTemperatureK)||!Positive(q.referenceDensityKgM3)||
      !(q.adiabaticIndex>1&&std::isfinite(q.adiabaticIndex))||
      !(q.polytropicIndex>=1&&std::isfinite(q.polytropicIndex))||
      (q.polytropicIndex!=1&&!(std::isfinite(q.blendBeginM)&&
          std::isfinite(q.blendEndM)&&q.blendEndM>q.blendBeginM)))
    return Return::Failure(StatusCode::InvalidConfiguration,
        "tube-isopoly species closure needs a positive reference and a regular thermal transition");

  const auto first=input.geometry(input.beginM);
  if(!first.ok())return Return::Failure(first.status.code,first.status.message);
  const double span=input.endM-input.beginM;
  const double scale=first.value.radiusM;
  // Solve in a normalized logarithmic coordinate. This gives the inner
  // acceleration region resolution without requiring a separate spatial
  // integrator for the long outer tail. The geometry and reported gradients
  // still use metres and include the coordinate Jacobian.
  const double extent=std::log1p(span/scale);
  if(!Positive(scale)||!Positive(extent))return Return::Failure(StatusCode::InvalidConfiguration,
      "tube-isopoly coordinate scale is invalid");
  const Problem p{input,scale,extent};
  auto previous=p.Sonic(0);
  if(!previous.ok())return Return::Failure(previous.status.code,previous.status.message);
  double previousX=0;
  std::vector<double> candidates;
  // A regular critical point must make both the momentum denominator and
  // numerator vanish. Sonic(x) enforces the first constraint and returns
  // the remaining residual for this bracket search. Scanning and bisection
  // only propose critical points. Full-support two-sided shooting below
  // decides whether a candidate is a physical wind solution.
  for(unsigned i=0;i<input.criticalScanIntervals;++i) {
    const double x=static_cast<double>(i+1)/input.criticalScanIntervals;
    const auto value=p.Sonic(x);
    if(!value.ok())return Return::Failure(value.status.code,value.status.message);
    if(value.value.second==0&&x<1)candidates.push_back(x);
    if(previous.value.second!=0&&value.value.second!=0&&
        std::signbit(previous.value.second)!=std::signbit(value.value.second)) {
      double lo=previousX, hi=x;
      for(int n=0;n<64&&hi-lo>roundoff;++n) {
        const double mid=lo+(hi-lo)/2;
        const auto probe=p.Sonic(mid);
        if(!probe.ok())return Return::Failure(probe.status.code,probe.status.message);
        if(probe.value.second==0) { lo=hi=mid; break; }
        if(std::signbit(probe.value.second)==std::signbit(previous.value.second))lo=mid;
        else hi=mid;
      }
      candidates.push_back(lo+(hi-lo)/2);
    }
    previous=value;
    previousX=x;
  }

  if(std::all_of(input.species.begin(),input.species.end(),[&](const IsoPolySpecies& q){
        return q.polytropicIndex==1||q.blendBeginM>=input.endM;
      })&&!candidates.empty()) {
    // The constant-temperature nozzle has a first integral with spatial
    // part H=c^2 ln(A)-Phi. At a regular sonic point the speed-dependent
    // part has its absolute minimum. A global solution requires H to have
    // an equally deep, UNIQUE interior minimum; numerical shooting alone
    // can step over a second, exactly tied minimum by roundoff/truncation.
    const double sound=std::exp(2*previous.value.first);
    const auto nozzle=[&](double x)->Result<double> {
      const auto g=p.Geometry(x);
      if(!g.ok())return Result<double>::Failure(g.status.code,g.status.message);
      const double areaTerm=sound*(std::log(g.value.areaM2)-std::log(first.value.areaM2));
      const double pressurePotential=g.value.effectivePotentialM2S2-first.value.effectivePotentialM2S2;
      return Result<double>::Success(areaTerm-pressurePotential);
    };
    std::vector<double> minima;
    minima.reserve(candidates.size());
    for(const double x:candidates) {
      const auto h=nozzle(x);
      if(!h.ok())return Return::Failure(h.status.code,h.status.message);
      minima.push_back(h.value);
    }
    const auto minimum=std::min_element(minima.begin(),minima.end());
    const double margin=8*input.relativeTolerance*sound+
        roundoff*(std::abs(first.value.effectivePotentialM2S2)+std::abs(*minimum));
    const auto end=nozzle(1);
    if(!end.ok())return Return::Failure(end.status.code,end.status.message);
    if(std::min(0.0,end.value)<=*minimum+margin||
        std::count_if(minima.begin(),minima.end(),[&](double h){return h<=*minimum+margin;})!=1)
      return Return::Failure(StatusCode::NoPhysicalSolution,
          "tube-isopoly isothermal nozzle lacks a unique interior global minimum");
    const double selected=candidates[static_cast<std::size_t>(minimum-minima.begin())];
    candidates.assign(1,selected);
  }

  IsoPolyTubeWind solved;
  solved.input_=input;
  solved.coordinateScaleM_=scale;
  solved.coordinateExtent_=extent;
  for(double x:candidates) {
    const auto value=p.Sonic(x);
    if(!value.ok())return Return::Failure(value.status.code,value.status.message);
    const auto slope=CriticalSlope(p,x,value.value.first);
    if(!slope.ok()&&slope.status.code==StatusCode::NoPhysicalSolution)continue;
    if(!slope.ok())return Return::Failure(slope.status.code,slope.status.message);
    const ShotNode sonic{x,value.value.first,slope.value};
    // Follow the same sonic slope in both directions. Reject candidates
    // that cannot stay on the subsonic inner and supersonic outer branches
    // over the complete guard domain. No radial extension from a partial
    // solution can provide valid thermodynamic or shock inputs.
    const auto inner=Shoot(p,sonic,-1);
    if(!inner.ok()&&inner.status.code==StatusCode::NoPhysicalSolution)continue;
    if(!inner.ok())return Return::Failure(inner.status.code,inner.status.message);
    const auto outer=Shoot(p,sonic,1);
    if(!outer.ok()&&outer.status.code==StatusCode::NoPhysicalSolution)continue;
    if(!outer.ok())return Return::Failure(outer.status.code,outer.status.message);
    if(!solved.nodes_.empty())return Return::Failure(StatusCode::NoPhysicalSolution,
        "tube-isopoly admits multiple full-support transonic branches");
    for(auto n=inner.value.rbegin();n!=inner.value.rend();++n)
      solved.nodes_.push_back({n->x,n->logSpeed,n->slope});
    solved.nodes_.push_back({sonic.x,sonic.logSpeed,sonic.slope});
    for(const auto& n:outer.value)solved.nodes_.push_back({n.x,n.logSpeed,n.slope});
    solved.criticalCoordinateM_=p.Coordinate(x);
  }
  if(solved.nodes_.empty())return Return::Failure(StatusCode::NoPhysicalSolution,
      "tube-isopoly has no regular transonic branch across the declared support");
  return Return::Success(std::move(solved));
}

Result<IsoPolyTubeSample> IsoPolyTubeWind::Evaluate(double s) const {
  using Return=Result<IsoPolyTubeSample>;
  if(!std::isfinite(s)||s<input_.beginM||s>input_.endM||nodes_.size()<2)
    return Return::Failure(StatusCode::OutOfDomain,"tube-isopoly query is outside solved support");
  const Problem p{input_,coordinateScaleM_,coordinateExtent_};
  const double x=s==input_.endM?1:p.X(s);
  const auto upper=std::upper_bound(nodes_.begin(),nodes_.end(),x,
      [](double key,const Node& node){return key<node.x;});
  const auto b=upper==nodes_.end()?upper-1:upper;
  const auto a=b-1;
  const double h=b->x-a->x, z=(x-a->x)/h;
  const double delta=(b->logSpeed-a->logSpeed)/h;
  // Differentiate the same Hermite polynomial used for the speed. The mass
  // flux and closure are applied to that speed at the exact query coordinate.
  // Computing derivatives from an unrelated pressure table would break the
  // discrete energy budget and obscure the convergence of shock inputs.
  const double v=a->logSpeed+h*z*(a->slope+
      z*(3*delta-2*a->slope-b->slope+z*(a->slope+b->slope-2*delta)));
  const double dvdx=a->slope+z*(2*(3*delta-2*a->slope-b->slope)+
      3*z*(a->slope+b->slope-2*delta));
  const auto g=p.GeometryAtCoordinate(s);
  if(!g.ok())return Return::Failure(g.status.code,g.status.message);
  const auto eq=p.At(x,g.value,v);
  if(!eq.ok())return Return::Failure(eq.status.code,eq.status.message);
  const double u=std::exp(v), rho=eq.value.density, dvds=dvdx/p.Jacobian(x);
  const double dLogDensity=-dvds-g.value.dLogAreaDsPerM;
  IsoPolyTubeSample value;
  value.coordinateM=s;
  value.wind={g.value.radiusM,g.value.areaM2,u,rho,eq.value.pressure,
      eq.value.soundSquared,g.value.effectivePotentialM2S2};
  value.speedDerivativePerS=u*dvds;
  value.densityDerivativeKgM4=rho*dLogDensity;
  value.pressureDerivativePaPerM=eq.value.soundSquared*value.densityDerivativeKgM4+
      eq.value.explicitPressureGradient;
  if(!std::isfinite(value.speedDerivativePerS)||!std::isfinite(value.densityDerivativeKgM4)||
      !std::isfinite(value.pressureDerivativePaPerM))return Return::Failure(
      StatusCode::NumericalFailure,"tube-isopoly dense gradients are unrepresentable");
  for(const auto& species:input_.species) {
    const auto q=Thermal(species,g.value.coordinateM,rho);
    value.temperatureK.push_back(q.temperature);
    value.partialPressurePa.push_back(q.pressure);
    // Net nonadiabatic energy per unit volume [W/m^3], including the effect
    // of any heat flux divergence and species energy exchange. This reports
    // the heat implied by a closure; the calculation does not supply it.
    const double heat=u*q.pressure/(species.adiabaticIndex-1)*
        (q.explicitLogGradient+(1+q.densityExponent-species.adiabaticIndex)*
        dLogDensity);
    if(!std::isfinite(heat))return Return::Failure(StatusCode::NumericalFailure,
        "tube-isopoly implied heat rate is unrepresentable");
    value.netHeatWPerM3.push_back(heat);
  }
  return Return::Success(std::move(value));
}

}  // namespace SEP::CoronalCME
