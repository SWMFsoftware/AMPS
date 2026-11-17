#include "sep_resonant_collision.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace SEP { namespace Transport { namespace ResonantCollision {
namespace {
const double C=299792458.0;
const double Pi=3.141592653589793238462643383279502884;
Status Error(const std::string& text) {
  return Status::Error(StatusCode::InvalidArgument,text);
}
bool Positive(double x){return std::isfinite(x)&&x>0.0;}
bool Nonnegative(double x){return std::isfinite(x)&&x>=0.0;}
int Sign(Branch branch){return branch==Branch::AlongB?1:-1;}
std::size_t Offset(Branch b,Helicity h) {
  return 2*static_cast<std::size_t>(b)+static_cast<std::size_t>(h);
}
Status CheckPitch(double q) {
  return std::isfinite(q)&&q>=-1.0&&q<=1.0?Status::Ok():
      Error("wave-frame pitch must lie in [-1,1]");
}
double ResonanceScale(const Spectrum& s,const Orbit& o) {
  return std::fabs(o.chargeC)*s.meanMagneticFieldT*o.waveGamma/
      o.momentumInWaveFrameKgMPerS;
}
double DimensionlessMass(const Orbit& o,double a,double b) {
  const double alpha=o.waveSpeedMPerS*o.momentumInWaveFrameKgMPerS/
      o.energyInWaveFrameJ;
  return (b-a)*(1.0+0.5*alpha*(a+b));
}
double Measure(const Orbit& o,double a,double b) {
  // J=Gamma^2*p'/c^2*(E'+u*p'*q), integrated without subtracting
  // two large cumulative masses at a narrow cell boundary.
  return o.waveGamma*o.waveGamma*o.momentumInWaveFrameKgMPerS/(C*C)*
      o.energyInWaveFrameJ*DimensionlessMass(o,a,b);
}
Status ValidateOrbit(const Orbit& o) {
  if((o.branch!=Branch::AlongB&&o.branch!=Branch::AgainstB)||
      !Positive(o.massKg)||!std::isfinite(o.chargeC)||o.chargeC==0.0||
      !std::isfinite(o.waveSpeedMPerS)||std::fabs(o.waveSpeedMPerS)>=C||
      !Positive(o.waveGamma)||!Positive(o.energyInWaveFrameJ)||
      !Nonnegative(o.momentumInWaveFrameKgMPerS)||!Positive(o.invariantJ))
    return Error("collision orbit is incomplete or outside relativistic domain");
  if(o.stationaryInWaveFrame!=(o.momentumInWaveFrameKgMPerS==0))
    return Error("collision orbit stationary-state flag is inconsistent");
  return Status::Ok();
}
std::vector<double> SplitPoints(const Spectrum& s,const Orbit& o,double a,double b) {
  std::vector<double> points{a,b};
  if(a<0&&b>0)points.push_back(0);
  const double scale=ResonanceScale(s,o);
  for(double k:s.kEdgesPerM) {
    const double q=scale/k;
    if(q>a&&q<b)points.push_back(q);
    if(-q>a&&-q<b)points.push_back(-q);
  }
  std::sort(points.begin(),points.end());
  points.erase(std::unique(points.begin(),points.end()),points.end());
  return points;
}
struct Resistance {Status status;bool blocked=false;double value=0;};
Resistance IntegratedResistance(const Spectrum& s,const Orbit& o,double a,double b) {
  Resistance result;result.status=Status::Ok();
  const auto points=SplitPoints(s,o,a,b);
  long double sum=0;
  for(std::size_t j=0;j+1<points.size();++j) {
    const double left=points[j],right=points[j+1];
    const double mid=left+0.5*(right-left);
    if(!(mid>left&&mid<right)) {
      result.status=Status::Error(StatusCode::UnresolvedCoefficient,
          "resonance mesh interval is below floating-point pitch resolution");
      return result;
    }
    const Resonance resonant=EvaluateResonance(s,o,mid);
    if(!resonant.status.ok()){result.status=resonant.status;return result;}
    if(!resonant.active) {
      result.blocked=true;result.value=std::numeric_limits<double>::infinity();
      return result;
    }
    if(left<0&&right>0) {
      result.status=Error("unpartitioned wave-frame resonance gap");return result;
    }
    // Integral |q|/(1-q^2) dq; log1p preserves precision near q=0.
    const double logLeft=std::log1p(-left*left);
    const double logRight=std::log1p(-right*right);
    const double integral=(mid>0?0.5:-0.5)*(logLeft-logRight);
    // The one-sided P_h integrates to the real circular-wave variance.
    // Random phase averaging contributes 1/2 to the force autocorrelation.
    const double denominator=0.5*Pi*o.chargeC*o.chargeC*resonant.magneticPowerT2M;
    if(!Positive(integral)||!Positive(denominator)) {
      result.status=Status::Error(StatusCode::UnresolvedCoefficient,
          "finite-bin resistance underflowed or is not positive");return result;
    }
    sum+=static_cast<long double>(integral)/denominator;
  }
  result.value=static_cast<double>(sum);
  if(!Positive(result.value))result.status=Status::Error(
      StatusCode::UnresolvedCoefficient,"finite-bin resistance overflowed");
  return result;
}
Status OpenUniform(RandomStream& random,double* value) {
  *value=random.UniformOpen01();
  return std::isfinite(*value)&&*value>0&&*value<1?Status::Ok():
      Error("collision random stream must provide an open-unit variate");
}
}

Status ValidateSpectrum(const Spectrum& s) {
  if(!s.polarizationResolved)return Status::Error(StatusCode::UnsupportedConfiguration,
      "resonant scattering requires four explicitly initialized branch/helicity spectra");
  if(!Positive(s.volumeM3)||!Positive(s.meanMagneticFieldT)||
      !Nonnegative(s.alfvenSpeedMPerS)||s.alfvenSpeedMPerS>=C||
      !Positive(s.vacuumPermeabilityHPerM)||!std::isfinite(s.epochS)||
      s.generation==0||s.identity.empty()||s.kEdgesPerM.size()<2)
    return Error("polarized spectrum lacks physical units, identity or positive bin geometry");
  for(std::size_t k=0;k<s.kEdgesPerM.size();++k)
    if(!Positive(s.kEdgesPerM[k])||(k&&s.kEdgesPerM[k]<=s.kEdgesPerM[k-1]))
      return Error("wavenumber edges must be finite, positive and strictly increasing");
  const std::size_t bins=s.kEdgesPerM.size()-1;
  for(const auto& branch:s.magneticEnergyJ) {
    if(branch.size()!=bins)return Error("all four polarized energy arrays must match the bin edges");
    for(double energy:branch)if(!Nonnegative(energy))
      return Error("polarized magnetic bin energies must be finite and nonnegative");
  }
  return Status::Ok();
}

Status ValidateDomain(const Spectrum& s) {
  const auto status=ValidateSpectrum(s);if(!status.ok())return status;
  if(s.domain.provenance.empty())return Status::Error(StatusCode::UnsupportedConfiguration,
      "collision domain requires explicit qualification or validation-study provenance");
  if(s.domain.authority==DomainAuthority::ValidationStudy)return Status::Ok();
  if(s.domain.authority!=DomainAuthority::ExplicitDeclaredLimits||
      !Positive(s.domain.maximumDeltaBOverB)||
      !Positive(s.domain.maximumKThermalIonGyroradius)||
      !Positive(s.domain.thermalIonGyroradiusM)||
      !Positive(s.domain.maximumAlfvenSpeedOverC)||
      !Positive(s.domain.ionInertialLengthM)||
      !Positive(s.domain.maximumKIonInertialLength)||
      !Positive(s.domain.ionCyclotronFrequencyRadPerS)||
      !Positive(s.domain.maximumWaveFrequencyOverIonCyclotron)||
      s.domain.plasmaCompositionIdentity.empty()||s.domain.thermalStateIdentity.empty())
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "quantitative weak-turbulence/MHD limits are not qualified; no numerical physical limits are assumed");
  long double energy=0;
  for(const auto& branch:s.magneticEnergyJ)for(double value:branch)energy+=value;
  const long double variance=2*static_cast<long double>(s.vacuumPermeabilityHPerM)*
      energy/s.volumeM3;
  const long double ratio=std::sqrt(variance)/s.meanMagneticFieldT;
  const long double kRho=static_cast<long double>(s.kEdgesPerM.back())*
      s.domain.thermalIonGyroradiusM;
  const long double kDi=static_cast<long double>(s.kEdgesPerM.back())*
      s.domain.ionInertialLengthM;
  const long double omegaRatio=static_cast<long double>(s.kEdgesPerM.back())*
      s.alfvenSpeedMPerS/s.domain.ionCyclotronFrequencyRadPerS;
  if(!std::isfinite(ratio)||!std::isfinite(kRho)||
      !std::isfinite(kDi)||!std::isfinite(omegaRatio)||
      ratio>s.domain.maximumDeltaBOverB||
      kRho>s.domain.maximumKThermalIonGyroradius||
      s.alfvenSpeedMPerS/C>s.domain.maximumAlfvenSpeedOverC||
      kDi>s.domain.maximumKIonInertialLength||
      omegaRatio>s.domain.maximumWaveFrequencyOverIonCyclotron)
    return Status::Error(StatusCode::OutOfDomain,
        "spectrum exceeds its explicitly qualified amplitude or MHD wavenumber domain");
  return Status::Ok();
}

Status RequireImplementedProductionClosure(bool exchange,bool localEvolution) {
  if(!exchange&&!localEvolution)return Status::Ok();
  return Status::Error(StatusCode::UnsupportedConfiguration,
      "approved polarized weak-Alfven production closure is unavailable: "
      "four initial branch/helicity spectra and qualified physical domain are required; "
      "the joint PIC collision adapter, atomic particle-wave commit and prescribed-background "
      "work ledger are not implemented. The available resonant kernel is validation-only; "
      "legacy two-branch aggregates are nonresonant diagnostics");
}

OrbitResult PrepareOrbit(const Spectrum& s,Branch branch,const Particle& particle) {
  OrbitResult result;result.status=ValidateDomain(s);if(!result.status.ok())return result;
  if((branch!=Branch::AlongB&&branch!=Branch::AgainstB)||
      !Positive(particle.massKg)||!std::isfinite(particle.chargeC)||particle.chargeC==0||
      !std::isfinite(particle.parallelMomentumKgMPerS)||
      !Nonnegative(particle.perpendicularMomentumKgMPerS)) {
    result.status=Error("collision particle requires signed charge, rest mass and finite momenta");
    return result;
  }
  Orbit o;o.branch=branch;o.massKg=particle.massKg;o.chargeC=particle.chargeC;
  o.waveSpeedMPerS=Sign(branch)*s.alfvenSpeedMPerS;
  const double beta=o.waveSpeedMPerS/C;
  o.waveGamma=1/std::sqrt((1-beta)*(1+beta));
  const double energy=std::hypot(particle.massKg*C,
      std::hypot(particle.parallelMomentumKgMPerS,particle.perpendicularMomentumKgMPerS))*C;
  const double parallelWave=o.waveGamma*(particle.parallelMomentumKgMPerS-
      o.waveSpeedMPerS*energy/(C*C));
  o.momentumInWaveFrameKgMPerS=std::hypot(parallelWave,particle.perpendicularMomentumKgMPerS);
  o.energyInWaveFrameJ=std::hypot(particle.massKg*C,o.momentumInWaveFrameKgMPerS)*C;
  o.invariantJ=o.energyInWaveFrameJ/o.waveGamma;
  o.stationaryInWaveFrame=o.momentumInWaveFrameKgMPerS==0;
  o.initialWavePitch=o.stationaryInWaveFrame?0:
      parallelWave/o.momentumInWaveFrameKgMPerS;
  result.status=ValidateOrbit(o);
  if(result.status.ok())result.status=CheckPitch(o.initialWavePitch);
  if(result.status.ok())result.orbit=o;
  return result;
}

ParticleResult ToPlasmaFrame(const Orbit& o,double pitch) {
  ParticleResult result;result.status=ValidateOrbit(o);if(!result.status.ok())return result;
  result.status=CheckPitch(pitch);if(!result.status.ok())return result;
  result.particle.massKg=o.massKg;result.particle.chargeC=o.chargeC;
  result.particle.parallelMomentumKgMPerS=o.waveGamma*(
      o.momentumInWaveFrameKgMPerS*pitch+o.waveSpeedMPerS*o.energyInWaveFrameJ/(C*C));
  result.particle.perpendicularMomentumKgMPerS=o.momentumInWaveFrameKgMPerS*
      std::sqrt((1-pitch)*(1+pitch));
  result.totalEnergyJ=o.waveGamma*(o.energyInWaveFrameJ+
      o.waveSpeedMPerS*o.momentumInWaveFrameKgMPerS*pitch);
  if(!std::isfinite(result.particle.parallelMomentumKgMPerS)||
      !Nonnegative(result.particle.perpendicularMomentumKgMPerS)||!Positive(result.totalEnergyJ))
    result.status=Error("Lorentz orbit transformation overflowed");
  return result;
}

Resonance EvaluateResonance(const Spectrum& s,const Orbit& o,double pitch) {
  Resonance result;result.status=ValidateSpectrum(s);if(!result.status.ok())return result;
  result.status=ValidateOrbit(o);if(!result.status.ok())return result;
  result.status=CheckPitch(pitch);if(!result.status.ok())return result;
  if(o.stationaryInWaveFrame){result.zeroReason="particle stationary in wave frame";return result;}
  if(pitch==0){result.zeroReason="finite-band wave-frame ninety-degree gap";
    result.kPerM=std::numeric_limits<double>::infinity();return result;}
  result.kPerM=ResonanceScale(s,o)/std::fabs(pitch);
  if(result.kPerM<s.kEdgesPerM.front()||result.kPerM>s.kEdgesPerM.back()) {
    result.zeroReason="resonance outside finite spectral support";return result;
  }
  result.bin=static_cast<std::size_t>(std::upper_bound(s.kEdgesPerM.begin(),
      s.kEdgesPerM.end(),result.kPerM)-s.kEdgesPerM.begin()-1);
  if(result.bin+1==s.kEdgesPerM.size())--result.bin;
  result.helicity=(o.chargeC>0)==(pitch>0)?Helicity::Right:Helicity::Left;
  const double energy=s.magneticEnergyJ[Offset(o.branch,result.helicity)][result.bin];
  if(energy==0){result.zeroReason="zero energy in the resonant helicity/bin";return result;}
  const double width=s.kEdgesPerM[result.bin+1]-s.kEdgesPerM[result.bin];
  result.magneticPowerT2M=2*s.vacuumPermeabilityHPerM*energy/(s.volumeM3*width);
  const auto plasma=ToPlasmaFrame(o,pitch);
  if(!plasma.status.ok()){result.status=plasma.status;return result;}
  // nu=pi*Omega^2*P_h/[B0^2*abs(v_parallel-u)], written without
  // subtracting almost equal velocities near the wave-frame resonance gap.
  result.frequencyPerS=Pi*o.chargeC*o.chargeC*C*C*o.waveGamma*
      result.magneticPowerT2M/(plasma.totalEnergyJ*
      o.momentumInWaveFrameKgMPerS*std::fabs(pitch));
  if(!Positive(result.magneticPowerT2M)||!Positive(result.frequencyPerS)) {
    result.status=Status::Error(StatusCode::UnresolvedCoefficient,
        "nonzero resonant power or rate underflowed/overflowed");return result;
  }
  result.active=true;
  return result;
}

GeneratorResult BuildGenerator(const Spectrum& s,const Orbit& o,std::size_t count) {
  GeneratorResult result;result.status=ValidateDomain(s);if(!result.status.ok())return result;
  result.status=ValidateOrbit(o);if(!result.status.ok())return result;
  if(count<2||count>10000000){result.status=Error("pitch mesh count is outside allocation bounds");return result;}
  Generator g;g.spectrum=s;g.orbit=o;
  g.withinDeclaredLimits=s.domain.authority==DomainAuthority::ExplicitDeclaredLimits;
  // Passing caller-supplied numerical bounds is not independent physical
  // qualification and cannot stand in for the unfinished production adapter.
  g.qualifiedForProduction=false;
  g.edges.reserve(count+2*s.kEdgesPerM.size()+2);
  g.edges.push_back(-1);g.edges.push_back(1);
  g.edges.push_back(0);
  if(!o.stationaryInWaveFrame) {
    const double scale=ResonanceScale(s,o);
    for(double k:s.kEdgesPerM) {
      const double edge=scale/k;
      if(edge>0&&edge<1){g.edges.push_back(-edge);g.edges.push_back(edge);}
    }
  }
  std::sort(g.edges.begin(),g.edges.end());
  g.edges.erase(std::unique(g.edges.begin(),g.edges.end()),g.edges.end());
  // Physical support interfaces are authoritative. A redundant optional
  // uniform-mesh point within roundoff of one is omitted; two distinct
  // physical interfaces are never merged or bridged.
  const auto physicalEdges=g.edges;
  for(std::size_t j=1;j<count;++j) {
    const double q=-1+2*static_cast<double>(j)/count;
    const auto next=std::lower_bound(physicalEdges.begin(),physicalEdges.end(),q);
    const double tolerance=8*std::numeric_limits<double>::epsilon();
    const bool nearRight=next!=physicalEdges.end()&&std::fabs(*next-q)<=tolerance;
    const bool nearLeft=next!=physicalEdges.begin()&&std::fabs(*(next-1)-q)<=tolerance;
    if(!nearRight&&!nearLeft)g.edges.push_back(q);
  }
  std::sort(g.edges.begin(),g.edges.end());
  const std::size_t cells=g.edges.size()-1;
  g.centers.resize(cells);g.measure.resize(cells);
  g.leftRatePerS.assign(cells,0);g.rightRatePerS.assign(cells,0);
  g.conductance.assign(cells-1,0);
  for(std::size_t j=0;j<cells;++j) {
    g.centers[j]=g.edges[j]+0.5*(g.edges[j+1]-g.edges[j]);
    // At p_wave=0 every collision is physically a no-op. Use a normalized
    // angular measure only to keep the inert generator well-formed.
    g.measure[j]=o.stationaryInWaveFrame?g.edges[j+1]-g.edges[j]:
        Measure(o,g.edges[j],g.edges[j+1]);
    if(!Positive(g.measure[j])||!(g.centers[j]>g.edges[j]&&g.centers[j]<g.edges[j+1])) {
      result.status=Status::Error(StatusCode::UnresolvedCoefficient,
          "phase-space pitch cell is below representable numerical resolution");return result;
    }
  }
  if(!o.stationaryInWaveFrame)for(std::size_t face=0;face+1<cells;++face) {
    const auto resistance=IntegratedResistance(s,o,g.centers[face],g.centers[face+1]);
    if(!resistance.status.ok()){result.status=resistance.status;return result;}
    if(resistance.blocked)continue;
    const double conductance=1/resistance.value;
    g.conductance[face]=conductance;
    g.rightRatePerS[face]=conductance/g.measure[face];
    g.leftRatePerS[face+1]=conductance/g.measure[face+1];
    if(!Positive(conductance)||!Positive(g.rightRatePerS[face])||
        !Positive(g.leftRatePerS[face+1])) {
      result.status=Status::Error(StatusCode::UnresolvedCoefficient,
          "finite-volume collision rate underflowed or overflowed");return result;
    }
  }
  result.generator=std::move(g);result.status=Status::Ok();return result;
}

Status ProbabilityDerivative(const Generator& g,const std::vector<double>& p,
                              std::vector<double>* derivative) {
  if(!derivative||p.size()!=g.centers.size()||g.measure.size()!=p.size()||
      g.conductance.size()+1!=p.size())return Error("probability vector does not match generator");
  for(double value:p)if(!Nonnegative(value))return Error("probability mass is not finite/nonnegative");
  derivative->assign(p.size(),0);
  for(std::size_t j=0;j+1<p.size();++j) {
    if(g.conductance[j]==0)continue;
    const double flux=g.conductance[j]*(p[j]/g.measure[j]-p[j+1]/g.measure[j+1]);
    if(!std::isfinite(flux))return Error("collision probability flux overflowed");
    (*derivative)[j]-=flux;(*derivative)[j+1]+=flux;
  }
  return Status::Ok();
}

std::size_t CellAt(const Generator& g,double q) {
  if(!CheckPitch(q).ok()||g.edges.size()<2)return std::numeric_limits<std::size_t>::max();
  if(q==1)return g.edges.size()-2;
  return static_cast<std::size_t>(std::upper_bound(g.edges.begin(),g.edges.end(),q)-g.edges.begin()-1);
}

ScalarResult TotalRate(const Generator& g,std::size_t cell) {
  ScalarResult result;
  if(cell>=g.centers.size()||g.leftRatePerS.size()!=g.centers.size()||
      g.rightRatePerS.size()!=g.centers.size()) {
    result.status=Error("collision cell is outside its generator");return result;
  }
  const double left=g.leftRatePerS[cell],right=g.rightRatePerS[cell],rate=left+right;
  if(!Nonnegative(left)||!Nonnegative(right)||!std::isfinite(rate)) {
    result.status=Error("collision generator contains an invalid event rate");return result;
  }
  result.value=rate;result.status=Status::Ok();return result;
}

Jump SampleJumpAtEvent(const Generator& g,std::size_t cell,RandomStream& random) {
  Jump result;result.targetCell=cell;
  const auto total=TotalRate(g,cell);result.status=total.status;
  if(!result.status.ok())return result;
  if(total.value==0){result.waitingTimeS=std::numeric_limits<double>::infinity();return result;}
  const double left=g.leftRatePerS[cell];
  double directionDraw,pitchDraw;
  result.status=OpenUniform(random,&directionDraw);if(!result.status.ok())return result;
  result.status=OpenUniform(random,&pitchDraw);if(!result.status.ok())return result;
  result.targetCell=directionDraw<left/total.value?cell-1:cell+1;
  if(result.targetCell>=g.centers.size()) {
    result.status=Error("collision jump has an invalid target or waiting time");return result;
  }
  const double a=g.edges[result.targetCell],b=g.edges[result.targetCell+1];
  const double alpha=g.orbit.waveSpeedMPerS*g.orbit.momentumInWaveFrameKgMPerS/
      g.orbit.energyInWaveFrameJ;
  const double densityAtLeft=1+alpha*a;
  const double mass=DimensionlessMass(g.orbit,a,b);
  const double discriminant=densityAtLeft*densityAtLeft+2*alpha*pitchDraw*mass;
  const double displacement=2*pitchDraw*mass/(densityAtLeft+std::sqrt(discriminant));
  result.targetWavePitch=a+displacement;
  if(!(discriminant>0)||!std::isfinite(result.targetWavePitch)||
      result.targetWavePitch<=a||result.targetWavePitch>=b) {
    result.status=Error("phase-space-weighted target pitch inversion failed");return result;
  }
  result.active=true;return result;
}

Jump SampleJump(const Generator& g,std::size_t cell,RandomStream& random) {
  Jump result;const auto rate=TotalRate(g,cell);result.status=rate.status;
  if(!result.status.ok())return result;
  if(rate.value==0){result.waitingTimeS=std::numeric_limits<double>::infinity();return result;}
  double draw;result.status=OpenUniform(random,&draw);if(!result.status.ok())return result;
  result=SampleJumpAtEvent(g,cell,random);
  if(!result.status.ok())return result;
  result.waitingTimeS=-std::log(draw)/rate.value;
  if(!Positive(result.waitingTimeS))result.status=Error("collision waiting time is not positive/finite");
  return result;
}

namespace {
// One deterministic implementation of the elastic wave-frame increments and
// resonance-bin partition serves generation and later transaction validation.
// It never chooses a jump: before/after pitches are supplied, and no RNG is
// accessed. Discrete neighbor eligibility remains a Generator responsibility.
Proposal ElasticProposal(const Spectrum& spectrum,const Orbit& orbit,double before,double after,double weight) {
  Proposal result;result.status=CheckPitch(before);if(!result.status.ok())return result;
  result.status=CheckPitch(after);if(!result.status.ok())return result;
  if(!Positive(weight)){result.status=Error("collision macro-particle weight must be positive");return result;}
  const auto first=ToPlasmaFrame(orbit,before),last=ToPlasmaFrame(orbit,after);
  if(!first.status.ok()){result.status=first.status;return result;}
  if(!last.status.ok()){result.status=last.status;return result;}
  result.before=first.particle;result.after=last.particle;
  result.beforeWavePitch=before;result.afterWavePitch=after;
  // Exact increments avoid subtracting two rest-mass-dominated energies.
  result.parallelMomentumChangeKgMPerS=orbit.waveGamma*
      orbit.momentumInWaveFrameKgMPerS*(after-before)*weight;
  result.particleEnergyChangeJ=orbit.waveSpeedMPerS*result.parallelMomentumChangeKgMPerS;
  result.spectrumIdentity=spectrum.identity;result.generation=spectrum.generation;
  result.epochS=spectrum.epochS;
  if(!std::isfinite(result.parallelMomentumChangeKgMPerS)||!std::isfinite(result.particleEnergyChangeJ)||
      (before!=after&&orbit.momentumInWaveFrameKgMPerS>0&&
       (result.parallelMomentumChangeKgMPerS==0||(orbit.waveSpeedMPerS!=0&&result.particleEnergyChangeJ==0)))) {
    result.status=Status::Error(StatusCode::UnresolvedCoefficient,"elastic proposal increment is outside floating-point resolution");return result;
  }
  if(before==after){result.status=Status::Ok();return result;}
  const auto points=SplitPoints(spectrum,orbit,std::min(before,after),std::max(before,after));
  for(std::size_t j=0;j+1<points.size();++j) {
    const double midpoint=points[j]+0.5*(points[j+1]-points[j]);
    const auto resonance=EvaluateResonance(spectrum,orbit,midpoint);
    if(!resonance.status.ok()){result.status=resonance.status;return result;}
    if(!resonance.active) {
      result.status=Status::Error(StatusCode::UnsupportedConfiguration,
          "collision proposal crosses a finite-band or empty-helicity resonance gap");
      result.exchange.clear();return result;
    }
    BinExchange exchange;exchange.branch=orbit.branch;exchange.helicity=resonance.helicity;
    exchange.bin=resonance.bin;
    const double fraction=(points[j+1]-points[j])/std::fabs(after-before);
    exchange.particleEnergyChangeJ=result.particleEnergyChangeJ*fraction;
    exchange.magneticEnergyChangeJ=-0.5*exchange.particleEnergyChangeJ;
    if(!std::isfinite(exchange.particleEnergyChangeJ)||!std::isfinite(exchange.magneticEnergyChangeJ)||
        (result.particleEnergyChangeJ!=0&&(exchange.particleEnergyChangeJ==0||exchange.magneticEnergyChangeJ==0))) {
      result.status=Status::Error(StatusCode::UnresolvedCoefficient,"resonance-bin exchange is outside floating-point resolution");
      result.exchange.clear();return result;
    }
    result.exchange.push_back(exchange);
  }
  result.status=Status::Ok();return result;
}

bool RoundoffMatches(double actual,double expected,long double scale) {
  // 512 binary64 eps covers the two Lorentz round trips and partition algebra;
  // scale has the quantity's own units. No max(1,SI value) absolute tolerance
  // can hide a forged tiny momentum or energy. Exact zeros require exact zero
  // when there is no nonzero physical reference scale.
  return std::isfinite(actual)&&std::isfinite(expected)&&scale>=0&&
      std::fabs(static_cast<long double>(actual)-expected)<=
      512*std::numeric_limits<double>::epsilon()*scale;
}
bool IncrementMatches(double actual,double expected) {
  return RoundoffMatches(actual,expected,std::max(std::fabs(static_cast<long double>(actual)),
      std::fabs(static_cast<long double>(expected))));
}
} // namespace

Proposal ProposeJump(const Generator& g,double before,double after,double weight) {
  Proposal rejected;rejected.status=CheckPitch(before);if(!rejected.status.ok())return rejected;
  rejected.status=CheckPitch(after);if(!rejected.status.ok())return rejected;
  const auto fromCell=CellAt(g,before),toCell=CellAt(g,after);
  if(fromCell>=g.centers.size()||toCell>=g.centers.size()||
      (fromCell>toCell?fromCell-toCell:toCell-fromCell)>1) {
    rejected.status=Error("collision proposal must join the same or neighboring pitch cells");return rejected;
  }
  return ElasticProposal(g.spectrum,g.orbit,before,after,weight);
}

Status ValidateProposal(const Spectrum& spectrum,Branch branch,const Proposal& proposal,double weight) {
  if(!proposal.status.ok())return Error("rejected collision proposal cannot enter a physical transaction");
  if(proposal.spectrumIdentity!=spectrum.identity||proposal.generation!=spectrum.generation||proposal.epochS!=spectrum.epochS)
    return Error("collision proposal spectrum identity/generation/epoch differs from the pinned state");
  const auto orbit=PrepareOrbit(spectrum,branch,proposal.before);if(!orbit.status.ok())return orbit.status;
  const auto expected=ElasticProposal(spectrum,orbit.orbit,proposal.beforeWavePitch,proposal.afterWavePitch,weight);
  if(!expected.status.ok())return expected.status;
  const long double momentumScale=orbit.orbit.waveGamma*(
      static_cast<long double>(orbit.orbit.momentumInWaveFrameKgMPerS)+
      std::fabs(static_cast<long double>(orbit.orbit.waveSpeedMPerS))*orbit.orbit.energyInWaveFrameJ/(C*C));
  const Particle* actual[2]={&proposal.before,&proposal.after};
  const Particle* canonical[2]={&expected.before,&expected.after};
  for(int i=0;i<2;++i)if(actual[i]->massKg!=canonical[i]->massKg||actual[i]->chargeC!=canonical[i]->chargeC||
      !RoundoffMatches(actual[i]->parallelMomentumKgMPerS,canonical[i]->parallelMomentumKgMPerS,momentumScale)||
      !RoundoffMatches(actual[i]->perpendicularMomentumKgMPerS,canonical[i]->perpendicularMomentumKgMPerS,momentumScale))
    return Status::Error(StatusCode::UnresolvedCoefficient,
        "proposal phase disagrees with its elastic wave-frame orbit beyond dimensionally scaled roundoff; phase was not moved");
  if(!IncrementMatches(proposal.parallelMomentumChangeKgMPerS,expected.parallelMomentumChangeKgMPerS)||
      !IncrementMatches(proposal.particleEnergyChangeJ,expected.particleEnergyChangeJ)||proposal.exchange.size()!=expected.exchange.size())
    return Error("collision proposal macro-weight/increments or resonance partition differ from the canonical elastic law");
  for(std::size_t j=0;j<expected.exchange.size();++j) {
    const auto& a=proposal.exchange[j];const auto& e=expected.exchange[j];
    if(a.branch!=e.branch||a.helicity!=e.helicity||a.bin!=e.bin||
        !IncrementMatches(a.particleEnergyChangeJ,e.particleEnergyChangeJ)||
        !IncrementMatches(a.magneticEnergyChangeJ,e.magneticEnergyChangeJ))
      return Error("collision proposal branch/helicity/bin debit differs from the canonical resonant path");
  }
  return Status::Ok();
}

} } }
