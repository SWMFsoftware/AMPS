#include "sep_resonant_collision_validation.h"
#include "sep_resonant_collision.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace SEP { namespace Testing { namespace {
namespace RC=Transport::ResonantCollision;
const double C=299792458.0;
const double Pi=3.14159265358979323846;
bool Relative(double a,double b,double tolerance=1e-11) {
  return std::isfinite(a)&&std::isfinite(b)&&
      std::fabs(a-b)<=tolerance*std::max(std::fabs(b),1e-300);
}
Result Finish(bool pass,const char* message) {
  Result result;result.status=pass?Status::Pass:Status::Fail;result.message=message;
  result.configuration.push_back("capability=neutral-collision-core;production_integration=false");
  result.configuration.push_back("physical_domain=validation-study-only;no-qualified-numerical-limits");
  return result;
}
Result Failure(const Transport::Status& status,const char* stage) {
  auto result=Finish(false,stage);result.message+=": "+status.message;return result;
}
RC::Spectrum Spectrum() {
  RC::Spectrum s;s.kEdgesPerM={1e-9,1e-8,1e-7,1e-6};
  s.magneticEnergyJ[0]={1e-8,2e-8,3e-8};
  s.magneticEnergyJ[1]={4e-8,5e-8,6e-8};
  s.magneticEnergyJ[2]={7e-8,8e-8,9e-8};
  s.magneticEnergyJ[3]={1e-7,1.1e-7,1.2e-7};
  s.polarizationResolved=true;s.volumeM3=1e6;s.meanMagneticFieldT=5e-9;
  s.alfvenSpeedMPerS=4e5;s.vacuumPermeabilityHPerM=1.256637062e-6;
  s.epochS=0;s.generation=1;s.identity="controlled-four-helicity:g1";
  s.domain.authority=RC::DomainAuthority::ValidationStudy;
  s.domain.provenance="manufactured polarized-slab mathematical verification";
  return s;
}
RC::Particle Particle(double speed=2e7,double mu=0.4) {
  RC::Particle p;p.massKg=1.67262192369e-27;p.chargeC=1.602176634e-19;
  const double momentum=p.massKg*speed/std::sqrt((1-speed/C)*(1+speed/C));
  p.parallelMomentumKgMPerS=momentum*mu;
  p.perpendicularMomentumKgMPerS=momentum*std::sqrt((1-mu)*(1+mu));
  return p;
}
Result FourSpectrumContract() {
  auto s=Spectrum();const auto p=Particle();
  auto orbit=RC::PrepareOrbit(s,RC::Branch::AlongB,p);
  if(!orbit.status.ok())return Failure(orbit.status,"orbit preparation");
  const auto resonance=RC::EvaluateResonance(s,orbit.orbit,0.4);
  if(!resonance.status.ok()||!resonance.active)return Failure(resonance.status,"active resonance");
  const auto plasma=RC::ToPlasmaFrame(orbit.orbit,0.4);
  const double energy=std::hypot(p.massKg*C,std::hypot(
      plasma.particle.parallelMomentumKgMPerS,plasma.particle.perpendicularMomentumKgMPerS))*C;
  const double omega=std::fabs(p.chargeC)*s.meanMagneticFieldT*C*C/energy;
  const double vParallel=plasma.particle.parallelMomentumKgMPerS*C*C/energy;
  const double k=omega/std::fabs(vParallel-s.alfvenSpeedMPerS);
  const auto bin=static_cast<std::size_t>(std::upper_bound(s.kEdgesPerM.begin(),s.kEdgesPerM.end(),k)-s.kEdgesPerM.begin()-1);
  const double power=2*s.vacuumPermeabilityHPerM*s.magneticEnergyJ[1][bin]/
      (s.volumeM3*(s.kEdgesPerM[bin+1]-s.kEdgesPerM[bin]));
  const double frequency=Pi*omega*omega*power/
      (s.meanMagneticFieldT*s.meanMagneticFieldT*std::fabs(vParallel-s.alfvenSpeedMPerS));
  bool pass=orbit.status.ok()&&resonance.status.ok()&&resonance.active&&
      resonance.helicity==RC::Helicity::Right&&Relative(resonance.kPerM,k)&&
      Relative(resonance.magneticPowerT2M,power)&&Relative(resonance.frequencyPerS,frequency);
  auto electron=p;electron.chargeC=-p.chargeC;
  const auto negative=RC::PrepareOrbit(s,RC::Branch::AlongB,electron);
  const auto flipped=RC::EvaluateResonance(s,negative.orbit,0.4);
  pass=pass&&flipped.status.ok()&&flipped.active&&flipped.helicity==RC::Helicity::Left&&
      !Relative(flipped.magneticPowerT2M,resonance.magneticPowerT2M);
  s.polarizationResolved=false;
  pass=pass&&!RC::PrepareOrbit(s,RC::Branch::AlongB,p).status.ok();
  s=Spectrum();s.domain.authority=RC::DomainAuthority::Unspecified;
  pass=pass&&!RC::PrepareOrbit(s,RC::Branch::AlongB,p).status.ok();
  s=Spectrum();const auto generator=RC::BuildGenerator(s,orbit.orbit,32);
  pass=pass&&generator.status.ok()&&!generator.generator.qualifiedForProduction;
  // Deliberately overbroad numerical bounds test only the capability guard:
  // passing declared numbers must never certify physical production validity.
  auto declared=s;declared.domain.authority=RC::DomainAuthority::ExplicitDeclaredLimits;
  declared.domain.maximumDeltaBOverB=100;
  declared.domain.maximumKThermalIonGyroradius=100;
  declared.domain.thermalIonGyroradiusM=1;
  declared.domain.maximumAlfvenSpeedOverC=1;
  declared.domain.ionInertialLengthM=1;
  declared.domain.maximumKIonInertialLength=100;
  declared.domain.ionCyclotronFrequencyRadPerS=1;
  declared.domain.maximumWaveFrequencyOverIonCyclotron=100;
  declared.domain.plasmaCompositionIdentity="manufactured-capability-check";
  declared.domain.thermalStateIdentity="manufactured-capability-check";
  declared.domain.provenance="intentionally overbroad numbers; not scientific qualification";
  const auto checked=RC::BuildGenerator(declared,orbit.orbit,16);
  pass=pass&&checked.status.ok()&&checked.generator.withinDeclaredLimits&&
      !checked.generator.qualifiedForProduction;
  declared.domain.maximumAlfvenSpeedOverC=0;
  pass=pass&&!RC::ValidateDomain(declared).ok();
  pass=pass&&RC::RequireImplementedProductionClosure(false,false).ok()&&
      !RC::RequireImplementedProductionClosure(true,false).ok()&&
      !RC::RequireImplementedProductionClosure(false,true).ok();
  return Finish(pass,"four explicitly polarized bin-energy spectra determine signed-charge resonance; unknown initialization/domain are rejected");
}

Result InvariantCoordinates() {
  auto s=Spectrum();
  // Test kinematics outside a production MHD qualification: v<v_A makes
  // plasma-pitch invariant coordinates fold, while wave-pitch stays unique.
  s.alfvenSpeedMPerS=0.6*C;const auto particle=Particle(0.2*C,0.7);
  const auto prepared=RC::PrepareOrbit(s,RC::Branch::AlongB,particle);
  if(!prepared.status.ok())return Failure(prepared.status,"slow-particle orbit preparation");
  bool pass=prepared.status.ok();
  const auto roundTrip=RC::ToPlasmaFrame(prepared.orbit,prepared.orbit.initialWavePitch);
  pass=pass&&roundTrip.status.ok()&&
      Relative(roundTrip.particle.parallelMomentumKgMPerS,particle.parallelMomentumKgMPerS)&&
      Relative(roundTrip.particle.perpendicularMomentumKgMPerS,particle.perpendicularMomentumKgMPerS);
  double maximumInvariantError=0;
  for(double pitch:{-1.0,-0.8,-0.1,0.0,0.3,0.9,1.0}) {
    const auto state=RC::ToPlasmaFrame(prepared.orbit,pitch);
    const double invariant=state.totalEnergyJ-s.alfvenSpeedMPerS*state.particle.parallelMomentumKgMPerS;
    const double error=std::fabs(invariant/prepared.orbit.invariantJ-1);
    maximumInvariantError=std::max(maximumInvariantError,error);
    pass=pass&&state.status.ok()&&error<1e-13;
  }
  // Independent determinant of (p,mu_plasma) with respect to (C,mu_wave).
  const double pitch=0.2,hPitch=1e-5;
  const double hInvariant=1e-5*(prepared.orbit.energyInWaveFrameJ-particle.massKg*C*C)/prepared.orbit.waveGamma;
  auto less=prepared.orbit,more=prepared.orbit;
  less.invariantJ-=hInvariant;more.invariantJ+=hInvariant;
  for(auto* orbit:{&less,&more}) {
    orbit->energyInWaveFrameJ=orbit->waveGamma*orbit->invariantJ;
    orbit->momentumInWaveFrameKgMPerS=std::sqrt(
        (orbit->energyInWaveFrameJ/C-particle.massKg*C)*
        (orbit->energyInWaveFrameJ/C+particle.massKg*C));
  }
  auto coordinates=[](const RC::ParticleResult& state) {
    const double p=std::hypot(state.particle.parallelMomentumKgMPerS,
        state.particle.perpendicularMomentumKgMPerS);
    return std::array<double,2>{{p,state.particle.parallelMomentumKgMPerS/p}};
  };
  const auto lowC=coordinates(RC::ToPlasmaFrame(less,pitch));
  const auto highC=coordinates(RC::ToPlasmaFrame(more,pitch));
  const auto lowQ=coordinates(RC::ToPlasmaFrame(prepared.orbit,pitch-hPitch));
  const auto highQ=coordinates(RC::ToPlasmaFrame(prepared.orbit,pitch+hPitch));
  const auto center=RC::ToPlasmaFrame(prepared.orbit,pitch);
  const auto centerCoordinates=coordinates(center);
  const double determinant=std::fabs(
      (highC[0]-lowC[0])*(highQ[1]-lowQ[1])-
      (highQ[0]-lowQ[0])*(highC[1]-lowC[1]))/(4*hInvariant*hPitch);
  const double numericalJacobian=centerCoordinates[0]*centerCoordinates[0]*determinant;
  const double exactJacobian=prepared.orbit.waveGamma*center.totalEnergyJ*
      prepared.orbit.momentumInWaveFrameKgMPerS/(C*C);
  const double jacobianError=std::fabs(numericalJacobian/exactJacobian-1);
  pass=pass&&jacobianError<1e-6;
  auto result=Finish(pass,"wave-pitch invariant coordinates remain unique below v_A and recover the independently differentiated relativistic phase-space Jacobian");
  result.metrics.push_back({"maximum_invariant_relative_error",maximumInvariantError,1e-13,"<","dimensionless"});
  result.metrics.push_back({"Jacobian_relative_error",jacobianError,1e-6,"<","dimensionless"});
  return result;
}

class CountingRandom final:public Transport::RandomStream {
 public:int draws=0;
  double UniformOpen01()override{++draws;return 0.5;}
  double Normal01()override{return 0;}
};
Result GeneratorContract() {
  const auto s=Spectrum();const auto orbit=RC::PrepareOrbit(s,RC::Branch::AlongB,Particle());
  if(!orbit.status.ok())return Failure(orbit.status,"generator orbit preparation");
  const auto built=RC::BuildGenerator(s,orbit.orbit,32);
  if(!built.status.ok())return Failure(built.status,"generator preparation");
  bool pass=built.status.ok();const auto& g=built.generator;
  const double total=std::accumulate(g.measure.begin(),g.measure.end(),0.0);
  std::vector<double> equilibrium(g.measure.size());
  for(std::size_t j=0;j<equilibrium.size();++j)equilibrium[j]=g.measure[j]/total;
  std::vector<double> derivative;
  pass=pass&&RC::ProbabilityDerivative(g,equilibrium,&derivative).ok();
  double residual=0;
  for(double value:derivative)residual=std::max(residual,std::fabs(value));
  pass=pass&&residual<1e-12;
  for(std::size_t j=0;j+1<g.measure.size();++j)
    pass=pass&&g.rightRatePerS[j]>=0&&g.leftRatePerS[j+1]>=0&&
        Relative(g.measure[j]*g.rightRatePerS[j],g.measure[j+1]*g.leftRatePerS[j+1]);
  CountingRandom random;
  const auto gap=RC::SampleJump(g,RC::CellAt(g,0),random);
  pass=pass&&gap.status.ok()&&!gap.active&&std::isinf(gap.waitingTimeS)&&random.draws==0;
  const auto crossed=RC::ProposeJump(g,-0.01,0.01,1);
  pass=pass&&!crossed.status.ok();
  const auto result=Finish(pass,"positive conservative rates satisfy exact phase-space detailed balance and finite-band gaps consume no collision draws");
  return result;
}

Result ExchangeProposal() {
  const auto s=Spectrum();const auto orbit=RC::PrepareOrbit(s,RC::Branch::AgainstB,Particle());
  if(!orbit.status.ok())return Failure(orbit.status,"exchange orbit preparation");
  const auto built=RC::BuildGenerator(s,orbit.orbit,64);const auto& g=built.generator;
  if(!built.status.ok())return Failure(built.status,"exchange generator preparation");
  const auto cell=RC::CellAt(g,0.45);CountingRandom random;
  const auto jump=RC::SampleJumpAtEvent(g,cell,random);
  const auto proposal=RC::ProposeJump(g,0.45,jump.targetWavePitch,7);
  double magnetic=0,particle=0;
  for(const auto& bin:proposal.exchange){magnetic+=bin.magneticEnergyChangeJ;particle+=bin.particleEnergyChangeJ;}
  bool pass=built.status.ok()&&jump.status.ok()&&jump.active&&random.draws==2&&
      proposal.status.ok()&&proposal.spectrumIdentity==s.identity&&proposal.generation==s.generation&&
      Relative(particle,proposal.particleEnergyChangeJ)&&
      Relative(2*magnetic,-proposal.particleEnergyChangeJ)&&
      Relative(proposal.particleEnergyChangeJ,
          orbit.orbit.waveSpeedMPerS*proposal.parallelMomentumChangeKgMPerS)&&
      proposal.exchange.size()>0;
  pass=pass&&RC::ValidateProposal(s,RC::Branch::AgainstB,proposal,7).ok();
  // These are wire-contract negative controls. The conservation identities
  // above and the independent Lorentz/Kubo tests remain the scientific
  // references; sharing validation/generation algebra is not a new oracle.
  auto tampered=proposal;tampered.after.perpendicularMomentumKgMPerS*=1.01;
  pass=pass&&!RC::ValidateProposal(s,RC::Branch::AgainstB,tampered,7).ok();
  tampered=proposal;tampered.exchange[0].magneticEnergyChangeJ*=1.01;
  pass=pass&&!RC::ValidateProposal(s,RC::Branch::AgainstB,tampered,7).ok();
  tampered=proposal;tampered.exchange[0].bin=s.kEdgesPerM.size();
  pass=pass&&!RC::ValidateProposal(s,RC::Branch::AgainstB,tampered,7).ok();
  tampered=proposal;tampered.exchange[0].helicity=tampered.exchange[0].helicity==RC::Helicity::Left?RC::Helicity::Right:RC::Helicity::Left;
  pass=pass&&!RC::ValidateProposal(s,RC::Branch::AgainstB,tampered,7).ok()&&
      !RC::ValidateProposal(s,RC::Branch::AgainstB,proposal,14).ok()&&
      !RC::ValidateProposal(s,RC::Branch::AlongB,proposal,7).ok();
  tampered=proposal;tampered.epochS+=1;
  pass=pass&&!RC::ValidateProposal(s,RC::Branch::AgainstB,tampered,7).ok()&&random.draws==2;
  return Finish(pass,"private neighboring-cell proposals obey exact wave-frame work and total-wave=2*magnetic energy closure without premature state publication");
}

double InitialMass(double left,double right,double alpha,double a,double b) {
  if(right<=left)return 0;
  const double frequency=Pi/(b-a);
  auto oscillatoryIntegral=[&](double x) {
    const double phase=frequency*(x-a);
    return std::sin(phase)/frequency+alpha*(x*std::sin(phase)/frequency+
        std::cos(phase)/(frequency*frequency));
  };
  return (right-left)*(1+0.5*alpha*(left+right))+
      0.25*(oscillatoryIntegral(right)-oscillatoryIntegral(left));
}

// Independent PDE reference: centered finite differences of (A f')'/J on
// nodal values, trapezoidal phase-space weights and Crank-Nicolson in time.
// It does not call BuildGenerator, its resistance integrals or event sampler.
std::vector<double> IndependentPde(double a,double b,double alpha,
    double coefficient,double finalTime,std::size_t intervals) {
  const std::size_t count=intervals+1;
  const double h=(b-a)/intervals,dt=finalTime/1000;
  std::vector<double> lower(count,0),upper(count,0),diagonal(count),upperRatio(count);
  std::vector<double> solution(count),rhs(count),forward(count);
  for(std::size_t i=0;i<count;++i) {
    const double x=a+i*h;
    const double weight=(1+alpha*x)*h*(i==0||i+1==count?0.5:1);
    if(i>0) {
      const double face=x-0.5*h;
      lower[i]=coefficient*(1-face*face)/face/(h*weight);
    }
    if(i+1<count) {
      const double face=x+0.5*h;
      upper[i]=coefficient*(1-face*face)/face/(h*weight);
    }
    solution[i]=1+0.25*std::cos(Pi*(x-a)/(b-a));
    const double raw=1+0.5*dt*(lower[i]+upper[i]);
    diagonal[i]=raw+(i>0?0.5*dt*lower[i]*upperRatio[i-1]:0);
    upperRatio[i]=-0.5*dt*upper[i]/diagonal[i];
  }
  for(int step=0;step<1000;++step) {
    for(std::size_t i=0;i<count;++i) {
      rhs[i]=(1-0.5*dt*(lower[i]+upper[i]))*solution[i];
      if(i>0)rhs[i]+=0.5*dt*lower[i]*solution[i-1];
      if(i+1<count)rhs[i]+=0.5*dt*upper[i]*solution[i+1];
      forward[i]=(rhs[i]+(i>0?0.5*dt*lower[i]*forward[i-1]:0))/diagonal[i];
    }
    solution.back()=forward.back();
    for(std::size_t i=count-1;i-->0;)
      solution[i]=forward[i]-upperRatio[i]*solution[i+1];
  }
  return solution;
}

// Exact semidiscrete CTMC evolution by positive uniformization. The Poisson
// series has mean <=2 per chunk and 64 terms; its omitted tail is <1e-70.
bool UniformizedEvolution(const RC::Generator& g,double finalTime,
                          std::vector<double>* probability) {
  double lambda=0;
  for(std::size_t i=0;i<g.centers.size();++i)
    lambda=std::max(lambda,g.leftRatePerS[i]+g.rightRatePerS[i]);
  if(lambda==0)return true;
  const auto chunks=static_cast<std::size_t>(std::ceil(lambda*finalTime/2));
  if(chunks==0||chunks>1000000)return false;
  const double mean=lambda*finalTime/chunks;
  std::vector<double> term,updated(probability->size()),sum(probability->size());
  for(std::size_t step=0;step<chunks;++step) {
    term=*probability;double poisson=std::exp(-mean);
    for(std::size_t i=0;i<term.size();++i)sum[i]=poisson*term[i];
    for(int order=1;order<64;++order) {
      for(std::size_t i=0;i<term.size();++i) {
        updated[i]=term[i]*(1-(g.leftRatePerS[i]+g.rightRatePerS[i])/lambda);
        if(i>0)updated[i]+=term[i-1]*g.rightRatePerS[i-1]/lambda;
        if(i+1<term.size())updated[i]+=term[i+1]*g.leftRatePerS[i+1]/lambda;
      }
      term.swap(updated);poisson*=mean/order;
      for(std::size_t i=0;i<term.size();++i)sum[i]+=poisson*term[i];
    }
    probability->swap(sum);
  }
  return true;
}

Result IndependentPdeConvergence() {
  auto spectrum=Spectrum();const auto prepared=RC::PrepareOrbit(spectrum,RC::Branch::AlongB,Particle());
  if(!prepared.status.ok())return Failure(prepared.status,"PDE orbit preparation");
  const auto& orbit=prepared.orbit;
  const double scale=std::fabs(orbit.chargeC)*spectrum.meanMagneticFieldT*
      orbit.waveGamma/orbit.momentumInWaveFrameKgMPerS;
  spectrum.kEdgesPerM={scale/0.8,scale/0.2};
  const double power=1e-12;
  for(auto& helicity:spectrum.magneticEnergyJ)helicity.assign(1,0);
  spectrum.magneticEnergyJ[1][0]=power*spectrum.volumeM3*
      (spectrum.kEdgesPerM[1]-spectrum.kEdgesPerM[0])/
      (2*spectrum.vacuumPermeabilityHPerM);
  const double a=scale/spectrum.kEdgesPerM.back(),b=scale/spectrum.kEdgesPerM.front();
  const double alpha=orbit.waveSpeedMPerS*orbit.momentumInWaveFrameKgMPerS/orbit.energyInWaveFrameJ;
  const double coefficient=0.5*Pi*orbit.chargeC*orbit.chargeC*power*C*C/
      (orbit.waveGamma*orbit.waveGamma*orbit.momentumInWaveFrameKgMPerS*orbit.energyInWaveFrameJ);
  const double finalTime=0.03/coefficient;
  const std::size_t referenceIntervals=1024;
  const auto reference=IndependentPde(a,b,alpha,coefficient,finalTime,referenceIntervals);
  const double normalization=InitialMass(a,b,alpha,a,b);
  std::vector<double> errors;
  double maximumMassError=0,minimumProbability=1;
  bool pass=true;
  for(std::size_t count:{40u,80u,160u}) {
    const auto built=RC::BuildGenerator(spectrum,orbit,count);
    if(!built.status.ok())return Failure(built.status,"PDE comparison generator preparation");
    const auto& g=built.generator;std::vector<double> probability(g.centers.size(),0);
    for(std::size_t i=0;i<probability.size();++i)
      probability[i]=InitialMass(std::max(a,g.edges[i]),
          std::min(b,g.edges[i+1]),alpha,a,b)/normalization;
    pass=pass&&UniformizedEvolution(g,finalTime,&probability);
    double error=0;
    for(std::size_t i=0;i<probability.size();++i) {
      minimumProbability=std::min(minimumProbability,probability[i]);
      if(g.centers[i]<=a||g.centers[i]>=b) {
        pass=pass&&probability[i]==0;continue;
      }
      const double x=(g.centers[i]-a)/(b-a)*referenceIntervals;
      const std::size_t node=std::min(referenceIntervals-1,static_cast<std::size_t>(x));
      const double referenceValue=reference[node]+(x-node)*(reference[node+1]-reference[node]);
      const double dimensionlessMass=(g.edges[i+1]-g.edges[i])*
          (1+0.5*alpha*(g.edges[i+1]+g.edges[i]));
      const double value=probability[i]*normalization/dimensionlessMass;
      error+=dimensionlessMass*std::fabs(value-referenceValue)/normalization;
    }
    errors.push_back(error);
    maximumMassError=std::max(maximumMassError,std::fabs(
        std::accumulate(probability.begin(),probability.end(),0.0)-1));
  }
  const double order=std::log(errors[1]/errors[2])/std::log(2.0);
  pass=pass&&errors[0]>errors[1]&&errors[1]>errors[2]&&order>1.5&&
      errors[2]<1e-3&&maximumMassError<1e-11&&minimumProbability>=0;
  auto result=Finish(pass,"positive finite-volume collision evolution converges to an independent centered-difference Crank-Nicolson FP solution and preserves the physical resonance gap");
  result.metrics.push_back({"fine_weighted_L1_error",errors[2],1e-3,"<","probability"});
  result.metrics.push_back({"observed_refinement_order",order,1.5,">","dimensionless"});
  result.metrics.push_back({"maximum_probability_mass_error",maximumMassError,1e-11,"<","probability"});
  result.metrics.push_back({"minimum_probability",minimumProbability,0,">=","probability"});
  return result;
}

struct ClockedState {
  RC::Particle particle=Particle();
  double opticalDepth=std::numeric_limits<double>::quiet_NaN();
  double representedEnergyChangeJ=0;
  std::size_t events=0,plusEvents=0,minusEvents=0;
};
Transport::Status AdvanceCombinedClock(const RC::Spectrum& spectrum,double duration,
    ClockedState* state,Transport::RandomStream& random,std::size_t stopAfterEvents=0) {
  double elapsed=0;
  while(elapsed<duration) {
    RC::Generator generators[2];std::size_t cells[2];double rates[2];
    for(int branch=0;branch<2;++branch) {
      const auto orbit=RC::PrepareOrbit(spectrum,branch==0?RC::Branch::AlongB:
          RC::Branch::AgainstB,state->particle);
      if(!orbit.status.ok())return orbit.status;
      const auto built=RC::BuildGenerator(spectrum,orbit.orbit,24);
      if(!built.status.ok())return built.status;
      generators[branch]=built.generator;
      cells[branch]=RC::CellAt(generators[branch],orbit.orbit.initialWavePitch);
      const auto rate=RC::TotalRate(generators[branch],cells[branch]);
      if(!rate.status.ok())return rate.status;
      rates[branch]=rate.value;
    }
    const double total=rates[0]+rates[1];
    if(total==0)return Transport::Status::Ok();
    if(std::isnan(state->opticalDepth))state->opticalDepth=-std::log(random.UniformOpen01());
    const double remaining=duration-elapsed;
    const double untilEvent=state->opticalDepth/total;
    if(untilEvent>remaining) {
      state->opticalDepth-=total*remaining;return Transport::Status::Ok();
    }
    elapsed+=untilEvent;
    const int branch=random.UniformOpen01()<rates[0]/total?0:1;
    const auto& generator=generators[branch];
    const auto jump=RC::SampleJumpAtEvent(generator,cells[branch],random);
    if(!jump.status.ok())return jump.status;
    if(!jump.active)return Transport::Status::Error(Transport::StatusCode::InvalidCoefficient,
        "combined event selected a zero-rate branch");
    const auto proposal=RC::ProposeJump(generator,generator.orbit.initialWavePitch,
        jump.targetWavePitch,3.0);
    if(!proposal.status.ok())return proposal.status;
    state->particle=proposal.after;
    state->representedEnergyChangeJ+=proposal.particleEnergyChangeJ;
    ++state->events;if(branch==0)++state->plusEvents;else++state->minusEvents;
    state->opticalDepth=-std::log(random.UniformOpen01());
    if(stopAfterEvents&&state->events>=stopAfterEvents)return Transport::Status::Ok();
    if(state->events>100000)return Transport::Status::Error(Transport::StatusCode::StepUnderflow,
        "bounded combined-clock fixture exceeded event limit");
  }
  return Transport::Status::Ok();
}
class ForcedBranchCoverage final:public Transport::RandomStream {
 public:std::size_t draws=0;
  double UniformOpen01()override {
    const auto index=draws++;
    if(index==0)return std::exp(-1e-6);
    switch((index-1)%4) {
      case 0:return ((index-1)/4)%2==0?0.001:0.999;
      case 1:return 0.75;
      case 2:return 0.5;
      default:return std::exp(-1e-6);
    }
  }
  double Normal01()override{return 0;}
};
Result CombinedHazardPartition() {
  const auto spectrum=Spectrum();ClockedState one,partitioned,restarted;
  const auto seed=TestSeed(83);
  Transport::KeyedRandomStream randomA(seed,7,5005,0),randomB(seed,7,5005,0),
      randomC(seed,7,5005,0);
  bool pass=AdvanceCombinedClock(spectrum,25,&one,randomA).ok()&&
      AdvanceCombinedClock(spectrum,3,&partitioned,randomB).ok();
  // Copy physical state and restore the complete saved stream at an actual
  // outer boundary; no scheduled future branch/pitch mark has been drawn.
  restarted=partitioned;const auto checkpoint=randomB.SaveState();
  pass=pass&&randomC.RestoreState(checkpoint).ok()&&
      AdvanceCombinedClock(spectrum,7,&partitioned,randomB).ok()&&
      AdvanceCombinedClock(spectrum,15,&partitioned,randomB).ok()&&
      AdvanceCombinedClock(spectrum,22,&restarted,randomC).ok();
  const auto a=randomA.SaveState(),b=randomB.SaveState(),c=randomC.SaveState();
  const auto untouched=Transport::KeyedRandomStream(seed,7,5005,0).SaveState();
  // Coverage uses a separate prescribed-variate prefix, never conditions the
  // Monte Carlo realization on an event-count outcome. Both physical branch
  // paths fire, with exactly one initial clock + four draws per event.
  ClockedState coverageState;ForcedBranchCoverage coverageRandom;
  const bool forcedCoverage=AdvanceCombinedClock(spectrum,1,&coverageState,
      coverageRandom,2).ok()&&coverageState.events==2&&
      coverageState.plusEvents==1&&coverageState.minusEvents==1&&coverageRandom.draws==9;
  // Event counts are observables, not universal lower bounds: a valid seed
  // may reach a finite-band collision gap early. Compare every accepted
  // event/history without conditioning the test on a favorable realization.
  pass=pass&&forcedCoverage&&a.integerState!=untouched.integerState&&
      one.events==partitioned.events&&one.events==restarted.events&&
      one.plusEvents==partitioned.plusEvents&&one.minusEvents==partitioned.minusEvents&&
      a.integerState==b.integerState&&a.integerState==c.integerState&&
      Relative(one.particle.parallelMomentumKgMPerS,
          partitioned.particle.parallelMomentumKgMPerS)&&
      Relative(one.particle.perpendicularMomentumKgMPerS,
          partitioned.particle.perpendicularMomentumKgMPerS)&&
      Relative(one.particle.parallelMomentumKgMPerS,
          restarted.particle.parallelMomentumKgMPerS)&&
      Relative(one.opticalDepth,partitioned.opticalDepth)&&
      Relative(one.opticalDepth,restarted.opticalDepth)&&
      Relative(one.representedEnergyChangeJ,partitioned.representedEnergyChangeJ);
  auto result=Finish(pass,"combined integrated hazard and at-event branch/pitch marks preserve physical state, exchange and complete RNG history across outer cuts and restored state");
  result.hasSeed=true;result.seed=seed;
  result.metrics.push_back({"physical_events",static_cast<double>(one.events),0,">=","count"});
  result.metrics.push_back({"plus_events",static_cast<double>(one.plusEvents),0,">=","count"});
  result.metrics.push_back({"minus_events",static_cast<double>(one.minusEvents),0,">=","count"});
  result.metrics.push_back({"generator_state_equal",a.integerState==b.integerState&&a.integerState==c.integerState?1.0:0.0,1,"==","boolean"});
  result.metrics.push_back({"forced_both_branch_paths",forcedCoverage?1.0:0.0,1,"==","boolean"});
  return result;
}

Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;
  descriptor.group="resonant-collision";descriptor.description=name;
  descriptor.initialization=InitializationLevel::None;
  descriptor.supportedBuildModes="native registry and identical dependency-light sanitizer callback";
  descriptor.runtime=RuntimeClass::Routine;
  descriptor.seedPolicy="manufactured fields and explicit deterministic streams";
  descriptor.seedCapability=std::string(id)=="RC06"?SeedCapability::Context:SeedCapability::Deterministic;
  descriptor.repeatSupported=true;
  descriptor.stateIsolation="private four-helicity spectra; validation-only domain; no production mutation";
  descriptor.callback=callback;return descriptor;
}
} // namespace
std::vector<Descriptor> ResonantCollisionDescriptors() {
  return {Make("RC01","Four-spectrum resonance and domain capabilities",FourSpectrumContract),
          Make("RC02","Relativistic invariant-coordinate Jacobian",InvariantCoordinates),
          Make("RC03","Conservative positive finite-band collision generator",GeneratorContract),
          Make("RC04","Exact private elastic exchange proposal",ExchangeProposal),
          Make("RC05","Independent polarized FP convergence",IndependentPdeConvergence),
          Make("RC06","Combined resonant-hazard partition continuity",CombinedHazardPartition)};
}
} }
