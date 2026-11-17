#include "sep_real_wave_normalization_validation.h"
#include "sep_resonant_collision.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_physics.h)

#include <array>
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace RC=Transport::ResonantCollision;
const double Pi=3.141592653589793238462643383279502884;
const double C=299792458.0;
using State=std::array<double,4>; // z,px,py,pz; SI, q=m=B0=1.
State Derivative(const State& x,double t,double k,double u,double amplitude,double phase,int handedness) {
  const double gammaMass=std::sqrt(1+(x[1]*x[1]+x[2]*x[2]+x[3]*x[3])/(C*C));
  const double vx=x[1]/gammaMass,vy=x[2]/gammaMass,vz=x[3]/gammaMass;
  const double angle=k*(x[0]-u*t)+phase;
  const double bx=amplitude*std::cos(angle),by=-handedness*amplitude*std::sin(angle);
  // Real circular field: deltaB=(a cos,-a sin,0), <deltaB^2>=a^2.
  // Ideal nondispersive E=-u*z_hat cross deltaB. These are actual Lorentz
  // forces; no QLT coefficient or resonance helper enters this integrator.
  return {{vz,u*by+vy-vz*by,-u*bx+vz*bx-vx,vx*by-vy*bx}};
}
void Rk4(State* state,double t,double dt,double k,double u,double a,double phase,int handedness) {
  const auto first=Derivative(*state,t,k,u,a,phase,handedness);State stage;
  for(int i=0;i<4;++i)stage[i]=(*state)[i]+0.5*dt*first[i];
  const auto second=Derivative(stage,t+0.5*dt,k,u,a,phase,handedness);
  for(int i=0;i<4;++i)stage[i]=(*state)[i]+0.5*dt*second[i];
  const auto third=Derivative(stage,t+0.5*dt,k,u,a,phase,handedness);
  for(int i=0;i<4;++i)stage[i]=(*state)[i]+dt*third[i];
  const auto fourth=Derivative(stage,t+dt,k,u,a,phase,handedness);
  for(int i=0;i<4;++i)(*state)[i]+=dt*(first[i]+2*second[i]+2*third[i]+fourth[i])/6;
}
struct Response {
  std::array<double,3> lorentzD{{0,0,0}},kuboD{{0,0,0}};
  double asymptoticD=0;
};
Response MeasureResponse(double u,double amplitude,double dt,int modes,
                         int handedness=1,double powerLawIndex=0) {
  const int phases=8;const double mu=0.5,v=1,omega=1;
  const double relative=v*mu-u,kCenter=omega/relative,halfWidth=1;
  const double dk=2*halfWidth/modes;
  const int endSteps=static_cast<int>(std::llround(80/dt));
  const int sampleSteps[3]={static_cast<int>(std::llround(20/dt)),
      static_cast<int>(std::llround(40/dt)),endSteps};
  const double times[3]={20,40,80};
  std::array<long double,3> sum{{0,0,0}},sum2{{0,0,0}},kubo{{0,0,0}};
  for(int mode=0;mode<modes;++mode) {
    const double k=kCenter-halfWidth+(mode+0.5)*dk;
    const double mismatch=k*relative-handedness*omega;
    // Uniform random k and this physical mode amplitude realize a normalized
    // one-sided k^-q variance spectrum. q=0 denotes the flat-band reference.
    const double modeVariance=powerLawIndex==0?amplitude*amplitude:
        amplitude*amplitude*2*halfWidth*(powerLawIndex-1)*std::pow(k,-powerLawIndex)/
        (std::pow(kCenter-halfWidth,1-powerLawIndex)-std::pow(kCenter+halfWidth,1-powerLawIndex));
    const double modeAmplitude=std::sqrt(modeVariance);
    const double angularForce=omega*std::sqrt(1-mu*mu)*(1-u*mu/v);
    for(int sample=0;sample<3;++sample) {
      const double argument=0.5*mismatch*times[sample];
      const double sinc=argument==0?1:std::sin(argument)/argument;
      // Direct finite-time correlation integral of the explicit real wave.
      kubo[sample]+=angularForce*angularForce*modeVariance*
          times[sample]*0.25*sinc*sinc/modes;
    }
    for(int phaseIndex=0;phaseIndex<phases;++phaseIndex) {
      const double phase=2*Pi*(phaseIndex+0.5)/phases;
      State state{{0,std::sqrt(1-mu*mu),0,mu}};int sample=0;
      for(int step=0;step<endSteps;++step) {
        Rk4(&state,step*dt,dt,k,u,modeAmplitude,phase,handedness);
        if(step+1==sampleSteps[sample]) {
          const double momentum=std::sqrt(state[1]*state[1]+state[2]*state[2]+state[3]*state[3]);
          const double change=state[3]/momentum-mu;
          sum[sample]+=change;sum2[sample]+=static_cast<long double>(change)*change;
          ++sample;if(sample==3)break;
        }
      }
    }
  }
  Response result;
  for(int sample=0;sample<3;++sample) {
    const long double mean=sum[sample]/(modes*phases);
    result.lorentzD[sample]=static_cast<double>(sum2[sample]/(modes*phases)-mean*mean)/(2*times[sample]);
    result.kuboD[sample]=static_cast<double>(kubo[sample]);
  }
  // P_h is ONE-SIDED REAL magnetic variance density, not total-wave energy:
  // integral_{kCenter-1}^{kCenter+1} P_h dk = <deltaB_h^2> = amplitude^2.
  const double power=powerLawIndex==0?amplitude*amplitude/(2*halfWidth):
      amplitude*amplitude*(powerLawIndex-1)*std::pow(kCenter,-powerLawIndex)/
      (std::pow(kCenter-halfWidth,1-powerLawIndex)-std::pow(kCenter+halfWidth,1-powerLawIndex));
  result.asymptoticD=handedness<0?0:0.5*Pi*omega*omega*(1-mu*mu)*
      std::pow(1-u*mu/v,2)*power/std::fabs(relative);
  return result;
}
Result LegacyBalancedNormalization() {
  namespace CP=Transport::CoefficientPhysics;
  const double amplitude=1e-6,variance=amplitude*amplitude,mu=0.5,v=1,q=5.0/3.0;
  // Equal independent real helicities are an explicitly manufactured field,
  // not an assumed initializer for the new polarized production state.
  const auto resonant=MeasureResponse(0,amplitude,0.01,128,1,q);
  const auto opposite=MeasureResponse(0,amplitude,0.01,128,-1,q);
  double maxLorentzError=0;
  for(int sample=0;sample<3;++sample) {
    const double measured=0.5*(resonant.lorentzD[sample]+opposite.lorentzD[sample]);
    const double kubo=0.5*(resonant.kuboD[sample]+opposite.kuboD[sample]);
    maxLorentzError=std::max(maxLorentzError,std::fabs(measured/kubo-1));
  }
  const double physicalJokipii=0.5*resonant.asymptoticD;
  CP::LocalInputView input;input.source="manufactured-real-balanced-helicities";
  input.representation="magnetic-variance";input.generation=1;input.checksum=1;
  input.heliocentricRadiusM=1;input.magneticFieldT=1;input.deltaB2T2=variance;
  input.deltaBPlus2T2=0.5*variance;input.deltaBMinus2T2=0.5*variance;
  CP::SpeciesProperties species;species.modelSpecies=0;species.name="controlled-q=m=1";
  species.restMassKg=1;species.signedChargeC=1;
  CP::SpectrumParameters band;band.referenceRadiusM=1;band.kMinAtReferencePerM=1;
  band.kMaxAtReferencePerM=3;band.spectralIndex=q;
  const auto jokipii=CP::EvaluateJokipiiSlab(input,band,species,v,mu);
  CP::FlorinskiyParameters smooth;smooth.spectralIndex=q;smooth.parallelCorrelationLengthM=1;
  // Independent analytic integral of the one-sided bend-over spectrum:
  // integral_0^infinity (1+(k/k0)^2)^(-q/2) dk
  // =k0*sqrt(pi)*Gamma((q-1)/2)/(2*Gamma(q/2)).
  const double k0=std::sqrt(Pi)*std::tgamma(q/2)/std::tgamma((q-1)/2);
  const double integral=k0*std::sqrt(Pi)*std::tgamma((q-1)/2)/(2*std::tgamma(q/2));
  const double totalPower=variance/integral*std::pow(1+std::pow(2/k0,2),-q/2);
  const double physicalFlorinskiy=0.25*Pi*(1-mu*mu)*totalPower/(v*mu);
  const auto florinskiy=CP::EvaluateFlorinskiySlab(input,smooth,species,v,mu);
  const double jokipiiRatio=jokipii.dMuMuPerS/physicalJokipii;
  const double florinskiyRatio=florinskiy.dMuMuPerS/physicalFlorinskiy;
  const bool pass=jokipii.status.ok()&&florinskiy.status.ok()&&maxLorentzError<5e-3&&
      std::fabs(jokipiiRatio-1)<1e-12&&std::fabs(florinskiyRatio-1)<1e-12;
  Result result;result.status=pass?Status::Pass:Status::Fail;
  result.message="legacy total magnetic variance maps to half the resonant one-helicity power for an explicitly balanced real slab field";
  result.configuration.push_back("manufactured_helicity=equal-independent-real-L-R;integral_positive_k_Ptotal=deltaBtotal_squared");
  result.configuration.push_back("reference=full-Lorentz-balanced-ensemble-plus-independent-positive-k-spectrum-integrals;no-production-polarization-assumed");
  result.metrics.push_back({"balanced_Lorentz_to_finite_Kubo_error",maxLorentzError,5e-3,"<","dimensionless"});
  result.metrics.push_back({"Jokipii_to_real_balanced_variance_ratio",jokipiiRatio,1,"==","dimensionless"});
  result.metrics.push_back({"Florinskiy_static_to_real_balanced_variance_ratio",florinskiyRatio,1,"==","dimensionless"});
  return result;
}
double CoreCoefficient(double u,double amplitude) {
  const double mu=0.5,v=1,relative=v*mu-u,kCenter=1/relative;
  RC::Spectrum spectrum;spectrum.kEdgesPerM={kCenter-1,kCenter+1};
  for(auto& component:spectrum.magneticEnergyJ)component.assign(1,0);
  spectrum.volumeM3=1;spectrum.meanMagneticFieldT=1;spectrum.alfvenSpeedMPerS=u;
  spectrum.vacuumPermeabilityHPerM=1.256637062e-6;
  spectrum.magneticEnergyJ[1][0]=amplitude*amplitude/(2*spectrum.vacuumPermeabilityHPerM);
  spectrum.polarizationResolved=true;spectrum.generation=1;spectrum.identity="real-circular-wave-calibration";
  spectrum.domain.authority=RC::DomainAuthority::ValidationStudy;
  spectrum.domain.provenance="independent finite-time Lorentz/Kubo normalization";
  RC::Particle particle;particle.massKg=1;particle.chargeC=1;
  particle.parallelMomentumKgMPerS=mu;particle.perpendicularMomentumKgMPerS=std::sqrt(1-mu*mu);
  const auto orbit=RC::PrepareOrbit(spectrum,RC::Branch::AlongB,particle);
  if(!orbit.status.ok())return std::numeric_limits<double>::quiet_NaN();
  const auto resonance=RC::EvaluateResonance(spectrum,orbit.orbit,orbit.orbit.initialWavePitch);
  if(!resonance.status.ok()||!resonance.active)return std::numeric_limits<double>::quiet_NaN();
  return 0.5*(1-mu*mu)*std::pow(1-u*mu/v,2)*resonance.frequencyPerS;
}
Result RealCircularWave() {
  const double amplitude=1e-6;bool pass=true;
  double maxLorentzError=0,maxCoreError=0,referenceRefinementError=0;
  Result result;
  for(double u:{0.0,0.1}) {
    const auto measured=MeasureResponse(u,amplitude,0.02,128);
    const auto refined=MeasureResponse(u,amplitude,0.01,128);
    for(int sample=0;sample<3;++sample) {
      const double error=std::fabs(refined.lorentzD[sample]/refined.kuboD[sample]-1);
      maxLorentzError=std::max(maxLorentzError,error);
      referenceRefinementError=std::max(referenceRefinementError,
          std::fabs(measured.lorentzD[sample]/refined.lorentzD[sample]-1));
    }
    const double coreRatio=CoreCoefficient(u,amplitude)/refined.asymptoticD;
    maxCoreError=std::max(maxCoreError,std::fabs(coreRatio-1));
    pass=pass&&refined.kuboD[0]<refined.kuboD[1]&&refined.kuboD[1]<refined.kuboD[2]&&
        std::fabs(refined.kuboD[2]/refined.asymptoticD-1)<0.03;
    result.metrics.push_back({u==0?"magnetostatic_core_to_real_wave_ratio":"moving_wave_core_to_real_wave_ratio",
        coreRatio,1,"==","dimensionless"});
    result.metrics.push_back({u==0?"magnetostatic_finite_time_to_limit":"moving_wave_finite_time_to_limit",
        refined.kuboD[2]/refined.asymptoticD,0.97,">","dimensionless"});
  }
  pass=pass&&maxLorentzError<5e-3&&referenceRefinementError<1e-3&&maxCoreError<1e-12;
  result.status=pass?Status::Pass:Status::Fail;
  result.message="explicit real circular-wave Lorentz orbits recover finite-time Kubo diffusion and the physical one-helicity magnetic-variance normalization";
  result.configuration.push_back("P_h_positive_k_integral=real_deltaB_h_squared;D=variance/(2*time)");
  result.configuration.push_back("real_field=deltaBx=a*cos(kz-ukt+phi),deltaBy=-a*sin(kz-ukt+phi);a=1e-6");
  result.configuration.push_back("reference=full-Lorentz-RK4;128-k-quadrature;8-uniform-phases-exact-at-quadratic-Kubo-order;T=20,40,80;dt=0.02,0.01");
  result.metrics.push_back({"maximum_Lorentz_to_finite_Kubo_error",maxLorentzError,5e-3,"<","dimensionless"});
  result.metrics.push_back({"Lorentz_time_refinement_error",referenceRefinementError,1e-3,"<","dimensionless"});
  result.metrics.push_back({"maximum_core_normalization_error",maxCoreError,1e-12,"<","dimensionless"});
  return result;
}
} // namespace
std::vector<Descriptor> RealWaveNormalizationDescriptors() {
  Descriptor descriptor;descriptor.id="RWN01";descriptor.name="Real-wave Lorentz/Kubo normalization";
  descriptor.group="resonant-collision";descriptor.description=descriptor.name;
  descriptor.initialization=InitializationLevel::None;descriptor.runtime=RuntimeClass::Extended;
  descriptor.supportedBuildModes="native registry and identical independent Lorentz sanitizer callback";
  descriptor.seedPolicy="deterministic uniform phase quadrature";
  descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;
  descriptor.stateIsolation="private prescribed real waves and independent Lorentz trajectories";
  descriptor.callback=RealCircularWave;
  Descriptor balanced=descriptor;balanced.id="RWN02";
  balanced.name="Legacy balanced real-wave variance normalization";balanced.description=balanced.name;
  balanced.callback=LegacyBalancedNormalization;return {descriptor,balanced};
}
} }
