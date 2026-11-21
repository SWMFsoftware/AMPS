#include "sep_maxwell_juttner_volume_validation.h"
#include "sep_maxwell_juttner_volume.h"
#include "sep_maxwell_juttner_flux.h"
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace V=Sources::MaxwellJuttnerVolume;
constexpr double C=Sources::MaxwellJuttner::SpeedOfLightMPerS;
constexpr double K=Sources::MaxwellJuttner::BoltzmannJPerK;
constexpr double Mass=1.67262192369e-27;
struct Moments {
  double sum=0,squares=0;unsigned count=0;
  void Add(double x){sum+=x;squares+=x*x;++count;}
  double Mean()const{return sum/count;}
  double Error()const{return std::sqrt(std::max(0.0,squares/count-Mean()*Mean())/count);}
};
V::Input Input(double theta,double beta=0) {
  V::Input input;input.restMassKg=Mass;input.properTemperatureK=theta*Mass*C*C/K;
  input.bulkVelocityMPerS[0]=beta*C;return input;
}
double MeanRestKineticGamma(double theta) {
  // Independent Simpson quadrature of the equilibrium momentum-volume law.
  // x=q/theta avoids cold Bessel underflow. At x=80 the exponential tail is
  // negligible compared with the declared Monte Carlo uncertainty; no source
  // sample is truncated by this reference integration interval.
  const unsigned intervals=40000;const double h=80.0/intervals;double norm=0,moment=0;
  for(unsigned i=0;i<=intervals;++i) {
    const double x=i*h,q=theta*x;
    const double density=(1+q)*std::sqrt(x*(2+q))*std::exp(-x);
    const double weight=i==0||i==intervals?1:(i%2?4:2);
    norm+=weight*density;moment+=weight*q*density;
  }
  return moment/norm;
}
Result Distribution() {
  bool pass=true;double largest=0;
  for(double theta:{1e-6,0.1,10.0}) {
    const auto input=Input(theta);Transport::KeyedRandomStream random(TestSeed(17),123,19,1);
    Moments energy,pressure;
    for(unsigned i=0;i<40000;++i) {
      V::Sample sample;const auto status=V::Draw(input,&random,&sample);
      if(!status.ok()){pass=false;break;}
      const double q=sample.restKineticGamma;energy.Add(q);
      // The ideal relativistic gas identity <p.v>/3=kT is independent of
      // the rejection envelope and of the quadrature energy reference.
      pressure.Add(q*(q+2)/(3*(1+q)*theta));
    }
    const double reference=MeanRestKineticGamma(theta);
    const double z=std::abs(energy.Mean()-reference)/energy.Error();largest=std::max(largest,z);
    pass=pass&&z<=6&&std::abs(pressure.Mean()-1)<=6*pressure.Error();
  }
  Result result;result.status=pass?Status::Pass:Status::Fail;
  result.message="volume Juttner sampling matches independent cold/relativistic quadrature and ideal-gas pressure";
  result.hasSeed=true;result.seed=TestSeed(17);result.metrics.push_back({"maximum_energy_z",largest,6,"<=","sigma"});
  result.metrics.push_back({"actual_seeded_observable",largest,0,"observed","dimensionless"});return result;
}
Result VolumeJacobian() {
  const double theta=0.2,beta=0.6;const auto input=Input(theta,beta);
  Transport::KeyedRandomStream random(TestSeed(41),123,21,1);Moments momentum,energy;
  bool pass=true;
  for(unsigned i=0;i<60000;++i) {
    V::Sample sample;const auto status=V::Draw(input,&random,&sample);
    if(!status.ok()){pass=false;break;}
    momentum.Add(sample.restMomentumKgMPerS[0]/(Mass*C));energy.Add(sample.simulationGamma);
  }
  // Lorentz transformation of an isotropic ideal-gas stress-energy tensor:
  // n_lab=gamma*n0, <p'_x>/(mc)=beta*theta, and
  // <gamma_lab>=gamma*(<gamma'>+beta^2*theta). The momentum test would fail
  // if the simulation-volume Jacobian were replaced by an isotropic rest draw.
  const double gamma=1.25,expectedEnergy=gamma*(1+MeanRestKineticGamma(theta)+beta*beta*theta);
  const double mz=std::abs(momentum.Mean()-beta*theta)/momentum.Error();
  const double ez=std::abs(energy.Mean()-expectedEnergy)/energy.Error();
  double actualGamma=0;pass=pass&&V::BulkGamma(input.bulkVelocityMPerS,&actualGamma).ok()&&
      std::abs(actualGamma-gamma)<2e-15&&mz<=6&&ez<=6;
  Result result;result.status=pass?Status::Pass:Status::Fail;
  result.message="simulation-volume thermal measure satisfies Lorentz number and stress-energy moments";
  result.hasSeed=true;result.seed=TestSeed(41);result.metrics.push_back({"momentum_z",mz,6,"<=","sigma"});
  result.metrics.push_back({"actual_seeded_observable",energy.Mean(),0,"observed","dimensionless"});
  result.metrics.push_back({"energy_z",ez,6,"<=","sigma"});return result;
}
Result FailureAtomicity() {
  bool pass=true;
  for(unsigned kind=0;kind<3;++kind) {
    auto input=Input(0.1);
    if(kind==0)input.properTemperatureK=-1;
    if(kind==1)input.restMassKg=0;
    if(kind==2)input.bulkVelocityMPerS[0]=C;
    Transport::KeyedRandomStream random(2,3,4,5);random.Normal01();auto expected=random;
    V::Sample output;output.restKineticGamma=19;output.restMomentumKgMPerS={{1,2,3}};
    pass=pass&&!V::Draw(input,&random,&output).ok()&&output.restKineticGamma==19&&output.restMomentumKgMPerS[1]==2;
    for(unsigned i=0;i<12;++i)pass=pass&&random.UniformOpen01()==expected.UniformOpen01()&&random.Normal01()==expected.Normal01();
  }
  auto cold=Input(0,0.6);Transport::KeyedRandomStream random(1,2,3,4);V::Sample output;
  pass=pass&&V::Draw(cold,&random,&output).ok()&&output.restKineticGamma==0&&output.simulationGamma==output.bulkGamma;
  Result result;result.status=pass?Status::Pass:Status::Fail;
  result.message="invalid thermal inputs preserve output and complete cached random state; cold limit is explicit";return result;
}
} // namespace
std::vector<Descriptor> MaxwellJuttnerVolumeDescriptors() {
  std::vector<Descriptor> result;const char* names[]={"Thermal volume equilibrium moments","Thermal simulation-volume Jacobian","Thermal volume failure atomicity"};
  const std::function<Result()> callbacks[]={Distribution,VolumeJacobian,FailureAtomicity};
  for(unsigned i=0;i<3;++i){Descriptor d;d.id="MJV0"+std::to_string(i+1);d.name=names[i];d.description=names[i];
    d.group="maxwell-juttner-volume";d.initialization=InitializationLevel::None;d.runtime=RuntimeClass::Routine;
    d.seedCapability=i<2?SeedCapability::Context:SeedCapability::Deterministic;
    d.repeatSupported=true;d.seedPolicy="TestSeed(default17/default41), private complete keyed streams and independent quadrature gates";
    d.supportedBuildModes="dependency-light/native identical production volume kernel; provider/native adapter require separate evidence";
    d.stateIsolation="private input, generator, moments and output";d.callback=callbacks[i];result.push_back(d);}
  return result;
}
} }
