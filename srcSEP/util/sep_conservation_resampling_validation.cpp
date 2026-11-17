#include "sep_conservation_resampling_validation.h"
#include "sep_conservative_resampling.h"
#include "sep_population_control.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <numeric>

namespace SEP { namespace Testing { namespace {
namespace C=PopulationControl::Conservative;
namespace K=PIC::ConservativeResampling;
const double c=299792458.0,mass=1.67262192369e-27,charge=1.602176634e-19;
Result Finish(bool ok,const char* message) {
  Result result;result.status=ok?Status::Pass:Status::Fail;result.message=message;
  result.metrics.push_back({"assertion_failures",ok?0.0:1.0,0,"<=","count"});
  result.configuration.push_back("evidence_level=conservative-weight-core;no-native-controller-or-MPI-claim");
  return result;
}
std::vector<C::Particle> Inputs() {
  const double beta[]={.05,.17,.31,.45,.58,.67,.78,.88,.93};
  const double mu[]={-.95,.72,-.43,.18,.84,-.77,.36,-.12,.59};
  const double s[]={.13,.73,.27,.91,.41,.05,.61,.84,.37};
  std::vector<C::Particle> particles;
  for(unsigned i=0;i<9;++i) {
    C::Particle p;p.identity=101+i;p.species=0;p.weight=1+.25*i;
    p.massKg=mass;p.chargeC=charge;p.localFraction=s[i];
    p.parallelMPerS=beta[i]*c*mu[i];
    p.normalMPerS=beta[i]*c*std::sqrt(1-mu[i]*mu[i]);
    particles.push_back(p);
  }
  return particles;
}
std::vector<long double> IndependentMoments(const std::vector<C::Particle>& particles,
                                           const std::vector<double>& weights) {
  std::vector<long double> moments(7,0.0L);
  for(std::size_t i=0;i<particles.size();++i) {
    const auto& p=particles[i];
    const long double vp=p.parallelMPerS,vn=p.normalMPerS;
    const long double speed=std::sqrt(vp*vp+vn*vn);
    const long double gamma=1/std::sqrt(1-speed*speed/(static_cast<long double>(c)*c));
    const long double mu=speed?vp/speed:0;
    // Independent energy expression; production uses the cancellation-free
    // gamma^2*m*v^2/(gamma+1), whereas this broad-speed fixture uses gamma-1.
    const long double rows[]={1.0L,p.chargeC,(gamma-1)*p.massKg*c*c,
        gamma*p.massKg*vp,mu,mu*mu,p.localFraction};
    for(unsigned row=0;row<7;++row)moments[row]+=weights[i]*rows[row];
  }
  return moments;
}
double MomentError(const std::vector<C::Particle>& particles,
                    const std::vector<double>& weights) {
  std::vector<double> original;
  for(const auto& p:particles)original.push_back(p.weight);
  const auto a=IndependentMoments(particles,original),b=IndependentMoments(particles,weights);
  double error=0;
  for(unsigned row=0;row<7;++row) {
    long double scale=0;
    for(std::size_t i=0;i<particles.size();++i) {
      std::vector<C::Particle> singleton(1,particles[i]);
      std::vector<double> weight(1,original[i]);
      scale+=std::fabs(IndependentMoments(singleton,weight)[row]);
    }
    error=std::max(error,scale?static_cast<double>(std::fabs(a[row]-b[row])/scale):
        (a[row]==b[row]?0.0:std::numeric_limits<double>::infinity()));
  }
  return error;
}
Result EndpointsConserve() {
  const auto particles=Inputs();const auto endpoints=C::Prepare(particles,1e-12);
  bool ok=endpoints.weights.status==K::Status::Success;
  double error=std::numeric_limits<double>::infinity();
  if(ok) {
    error=std::max(MomentError(particles,endpoints.weights.positive),
                   MomentError(particles,endpoints.weights.negative));
    ok=error<=1e-12;
    for(const auto* weights:{&endpoints.weights.positive,&endpoints.weights.negative}) {
      ok=ok&&std::count(weights->begin(),weights->end(),0.0)>0;
      for(double weight:*weights)ok=ok&&std::isfinite(weight)&&weight>=0;
    }
  }
  Result result=Finish(ok,"both positive endpoints conserve independently evaluated relativistic SI and pitch/spatial moments");
  result.metrics.push_back({"independent_maximum_relative_residual",error,1e-12,"<=","dimensionless"});
  return result;
}
Result AnalyticNullspace() {
  std::vector<double> weights(7,1.0);
  std::vector<std::vector<double>> rows(6,std::vector<double>(7));
  for(unsigned j=0;j<6;++j)for(unsigned i=0;i<7;++i)rows[j][i]=std::pow(static_cast<double>(i),j);
  const auto endpoints=K::Prepare(weights,rows,1e-12);
  const int finiteDifference[]={1,-6,15,-20,15,-6,1};
  bool ok=endpoints.status==K::Status::Success&&endpoints.rank==6;
  double error=std::numeric_limits<double>::infinity();
  if(ok) {
    error=std::fabs(endpoints.positiveProbability-4.0/7.0);
    for(unsigned i=0;i<7;++i) {
      error=std::max(error,std::fabs(endpoints.positive[i]-(1+finiteDifference[i]/20.0)));
      error=std::max(error,std::fabs(endpoints.negative[i]-(1-finiteDifference[i]/15.0)));
    }
    ok=error<1e-11;
  }
  Result result=Finish(ok,"pivoted QR matches the independent sixth finite-difference null vector and exact endpoint probabilities");
  result.metrics.push_back({"analytic_endpoint_absolute_error",error,1e-11,"<","dimensionless"});
  return result;
}
Result RankAndSparse() {
  auto particles=Inputs();particles.resize(2);
  bool ok=C::Prepare(particles,1e-12).weights.status==K::Status::NoNullDirection;
  particles=Inputs();
  for(auto& p:particles) {
    p.parallelMPerS=1e7;p.normalMPerS=2e7;p.localFraction=.4;
  }
  const auto identical=C::Prepare(particles,1e-12);
  ok=ok&&identical.weights.status==K::Status::Success&&identical.weights.rank==1;
  for(std::size_t i=0;i<particles.size();++i)particles[i].localFraction=.1+.08*i;
  const auto spatial=C::Prepare(particles,1e-12);
  ok=ok&&spatial.weights.status==K::Status::Success&&spatial.weights.rank==2;
  return Finish(ok,"full-rank sparse candidates explicitly refuse reduction; cold/duplicate constraints retain their actual numerical rank");
}
Result UnbiasedDistribution() {
  const auto particles=Inputs();const auto endpoints=C::Prepare(particles,1e-12);
  if(endpoints.weights.status!=K::Status::Success)return Finish(false,endpoints.weights.message.c_str());
  const double probability=endpoints.weights.positiveProbability;
  bool ok=true;double maximumError=0;
  for(std::size_t i=0;i<particles.size();++i) {
    const long double expected=probability*static_cast<long double>(endpoints.weights.positive[i])+
        (1.0L-probability)*endpoints.weights.negative[i];
    maximumError=std::max(maximumError,static_cast<double>(std::fabs(expected-particles[i].weight)/particles[i].weight));
  }
  ok=maximumError<=1e-12;
  std::vector<double> selected{71.0};
  ok=ok&&!C::Select(endpoints,0.0,&selected).ok()&&selected==std::vector<double>{71.0};
  ok=ok&&C::Select(endpoints,.5*probability,&selected).ok()&&selected==endpoints.weights.positive;
  ok=ok&&C::Select(endpoints,.5*(1+probability),&selected).ok()&&selected==endpoints.weights.negative;
  const unsigned trials=20000;unsigned positive=0;double firstDraw=0.0;
  Transport::KeyedRandomStream random(TestSeed(0),77,UINT64_C(0x4d45524745),0);
  for(unsigned i=0;i<trials;++i) {
    const double actual=random.UniformOpen01();if(i==0)firstDraw=actual;
    positive+=actual<probability;
  }
  const double sigma=std::sqrt(trials*probability*(1-probability));
  const double z=std::fabs(positive-trials*probability)/sigma;
  ok=ok&&z<=5.0;
  Result result=Finish(ok,"conditional expected weight of every phase point is unchanged; a fixed five-sigma branch-frequency gate passes");
  result.hasSeed=true;result.seed=TestSeed(0);
  result.metrics.push_back({"maximum_conditional_weight_error",maximumError,1e-12,"<=","dimensionless"});
  result.metrics.push_back({"branch_frequency_z_score",z,5.0,"<=","sigma"});
  result.metrics.push_back({"actual_seeded_observable",firstDraw,0.0,"observed","dimensionless"});
  return result;
}
Result SiNegativeControls() {
  PopulationControl::ParticleRecord p;p.stableId=1;p.species=0;p.weight=1;
  p.chargeC=charge;p.kineticEnergyJ=1e-14;p.parallelMomentumKgMPerS=1e-19;p.mu=.2;
  std::vector<PopulationControl::ParticleRecord> before{p,p};before[1].stableId=2;
  auto after=before;
  for(auto& record:after) {record.kineticEnergyJ*=.5;record.parallelMomentumKgMPerS*=.5;}
  const auto lost=PopulationControl::CompareInvariants(before,after,1e-12);
  bool ok=!lost.status.ok()&&lost.maximumRelativeError>.49&&
      std::find(lost.failedInvariants.begin(),lost.failedInvariants.end(),"energy")!=lost.failedInvariants.end()&&
      std::find(lost.failedInvariants.begin(),lost.failedInvariants.end(),"parallel-momentum")!=lost.failedInvariants.end();
  before[1].parallelMomentumKgMPerS=-before[0].parallelMomentumKgMPerS;
  after=before;after[0].parallelMomentumKgMPerS*=1.1;
  const auto cancellation=PopulationControl::CompareInvariants(before,after,1e-12);
  ok=ok&&!cancellation.status.ok()&&cancellation.maximumRelativeError>.049;
  return Finish(ok,"SI-scale negative controls reject half-energy/momentum loss and a perturbation of zero net signed momentum");
}
Result DomainRejection() {
  const auto input=Inputs();bool ok=true;
  for(unsigned defect=0;defect<8;++defect) {
    auto particles=input;
    if(defect==0)particles[0].identity=particles[1].identity;
    if(defect==1)particles[0].weight=-1;
    if(defect==2)particles[0].normalMPerS=-1;
    if(defect==3)particles[0].parallelMPerS=c;
    if(defect==4)particles[0].localFraction=1.1;
    if(defect==5)particles[0].massKg=2*mass;
    if(defect==6)particles[0].chargeC=std::numeric_limits<double>::quiet_NaN();
    if(defect==7)particles[0].weight=std::numeric_limits<double>::infinity();
    const auto endpoints=C::Prepare(particles,1e-12);
    ok=ok&&endpoints.weights.status==K::Status::InvalidInput&&endpoints.weights.positive.empty();
  }
  const auto strict=C::Prepare(input,0.0);
  ok=ok&&strict.weights.status==K::Status::NumericalRejection;
  return Finish(ok,"invalid domains and an unmet exact-zero residual budget return explicit rejection without replacement physics or weights");
}
Result StressDiagnostics() {
  const auto input=Inputs();const auto endpoints=C::Prepare(input,1e-12);
  if(endpoints.weights.status!=K::Status::Success)return Finish(false,endpoints.weights.message.c_str());
  const double p=endpoints.weights.positiveProbability;
  const double before[]={endpoints.before.parallelCurrentCMPerS,
      endpoints.before.parallelStressJ,endpoints.before.perpendicularStressJ};
  const double positive[]={endpoints.positive.parallelCurrentCMPerS,
      endpoints.positive.parallelStressJ,endpoints.positive.perpendicularStressJ};
  const double negative[]={endpoints.negative.parallelCurrentCMPerS,
      endpoints.negative.parallelStressJ,endpoints.negative.perpendicularStressJ};
  bool ok=true,changed=false;
  for(unsigned i=0;i<3;++i) {
    const double scale=std::max(std::fabs(before[i]),std::max(std::fabs(positive[i]),std::fabs(negative[i])));
    ok=ok&&std::isfinite(scale)&&scale>0&&
        std::fabs(p*positive[i]+(1-p)*negative[i]-before[i])<=1e-12*scale;
    changed=changed||std::fabs(positive[i]-before[i])>1e-10*scale;
  }
  return Finish(ok&&changed,"unpreserved current/stress changes are explicit and average to their original moments; mu-squared is not mislabelled stress conservation");
}
Result RepeatedRefinement() {
  auto particles=Inputs();const auto original=particles;
  unsigned passes=0;bool ok=true;
  for(;passes<9;++passes) {
    const auto endpoints=C::Prepare(particles,1e-12);
    if(endpoints.weights.status==K::Status::NoNullDirection)break;
    if(endpoints.weights.status!=K::Status::Success) {ok=false;break;}
    std::vector<C::Particle> next;
    for(std::size_t i=0;i<particles.size();++i)if(endpoints.weights.positive[i]>0) {
      auto survivor=particles[i];survivor.weight=endpoints.weights.positive[i];next.push_back(survivor);
    }
    ok=ok&&next.size()<particles.size();particles.swap(next);
  }
  std::vector<double> expanded(original.size(),0.0);
  for(std::size_t i=0;i<original.size();++i)for(const auto& p:particles)
    if(original[i].identity==p.identity)expanded[i]=p.weight;
  const double error=MomentError(original,expanded);
  ok=ok&&passes>0&&particles.size()<=6&&error<=4e-12;
  Result result=Finish(ok,"repeated deletion stops at constraint rank without changing surviving phase points or relaxing the original SI moments");
  result.metrics.push_back({"repeated_maximum_relative_residual",error,4e-12,"<=","dimensionless"});
  return result;
}
Result ThousandParticleStress() {
  const auto start=std::chrono::steady_clock::now();
  std::vector<C::Particle> particles;
  for(unsigned i=0;i<1000;++i) {
    C::Particle p;p.identity=i+1;p.species=0;p.massKg=mass;p.chargeC=charge;
    p.weight=1+std::fmod(.137*i,3.0);p.localFraction=std::fmod(.6180339887498949*(i+.5),1.0);
    const double speed=c*(.02+.95*(i+.5)/1000),mu=.98*std::sin(.73*i);
    p.parallelMPerS=speed*mu;p.normalMPerS=speed*std::sqrt(1-mu*mu);particles.push_back(p);
  }
  const auto initial=particles;
  std::vector<std::vector<double>> rows;std::vector<double> before;
  bool ok=C::EvaluateConstraints(initial,c,&rows,&before).ok();
  std::vector<double> current=before;
  Transport::KeyedRandomStream random(TestSeed(0),1000,UINT64_C(0x4d45524745),0);
  unsigned passes=0;double firstDraw=0.0;
  while(ok&&particles.size()>800&&passes<200) {
    const auto endpoints=C::Prepare(particles,1e-12);
    if(endpoints.weights.status!=K::Status::Success) {ok=false;break;}
    for(const auto* endpoint:{&endpoints.weights.positive,&endpoints.weights.negative}) {
      auto candidate=current;
      for(std::size_t i=0;i<particles.size();++i)candidate[particles[i].identity-1]=(*endpoint)[i];
      ok=ok&&C::MaximumRelativeResidual(rows,before,candidate)<=1e-12;
    }
    if(!ok)break;
    std::vector<double> selected;
    const double actual=random.UniformOpen01();if(passes==0)firstDraw=actual;
    ok=C::Select(endpoints,actual,&selected).ok();
    std::vector<C::Particle> survivors;
    for(std::size_t i=0;i<particles.size();++i) {
      current[particles[i].identity-1]=selected[i];
      if(selected[i]>0) {auto p=particles[i];p.weight=selected[i];survivors.push_back(p);}
    }
    particles.swap(survivors);++passes;
  }
  const double error=C::MaximumRelativeResidual(rows,before,current);
  const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  ok=ok&&particles.size()<=800&&passes<=200&&error<=1e-12;
  Result result=Finish(ok,"bounded 1000-to-800 stress retains phase-point identities and original SI moments; runtime is evidence, not a portable performance claim");
  result.hasSeed=true;result.seed=TestSeed(0);
  result.configuration.push_back("stress_initial=1000;target=800;maximum_deletions=200;observed_runtime_s="+std::to_string(seconds));
  result.metrics.push_back({"stress_maximum_relative_residual",error,1e-12,"<=","dimensionless"});
  result.metrics.push_back({"actual_seeded_observable",firstDraw,0.0,"observed","dimensionless"});
  return result;
}
Result PredictableStopping() {
  // Exact two-stage probability tree. The negative control selectively rolls
  // back a stopped branch and thereby biases the original weights. Keeping the
  // last accepted prefix instead gives a bounded stopped weight martingale.
  const auto first=K::Prepare({1,1},{{1,1}},1e-12);
  if(first.status!=K::Status::Success)return Finish(false,"first stopping-tree endpoint failed");
  const auto second=K::Prepare({first.positive[0],1},{{1,1}},1e-12);
  if(second.status!=K::Status::Success)return Finish(false,"second stopping-tree endpoint failed");
  const long double p=first.positiveProbability,q=second.positiveProbability;
  const long double stopped[]={
      p*(q*second.positive[0]+(1-q)*second.negative[0])+(1-p)*first.negative[0],
      p*first.positive[1]+(1-p)*first.negative[1],
      p*(q*second.positive[1]+(1-q)*second.negative[1])+(1-p)};
  const long double rolledBack[]={
      p*(q*second.positive[0]+(1-q)*second.negative[0])+(1-p),
      p*first.positive[1]+(1-p),
      p*(q*second.positive[1]+(1-q)*second.negative[1])+(1-p)};
  bool ok=true;double bias=0;
  for(unsigned i=0;i<3;++i) {
    ok=ok&&std::fabs(stopped[i]-1.0L)<=1e-12;
    bias=std::max(bias,static_cast<double>(std::fabs(rolledBack[i]-1.0L)));
  }
  ok=ok&&bias>.49;
  return Finish(ok,"predictable stopping preserves each expected weight; outcome-dependent whole-pass rollback negative control produces 50 percent bias");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.group="conservative-resampling";
  descriptor.description=name;descriptor.callback=callback;
  descriptor.initialization=InitializationLevel::None;descriptor.runtime=RuntimeClass::Routine;
  descriptor.supportedBuildModes="native registry and strict dependency-light C++11 sanitizer launcher";
  const bool stochastic=descriptor.id=="CRS04"||descriptor.id=="CRS09";
  descriptor.seedCapability=stochastic?SeedCapability::Context:SeedCapability::Deterministic;
  descriptor.repeatSupported=true;
  descriptor.seedPolicy=stochastic?"TestSeed schema1; independent population stream default0":
      "deterministic endpoints; no RNG draws";
  descriptor.stateIsolation="private weights/constraints only; no PIC or transport-stream mutation";
  return descriptor;
}
} // namespace
std::vector<Descriptor> ConservationResamplingDescriptors() {
  return {Make("CRS01","Conservative relativistic SI endpoints",EndpointsConserve),
          Make("CRS02","Independent analytic nullspace endpoints",AnalyticNullspace),
          Make("CRS03","Constraint rank and sparse no-reduction",RankAndSparse),
          Make("CRS04","Conditional distribution preservation",UnbiasedDistribution),
          Make("CRS05","SI conservation negative controls",SiNegativeControls),
          Make("CRS06","Atomic numerical/domain rejection",DomainRejection),
          Make("CRS07","Unpreserved stress/current diagnostics",StressDiagnostics),
          Make("CRS08","Repeated conservative deletion refinement",RepeatedRefinement),
          Make("CRS09","Bounded thousand-particle conservation stress",ThousandParticleStress),
          Make("CRS10","Predictable numerical stopping unbiasedness",PredictableStopping)};
}
}}
