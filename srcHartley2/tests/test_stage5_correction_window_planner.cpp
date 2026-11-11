#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "srcHartley2/CorrectionWindowPlanner.h"

namespace Window=Comet::CorrectionWindowPlanner;

namespace {

void Require(bool condition,const char *message) {
  if (!condition) throw std::runtime_error(message);
}

bool Close(long double left,long double right,long double tolerance=1.0e-13L) {
  return std::fabs(left-right)<=tolerance*
      std::max(1.0L,std::max(std::fabs(left),std::fabs(right)));
}

template<class Exception,class Function>
void RequireThrows(Function function,const char *message) {
  bool threw=false;
  try { function(); }
  catch (const Exception&) { threw=true; }
  Require(threw,message);
}

long double Sum(const std::vector<double> &values) {
  long double result=0.0L;
  for (double value:values) result+=value;
  return result;
}

double ProductionOrderSum(const std::vector<double> &values) {
  // Mirror the runtime assignment order in binary64.  The final remainder
  // must restore the incoming carrier exactly, not only within a tolerance.
  double result=0.0;
  for (double value:values) result+=value;
  return result;
}

void TestPolicyValidationAndWindowBoundaries() {
  const Window::cPolicy enabled(0.25,1.0,4.0,32);
  const Window::cPolicy disabled;

  Window::cParticlePlan plan=Window::BuildParticlePlan(
      Window::cRequest(0,0.01),disabled);
  Require(plan.Action==Window::eAction::Keep &&
      plan.OutputCorrections==std::vector<double>{0.01},
      "disabled correction window changed a particle");

  for (double correction:{0.25,1.0,4.0}) {
    plan=Window::BuildParticlePlan(Window::cRequest(0,correction),enabled);
    Require(plan.Action==Window::eAction::Keep &&
        plan.OutputCorrections.size()==1 &&
        plan.OutputCorrections[0]==correction,
        "inclusive correction-window boundary was not retained");
  }

  RequireThrows<std::invalid_argument>([]() {
    Window::ValidatePolicy(Window::cPolicy(0.0,1.0,4.0,4));
  },"zero lower correction was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Window::ValidatePolicy(Window::cPolicy(1.0,1.0,4.0,4));
  },"degenerate lower/target correction was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Window::ValidatePolicy(Window::cPolicy(0.25,4.0,4.0,4));
  },"degenerate target/upper correction was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Window::ValidatePolicy(Window::cPolicy(0.25,1.0,4.0,1));
  },"enabled one-particle split limit was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Window::ValidatePolicy(Window::cPolicy(0.9,1.0,1.1,8));
  },"window incompatible with ceil splitting was accepted");
  RequireThrows<std::invalid_argument>([&]() {
    Window::BuildParticlePlan(Window::cRequest(0,
        std::numeric_limits<double>::infinity()),enabled);
  },"infinite particle correction was accepted");
}

void TestRouletteMeanAndVariance() {
  const Window::cPolicy policy(0.25,1.0,4.0,32);
  const double correction=0.125;
  const Window::cParticlePlan plan=Window::BuildParticlePlan(
      Window::cRequest(0,correction),policy);
  Require(plan.Action==Window::eAction::Roulette,
      "sub-window correction did not select roulette");
  Require(plan.OutputCorrections.size()==1 &&
      plan.OutputCorrections[0]==policy.TargetCorrection,
      "roulette survivor correction is not the target");
  Require(plan.SurvivalProbability==correction/policy.TargetCorrection,
      "roulette survival probability is not c/c_target");

  constexpr long trials=400000;
  std::mt19937_64 generator(32452843);
  std::uniform_real_distribution<double> uniform(0.0,1.0);
  long double first=0.0L,second=0.0L;
  for (long trial=0;trial<trials;trial++) {
    const double represented=uniform(generator)<plan.SurvivalProbability ?
        plan.OutputCorrections[0] : 0.0;
    first+=represented;
    second+=represented*represented;
  }
  const long double mean=first/trials;
  const long double variance=second/trials-mean*mean;
  const long double expectedVariance=correction*
      (policy.TargetCorrection-correction);
  const long double meanSigma=std::sqrt(expectedVariance/trials);
  const long double varianceTolerance=6.0L*expectedVariance/
      std::sqrt(static_cast<long double>(trials));
  Require(std::fabs(mean-correction)<=6.0L*meanSigma,
      "roulette mean is inconsistent with the input correction");
  Require(std::fabs(variance-expectedVariance)<=varianceTolerance,
      "roulette variance is inconsistent with Bernoulli theory");
}

void TestSplitSumWindowAndLimits() {
  const Window::cPolicy policy(0.25,1.0,4.0,32);
  Window::cParticlePlan plan=Window::BuildParticlePlan(
      Window::cRequest(0,10.0),policy);
  Require(plan.Action==Window::eAction::Split &&
      plan.OutputCorrections.size()==10,
      "correction ten did not split to the target count");
  Require(Close(Sum(plan.OutputCorrections),10.0L),
      "split corrections do not sum to the input correction");
  Require(ProductionOrderSum(plan.OutputCorrections)==10.0,
      "split corrections are not exact in production order");
  for (double output:plan.OutputCorrections)
    Require(output>=policy.LowerCorrection &&
        output<=policy.UpperCorrection,
        "split correction is outside the configured window");

  // A non-integer ratio uses ceil(c/target), not a capacity-dependent count.
  plan=Window::BuildParticlePlan(Window::cRequest(0,10.25),policy);
  Require(plan.OutputCorrections.size()==11 &&
      Close(Sum(plan.OutputCorrections),10.25L),
      "ceil split did not conserve correction");
  Require(ProductionOrderSum(plan.OutputCorrections)==10.25,
      "ceil split is not exact in production order");
  for (double output:plan.OutputCorrections)
    Require(output>=policy.LowerCorrection &&
        output<=policy.UpperCorrection,
        "bounded split output is outside its window");

  RequireThrows<std::length_error>([&]() {
    Window::BuildParticlePlan(Window::cRequest(0,100.0),policy);
  },"unachievable split was not rejected by the hard limit");
}

void TestBatchCapacitySpeciesAndRepeatedProducts() {
  const std::vector<Window::cPolicy> policies={
      Window::cPolicy(0.25,1.0,4.0,32),
      Window::cPolicy(0.5,2.0,8.0,16),
      Window::cPolicy()
  };
  const std::vector<Window::cRequest> requests={
      Window::cRequest(0,10.0),
      Window::cRequest(0,0.125),
      Window::cRequest(1,20.0),
      Window::cRequest(2,64.0)
  };

  const Window::cBatchPlan accepted=Window::BuildPlan(requests,policies,18);
  Require(accepted.ParticlePlans.size()==requests.size(),
      "multi-species batch lost a request");
  Require(accepted.RequiredAdditionalParticles==18,
      "multi-species batch clone capacity is not exact");
  Require(accepted.ParticlePlans[3].Action==Window::eAction::Keep,
      "disabled species was modified");

  // Planning is pure.  A capacity failure must leave the input vector
  // byte-for-byte unchanged, which models runtime preflight before mutation.
  const std::vector<Window::cRequest> saved=requests;
  RequireThrows<std::length_error>([&]() {
    Window::BuildPlan(requests,policies,17);
  },"insufficient aggregate capacity did not reject the batch");
  Require(requests.size()==saved.size(),"failed batch changed request count");
  for (std::size_t index=0;index<requests.size();index++)
    Require(requests[index].Species==saved[index].Species &&
        requests[index].Correction==saved[index].Correction,
        "failed batch mutated a request");

  // Two identical concrete products must remain two independent requests.
  // Their corrections are derived from the same parent carrier divided by
  // the same daughter base carrier; aggregating them would erase
  // stoichiometry and under-count the hard transaction size.
  const double parentCarrier=8.0,daughterBaseCarrier=2.0;
  const std::vector<Window::cRequest> repeatedProducts={
      Window::cRequest(0,parentCarrier/daughterBaseCarrier),
      Window::cRequest(0,parentCarrier/daughterBaseCarrier)
  };
  const Window::cBatchPlan repeated=Window::BuildPlan(
      repeatedProducts,policies,0);
  Require(repeated.ParticlePlans.size()==2 &&
      repeated.RequiredAdditionalParticles==0,
      "repeated stoichiometric products were aggregated or expanded");
  Require(Close(Sum(repeated.ParticlePlans[0].OutputCorrections),4.0L) &&
      Close(Sum(repeated.ParticlePlans[1].OutputCorrections),4.0L),
      "repeated product physical carriers changed");

  const Window::cBatchPlan empty=Window::BuildPlan({},policies,0);
  Require(empty.ParticlePlans.empty() &&
      empty.RequiredAdditionalParticles==0,
      "empty species population created a correction-window particle");

  RequireThrows<std::out_of_range>([&]() {
    Window::BuildPlan({Window::cRequest(3,1.0)},policies,0);
  },"request without a species policy was accepted");
}

} // namespace

int main() {
  try {
    TestPolicyValidationAndWindowBoundaries();
    TestRouletteMeanAndVariance();
    TestSplitSumWindowAndLimits();
    TestBatchCapacitySpeciesAndRepeatedProducts();
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,"Stage-5 correction-window planner failed: %s\n",
        exception.what());
    return 1;
  }

  std::printf("PASS Stage-5 correction-window planner policy=1 roulette=1 "
      "split=1 capacity=1 multi-species=1 repeated-products=1 empty=1\n");
  return 0;
}
