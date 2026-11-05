#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

#include "srcHartley2/ChannelImportancePlanner.h"

namespace Importance=Comet::ChannelImportancePlanner;

namespace {

void Require(bool condition,const char *message) {
  if (!condition) throw std::runtime_error(message);
}

bool Close(double left,double right,double tolerance=2.0e-13) {
  return std::fabs(left-right)<=tolerance*
      std::max({1.0,std::fabs(left),std::fabs(right)});
}

template<class Exception,class Function>
void RequireThrows(Function function,const char *message) {
  bool threw=false;
  try { function(); }
  catch (const Exception&) { threw=true; }
  Require(threw,message);
}

double Sum(const std::vector<double> &values) {
  double result=0.0;
  for (double value:values) result+=value;
  return result;
}

std::size_t HistoricalSelector(const std::vector<double> &rates,
    double total,double unitUniform) {
  double selector=unitUniform*total;
  std::size_t selected=rates.size()-1;
  for (std::size_t channel=0;channel<rates.size();channel++) {
    if (!(rates[channel]>0.0)) continue;
    selector-=rates[channel];
    if (selector<=0.0) {
      selected=channel;
      break;
    }
  }
  return selected;
}

void TestAnalogIdentityAndBoundaries() {
  const std::vector<double> rates={2.0,3.0,5.0};
  const Importance::cPlan plan=Importance::BuildPlan(rates);
  Require(!plan.ImportanceSamplingEnabled,
      "empty q* policy did not preserve analog mode");
  Require(Sum(plan.PhysicalProbability)==1.0 &&
      Sum(plan.BiasedProbability)==1.0,
      "analog probabilities do not have an exact final remainder");
  for (double likelihood:plan.LikelihoodRatio)
    Require(likelihood==1.0,"analog likelihood is not exactly one");

  const std::vector<double> uniforms={
      0.0,
      std::nextafter(0.2,0.0),0.2,std::nextafter(0.2,1.0),
      std::nextafter(0.5,0.0),0.5,std::nextafter(0.5,1.0),
      std::nextafter(1.0,0.0)};
  for (double uniform:uniforms) {
    const std::size_t historical=HistoricalSelector(rates,10.0,uniform);
    Require(Importance::SelectAnalogRateSubtract(rates,10.0,uniform)==
        historical,"analog rate-subtraction identity changed");
    Require(Importance::SelectBiasedProbabilitySubtract(plan,uniform)==
        historical,"identity q* selection changed a channel boundary");
  }
}

void TestTwoChannelMeanVarianceAndEss() {
  const Importance::cPlan plan=Importance::BuildPlan(
      {99.0,1.0},{0.5,0.5});
  Require(plan.ImportanceSamplingEnabled,
      "two-channel bias was not detected");
  Require(Close(plan.PhysicalProbability[0],0.99) &&
      Close(plan.PhysicalProbability[1],0.01),
      "two-channel physical probabilities are incorrect");
  Require(Close(plan.LikelihoodRatio[0],1.98) &&
      Close(plan.LikelihoodRatio[1],0.02),
      "two-channel likelihoods are incorrect");

  const double rareVariance=0.01*0.01*(1.0/0.5-1.0);
  Require(Close(plan.WeightedChannelVariance[1],rareVariance,1.0e-12),
      "rare-channel analytic variance is incorrect");
  Require(plan.WeightedChannelVariance[1]<0.01*0.99,
      "rare-channel bias did not reduce variance");

  // Midpoint-stratified uniforms produce exactly half the selections in each
  // q*=0.5 interval.  This tests the implemented selector, not only formulas.
  const std::size_t samples=100000;
  long double rareWeight=0.0L,weightSum=0.0L,weightSquares=0.0L;
  for (std::size_t sample=0;sample<samples;sample++) {
    const double uniform=(static_cast<double>(sample)+0.5)/samples;
    const std::size_t channel=
        Importance::SelectBiasedProbabilitySubtract(plan,uniform);
    const double weight=plan.LikelihoodRatio[channel];
    if (channel==1) rareWeight+=weight;
    weightSum+=weight;
    weightSquares+=weight*weight;
  }
  Require(Close(static_cast<double>(rareWeight/samples),0.01,1.0e-12),
      "biased rare-channel weighted mean is not physical");
  Require(Close(static_cast<double>(weightSum/samples),1.0,1.0e-12),
      "mean likelihood is not one");
  const double expectedEss=Importance::ExpectedEffectiveSampleSize(
      plan,samples);
  const double sampleEss=Importance::SampleEffectiveSampleSize(
      static_cast<double>(weightSum),static_cast<double>(weightSquares));
  Require(Close(expectedEss,sampleEss,1.0e-12) &&
      sampleEss>0.0 && sampleEss<=samples,
      "channel effective sample size is incorrect");
}

void TestManyChannelAndFinalRemainder() {
  const Importance::cPlan plan=Importance::BuildPlan(
      {50.0,30.0,15.0,5.0},{0.25,0.25,0.25,0.25});
  Require(Sum(plan.PhysicalProbability)==1.0 &&
      Sum(plan.BiasedProbability)==1.0,
      "many-channel probability remainder is not exact");
  double expectedLikelihood=0.0;
  for (std::size_t channel=0;channel<4;channel++)
    expectedLikelihood+=plan.BiasedProbability[channel]*
        plan.LikelihoodRatio[channel];
  Require(Close(expectedLikelihood,1.0,1.0e-14),
      "many-channel mean likelihood is not one");
  Require(Importance::SelectBiasedProbabilitySubtract(plan,0.0)==0 &&
      Importance::SelectBiasedProbabilitySubtract(plan,0.25)==0 &&
      Importance::SelectBiasedProbabilitySubtract(
          plan,std::nextafter(0.25,1.0))==1 &&
      Importance::SelectBiasedProbabilitySubtract(
          plan,std::nextafter(1.0,0.0))==3,
      "many-channel interval boundaries are incorrect");

  // Decimal inputs do not generally add to one exactly in binary64.  The
  // final active interval must nevertheless close the family exactly.
  const Importance::cPlan decimal=Importance::BuildPlan(
      {1.0,1.0,1.0,1.0},{0.1,0.2,0.3,0.4});
  Require(Sum(decimal.BiasedProbability)==1.0,
      "decimal q* family did not receive an exact final remainder");
}

void TestDominantTinyAndDisabledChannels() {
  const Importance::cPlan tiny=Importance::BuildPlan(
      {1.0,1.0e-12},{0.5,0.5});
  Require(tiny.LikelihoodRatio[0]>1.0 &&
      tiny.LikelihoodRatio[1]>0.0 &&
      std::isfinite(tiny.LikelihoodRatio[1]) &&
      tiny.WeightedChannelVariance[1]>=0.0,
      "dominant/tiny channel plan is invalid");

  const Importance::cPlan disabled=Importance::BuildPlan(
      {1.0,0.0,3.0},{0.25,0.0,0.75});
  Require(disabled.PhysicalProbability[1]==0.0 &&
      disabled.BiasedProbability[1]==0.0 &&
      disabled.LikelihoodRatio[1]==0.0 &&
      disabled.WeightedChannelVariance[1]==0.0,
      "disabled channel acquired statistical support");
  for (std::size_t sample=0;sample<1000;sample++) {
    const double uniform=(sample+0.5)/1000.0;
    Require(Importance::SelectBiasedProbabilitySubtract(
        disabled,uniform)!=1,"disabled channel was selected");
  }
}

void TestInvalidInputs() {
  const double nan=std::numeric_limits<double>::quiet_NaN();
  RequireThrows<std::invalid_argument>([]() {
    Importance::BuildPlan({});
  },"empty family was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::BuildPlan({0.0,0.0});
  },"all-disabled family was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::BuildPlan({1.0,-1.0});
  },"negative physical rate was accepted");
  RequireThrows<std::invalid_argument>([nan]() {
    Importance::BuildPlan({1.0,nan});
  },"non-finite physical rate was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::BuildPlan({1.0,1.0},{1.0});
  },"q* size mismatch was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::BuildPlan({1.0,1.0},{0.4,0.5});
  },"nonnormalized q* was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::BuildPlan({1.0,1.0},{1.0,0.0});
  },"positive channel with q*=0 was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::BuildPlan({1.0,0.0},{0.5,0.5});
  },"disabled channel with q*>0 was accepted");
  RequireThrows<std::invalid_argument>([nan]() {
    Importance::BuildPlan({1.0,1.0},{nan,0.5});
  },"non-finite q* was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::SelectAnalogRateSubtract({1.0,1.0},2.0,1.0);
  },"unit endpoint was accepted as a uniform variate");
  RequireThrows<std::invalid_argument>([]() {
    Importance::SelectAnalogRateSubtract({1.0,1.0},3.0,0.5);
  },"inconsistent analog total was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Importance::SampleEffectiveSampleSize(0.0,1.0);
  },"zero sample weight was accepted by ESS");
  RequireThrows<std::invalid_argument>([]() {
    const Importance::cPlan plan=Importance::BuildPlan({1.0,1.0});
    Importance::ExpectedEffectiveSampleSize(plan,0);
  },"zero expected-ESS sample count was accepted");
}

} // namespace

int main() {
  try {
    TestAnalogIdentityAndBoundaries();
    TestTwoChannelMeanVarianceAndEss();
    TestManyChannelAndFinalRemainder();
    TestDominantTinyAndDisabledChannels();
    TestInvalidInputs();
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,"Stage-7 channel-importance planner failed: %s\n",
        exception.what());
    return 1;
  }
  std::printf("PASS Stage-7 channel-importance planner analog-identity=1 "
      "two-channel=1 many-channel=1 dominant-tiny=1 disabled=1 "
      "boundaries=1 weighted-mean=1 variance=1 ess=1 errors=1\n");
  return 0;
}
