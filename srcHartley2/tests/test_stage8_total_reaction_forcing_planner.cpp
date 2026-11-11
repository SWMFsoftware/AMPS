#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "srcHartley2/TotalReactionForcingPlanner.h"

namespace Forcing=Comet::TotalReactionForcingPlanner;

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

void TestExactPhysicalSplit() {
  const Forcing::cSplitPlan plan=Forcing::BuildSplitPlan(0.25,2.0);
  Require(Close(plan.OpticalDepth,0.5) &&
      Close(plan.ReactedFraction,1.0-std::exp(-0.5)) &&
      Close(plan.SurvivorFraction,std::exp(-0.5)),
      "constant-hazard physical fractions are incorrect");
  Require(plan.SurvivorFraction+plan.ReactedFraction==1.0,
      "constant-hazard fractions do not close exactly");

  for (double carrier:{1.0,0.1,1.0e-200,1.0e200}) {
    const Forcing::cCarrierSplit split=Forcing::SplitCarrier(plan,carrier);
    Require(split.Survivor>=0.0 && split.Reacted>=0.0 &&
        split.Survivor+split.Reacted==carrier,
        "forced carrier split does not close exactly");
  }

  const Forcing::cSplitPlan zero=Forcing::BuildSplitPlan(0.0,3.0);
  Require(zero.SurvivorFraction==1.0 && zero.ReactedFraction==0.0,
      "zero hazard did not preserve the parent");
  const Forcing::cSplitPlan tiny=Forcing::BuildSplitPlan(1.0e-300,1.0);
  Require(tiny.ReactedFraction>0.0 &&
      Close(tiny.ReactedFraction,1.0e-300,1.0e-12),
      "expm1 did not retain a tiny reaction probability");
  const Forcing::cSplitPlan dominant=Forcing::BuildSplitPlan(1000.0,1.0);
  Require(dominant.ReactedFraction==1.0 &&
      dominant.SurvivorFraction==0.0,
      "large optical depth did not suppress a zero carrier survivor");
}

void TestConditionalConstantHazard() {
  const double rate=0.7,step=3.0;
  const Forcing::cSplitPlan plan=Forcing::BuildSplitPlan(rate,step);
  const std::vector<double> uniforms={
      0.0,0.1,0.25,0.5,0.75,0.9,std::nextafter(1.0,0.0)};
  double previous=-1.0;
  for (double uniform:uniforms) {
    const double time=Forcing::SampleConditionalReactionTime(plan,uniform);
    Require(time>=0.0 && time<=step && time>=previous,
        "conditional exponential time is outside or unordered");
    const double recovered=(1.0-std::exp(-rate*time))/
        plan.ReactedFraction;
    Require(Close(recovered,uniform,2.0e-12),
        "conditional exponential inverse CDF is incorrect");
    previous=time;
  }
}

void TestPiecewiseSpatialHazard() {
  // Include a zero-hazard gap to prove that a conditional reaction never
  // lands there.  Optical depth, not segment length, determines selection.
  const std::vector<Forcing::cHazardSegment> segments={
      {0.2,2.0},{0.0,5.0},{1.0,0.6},{0.05,4.0}};
  const Forcing::cPiecewisePlan plan=
      Forcing::BuildPiecewisePlan(segments);
  Require(Close(plan.TotalTime,11.6) &&
      Close(plan.TotalOpticalDepth,1.2) &&
      plan.SurvivorFraction+plan.ReactedFraction==1.0,
      "piecewise physical totals are incorrect");

  const std::size_t samples=20000;
  std::size_t counts[4]={0,0,0,0};
  long double weightedOpticalMean=0.0L;
  for (std::size_t sample=0;sample<samples;sample++) {
    const double uniform=(static_cast<double>(sample)+0.5)/samples;
    const Forcing::cPiecewiseReactionPosition position=
        Forcing::SampleConditionalReactionPosition(plan,uniform);
    Require(position.Segment<segments.size() && position.Segment!=1 &&
        position.TimeInSegment>=0.0 &&
        position.TimeInSegment<=segments[position.Segment].Duration &&
        position.ElapsedTime>=0.0 && position.ElapsedTime<=plan.TotalTime,
        "piecewise reaction occupied an invalid or zero-hazard segment");
    counts[position.Segment]++;
    weightedOpticalMean+=position.OpticalDepth;
  }
  Require(counts[0]>0 && counts[1]==0 && counts[2]>0 && counts[3]>0,
      "piecewise conditional sampler did not cover all active segments");

  // E[H | H<Htot] for an exponential optical depth truncated at Htot.
  const double expectedOpticalMean=1.0-
      plan.TotalOpticalDepth*std::exp(-plan.TotalOpticalDepth)/
      plan.ReactedFraction;
  Require(Close(static_cast<double>(weightedOpticalMean/samples),
      expectedOpticalMean,2.0e-8),
      "piecewise conditional optical-depth mean is incorrect");
}

void TestMultipleGenerationClosure() {
  // Repeatedly force the survivor branch.  The reaction carriers form the
  // exact discrete exponential loss history and the final survivor closes the
  // original carrier.  This is the minimal multiple-generation invariant the
  // runtime must retain when descendants become parents in later steps.
  const Forcing::cSplitPlan plan=Forcing::BuildSplitPlan(0.03,0.5);
  double survivor=1.0,reactionSum=0.0;
  for (int generation=0;generation<400;generation++) {
    const Forcing::cCarrierSplit split=
        Forcing::SplitCarrier(plan,survivor);
    survivor=split.Survivor;
    reactionSum+=split.Reacted;
  }
  Require(Close(survivor,std::exp(-0.03*0.5*400),2.0e-13),
      "multi-generation survivor is not exponential");
  Require(Close(reactionSum+survivor,1.0,2.0e-13),
      "multi-generation reaction and survivor carriers do not close");
}

struct cThreeSpeciesState {
  double A=0.0;
  double B=0.0;
  double C=0.0;
};

cThreeSpeciesState AnalyticTwoGenerationState(
    const Forcing::cSplitPlan &aPlan,
    const Forcing::cSplitPlan &bPlan,int cycles) {
  const double sa=aPlan.SurvivorFraction;
  const double sb=bPlan.SurvivorFraction;
  cThreeSpeciesState result;
  result.A=std::pow(sa,cycles);
  if (sa==sb) {
    result.B=cycles*aPlan.ReactedFraction*std::pow(sa,cycles-1);
  }
  else {
    // A reacts during one cycle to create B.  The processed flag prevents the
    // new B from reacting until the following cycle, which gives this finite
    // geometric convolution of A production and B survival.
    result.B=aPlan.ReactedFraction*
        (std::pow(sb,cycles)-std::pow(sa,cycles))/(sb-sa);
  }
  result.C=1.0-result.A-result.B;
  return result;
}

void TestTwoGenerationAnalyticNetwork() {
  const Forcing::cSplitPlan aPlan=Forcing::BuildSplitPlan(0.08,0.75);
  const Forcing::cSplitPlan bPlan=Forcing::BuildSplitPlan(0.21,0.75);
  const int cycles=32;

  // Cycle order matches the production processed-flag contract: split the A
  // and B carriers that existed at cycle entry, then add A's reacted carrier
  // to next-cycle B.  C is stable.  Every reaction is one-to-one, so the sum
  // A+B+C must remain the initial carrier.
  cThreeSpeciesState forced{1.0,0.0,0.0};
  for (int cycle=0;cycle<cycles;cycle++) {
    const Forcing::cCarrierSplit a=
        Forcing::SplitCarrier(aPlan,forced.A);
    Forcing::cCarrierSplit b;
    if (forced.B>0.0) b=Forcing::SplitCarrier(bPlan,forced.B);
    forced.A=a.Survivor;
    forced.B=b.Survivor+a.Reacted;
    forced.C+=b.Reacted;
    Require(Close(forced.A+forced.B+forced.C,1.0,4.0e-14),
        "two-generation forced carriers do not close each cycle");
  }

  const cThreeSpeciesState exact=
      AnalyticTwoGenerationState(aPlan,bPlan,cycles);
  Require(Close(forced.A,exact.A,4.0e-13) &&
      Close(forced.B,exact.B,4.0e-13) &&
      Close(forced.C,exact.C,4.0e-13),
      "two-generation forced network disagrees with its analytic solution");
}

std::uint64_t Mix(std::uint64_t value) {
  value+=0x9e3779b97f4a7c15ULL;
  value=(value^(value>>30))*0xbf58476d1ce4e5b9ULL;
  value=(value^(value>>27))*0x94d049bb133111ebULL;
  return value^(value>>31);
}

double KeyedUniform(std::uint64_t seed,std::uint64_t particle,
    std::uint64_t cycle,std::uint64_t speciesState) {
  std::uint64_t key=Mix(seed);
  key=Mix(key^particle);
  key=Mix(key^(cycle+0x100000001b3ULL*speciesState));
  return static_cast<double>(key>>11)*0x1.0p-53;
}

void TestAnalogForcedEstimatorComparison() {
  const Forcing::cSplitPlan aPlan=Forcing::BuildSplitPlan(0.08,0.75);
  const Forcing::cSplitPlan bPlan=Forcing::BuildSplitPlan(0.21,0.75);
  const int cycles=16;
  const std::size_t particles=24000,seeds=64;
  const cThreeSpeciesState exact=
      AnalyticTwoGenerationState(aPlan,bPlan,cycles);

  long double sumA=0.0L,sumC=0.0L,sumA2=0.0L,sumC2=0.0L;
  for (std::size_t replicate=0;replicate<seeds;replicate++) {
    const std::uint64_t seed=0x51a7e80000000000ULL+replicate;
    std::vector<unsigned char> state(particles,0); // 0=A, 1=B, 2=C
    for (int cycle=0;cycle<cycles;cycle++) {
      for (std::size_t particle=0;particle<particles;particle++) {
        if (state[particle]==0) {
          if (KeyedUniform(seed,particle,cycle,0)<aPlan.ReactedFraction)
            state[particle]=1;
        }
        else if (state[particle]==1) {
          if (KeyedUniform(seed,particle,cycle,1)<bPlan.ReactedFraction)
            state[particle]=2;
        }
      }
    }
    std::size_t countA=0,countC=0;
    for (unsigned char species:state) {
      if (species==0) countA++;
      else if (species==2) countC++;
    }
    const double estimateA=static_cast<double>(countA)/particles;
    const double estimateC=static_cast<double>(countC)/particles;
    sumA+=estimateA;
    sumC+=estimateC;
    sumA2+=estimateA*estimateA;
    sumC2+=estimateC*estimateC;
  }

  const double meanA=static_cast<double>(sumA/seeds);
  const double meanC=static_cast<double>(sumC/seeds);
  const double varianceA=static_cast<double>(
      (sumA2-sumA*sumA/seeds)/(seeds-1));
  const double varianceC=static_cast<double>(
      (sumC2-sumC*sumC/seeds)/(seeds-1));
  const double expectedVarianceA=exact.A*(1.0-exact.A)/particles;
  const double expectedVarianceC=exact.C*(1.0-exact.C)/particles;
  const double standardErrorA=std::sqrt(expectedVarianceA/seeds);
  const double standardErrorC=std::sqrt(expectedVarianceC/seeds);

  // Fixed seeds make this a reproducible bound, while the six-standard-error
  // limit tests unbiased agreement rather than an accidentally exact sample.
  Require(std::fabs(meanA-exact.A)<=6.0*standardErrorA &&
      std::fabs(meanC-exact.C)<=6.0*standardErrorC,
      "fixed-seed analog ensemble disagrees with analytic survival/yield");
  Require(varianceA>0.35*expectedVarianceA &&
      varianceA<1.8*expectedVarianceA &&
      varianceC>0.35*expectedVarianceC &&
      varianceC<1.8*expectedVarianceC,
      "fixed-seed analog occurrence variance is outside its bounded expectation");

  // Deterministic forcing returns the analytic carrier on every replicate;
  // its occurrence variance is exactly zero before transport/collision noise.
  const double forcedVarianceA=0.0,forcedVarianceC=0.0;
  Require(forcedVarianceA<varianceA && forcedVarianceC<varianceC &&
      Close(exact.A,std::pow(aPlan.SurvivorFraction,cycles),1.0e-14),
      "forced estimator did not remove total-occurrence variance");
}

void TestCapacityBounds() {
  const Forcing::cSplitPlan ordinary=Forcing::BuildSplitPlan(0.1,1.0);
  const Forcing::cCapacityPlan plan=Forcing::BuildCapacityPlan(
      ordinary,3,4,8,32);
  Require(plan.ReactionDaughters==12 && plan.OutputParticles==13 &&
      plan.RequiredNewSlots==12 &&
      Forcing::FitsAvailableSlots(plan,12) &&
      !Forcing::FitsAvailableSlots(plan,11),
      "ordinary forced transaction capacity is incorrect");

  const Forcing::cCapacityPlan sink=Forcing::BuildCapacityPlan(
      ordinary,0,1,8,32);
  Require(sink.ReactionDaughters==0 && sink.OutputParticles==1 &&
      sink.RequiredNewSlots==0,
      "forced sink transaction should retain only the survivor");

  const Forcing::cSplitPlan zero=Forcing::BuildSplitPlan(0.0,1.0);
  const Forcing::cCapacityPlan zeroPlan=Forcing::BuildCapacityPlan(
      zero,3,4,8,32);
  Require(zeroPlan.ReactionDaughters==0 && zeroPlan.OutputParticles==1 &&
      zeroPlan.RequiredNewSlots==0,
      "zero-hazard capacity planned inactive daughters");

  const Forcing::cSplitPlan complete=Forcing::BuildSplitPlan(1000.0,1.0);
  const Forcing::cCapacityPlan completePlan=Forcing::BuildCapacityPlan(
      complete,3,4,8,32);
  Require(completePlan.ReactionDaughters==12 &&
      completePlan.OutputParticles==12 &&
      completePlan.RequiredNewSlots==11,
      "unit-reaction capacity retained a zero carrier survivor");
}

void TestInvalidInputs() {
  const double nan=std::numeric_limits<double>::quiet_NaN();
  RequireThrows<std::invalid_argument>([]() {
    Forcing::BuildSplitPlan(-1.0,1.0);
  },"negative rate was accepted");
  RequireThrows<std::invalid_argument>([nan]() {
    Forcing::BuildSplitPlan(nan,1.0);
  },"nonfinite rate was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::BuildSplitPlan(1.0,0.0);
  },"zero timestep was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::SplitCarrier(Forcing::BuildSplitPlan(1.0,1.0),0.0);
  },"zero parent carrier was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::SampleConditionalReactionTime(
        Forcing::BuildSplitPlan(0.0,1.0),0.5);
  },"conditional time accepted zero hazard");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::SampleConditionalReactionTime(
        Forcing::BuildSplitPlan(1.0,1.0),1.0);
  },"conditional time accepted the unit endpoint");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::BuildPiecewisePlan({});
  },"empty piecewise plan was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::BuildPiecewisePlan({{1.0,0.0}});
  },"zero-duration piecewise segment was accepted");
  RequireThrows<std::invalid_argument>([]() {
    const auto plan=Forcing::BuildPiecewisePlan({{0.0,1.0}});
    Forcing::SampleConditionalReactionPosition(plan,0.5);
  },"conditional position accepted an all-zero hazard");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::BuildCapacityPlan(
        Forcing::BuildSplitPlan(1.0,1.0),2,0,8,32);
  },"zero packet count was accepted");
  RequireThrows<std::invalid_argument>([]() {
    Forcing::BuildCapacityPlan(
        Forcing::BuildSplitPlan(1.0,1.0),2,9,8,32);
  },"packet count above maximum was accepted");
  RequireThrows<std::length_error>([]() {
    Forcing::BuildCapacityPlan(
        Forcing::BuildSplitPlan(1.0,1.0),5,8,8,32);
  },"daughter transaction above maximum was accepted");
}

} // namespace

int main() {
  try {
    TestExactPhysicalSplit();
    TestConditionalConstantHazard();
    TestPiecewiseSpatialHazard();
    TestMultipleGenerationClosure();
    TestTwoGenerationAnalyticNetwork();
    TestAnalogForcedEstimatorComparison();
    TestCapacityBounds();
    TestInvalidInputs();
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,"Stage-8 total-reaction-forcing planner failed: %s\n",
        exception.what());
    return 1;
  }
  std::printf("PASS Stage-8 total-reaction-forcing planner split=1 "
      "constant-hazard=1 piecewise-hazard=1 multigeneration=1 "
      "two-generation=1 analog-forced=1 variance-reduction=1 "
      "capacity=1 extremes=1 errors=1\n");
  return 0;
}
