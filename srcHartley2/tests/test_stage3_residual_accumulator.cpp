#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "src/models/chemical_solver/SystematicResidualAccumulator.h"

namespace Accumulator=ChemicalSolver::SystematicResidualAccumulator;

namespace {

bool Close(double a,double b,double tolerance) {
  return std::fabs(a-b)<=tolerance;
}

bool CheckCase(double mean,long int events,const char *name) {
  constexpr long int trials=20000;
  const double sum=events*mean;
  const double fraction=sum-std::floor(sum);
  const double systematicVariance=fraction*(1.0-fraction);
  const double eventFraction=mean-std::floor(mean);
  const double bernoulliVariance=events*eventFraction*(1.0-eventFraction);

  double systematicFirst=0.0,systematicSecond=0.0;
  for (long int trial=0;trial<trials;trial++) {
    double phase=(trial+0.5)/trials;
    long int total=0;
    for (long int event=0;event<events;event++) {
      const Accumulator::cPlan plan=Accumulator::Plan(phase,mean);
      if (plan.Multiplicity<0 || !Accumulator::IsValidPhase(plan.NextPhase))
        return false;
      phase=plan.NextPhase;
      total+=plan.Multiplicity;
      if (!(std::fabs(total-(event+1)*mean)<1.0)) return false;
    }
    systematicFirst+=total;
    systematicSecond+=static_cast<double>(total)*total;
  }
  const double systematicMean=systematicFirst/trials;
  const double measuredSystematicVariance=systematicSecond/trials-
      systematicMean*systematicMean;

  std::mt19937_64 generator(982451653);
  std::uniform_real_distribution<double> uniform(0.0,1.0);
  double bernoulliFirst=0.0,bernoulliSecond=0.0;
  for (long int trial=0;trial<trials;trial++) {
    long int total=0;
    for (long int event=0;event<events;event++) {
      total+=static_cast<long int>(std::floor(mean));
      if (uniform(generator)<eventFraction) total++;
    }
    bernoulliFirst+=total;
    bernoulliSecond+=static_cast<double>(total)*total;
  }
  const double bernoulliMean=bernoulliFirst/trials;
  const double measuredBernoulliVariance=bernoulliSecond/trials-
      bernoulliMean*bernoulliMean;

  const bool pass=Close(systematicMean,sum,1.0e-12) &&
      Close(measuredSystematicVariance,systematicVariance,1.0e-12) &&
      Close(bernoulliMean,sum,0.04*std::sqrt(bernoulliVariance+1.0)) &&
      Close(measuredBernoulliVariance,bernoulliVariance,
          0.06*(bernoulliVariance+1.0));
  std::cout << (pass ? "PASS" : "FAIL") << " stage3-case=" << name
      << " mu=" << mean << " events=" << events
      << " systematic-mean=" << systematicMean
      << " systematic-variance=" << measuredSystematicVariance
      << " bernoulli-variance=" << measuredBernoulliVariance << '\n';
  return pass;
}

} // namespace

int main() {
  bool pass=true;
  pass=CheckCase(0.25,37,"mu-lt1") && pass;
  pass=CheckCase(1.0,37,"mu-integer") && pass;
  pass=CheckCase(2.75,37,"mu-gt1") && pass;

  // Planning is side-effect free: a failed physical transaction can discard
  // this value, and retrying from the same phase produces an identical plan.
  const double originalPhase=0.83;
  const Accumulator::cPlan first=Accumulator::Plan(originalPhase,0.42);
  const Accumulator::cPlan retry=Accumulator::Plan(originalPhase,0.42);
  pass=first.Multiplicity==retry.Multiplicity &&
      std::memcmp(&first.NextPhase,&retry.NextPhase,sizeof(double))==0 && pass;

  // Concrete products with the same packet phase and carrier remain exactly
  // synchronized. Different means use distinct state slots and diverge.
  double phase0=Accumulator::UnitPhase(UINT64_C(1234567));
  double phase1=phase0;
  for (int event=0;event<100;event++) {
    const Accumulator::cPlan p0=Accumulator::Plan(phase0,0.37);
    const Accumulator::cPlan p1=Accumulator::Plan(phase1,0.37);
    pass=p0.Multiplicity==p1.Multiplicity &&
        p0.NextPhase==p1.NextPhase && pass;
    phase0=p0.NextPhase;
    phase1=p1.NextPhase;
  }

  // A raw state-vector round trip must continue bitwise identically. This
  // complements the integrated core restart test that serializes block data.
  std::vector<double> saved={1.0,7.0,phase0,1.0,0.125,0.0};
  std::vector<unsigned char> bytes(saved.size()*sizeof(double));
  std::memcpy(bytes.data(),saved.data(),bytes.size());
  std::vector<double> loaded(saved.size(),0.0);
  std::memcpy(loaded.data(),bytes.data(),bytes.size());
  pass=std::memcmp(saved.data(),loaded.data(),bytes.size())==0 && pass;
  const Accumulator::cPlan continued0=Accumulator::Plan(saved[2],2.75);
  const Accumulator::cPlan continued1=Accumulator::Plan(loaded[2],2.75);
  pass=continued0.Multiplicity==continued1.Multiplicity &&
      continued0.NextPhase==continued1.NextPhase && pass;

  pass=Accumulator::Plan(-0.1,1.0).Multiplicity<0 && pass;
  pass=Accumulator::Plan(0.1,-1.0).Multiplicity<0 && pass;
  pass=Accumulator::Plan(0.1,
      std::numeric_limits<double>::infinity()).Multiplicity<0 && pass;
  pass=Accumulator::Plan(0.1,static_cast<double>(
      std::numeric_limits<long int>::max())).Multiplicity<0 && pass;

  std::cout << (pass ? "PASS" : "FAIL")
      << " Stage-3 pure residual accumulator plan-commit=1 packet=1 "
      << "round-trip=1 invalid-input=1\n";
  return pass ? 0 : 1;
}
