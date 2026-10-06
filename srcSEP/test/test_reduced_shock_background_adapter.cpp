#include "reduced_shock_background_adapter.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

int passed=0;
int failed=0;

void Check(bool condition,const char* id,const char* explanation) {
  if(condition) {
    ++passed;
    std::cout<<id<<" PASS: "<<explanation<<'\n';
  }
  else {
    ++failed;
    std::cerr<<id<<" FAIL: "<<explanation<<'\n';
  }
}

bool FinitePositive(double value) {
  return std::isfinite(value)&&value>0.0;
}

} // namespace

int main(int argc,char** argv) {
  if(argc!=2) {
    std::cerr<<"usage: test_reduced_shock_background_adapter EVENT\n";
    return 2;
  }

  std::string error;
  Check(SEP::ReducedShock::Configure(argv[1],&error),"RSHSEP01",
        "the production adapter resolves the event and relative asset");
  if(!SEP::ReducedShock::Enabled()) {
    std::cerr<<"SUMMARY PASS="<<passed<<" FAIL="<<failed
             <<" SKIP=0 ERROR=0\n";
    return 1;
  }
  Check(!SEP::ReducedShock::EventIdentity().empty(),"RSHSEP02",
        "the shared normalized physics fingerprint is exposed");

  Check(SEP::ReducedShock::Prepare(0.0,&error),"RSHSEP03",
        "event epoch zero commits as the first immutable generation");
  const SEP::ReducedShock::EpochMetadata first=SEP::ReducedShock::Metadata();
  Check(first.epochS==0.0&&first.generation==1&&
        first.validUntilS==SEP::ReducedShock::BackgroundCadenceS(),
        "RSHSEP04","epoch validity follows the declared background cadence");

  // This position lies on the example's +Z native field line at the maintained
  // 20-Rsun inner boundary.  The assertions deliberately cover primitives,
  // not hard-coded profile values, so the canonical provider remains the sole
  // owner of the coronal/Parker equations and unit conversions.
  SEP::ReducedShock::AmbientSample sample;
  const std::array<double,3> position{{0.0,0.0,1.3914e10}};
  const bool queried=SEP::ReducedShock::EvaluateAmbient(position,&sample,&error);
  Check(queried,"RSHSEP05","a native-line HCI position returns ambient state");
  Check(queried&&FinitePositive(sample.numberDensityM3)&&
        FinitePositive(sample.pressurePa)&&
        FinitePositive(sample.protonTemperatureK)&&
        std::isfinite(sample.magneticFieldT[2])&&
        std::isfinite(sample.velocityMPerS[2]),"RSHSEP06",
        "ambient plasma and IMF primitives are finite and physical");

  // One second is neither an application epoch nor a legal provider cadence
  // tick for this fixture.  Its rejection must leave the committed generation
  // intact; otherwise native vertices could acquire a label for an unrealized
  // provider state.
  Check(!SEP::ReducedShock::Prepare(1.0,&error),"RSHSEP07",
        "an off-cadence candidate is rejected transactionally");
  Check(SEP::ReducedShock::Metadata().epochS==first.epochS&&
        SEP::ReducedShock::Metadata().generation==first.generation,
        "RSHSEP08","failed preparation preserves committed metadata");

  const double cadence=SEP::ReducedShock::BackgroundCadenceS();
  const double initialApex=SEP::ReducedShock::CurrentFrontSummary().apexRadiusM;
  Check(SEP::ReducedShock::Prepare(cadence,&error),"RSHSEP09",
        "the next cadence-aligned epoch commits");
  const SEP::ReducedShock::FrontSummary next=
      SEP::ReducedShock::CurrentFrontSummary();
  Check(SEP::ReducedShock::Metadata().generation==2&&
        next.apexRadiusM>initialApex,"RSHSEP10",
        "generation advances with the prescribed outward front");

  std::cout<<"SUMMARY PASS="<<passed<<" FAIL="<<failed
           <<" SKIP=0 ERROR=0\n";
  return failed==0?0:1;
}
