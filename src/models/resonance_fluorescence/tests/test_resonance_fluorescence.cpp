#include "resonance_fluorescence.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

namespace RF=AMPS::ResonanceFluorescence;

namespace {
constexpr double AU=1.495978707e11;

bool Relative(double value,double expected,double tolerance) {
  return std::isfinite(value) && std::isfinite(expected) && expected!=0.0 &&
      std::fabs(value/expected-1.0)<=tolerance;
}

bool Check(bool condition,const char *message) {
  if (!condition) std::fprintf(stderr,"FAIL: %s\n",message);
  return condition;
}
} // namespace

int main() {
  bool pass=true;
  const char *error=nullptr;
  pass=Check(RF::ValidateReferenceData(&error),
      error==nullptr ? "compiled reference validation" : error) && pass;

  RF::Model model=RF::Model::HydrogenLymanAlphaIKL2020;
  pass=Check(RF::ModelFromName("H_LYA_IKL2020",&model) &&
      model==RF::Model::HydrogenLymanAlphaIKL2020,
      "strict H model-name conversion") && pass;
  pass=Check(!RF::ModelFromName("H",&model),"unknown name rejection") && pass;

  // These six independent scalars are the reviewed 3I/ATLAS case metadata.
  // The old calculation used slightly rounded constants, so 1e-4 tests both
  // the complete line-profile shape and the physical normalization.
  const double distanceAU[]={1.75,1.66,1.54,1.36,1.36,1.37};
  const double velocityKMS[]={-39.7,-36.2,-29.6,-5.8,0.4,7.1};
  const double fluxCM2[]={5.1040518742e11,5.2756309969e11,
      5.4427133660e11,5.2444958674e11,5.4475114410e11,
      5.4492645419e11};
  const double expectedH[]={9.6590204206e-4,1.0721749589e-3,
      1.2024559598e-3,1.3079566347e-3,1.3707810457e-3,
      1.3575134628e-3};
  const double expectedOHEnergyErgS[]={1.16e-15,1.21e-15,1.29e-15,
      9.26e-16,8.12e-16,1.44e-15};

  for (int i=0;i<6;i++) {
    const RF::Conditions conditions={velocityKMS[i]*1.0e3,
        distanceAU[i]*AU,fluxCM2[i]*1.0e4};
    const RF::Result h=RF::Evaluate(
        RF::Model::HydrogenLymanAlphaIKL2020,conditions);
    pass=Check(h.status==RF::Status::Ok &&
        Relative(h.photonRateS1,expectedH[i],1.0e-4),
        "IKL H case-center regression") && pass;

    const RF::Result oh=RF::Evaluate(
        RF::Model::HydroxylAX00SchleicherAhearn1988Unquenched,conditions);
    pass=Check(oh.status==RF::Status::Ok &&
        Relative(oh.emittedPowerW,expectedOHEnergyErgS[i]*1.0e-7,5.0e-3) &&
        std::isnan(oh.photonRateS1),"Lowell OH case-center regression") && pass;
  }

  RF::Conditions conditions={0.0,AU,5.2e15};
  const RF::Result h1=RF::Evaluate(
      RF::Model::HydrogenLymanAlphaIKL2020,conditions);
  conditions.heliocentricDistanceM=2.0*AU;
  const RF::Result h2=RF::Evaluate(
      RF::Model::HydrogenLymanAlphaIKL2020,conditions);
  pass=Check(h1.status==RF::Status::Ok && h2.status==RF::Status::Ok &&
      Relative(h2.photonRateS1,h1.photonRateS1/4.0,2.0e-15),
      "inverse-square scaling") && pass;

  conditions.heliocentricDistanceM=AU;
  conditions.radialVelocityMPS=200001.0;
  pass=Check(RF::Evaluate(RF::Model::HydrogenLymanAlphaIKL2020,conditions).status==
      RF::Status::RadialVelocityOutOfRange,"H extrapolation rejection") && pass;
  conditions.radialVelocityMPS=60001.0;
  pass=Check(RF::Evaluate(
      RF::Model::HydroxylAX00SchleicherAhearn1988Unquenched,conditions).status==
      RF::Status::RadialVelocityOutOfRange,"OH extrapolation rejection") && pass;
  conditions.radialVelocityMPS=0.0;
  conditions.solarLymanAlphaFlux1AUPhotonM2S=
      std::numeric_limits<double>::quiet_NaN();
  pass=Check(RF::Evaluate(RF::Model::HydrogenLymanAlphaIKL2020,conditions).status==
      RF::Status::MissingSolarLymanAlphaFlux,"missing H irradiance rejection") && pass;

  if (!pass) return 1;
  std::puts("PASS resonance-fluorescence core evaluator");
  return 0;
}
