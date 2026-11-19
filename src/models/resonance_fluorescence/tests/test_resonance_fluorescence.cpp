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
        distanceAU[i]*AU,fluxCM2[i]*1.0e4,0.0,0.0};
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

  RF::Conditions conditions={0.0,AU,5.2e15,0.0,0.0};
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

  // ------------------------------------------------------------ O I 130.4 nm
  const char *oxygenNames[]={"O_I_1304_KVB2022_DIRECT_J2",
      "O_I_1304_SUMER2001_DIRECT_J2","O_I_1304_SUMER2001_LYB_BOWEN_J2"};
  const RF::Model oxygenModels[]={RF::Model::OxygenI1304DirectKVB2022GroundJ2,
      RF::Model::OxygenI1304DirectSUMER2001GroundJ2,
      RF::Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2};
  for (int i=0;i<3;i++) {
    RF::Model parsed=RF::Model::HydrogenLymanAlphaIKL2020;
    const RF::ModelMetadata *metadata=RF::GetMetadata(oxygenModels[i]);
    pass=Check(RF::ModelFromName(oxygenNames[i],&parsed) &&
        parsed==oxygenModels[i] && metadata!=nullptr &&
        std::strcmp(metadata->species,"O")==0 &&
        std::strcmp(RF::ModelName(oxygenModels[i]),oxygenNames[i])==0 &&
        metadata->photonRateAvailable && !metadata->requiresSolarLymanAlphaFlux,
        "strict O model-name conversion and metadata") && pass;
  }

  // Killen et al. (2022) 130.217-nm values at 0.352 au taken verbatim from
  // reference/g_factors/o_i_1304_kvb2022/oI_1302_gvalues.txt.  At the
  // reference triplet flux implied by their three tables (8.41157e13
  // photons m^-2 s^-1; SORCE SOLSTICE 2009-03-01 gives 8.1e13 for the same
  // window) the 1-AU rate is the tabulated value times 0.352^2.
  const double kvbReferenceFlux=8.4115654293078703e13;
  const double kvbScale=0.352*0.352;
  RF::Conditions oxygen={0.0,AU,0.0,kvbReferenceFlux,0.0};
  RF::Result o=RF::Evaluate(RF::Model::OxygenI1304DirectKVB2022GroundJ2,oxygen);
  pass=Check(o.status==RF::Status::Ok &&
      Relative(o.photonRateS1,6.572680e-05*kvbScale,1.0e-12) &&
      Relative(o.emittedPowerW,o.photonRateS1*1.5239488142178161e-18,1.0e-12) &&
      std::isnan(o.radiationPressureRatio),"KVB O I line-center value") && pass;
  oxygen.radialVelocityMPS=-39.7e3;
  oxygen.heliocentricDistanceM=1.75*AU;
  oxygen.solarOxygen1304TripletFlux1AUPhotonM2S=1.214*kvbReferenceFlux;
  o=RF::Evaluate(RF::Model::OxygenI1304DirectKVB2022GroundJ2,oxygen);
  pass=Check(o.status==RF::Status::Ok &&
      Relative(o.photonRateS1,2.196117e-07*kvbScale*1.214/(1.75*1.75),1.0e-12),
      "KVB O I flux and inverse-square scaling") && pass;
  oxygen.radialVelocityMPS=-48.0e3;
  pass=Check(RF::Evaluate(RF::Model::OxygenI1304DirectKVB2022GroundJ2,
      oxygen).status==RF::Status::Ok,"KVB O I domain endpoint") && pass;
  oxygen.radialVelocityMPS=-48.001e3;
  pass=Check(RF::Evaluate(RF::Model::OxygenI1304DirectKVB2022GroundJ2,
      oxygen).status==RF::Status::RadialVelocityOutOfRange,
      "KVB O I truncated-profile rejection") && pass;
  oxygen.radialVelocityMPS=48.001e3;
  pass=Check(RF::Evaluate(RF::Model::OxygenI1304DirectKVB2022GroundJ2,
      oxygen).status==RF::Status::RadialVelocityOutOfRange,
      "KVB O I positive-domain rejection") && pass;
  oxygen.radialVelocityMPS=0.0;
  oxygen.solarOxygen1304TripletFlux1AUPhotonM2S=0.0;
  pass=Check(RF::Evaluate(RF::Model::OxygenI1304DirectKVB2022GroundJ2,
      oxygen).status==RF::Status::MissingSolarOxygen1304Flux &&
      RF::Evaluate(RF::Model::OxygenI1304DirectSUMER2001GroundJ2,
      oxygen).status==RF::Status::MissingSolarOxygen1304Flux,
      "missing O I triplet irradiance rejection") && pass;

  // SUMER-atlas tables against an independent NumPy prototype that
  // interpolates the BASS2000 atlas directly (1 AU; triplet flux 1.214 and
  // 1.267 times the KVB reference; Lyman-beta 4.842e13 and 4.961e13).
  const double sumerVelocityKMS[]={-39.7,0.4};
  const double sumerTripletRatio[]={1.214,1.267};
  const double sumerLymanBetaM2[]={4.842e13,4.961e13};
  const double expectedSumerDirect[]={1.749e-6,1.340e-5};
  const double expectedBowen[]={3.909e-7,4.248e-7};
  for (int i=0;i<2;i++) {
    const RF::Conditions sumer={sumerVelocityKMS[i]*1.0e3,AU,0.0,
        sumerTripletRatio[i]*kvbReferenceFlux,sumerLymanBetaM2[i]};
    const RF::Result direct=RF::Evaluate(
        RF::Model::OxygenI1304DirectSUMER2001GroundJ2,sumer);
    const RF::Result bowen=RF::Evaluate(
        RF::Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2,sumer);
    pass=Check(direct.status==RF::Status::Ok &&
        Relative(direct.photonRateS1,expectedSumerDirect[i],3.0e-3),
        "SUMER O I direct prototype regression") && pass;
    pass=Check(bowen.status==RF::Status::Ok &&
        Relative(bowen.photonRateS1,expectedBowen[i],3.0e-3),
        "SUMER O I Bowen prototype regression") && pass;
  }
  RF::Conditions bowenOnly={-99.0e3,2.0*AU,0.0,
      std::numeric_limits<double>::quiet_NaN(),4.9e13};
  const RF::Result far=RF::Evaluate(
      RF::Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2,bowenOnly);
  bowenOnly.heliocentricDistanceM=AU;
  const RF::Result near=RF::Evaluate(
      RF::Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2,bowenOnly);
  pass=Check(far.status==RF::Status::Ok && near.status==RF::Status::Ok &&
      Relative(far.photonRateS1,near.photonRateS1/4.0,2.0e-15),
      "Bowen ignores the triplet flux and scales by inverse square") && pass;
  bowenOnly.radialVelocityMPS=100.001e3;
  pass=Check(RF::Evaluate(RF::Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2,
      bowenOnly).status==RF::Status::RadialVelocityOutOfRange,
      "SUMER atlas window rejection") && pass;
  bowenOnly.radialVelocityMPS=0.0;
  bowenOnly.solarLymanBetaFlux1AUPhotonM2S=-1.0;
  pass=Check(RF::Evaluate(RF::Model::OxygenI1304LymanBetaBowenSUMER2001GroundJ2,
      bowenOnly).status==RF::Status::MissingSolarLymanBetaFlux,
      "missing Lyman-beta irradiance rejection") && pass;

  if (!pass) return 1;
  std::puts("PASS resonance-fluorescence core evaluator");
  return 0;
}
