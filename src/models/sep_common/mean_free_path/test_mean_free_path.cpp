#include "mean_free_path.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace MFP = SEP::MeanFreePath;

namespace {

constexpr double ProtonMassKg = 1.67262192369e-27;
constexpr double ElectronMassKg = 9.1093837015e-31;
constexpr double Pi = 3.14159265358979323846;

struct Tests {
  int pass = 0;
  int fail = 0;
  void Check(bool condition, const std::string& name) {
    if (condition) { ++pass; std::cout << "PASS  " << name << '\n'; }
    else { ++fail; std::cerr << "FAIL  " << name << '\n'; }
  }
  void Near(double actual, double expected, double relative,
            const std::string& name) {
    const double scale = std::max(std::abs(expected), 1.0e-300);
    Check(std::isfinite(actual) && std::abs(actual - expected) <= relative * scale,
          name + " actual=" + std::to_string(actual) +
          " expected=" + std::to_string(expected));
  }
};

MFP::ParticleState AtKineticEnergy(double massKg, double chargeC,
                                   double kineticEV,
                                   const std::string& species) {
  const double restJ = massKg * MFP::SpeedOfLightMPerS * MFP::SpeedOfLightMPerS;
  const double kineticJ = kineticEV * MFP::ElectronVoltJ;
  const double pcJ = std::sqrt(kineticJ * (kineticJ + 2.0 * restJ));
  MFP::ParticleState p;
  p.massKg = massKg;
  p.chargeC = chargeC;
  p.momentumKgMPerS = pcJ / MFP::SpeedOfLightMPerS;
  p.speciesId = species;
  return p;
}

MFP::ParticleState AtRigidity(double massKg, double chargeC, double rigidityV,
                              const std::string& species) {
  MFP::ParticleState p;
  p.massKg = massKg;
  p.chargeC = chargeC;
  p.momentumKgMPerS = rigidityV * std::abs(chargeC) /
                      MFP::SpeedOfLightMPerS;
  p.speciesId = species;
  return p;
}

MFP::Configuration Configure(const std::string& id,
                             std::initializer_list<MFP::InputParameter> p,
                             Tests* tests) {
  MFP::Configuration c;
  const MFP::Status status = MFP::BuildConfiguration(id, p, &c);
  tests->Check(status.ok(), "configure " + id +
               (status.ok() ? "" : ": " + status.detail));
  return c;
}

void TestRegistryAndParser(Tests* t) {
  t->Check(MFP::ModelRegistry().size() >= 50,
           "registry contains the SEP/GCR inventory");
  const MFP::ModelDescriptor* blocked = MFP::FindModel("shock-parasol");
  t->Check(blocked && !blocked->executable &&
           blocked->declaredState == MFP::RuntimeState::RequiresSourceOrCodeAudit,
           "PARASOL remains a typed U-6 source gate");

  MFP::Configuration invalid;
  MFP::Status status = MFP::BuildConfiguration("SHOCK-PARASOL", {}, &invalid);
  t->Check(status.code == MFP::StatusCode::RequiresSourceOrCodeAudit,
           "blocked model reports source-audit status");
  status = MFP::BuildConfiguration("SEP-PATH09",
      {{"lambda0_m", "1"}, {"lambda_kind", "unspecified"},
       {"momentum_variable", "rigidity_v"}, {"x0_SI", "1"},
       {"momentum_exponent", "1/3"}, {"radius0_m", "1"},
       {"radial_exponent", "2/3"}}, &invalid);
  t->Check(status.code == MFP::StatusCode::RequiresUserDecision,
           "lambda_unspecified is refused (U-11)");
  status = MFP::BuildConfiguration("SEP-MFLAMPA25",
      {{"lambda0_m", "1 AU"}, {"lambda_kind", "parallel"},
       {"momentum_variable", "momentum_pc_ev"}, {"x0_SI", "1e9"},
       {"momentum_exponent", "1/3"}, {"radius0_m", "1"},
       {"radial_exponent", "1"}}, &invalid);
  t->Check(status.code == MFP::StatusCode::InvalidConfiguration,
           "SI parser rejects unit suffixes");

  MFP::PublishedScalar scalar;
  status = MFP::ParsePublishedScalar("1/3", true, &scalar);
  t->Check(status.ok() && scalar.kind == MFP::PublishedScalarKind::ExactRational,
           "dimensionless exact rational accepted");
  status = MFP::ParsePublishedScalar("1/3", false, &scalar);
  t->Check(!status.ok(), "dimensional exact rational rejected");
  status = MFP::ParsePublishedScalar("1.5/2", true, &scalar);
  t->Check(!status.ok(), "non-integer rational syntax rejected");
  status = MFP::ParsePublishedScalar("0.3 AU", false, &scalar);
  t->Check(status.ok() && scalar.kind == MFP::PublishedScalarKind::NumberWithUnit &&
           scalar.unit == "AU", "published unit retained rather than guessed");
  status = MFP::ParsePublishedScalar("0.1 to 0.4 AU", false, &scalar);
  t->Check(status.code == MFP::StatusCode::RequiresUserDecision,
           "published range is not replaced by midpoint");
}

void TestKinematicsAndGeometry(Tests* t) {
  MFP::Kinematics k;
  MFP::ParticleState proton = AtKineticEnergy(
      ProtonMassKg, MFP::ElementaryChargeC, 1.0e6, "proton");
  MFP::Status status = MFP::ComputeKinematics(proton, &k);
  t->Check(status.ok(), "1 MeV proton kinematics status");
  t->Near(k.rigidityV / 1.0e6, 43.3306378480744, 1.0e-12,
          "F-KIN-01 proton rigidity");
  t->Near(k.beta, 0.0461321467913892, 1.0e-12,
          "F-KIN-01 proton beta");

  MFP::ParticleState electron = AtKineticEnergy(
      ElectronMassKg, -MFP::ElementaryChargeC, 0.094e6, "electron");
  status = MFP::ComputeKinematics(electron, &k);
  t->Near(k.rigidityV / 1.0e6, 0.323888565094970, 1.0e-12,
          "F-KIN-01 electron rigidity");

  MFP::LambdaValue parallel{2.0, MFP::LambdaKind::Parallel, std::nullopt};
  MFP::LambdaValue radial;
  status = MFP::ParallelToRadialSEP(parallel, Pi / 4.0, &radial);
  t->Check(status.ok(), "parallel to radial SEP status");
  t->Near(radial.metres, 1.0, 1.0e-14, "F-SEP-02 cos-squared projection");
  MFP::LambdaValue back;
  status = MFP::RadialSEPToParallel(radial, Pi / 4.0, &back);
  t->Near(back.metres, 2.0, 1.0e-14, "radial SEP round trip");
  status = MFP::RadialSEPToParallel(radial, Pi / 2.0, &back);
  t->Check(status.code == MFP::StatusCode::OutsideModelDomain,
           "singular radial conversion rejected");
  status = MFP::ParallelToRadialTensor(parallel, 1.0, 3.0, Pi / 4.0, &radial);
  t->Near(radial.metres, 1.5, 1.0e-14,
          "Eq. (3) retains perpendicular tensor term");
}

void TestSEPAndChen(Tests* t) {
  MFP::Configuration c = Configure("SEP-MFLAMPA25",
      {{"lambda0_m", "44879361210"}, {"lambda_kind", "parallel"},
       {"momentum_variable", "momentum_pc_ev"}, {"x0_SI", "1e9"},
       {"momentum_exponent", "1/3"}, {"radius0_m", "149597870700"},
       {"radial_exponent", "1"}}, t);
  MFP::ParticleState proton = AtKineticEnergy(
      ProtonMassKg, MFP::ElementaryChargeC, 10.0e6, "proton");
  MFP::LocalState local; local.radiusM = MFP::AstronomicalUnitM;
  MFP::Result result = MFP::Evaluate(proton, local, c);
  t->Check(result.status.ok() && result.lambda.has_value(),
           "M-FLAMPA25 evaluation status");
  t->Near(result.lambda->metres / MFP::AstronomicalUnitM,
          0.154786263973336, 2.0e-12, "F-SEP-03 power law");
  t->Near(*result.kappaParallelM2PerS, 3.3516433607065e17, 2.0e-12,
          "F-SEP-03 kappa=v lambda/3");

  c = Configure("SEP-CHEN24",
      {{"output_quantity", "lambda_parallel"}, {"species_id", "proton"},
       {"domain_policy", "error"}}, t);
  proton = AtKineticEnergy(ProtonMassKg, MFP::ElementaryChargeC, 1.0e6, "proton");
  local.radiusM = 0.5 * MFP::AstronomicalUnitM;
  result = MFP::Evaluate(proton, local, c);
  t->Check(result.status.ok() && result.provenance.derived,
           "confirmed Chen mapping returns explicitly derived lambda");
  t->Near(*result.kappaParallelM2PerS * 1.0e4, 3.09346072622093e20,
          2.0e-12, "F-SEP-06 Chen kappa");
  t->Near(result.lambda->metres / MFP::AstronomicalUnitM,
          0.0448555391543603, 2.0e-12, "F-SEP-06 Chen proton lambda");
  MFP::ParticleState mismatched = proton; mismatched.speciesId = "ion";
  result = MFP::Evaluate(mismatched, local, c);
  t->Check(result.status.code == MFP::StatusCode::InvalidInput,
           "Chen derived lambda rejects species mismatch");
}

void TestPitchAngleAndFocusing(Tests* t) {
  MFP::LambdaValue lambda{1.0, MFP::LambdaKind::Parallel, std::nullopt};
  double amplitude = 0.0;
  MFP::Status status = MFP::PitchAmplitudeFromLambda(
      "PA-QFORM", 5.0 / 3.0, 0.01, 0.0,
      MFP::OperatorConvention::Standard, lambda, 1.0, &amplitude);
  t->Check(status.ok(), "q-form amplitude status");
  t->Near(amplitude, 1.60199462084936, 2.0e-11,
          "F-PA-01 q-form normalized amplitude");
  MFP::LambdaValue recovered;
  status = MFP::LambdaFromPitchAmplitude("PA-QFORM", 5.0 / 3.0, 0.01, 0.0,
      MFP::OperatorConvention::Standard, amplitude, 1.0, &recovered);
  t->Near(recovered.metres, 1.0, 2.0e-12, "q-form amplitude round trip");

  status = MFP::PitchAmplitudeFromLambda("PA-EPS", 0.0, 0.048, 0.0,
      MFP::OperatorConvention::Standard, lambda, 1.0, &amplitude);
  t->Near(amplitude, 4.597261773672, 2.0e-11,
          "F-PA-02 epsilon exact phi");

  MFP::Configuration c = Configure("PA-EPREM",
      {{"lambda_parallel_m", "1"}, {"operator_convention", "half_d"},
       {"lambda_interpretation", "printed_parameter"}}, t);
  MFP::ParticleState proton = AtKineticEnergy(
      ProtonMassKg, MFP::ElementaryChargeC, 1.0e6, "proton");
  MFP::Result result = MFP::Evaluate(proton, {}, c);
  t->Near(result.lambda->metres, 2.0, 1.0e-13,
          "F-PA-04 EPREM printed parameter gives 2 lambda transport MFP");

  c = Configure("PA-DROGE-VA",
      {{"lambda_parallel_m", "1"}, {"operator_convention", "standard"},
       {"q", "5/3"}, {"alfven_to_particle_speed", "0.1"}}, t);
  result = MFP::Evaluate(proton, {}, c);
  t->Near(result.lambda->metres, 0.5985115566934, 3.0e-11,
          "F-PA-03 Droge effective/nominal lambda");

  double ratio = 0.0;
  status = MFP::HeWanFocusingRatio(1.0e-6, &ratio);
  t->Near(ratio, 0.9999999999996, 2.0e-15,
          "F-NUM-01 stable He-Wan small-x series");
  status = MFP::HeWanFocusingRatio(3.0, &ratio);
  t->Near(ratio, 0.222771694034808, 2.0e-13,
          "F-PA-06 He-Wan direct branch");
}

void TestTurbulenceShockAndGCR(Tests* t) {
  MFP::ParticleState proton = AtRigidity(
      ProtonMassKg, MFP::ElementaryChargeC, 1.0e6, "proton");
  MFP::LocalState local;
  local.meanFieldMagnitudeT = 4.12e-9;
  local.slabVarianceT2 = 4.12e-9 * 4.12e-9;  // fixture reports lambda*(dB/B)^2
  local.kMinPerM = 1.0e-10;
  MFP::Configuration c = Configure("QLT-TS03-P",
      {{"inertial_index", "5/3"}}, t);
  MFP::Result result = MFP::Evaluate(proton, local, c);
  // At 1 MV Eq. (21) includes both asymptotes.  Compare against an independent
  // literal evaluation of the published equation rather than a generated value.
  const double rL = 1.0e6 / (MFP::SpeedOfLightMPerS * 4.12e-9);
  const double R = rL * 1.0e-10;
  const double s = 5.0 / 3.0;
  const double expected = 3.0 * s * rL * rL * 1.0e-10 /
      (4.0 * Pi * (s - 1.0)) *
      (1.0 + 8.0 / ((2.0 - s) * (4.0 - s) * std::pow(R, s)));
  t->Near(result.lambda->metres, expected, 2.0e-13,
          "F-QLT-01 TS2003 Eq. (21)");

  proton = AtKineticEnergy(ProtonMassKg, MFP::ElementaryChargeC, 1.0e6, "proton");
  local = MFP::LocalState{};
  local.meanFieldMagnitudeT = 5.0e-9;
  local.shockSide = MFP::LocalState::ShockSide::Upstream;
  c = Configure("SHOCK-BOHM", {}, t);
  result = MFP::Evaluate(proton, local, c);
  t->Near(result.lambda->metres, 28907090.0163035, 2.0e-12,
          "F-SH-01 Bohm gyroradius");

  proton = AtRigidity(ProtonMassKg, MFP::ElementaryChargeC, 1.0e9, "proton");
  local = MFP::LocalState{};
  local.meanFieldMagnitudeT = 1.0e-9;
  c = Configure("GCR-NWU14",
      {{"K0_m2_per_s", "1e18"}, {"field_reference_T", "1e-9"},
       {"field_normalization", "magnitude"}, {"rigidity0_V", "1e9"},
       {"break_rigidity_V", "4e9"}, {"low_slope", "0.56"},
       {"high_slope", "1.95"}, {"smoothness", "3"}}, t);
  result = MFP::Evaluate(proton, local, c);
  MFP::Kinematics kin; MFP::ComputeKinematics(proton, &kin);
  t->Near(*result.kappaParallelM2PerS, 1.0e18 * kin.beta, 2.0e-13,
          "F-GCR-01 NWU G(1 GV)=1");

  c = Configure("GCR-DUAN25",
      {{"K0_m2_per_s", "1"}, {"equatorial_field_T", "1e-9"},
       {"field_normalization", "magnitude"}, {"break_rigidity_V", "4.3e9"},
       {"a", "0.8"}, {"b", "1.7"}, {"c", "2.2"},
       {"formula_variant", "duan25_printed_d33"}}, t);
  proton = AtRigidity(ProtonMassKg, MFP::ElementaryChargeC, 4.3e9, "proton");
  result = MFP::Evaluate(proton, local, c);
  MFP::ComputeKinematics(proton, &kin);
  t->Near(*result.kappaParallelM2PerS / kin.beta, 4.59479341998814,
          2.0e-13, "F-GCR-08 Duan printed value at break");

  c = Configure("GCR-HELMOD17",
      {{"K0_AU2_per_s_numeric", "0.0003059"}, {"g_low", "0.3"},
       {"normalization_convention", "d21_numeric_use"}}, t);
  local = MFP::LocalState{};
  local.radiusM = MFP::AstronomicalUnitM;
  proton = AtRigidity(ProtonMassKg, MFP::ElementaryChargeC, 1.0e9, "proton");
  result = MFP::Evaluate(proton, local, c);
  t->Check(result.status.ok(), "HelMod does not invent or require a field factor");
  t->Near(*result.kappaParallelM2PerS / (MFP::AstronomicalUnitM *
          MFP::AstronomicalUnitM), 0.000193335542790037, 2.0e-12,
          "F-GCR-05 HelMod numerical-use convention");
}

void TestManagerBatchAndData(const std::string& bundle, Tests* t) {
  MFP::Configuration good = Configure("SHOCK-BOHM", {}, t);
  MFP::Status status = MFP::SetActiveConfiguration(good);
  t->Check(status.ok() && MFP::ActiveModelFunction != nullptr,
           "active dispatch pointer installed");
  const std::string fingerprint = MFP::GetActiveConfiguration().fingerprint;
  status = MFP::ConfigureActiveModel("QLT-ZANK98", {});
  t->Check(!status.ok() && MFP::GetActiveConfiguration().fingerprint == fingerprint,
           "failed reconfiguration is transactional");

  std::vector<MFP::Result> batch;
  MFP::ParticleState p = AtKineticEnergy(
      ProtonMassKg, MFP::ElementaryChargeC, 1.0e6, "proton");
  MFP::LocalState local; local.meanFieldMagnitudeT = 5.0e-9;
  local.shockSide = MFP::LocalState::ShockSide::Upstream;
  status = MFP::EvaluateBatch({p, p}, {local}, good, &batch);
  t->Check(status.code == MFP::StatusCode::InvalidInput && batch.empty(),
           "batch shape error writes no partial output");
  status = MFP::EvaluateBatch({p, p}, {local, local}, good, &batch);
  t->Check(status.ok() && batch.size() == 2 && batch[0].status.ok(),
           "equal-shape batch evaluation");

  MFP::SourceRecord record;
  status = MFP::LoadSourceRecord(bundle, "parameters/model_registry.json",
                                 "SEP-CHEN24", &record);
  t->Check(status.ok() && record.rawJson.find("Chen2024") != std::string::npos &&
           record.manifestSha256.size() == 64,
           "read-only source-key loader retains raw JSON and manifest SHA-256");
}

}  // namespace

int main(int argc, char** argv) {
  std::string bundle = "MEAN_FREE_PATH_MODEL_DATA";
  if (argc == 3 && std::string(argv[1]) == "--bundle") bundle = argv[2];
  else if (argc != 1) {
    std::cerr << "usage: test_mean_free_path [--bundle DIRECTORY]\n";
    return 2;
  }
  Tests tests;
  TestRegistryAndParser(&tests);
  TestKinematicsAndGeometry(&tests);
  TestSEPAndChen(&tests);
  TestPitchAngleAndFocusing(&tests);
  TestTurbulenceShockAndGCR(&tests);
  TestManagerBatchAndData(bundle, &tests);
  std::cout << "SUMMARY pass=" << tests.pass << " fail=" << tests.fail << '\n';
  return tests.fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
