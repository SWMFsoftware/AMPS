#include "parallel_diffusion.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace PD = SEP::ParallelDiffusion;

namespace {

// Tests intentionally use fixed mathematical inputs and independent analytic
// identities rather than expected coefficients saved from the selected model
// evaluator.  The central kinematic helper is first checked against a
// separately constructed fixture before later coefficient tests consume it.
// The proton constants below belong to the specification's Section 15.2
// regression fixture: c and e are exact SI defining constants, while this
// proton mass is the documented CODATA 2018 value.  They are not runtime
// species defaults or a model calibration.
constexpr double ProtonMassKg = 1.67262192369e-27;
constexpr double ElementaryChargeC = 1.602176634e-19;

struct TestRecord {
  std::string id;
  std::string model;
  bool passed = false;
  std::string fixture;
  std::string detail;
};

bool Relative(double actual, double expected, double tolerance) {
  // The max(|expected|,1) scale provides relative comparison for the order-one
  // and larger quantities used here while retaining an absolute floor near
  // zero.  Exact-zero behavior is checked through status/option presence or
  // direct identities rather than using this helper as a physical tolerance.
  return std::fabs(actual - expected) <=
      tolerance * std::max(std::fabs(expected), 1.0);
}

PD::ParticleState ParticleAtRigidity(double rigidityV, double massKg,
                                     double chargeC) {
  // Invert R=pc/|q| analytically.  This fixture helper does not invoke the
  // production kinematics routine, so equal-rigidity species tests are not
  // circular checks of the implementation under test.
  PD::ParticleState particle;
  particle.massKg = massKg;
  particle.chargeC = chargeC;
  particle.momentumKgMPerS =
      rigidityV * std::fabs(chargeC) / PD::SpeedOfLightMPerS;
  return particle;
}

PD::ParticleState TenMeVProton() {
  // This test-side conversion starts from T and uses pc=sqrt(T(T+2mc^2)),
  // independently of the production momentum-to-kinematics implementation.
  const long double electronVoltJ = 1.602176634e-19L;
  const long double kineticJ = 10.0e6L * electronVoltJ;
  const long double massKg = 1.67262192369e-27L;
  const long double c = 299792458.0L;
  const long double pc =
      std::sqrt(kineticJ * (kineticJ + 2.0L * massKg * c * c));
  PD::ParticleState particle;
  particle.massKg = static_cast<double>(massKg);
  particle.chargeC = ElementaryChargeC;
  particle.momentumKgMPerS = static_cast<double>(pc / c);
  particle.nucleonCount = 1.0;
  return particle;
}

void Add(std::vector<TestRecord>* tests, const std::string& id,
         const std::string& model, bool passed, const std::string& fixture,
         const std::string& detail) {
  tests->push_back(TestRecord{id, model, passed, fixture, detail});
}

std::string JsonEscape(const std::string& text) {
  std::string escaped;
  for (char c : text) {
    if (c == '"' || c == '\\') escaped.push_back('\\');
    if (c == '\n') escaped += "\\n";
    else escaped.push_back(c);
  }
  return escaped;
}

bool TestKinematics(std::vector<TestRecord>* tests) {
  PD::ParticleKinematics value;
  const PD::Status status =
      PD::ComputeParticleKinematics(TenMeVProton(), &value);
  const bool pass = status.ok() &&
      Relative(value.gamma, 1.01065788925, 5.0e-12) &&
      Relative(value.beta, 0.144844004126, 5.0e-12) &&
      Relative(value.speedMPerS, 4.34231400235e7, 5.0e-12) &&
      Relative(value.rigidityV, 1.37351526250e8, 5.0e-12);
  Add(tests, "PD01-UNITS", "core", pass,
      "Section 15.2 printed 10 MeV proton fixture; p computed independently from T",
      pass ? "relativistic speed, gamma, beta, and rigidity agree"
           : "relativistic conversion differs from the Section 15.2 fixture");

  PD::ParticleState neutral = TenMeVProton();
  neutral.chargeC = 0.0;
  const PD::Status neutralStatus =
      PD::ComputeParticleKinematics(neutral, &value);
  const bool invalidPass =
      neutralStatus.code == PD::StatusCode::InvalidParticle;
  Add(tests, "PD01-INVALID-PARTICLE", "core", invalidPass,
      "Section 14.2 charged-particle domain",
      invalidPass ? "neutral particle is rejected before coefficient evaluation"
                  : "neutral particle entered the charged-particle API");
  return pass && invalidPass;
}

bool TestConstants(std::vector<TestRecord>* tests) {
  // These tests independently apply Equation (4) in both directions.  The
  // negative control ensures constant_lambda and constant_kappa have not been
  // collapsed into one superficially similar configuration.
  const PD::ParticleState particle = TenMeVProton();
  PD::ParticleKinematics kinematics;
  PD::ComputeParticleKinematics(particle, &kinematics);
  PD::LocalState local;

  PD::ModelConfiguration lambda;
  lambda.model = PD::ModelId::ConstantLambda;
  lambda.constantLambda.lambdaParallelM = 0.1 * 149597870700.0;
  const PD::ParallelResult lambdaResult = PD::Evaluate(particle, local, lambda);
  const double expectedKappa =
      kinematics.speedMPerS * lambda.constantLambda.lambdaParallelM / 3.0;
  const bool lambdaPass = lambdaResult.status.ok() &&
      lambdaResult.lambdaParallelM.has_value() &&
      lambdaResult.kappaParallelM2PerS.has_value() &&
      Relative(*lambdaResult.lambdaParallelM,
               lambda.constantLambda.lambdaParallelM, 1.0e-15) &&
      Relative(*lambdaResult.kappaParallelM2PerS, expectedKappa, 1.0e-15) &&
      lambdaResult.gradKappaParallelMPerS.has_value();
  Add(tests, "PD02-CONSTANT-LAMBDA", "constant_lambda", lambdaPass,
      "Equation (11) direct identity",
      lambdaPass ? "lambda is constant and kappa=v*lambda/3"
                 : "constant-lambda identity failed");

  PD::ModelConfiguration kappa;
  kappa.model = PD::ModelId::ConstantKappa;
  kappa.constantKappa.kappaParallelM2PerS = 7.5e17;
  const PD::ParallelResult kappaResult = PD::Evaluate(particle, local, kappa);
  const double expectedLambda =
      3.0 * kappa.constantKappa.kappaParallelM2PerS /
      kinematics.speedMPerS;
  const bool kappaPass = kappaResult.status.ok() &&
      Relative(*kappaResult.kappaParallelM2PerS,
               kappa.constantKappa.kappaParallelM2PerS, 1.0e-15) &&
      Relative(*kappaResult.lambdaParallelM, expectedLambda, 1.0e-15) &&
      !Relative(*kappaResult.kappaParallelM2PerS, expectedKappa, 1.0e-6);
  Add(tests, "PD02-CONSTANT-KAPPA", "constant_kappa", kappaPass,
      "Equation (12) direct identity and distinct-model negative control",
      kappaPass ? "kappa is constant and lambda=3*kappa/v"
                : "constant-kappa identity or distinct-model control failed");
  return lambdaPass && kappaPass;
}

bool TestPowerLaw(std::vector<TestRecord>* tests) {
  // Ratios are chosen so the expected products can be evaluated directly.
  // This exercises normalization and separability without importing an
  // observationally calibrated value or duplicating the production log path.
  const double rigidity0 = 2.0e8;
  PD::ParticleState particle =
      ParticleAtRigidity(rigidity0, ProtonMassKg, ElementaryChargeC);
  PD::LocalState local;
  local.positionM = {{149597870700.0, 0.0, 0.0}};

  PD::ModelConfiguration configuration;
  configuration.model = PD::ModelId::PowerLawLambda;
  PD::PowerLawLambdaParameters& p = configuration.powerLawLambda;
  p.lambda0M = 1.3e10;
  p.independentVariable = PD::IndependentVariable::Rigidity;
  p.independentReferenceSI = rigidity0;
  p.independentExponent = 1.0 / 3.0;
  p.useRadialFactor = true;
  p.radius0M = local.positionM[0];
  p.radialExponent = 0.5;
  const PD::ParallelResult atReference =
      PD::Evaluate(particle, local, configuration);
  PD::ParticleKinematics kinematics;
  PD::ComputeParticleKinematics(particle, &kinematics);
  const bool referencePass = atReference.status.ok() &&
      Relative(*atReference.lambdaParallelM, p.lambda0M, 1.0e-14) &&
      Relative(*atReference.dLnKappaDLnRigidity,
               p.independentExponent +
                   1.0 / (kinematics.gamma * kinematics.gamma),
               1.0e-14);
  Add(tests, "PD02-POWER-REFERENCE", "power_law_lambda", referencePass,
      "Equations (13),(16), ratios equal unity",
      referencePass ? "reference normalization and beta slope agree"
                    : "power-law reference normalization or slope failed");

  const double scale = 4.0;
  particle = ParticleAtRigidity(scale * rigidity0, ProtonMassKg,
                                ElementaryChargeC);
  local.positionM[0] = 9.0 * p.radius0M;
  const PD::ParallelResult scaled = PD::Evaluate(particle, local, configuration);
  const double expectedLambda = p.lambda0M *
      std::pow(scale, p.independentExponent) *
      std::pow(9.0, p.radialExponent);
  const bool scalePass = scaled.status.ok() &&
      Relative(*scaled.lambdaParallelM, expectedLambda, 1.0e-14);
  Add(tests, "PD02-POWER-SCALING", "power_law_lambda", scalePass,
      "Equation (13) independently multiplied rigidity and radial ratios",
      scalePass ? "separable rigidity/radial scaling agrees"
                : "separable power-law scaling failed");

  PD::ModelConfiguration perNucleon = configuration;
  perNucleon.powerLawLambda.useRadialFactor = false;
  perNucleon.powerLawLambda.independentVariable =
      PD::IndependentVariable::EnergyPerNucleon;
  perNucleon.powerLawLambda.independentReferenceSI = 1.0e-12;
  particle.nucleonCount.reset();
  const PD::ParallelResult missing = PD::Evaluate(particle, local, perNucleon);
  const bool missingPass =
      missing.status.code == PD::StatusCode::MissingInput;
  Add(tests, "PD02-POWER-SPECIES", "power_law_lambda", missingPass,
      "Section 2.1 no-inferred-mass-number rule",
      missingPass ? "energy per nucleon fails without nucleon count"
                  : "energy per nucleon silently inferred species data");

  particle = TenMeVProton();
  PD::ParticleKinematics variables;
  PD::ComputeParticleKinematics(particle, &variables);
  PD::ModelConfiguration variant = configuration;
  PD::PowerLawLambdaParameters& variantParameters = variant.powerLawLambda;
  variantParameters.useRadialFactor = false;
  variantParameters.lambda0M = 6.2e9;
  variantParameters.independentExponent = -0.4;
  bool variantsPass = true;
  const PD::IndependentVariable choices[] = {
      PD::IndependentVariable::Rigidity,
      PD::IndependentVariable::TotalKineticEnergy,
      PD::IndependentVariable::EnergyPerNucleon,
      PD::IndependentVariable::Speed};
  const double references[] = {
      variables.rigidityV, variables.kineticEnergyJ,
      variables.kineticEnergyJ, variables.speedMPerS};
  for (std::size_t i = 0; i < 4; ++i) {
    variantParameters.independentVariable = choices[i];
    variantParameters.independentReferenceSI = references[i];
    const PD::ParallelResult evaluated = PD::Evaluate(particle, local, variant);
    variantsPass = variantsPass && evaluated.status.ok() &&
        Relative(*evaluated.lambdaParallelM, variantParameters.lambda0M,
                 2.0e-14);
  }
  Add(tests, "PD02-POWER-VARIANTS", "power_law_lambda", variantsPass,
      "Equations (13)--(15), all four variables at their reference state",
      variantsPass ? "rigidity, total-energy, per-nucleon, and speed variants agree"
                   : "one or more explicit independent-variable variants failed");

  variantParameters.independentVariable = PD::IndependentVariable::Rigidity;
  variantParameters.independentReferenceSI = variables.rigidityV;
  variantParameters.useFieldFactor = true;
  variantParameters.fieldExponent = 0.0;
  variantParameters.fieldReferenceT = 0.0;
  const PD::ParallelResult zeroFieldExponent =
      PD::Evaluate(particle, PD::LocalState(), variant);
  const bool disabledFieldPass = zeroFieldExponent.status.ok();
  Add(tests, "PD02-POWER-DISABLED-FIELD", "power_law_lambda",
      disabledFieldPass, "Equation (13) eta=0 exact unity factor",
      disabledFieldPass ? "zero field exponent introduces no B requirement"
                        : "disabled field factor incorrectly required B");

  variantParameters.useFieldFactor = false;
  variantParameters.useTimeFactor = true;
  const PD::ParallelResult missingTime =
      PD::Evaluate(particle, PD::LocalState(), variant);
  const bool optionalInputPass =
      missingTime.status.code == PD::StatusCode::MissingInput;
  Add(tests, "PD02-POWER-OPTIONAL-INPUT", "power_law_lambda",
      optionalInputPass, "Equation (13) enabled time-factor contract",
      optionalInputPass ? "enabled factor is required and never defaulted"
                        : "enabled time factor was silently defaulted");
  return referencePass && scalePass && missingPass && variantsPass &&
      disabledFieldPass && optionalInputPass;
}

bool TestBrokenRigidity(std::vector<TestRecord>* tests) {
  // The reference, equal-slope, asymptotic, derivative, and equal-rigidity
  // species checks are independent consequences of Equations (17)--(19).
  // Together they test normalization, transition shape, beta ownership, and
  // the distinction between lambda and kappa rather than one sample value.
  const double rigidity0 = 1.0e9;
  PD::ParticleState particle =
      ParticleAtRigidity(rigidity0, ProtonMassKg, ElementaryChargeC);
  PD::LocalState local;
  PD::ModelConfiguration configuration;
  configuration.model = PD::ModelId::BrokenRigidityKappa;
  PD::BrokenRigidityKappaParameters& p = configuration.brokenRigidityKappa;
  p.kStarM2PerS = 2.0e18;
  p.rigidity0V = rigidity0;
  p.breakRigidityV = 3.0e9;
  p.lowSlope = 0.25;
  p.highSlope = 1.75;
  p.smoothness = 2.5;
  const PD::ParallelResult reference = PD::Evaluate(particle, local, configuration);
  PD::ParticleKinematics kinematics;
  PD::ComputeParticleKinematics(particle, &kinematics);
  const bool referencePass = reference.status.ok() &&
      Relative(*reference.kappaParallelM2PerS,
               p.kStarM2PerS * kinematics.beta, 2.0e-15) &&
      Relative(*reference.lambdaParallelM,
               3.0 * p.kStarM2PerS / PD::SpeedOfLightMPerS, 2.0e-15);
  Add(tests, "PD02-BROKEN-REFERENCE", "broken_rigidity_kappa", referencePass,
      "Equations (17),(18), H(R0)=1",
      referencePass ? "speed-factored K_star normalization agrees"
                    : "broken-law reference normalization failed");

  PD::ModelConfiguration single = configuration;
  single.brokenRigidityKappa.highSlope =
      single.brokenRigidityKappa.lowSlope;
  const double testRigidity = 7.0e8;
  particle = ParticleAtRigidity(testRigidity, ProtonMassKg, ElementaryChargeC);
  const PD::ParallelResult singleResult = PD::Evaluate(particle, local, single);
  const double expectedLambda = 3.0 * p.kStarM2PerS /
      PD::SpeedOfLightMPerS *
      std::pow(testRigidity / rigidity0, p.lowSlope);
  const bool singlePass = singleResult.status.ok() &&
      Relative(*singleResult.lambdaParallelM, expectedLambda, 2.0e-14);
  Add(tests, "PD02-BROKEN-SINGLE", "broken_rigidity_kappa", singlePass,
      "Equation (19) exact a=b limit independent of break",
      singlePass ? "a=b reduces to one power law"
                 : "a=b broken-law limit failed");

  const double h = 1.0e-3;
  const double centerRigidity = 2.2e9;
  const PD::ParallelResult coarseMinus = PD::Evaluate(
      ParticleAtRigidity(centerRigidity * std::exp(-h), ProtonMassKg,
                         ElementaryChargeC), local, configuration);
  const PD::ParallelResult center = PD::Evaluate(
      ParticleAtRigidity(centerRigidity, ProtonMassKg, ElementaryChargeC),
      local, configuration);
  const PD::ParallelResult coarsePlus = PD::Evaluate(
      ParticleAtRigidity(centerRigidity * std::exp(h), ProtonMassKg,
                         ElementaryChargeC), local, configuration);
  const PD::ParallelResult fineMinus = PD::Evaluate(
      ParticleAtRigidity(centerRigidity * std::exp(-0.5 * h), ProtonMassKg,
                         ElementaryChargeC), local, configuration);
  const PD::ParallelResult finePlus = PD::Evaluate(
      ParticleAtRigidity(centerRigidity * std::exp(0.5 * h), ProtonMassKg,
                         ElementaryChargeC), local, configuration);
  const double coarseDifference =
      (std::log(*coarsePlus.kappaParallelM2PerS) -
       std::log(*coarseMinus.kappaParallelM2PerS)) / (2.0 * h);
  const double fineDifference =
      (std::log(*finePlus.kappaParallelM2PerS) -
       std::log(*fineMinus.kappaParallelM2PerS)) / h;
  // Symmetric differences have O(h^2) error. Richardson cancellation gives
  // an independent O(h^4) reference without copying Equation (19).
  const double finiteDifference =
      (4.0 * fineDifference - coarseDifference) / 3.0;
  const bool derivativePass = coarseMinus.status.ok() && center.status.ok() &&
      coarsePlus.status.ok() && fineMinus.status.ok() && finePlus.status.ok() &&
      Relative(*center.dLnKappaDLnRigidity, finiteDifference, 2.0e-11);
  Add(tests, "PD02-BROKEN-DERIVATIVE", "broken_rigidity_kappa",
      derivativePass, "central log-rigidity finite difference independent of analytic derivative",
      derivativePass ? "Equation (19), including beta, agrees under refinement"
                     : "broken-law logarithmic derivative failed");

  const PD::ParallelResult low = PD::Evaluate(
      ParticleAtRigidity(p.breakRigidityV * 1.0e-6, ProtonMassKg,
                         ElementaryChargeC), local, configuration);
  const PD::ParallelResult high = PD::Evaluate(
      ParticleAtRigidity(p.breakRigidityV * 1.0e6, ProtonMassKg,
                         ElementaryChargeC), local, configuration);
  const bool asymptotePass = low.status.ok() && high.status.ok() &&
      Relative(*low.dLnLambdaDLnRigidity, p.lowSlope, 2.0e-14) &&
      Relative(*high.dLnLambdaDLnRigidity, p.highSlope, 2.0e-14);
  Add(tests, "PD02-BROKEN-ASYMPTOTES", "broken_rigidity_kappa",
      asymptotePass, "Equation (19) at R/Rb=1e-6 and 1e6",
      asymptotePass ? "low/high logarithmic slopes reach a and b"
                    : "broken-law asymptotic slopes failed");

  const PD::ParticleState proton =
      ParticleAtRigidity(centerRigidity, ProtonMassKg, ElementaryChargeC);
  const PD::ParticleState alpha =
      ParticleAtRigidity(centerRigidity, 4.0 * ProtonMassKg,
                         2.0 * ElementaryChargeC);
  const PD::ParallelResult protonResult = PD::Evaluate(proton, local, configuration);
  const PD::ParallelResult alphaResult = PD::Evaluate(alpha, local, configuration);
  const bool speciesPass = protonResult.status.ok() && alphaResult.status.ok() &&
      Relative(*protonResult.lambdaParallelM,
               *alphaResult.lambdaParallelM, 2.0e-15) &&
      !Relative(*protonResult.kappaParallelM2PerS,
                *alphaResult.kappaParallelM2PerS, 1.0e-6);
  Add(tests, "PD02-BROKEN-SPECIES", "broken_rigidity_kappa", speciesPass,
      "Equation (18) equal-rigidity proton/alpha identity",
      speciesPass ? "lambda agrees and speed-dependent kappa differs"
                  : "equal-rigidity species behavior failed");
  return referencePass && singlePass && derivativePass && asymptotePass &&
      speciesPass;
}

bool TestBohm(std::vector<TestRecord>* tests) {
  // Evaluate r_L=p/(|q|B) directly from the particle fixture for both field
  // conventions.  The missing-field case proves field validation remains
  // selective instead of installing an undocumented B value.
  const PD::ParticleState particle = TenMeVProton();
  PD::LocalState local;
  local.meanFieldT = std::array<double, 3>{{0.0, 0.0, 5.0e-9}};
  PD::ModelConfiguration configuration;
  configuration.model = PD::ModelId::Bohm;
  configuration.bohm.etaB = 2.5;
  configuration.bohm.fieldDefinition = PD::BohmFieldDefinition::MeanField;
  const PD::ParallelResult result = PD::Evaluate(particle, local, configuration);
  const double expectedLambda = 2.5 * particle.momentumKgMPerS /
      (std::fabs(particle.chargeC) * 5.0e-9);
  const bool pass = result.status.ok() &&
      Relative(*result.lambdaParallelM, expectedLambda, 2.0e-15);
  Add(tests, "PD02-BOHM", "bohm", pass,
      "Equation (61) direct gyroradius identity",
      pass ? "configured mean-field Bohm coefficient agrees"
           : "Bohm coefficient failed");

  PD::LocalState missing;
  const PD::ParallelResult rejected =
      PD::Evaluate(particle, missing, configuration);
  const bool missingPass =
      rejected.status.code == PD::StatusCode::MissingInput;
  Add(tests, "PD02-SELECTIVE-FIELD", "bohm", missingPass,
      "Section 14.2 selective input validation",
      missingPass ? "Bohm requires B while constant models do not"
                  : "Bohm accepted an absent configured field");

  configuration.bohm.fieldDefinition =
      PD::BohmFieldDefinition::EffectiveField;
  PD::LocalState effective;
  effective.effectiveFieldMagnitudeT = 8.0e-9;
  const PD::ParallelResult effectiveResult =
      PD::Evaluate(particle, effective, configuration);
  const double expectedEffective = 2.5 * particle.momentumKgMPerS /
      (std::fabs(particle.chargeC) * 8.0e-9);
  const bool effectivePass = effectiveResult.status.ok() &&
      Relative(*effectiveResult.lambdaParallelM, expectedEffective, 2.0e-15);
  Add(tests, "PD02-BOHM-FIELD-DEFINITION", "bohm", effectivePass,
      "Equation (61) explicit effective-field variant",
      effectivePass ? "effective field is selected without relabeling mean B"
                    : "effective-field Bohm selection failed");
  return pass && missingPass && effectivePass;
}

bool TestParserAndDispatch(std::vector<TestRecord>* tests) {
  // Parser tests treat selection as a transaction.  The function pointer and
  // configuration must remain matched after a later candidate fails, because
  // a mixed pair could evaluate valid arithmetic for the wrong physical law.
  const std::vector<PD::InputParameter> lambdaParameters = {
      {"lambda_parallel_m", "4.5e9"}};
  const PD::Status selected =
      PD::ConfigureActiveModel("constant_lambda", lambdaParameters);
  PD::ModelFunction pointer = PD::ActiveModelFunction;
  const PD::ParallelResult active =
      PD::EvaluateActive(TenMeVProton(), PD::LocalState());
  const std::vector<PD::InputParameter> invalidParameters = {
      {"kappa_parallel_m2_per_s", "-1"}};
  const PD::Status invalid =
      PD::ConfigureActiveModel("constant_kappa", invalidParameters);
  const PD::ParallelResult afterFailure =
      PD::EvaluateActive(TenMeVProton(), PD::LocalState());
  const bool transactionPass = selected.ok() && active.status.ok() &&
      !invalid.ok() && pointer == PD::ActiveModelFunction &&
      afterFailure.provenance.evaluatedModelId == "constant_lambda" &&
      Relative(*afterFailure.lambdaParallelM, 4.5e9, 1.0e-15);
  Add(tests, "PD01-DISPATCH", "registry", transactionPass,
      "parser-neutral configuration and transactional function-pointer selection",
      transactionPass ? "validated selection updates pointer; rejected input preserves it"
                      : "active dispatch transaction failed");

  PD::ModelConfiguration polynomial;
  const PD::Status polynomialStatus = PD::BuildConfiguration(
      "nlgce_f_2014", std::vector<PD::InputParameter>(), &polynomial);
  const bool polynomialPass = polynomialStatus.ok() &&
      polynomial.model == PD::ModelId::NlgceF2014;
  Add(tests, "PD05-PARSER", "nlgce_f_2014", polynomialPass,
      "fixed published coefficient-set schema",
      polynomialPass ? "published NLGCE-F backend is selectable"
                     : "implemented NLGCE-F backend was not selectable");

  PD::ModelConfiguration duplicate;
  const PD::Status duplicateStatus = PD::BuildConfiguration(
      "constant_lambda",
      {{"lambda_parallel_m", "1"}, {"lambda_parallel_m", "2"}},
      &duplicate);
  const bool duplicatePass =
      duplicateStatus.code == PD::StatusCode::InvalidConfiguration;
  Add(tests, "PD01-SCHEMA", "registry", duplicatePass,
      "duplicate-key fail-closed schema",
      duplicatePass ? "duplicate parameter rejected"
                    : "duplicate parameter silently replaced");

  PD::ModelConfiguration parsedPower;
  const PD::Status parsedPowerStatus = PD::BuildConfiguration(
      "power_law_lambda",
      {{"lambda0_m", "9e9"},
       {"independent_variable", "rigidity"},
       {"rigidity0_V", "2e8"},
       {"independent_exponent", "0.5"},
       {"radius0_m", "1e11"},
       {"radial_exponent", "1.0"},
       {"use_time_factor", "false"},
       {"use_region_factor", "true"}},
      &parsedPower);
  PD::LocalState parsedLocal;
  parsedLocal.positionM = {{2.0e11, 0.0, 0.0}};
  parsedLocal.regionFactor = 0.25;
  const PD::ParallelResult parsedResult = parsedPowerStatus.ok()
      ? PD::Evaluate(ParticleAtRigidity(8.0e8, ProtonMassKg,
                                       ElementaryChargeC),
                     parsedLocal, parsedPower)
      : PD::ParallelResult();
  // 9e9 * (8e8/2e8)^0.5 * (2e11/1e11)^1 * 0.25 = 9e9.
  const bool parserPass = parsedPowerStatus.ok() && parsedResult.status.ok() &&
      Relative(*parsedResult.lambdaParallelM, 9.0e9, 2.0e-14);
  Add(tests, "PD01-PARSER-BRIDGE", "power_law_lambda", parserPass,
      "complete string-to-typed schema with radial and region factors",
      parserPass ? "parser bridge preserves every supplied factor"
                 : "parser bridge lost or misinterpreted a model parameter");

  const std::vector<PD::ModelDescriptor>& registry = PD::ModelRegistry();
  std::size_t implemented = 0;
  bool unique = registry.size() == 16;
  for (std::size_t i = 0; i < registry.size(); ++i) {
    if (registry[i].implemented) ++implemented;
    for (std::size_t j = i + 1; j < registry.size(); ++j)
      unique = unique && std::string(registry[i].stableId) != registry[j].stableId;
  }
  PD::ModelConfiguration first;
  first.model = PD::ModelId::ConstantLambda;
  first.constantLambda.lambdaParallelM = 1.0;
  PD::ModelConfiguration second = first;
  second.constantLambda.lambdaParallelM = 2.0;
  const std::string firstFingerprint = PD::ConfigurationFingerprint(first);
  const bool shaShape = firstFingerprint.size() == 64 &&
      firstFingerprint.find_first_not_of("0123456789abcdef") ==
          std::string::npos &&
      firstFingerprint ==
          "a6b8179215917864f8e8913fd66f626dd4583ed8de2c67d6f1e100ce065a14ba";
  const bool registryPass = unique && implemented == 16 && shaShape &&
      firstFingerprint != PD::ConfigurationFingerprint(second);
  Add(tests, "PD01-REGISTRY-PROVENANCE", "registry", registryPass,
      "Section 3 complete stable-ID inventory and configuration identity",
      registryPass ? "16 unique IDs, all first-release backends implemented"
                   : "registry inventory or fingerprint sensitivity failed");
  return transactionPass && polynomialPass && duplicatePass && parserPass &&
      registryPass;
}

std::vector<PD::SectionLine> Section(
    const std::vector<std::string>& text, std::size_t firstLine) {
  // Number the fixture lines exactly as a host would number file lines, so
  // diagnostics can be checked against the line that contains the defect.
  std::vector<PD::SectionLine> lines;
  for (std::size_t i = 0; i < text.size(); ++i)
    lines.push_back(PD::SectionLine{firstLine + i, text[i]});
  return lines;
}

bool Contains(const std::string& text, const std::string& fragment) {
  return text.find(fragment) != std::string::npos;
}

double SpeedFromMomentumIndependent(const PD::ParticleState& particle) {
  // v = p c^2 / E with E = sqrt((pc)^2 + (mc^2)^2).  Evaluated here in long
  // double, separately from the production kinematics helper.
  const long double c = PD::SpeedOfLightMPerS;
  const long double pc = static_cast<long double>(particle.momentumKgMPerS) * c;
  const long double rest = static_cast<long double>(particle.massKg) * c * c;
  return static_cast<double>(pc * c / std::sqrt(pc * pc + rest * rest));
}

bool TestUnconfiguredBoundDispatch(std::vector<TestRecord>* tests) {
  // Must run before any test installs a model: the mover-facing pointer has
  // to fail closed, with no coefficient, until a section has been accepted.
  const PD::ParallelResult result =
      PD::ActiveParallelDiffusion(TenMeVProton(), PD::LocalState());
  const bool pass = !PD::HasActiveConfiguration() &&
      result.status.code == PD::StatusCode::InvalidConfiguration &&
      !result.kappaParallelM2PerS.has_value() &&
      !result.lambdaParallelM.has_value();
  Add(tests, "PD13-UNCONFIGURED", "registry", pass,
      "mover-facing pointer before any configuration",
      pass ? "unconfigured pointer returns InvalidConfiguration and no value"
           : "unconfigured pointer produced a coefficient or wrong status");
  return pass;
}

bool TestSectionParser(std::vector<TestRecord>* tests) {
  // Grammar: '#' and '!' comments, blank lines, CRLF endings, blanks around
  // '=', case-insensitive model ID, case-sensitive parameter keys.  The
  // expected configuration is built field by field, not by the parser.
  const std::vector<PD::SectionLine> valid = Section({
      "# srcSEP3D-style comment",
      "",
      "  model   =   CONSTANT_KAPPA   ! Fortran-style comment",
      "kappa_parallel_m2_per_s = 1.25e18 # trailing comment\r",
      "   ! comment-only line"}, 40);
  PD::ParsedSection parsed;
  const PD::Status status = PD::ParseSection(valid, &parsed);
  PD::ModelConfiguration expected;
  expected.model = PD::ModelId::ConstantKappa;
  expected.constantKappa.kappaParallelM2PerS = 1.25e18;
  const bool grammarPass = status.ok() &&
      parsed.modelId == "CONSTANT_KAPPA" && parsed.modelLineNumber == 42 &&
      parsed.parameters.size() == 1 &&
      parsed.parameters[0].name == "kappa_parallel_m2_per_s" &&
      parsed.parameters[0].value == "1.25e18" &&
      parsed.parameters[0].lineNumber == 43 &&
      parsed.configuration.model == PD::ModelId::ConstantKappa &&
      parsed.configuration.constantKappa.kappaParallelM2PerS == 1.25e18 &&
      parsed.configurationFingerprint ==
          PD::ConfigurationFingerprint(expected) &&
      !PD::HasActiveConfiguration();
  Add(tests, "PD13-SECTION-GRAMMAR", "constant_kappa", grammarPass,
      "comments, blanks, CRLF, case rules, and line provenance",
      grammarPass ? "section parsed to the independently built configuration"
                  : "section grammar or provenance failed: " + status.detail);

  // Each defect must fail with its typed status, name the defective line,
  // and leave the caller's output untouched (transactional parse).
  struct ErrorCase {
    const char* name;
    std::vector<std::string> text;
    PD::StatusCode code;
    const char* fragment;
  };
  const std::vector<ErrorCase> cases = {
      {"missing model", {"kappa_parallel_m2_per_s = 1e18"},
       PD::StatusCode::MissingInput, "model = <id>"},
      {"no separator", {"model = constant_kappa", "kappa_parallel_m2_per_s"},
       PD::StatusCode::InvalidConfiguration, "line 11:"},
      {"empty value", {"model = constant_kappa", "kappa_parallel_m2_per_s ="},
       PD::StatusCode::InvalidConfiguration, "line 11:"},
      {"duplicate model", {"model = constant_kappa", "model = constant_lambda"},
       PD::StatusCode::InvalidConfiguration, "first assigned at line 10"},
      {"duplicate parameter",
       {"model = constant_kappa", "kappa_parallel_m2_per_s = 1e18",
        "kappa_parallel_m2_per_s = 2e18"},
       PD::StatusCode::InvalidConfiguration, "line 12:"},
      {"other model's key",
       {"model = constant_kappa", "kappa_parallel_m2_per_s = 1e18",
        "lambda_parallel_m = 1e9"},
       PD::StatusCode::InvalidConfiguration, "line 12: unknown parameter"},
      {"key case", {"model = constant_kappa", "KAPPA_PARALLEL_M2_PER_S = 1e18"},
       PD::StatusCode::MissingInput, "line 10:"},
      {"unknown model", {"! header", "model = no_such_model"},
       PD::StatusCode::UnsupportedModel, "line 11:"},
      {"missing required key", {"model = constant_lambda"},
       PD::StatusCode::MissingInput, "line 10: missing parameter"},
      {"malformed value", {"model = constant_lambda", "lambda_parallel_m = 1e9m"},
       PD::StatusCode::InvalidConfiguration, "line 11:"},
      {"comment hides value", {"model = constant_lambda", "lambda_parallel_m = !1e9"},
       PD::StatusCode::InvalidConfiguration, "line 11:"},
      {"negative value", {"model = constant_kappa", "kappa_parallel_m2_per_s = -1"},
       PD::StatusCode::InvalidConfiguration, "line 11:"}};
  bool errorsPass = true;
  std::string failedCase;
  for (const ErrorCase& item : cases) {
    PD::ParsedSection untouched = parsed;
    const PD::Status rejected =
        PD::ParseSection(Section(item.text, 10), &untouched);
    const bool ok = rejected.code == item.code &&
        Contains(rejected.detail, item.fragment) &&
        untouched.configurationFingerprint == parsed.configurationFingerprint;
    if (!ok && failedCase.empty())
      failedCase = std::string(item.name) + " -> " + rejected.detail;
    errorsPass = errorsPass && ok;
  }
  Add(tests, "PD13-SECTION-ERRORS", "registry", errorsPass,
      "typed, line-attributed, transactional section rejection",
      errorsPass ? "12 malformed sections rejected with their line numbers"
                 : "section defect not rejected as specified: " + failedCase);
  return grammarPass && errorsPass;
}

bool TestBoundDispatch(std::vector<TestRecord>* tests) {
  // Registry: one distinct bound evaluator per stable model identity.
  bool registryPass = true;
  std::vector<PD::ActiveModelEvaluator> seen;
  for (const PD::ModelDescriptor& descriptor : PD::ModelRegistry()) {
    const PD::ActiveModelEvaluator bound =
        PD::BoundFunctionForModel(descriptor.id);
    for (const PD::ActiveModelEvaluator previous : seen)
      registryPass = registryPass && previous != bound;
    registryPass = registryPass && bound != nullptr;
    seen.push_back(bound);
  }
  registryPass = registryPass && seen.size() == 16;
  Add(tests, "PD13-BOUND-REGISTRY", "registry", registryPass,
      "16 distinct model-specific mover-facing evaluators",
      registryPass ? "every stable ID has its own bound evaluator"
                   : "bound evaluator missing or shared between models");

  // Selection through a section installs the model-specific pointer; the
  // value is checked against kappa=v*lambda/3 with v computed independently.
  const PD::ParticleState proton = TenMeVProton();
  const double speed = SpeedFromMomentumIndependent(proton);
  PD::ParsedSection installed;
  const PD::Status lambdaStatus = PD::ConfigureActiveModelFromSection(
      Section({"model = constant_lambda", "lambda_parallel_m = 4.5e9"}, 1),
      &installed);
  const PD::ActiveModelEvaluator lambdaPointer = PD::ActiveParallelDiffusion;
  const PD::ParallelResult lambdaResult =
      PD::ActiveParallelDiffusion(proton, PD::LocalState());
  const PD::ParallelResult viaEvaluateActive =
      PD::EvaluateActive(proton, PD::LocalState());
  const bool lambdaPass = lambdaStatus.ok() && PD::HasActiveConfiguration() &&
      lambdaPointer == PD::BoundFunctionForModel(PD::ModelId::ConstantLambda) &&
      lambdaResult.status.ok() &&
      Relative(*lambdaResult.lambdaParallelM, 4.5e9, 1.0e-15) &&
      Relative(*lambdaResult.kappaParallelM2PerS, speed * 4.5e9 / 3.0,
               1.0e-14) &&
      lambdaResult.provenance.evaluatedModelId == "constant_lambda" &&
      lambdaResult.provenance.configurationFingerprint ==
          installed.configurationFingerprint &&
      viaEvaluateActive.status.ok() &&
      *viaEvaluateActive.kappaParallelM2PerS ==
          *lambdaResult.kappaParallelM2PerS;

  // Switching models changes the pointer.  Expected value, as in
  // PD01-PARSER-BRIDGE: 9e9*(8e8/2e8)^0.5*(2e11/1e11)^1*0.25 = 9e9 m.
  const PD::Status powerStatus = PD::ConfigureActiveModelFromSection(
      Section({"model = power_law_lambda", "lambda0_m = 9e9",
               "independent_variable = rigidity", "rigidity0_V = 2e8",
               "independent_exponent = 0.5", "radius0_m = 1e11",
               "radial_exponent = 1.0", "use_region_factor = true"}, 1),
      nullptr);
  PD::LocalState powerLocal;
  powerLocal.positionM = {{2.0e11, 0.0, 0.0}};
  powerLocal.regionFactor = 0.25;
  const PD::ParallelResult powerResult = PD::ActiveParallelDiffusion(
      ParticleAtRigidity(8.0e8, ProtonMassKg, ElementaryChargeC), powerLocal);
  const bool switchPass = powerStatus.ok() &&
      PD::ActiveParallelDiffusion ==
          PD::BoundFunctionForModel(PD::ModelId::PowerLawLambda) &&
      PD::ActiveParallelDiffusion != lambdaPointer &&
      powerResult.status.ok() &&
      Relative(*powerResult.lambdaParallelM, 9.0e9, 2.0e-14);

  // A rejected section must leave pointer, parameters, and output unchanged.
  const PD::ActiveModelEvaluator beforeFailure = PD::ActiveParallelDiffusion;
  PD::ParsedSection unchanged = installed;
  const PD::Status rejected = PD::ConfigureActiveModelFromSection(
      Section({"model = constant_kappa", "kappa_parallel_m2_per_s = -1"}, 1),
      &unchanged);
  const PD::ParallelResult afterFailure = PD::ActiveParallelDiffusion(
      ParticleAtRigidity(8.0e8, ProtonMassKg, ElementaryChargeC), powerLocal);
  const bool transactionPass = !rejected.ok() &&
      PD::ActiveParallelDiffusion == beforeFailure &&
      unchanged.configurationFingerprint == installed.configurationFingerprint &&
      afterFailure.status.ok() &&
      afterFailure.provenance.evaluatedModelId == "power_law_lambda" &&
      Relative(*afterFailure.lambdaParallelM, 9.0e9, 2.0e-14);
  const bool dispatchPass = lambdaPass && switchPass && transactionPass;
  Add(tests, "PD13-BOUND-DISPATCH", "registry", dispatchPass,
      "section-selected model-specific pointer, model switch, failed update",
      dispatchPass ? "pointer tracks the accepted section; rejected section preserves it"
                   : std::string("bound dispatch failed (lambda=") +
                         (lambdaPass ? "ok" : "bad") + ", switch=" +
                         (switchPass ? "ok" : "bad") + ", transaction=" +
                         (transactionPass ? "ok" : "bad") + ")");
  return registryPass && dispatchPass;
}

bool TestHostInputAvailability(std::vector<TestRecord>* tests) {
  // The gate inspects only the selected model's declared needs; configurations
  // are constructed directly because the gate does not re-validate values.
  const PD::HostInputAvailability none;
  PD::HostInputAvailability all;
  all.nucleonCount = true;
  all.turbulenceDecomposition = true;
  all.timeFactor = true;
  all.regionFactor = true;
  all.radialFactor = true;
  all.effectiveFieldMagnitude = true;

  bool pass = true;
  std::string failed;
  const auto expect = [&](const char* name, const PD::ModelConfiguration& c,
                          bool acceptedWithoutInputs) {
    const PD::Status withNone = PD::CheckHostInputAvailability(c, none);
    const PD::Status withAll = PD::CheckHostInputAvailability(c, all);
    const bool ok = withAll.ok() &&
        (acceptedWithoutInputs
             ? withNone.ok()
             : withNone.code == PD::StatusCode::InvalidConfiguration);
    if (!ok && failed.empty()) failed = name;
    pass = pass && ok;
  };

  PD::ModelConfiguration c;
  c.model = PD::ModelId::ConstantLambda;
  expect("constant_lambda", c, true);
  c.model = PD::ModelId::ConstantKappa;
  expect("constant_kappa", c, true);
  c.model = PD::ModelId::PrescribedLambdaMuShape;
  expect("prescribed_lambda_mu_shape", c, true);
  for (const PD::ModelId turbulent :
       {PD::ModelId::QltSlabSpectrum, PD::ModelId::QltSlabInertial,
        PD::ModelId::BroadenedSlab, PD::ModelId::NlpaGivenPerp,
        PD::ModelId::NlgcE, PD::ModelId::NlgceN, PD::ModelId::NlgceF2014,
        PD::ModelId::TurbulenceAdapter, PD::ModelId::WaveSpectrumAdapter}) {
    c.model = turbulent;
    expect(PD::ModelName(turbulent), c, false);
  }
  c = PD::ModelConfiguration();
  c.model = PD::ModelId::PowerLawLambda;
  expect("power_law_lambda rigidity", c, true);
  c.powerLawLambda.independentVariable = PD::IndependentVariable::EnergyPerNucleon;
  expect("power_law_lambda energy_per_nucleon", c, false);
  c.powerLawLambda.independentVariable = PD::IndependentVariable::Rigidity;
  c.powerLawLambda.useTimeFactor = true;
  expect("power_law_lambda time_factor", c, false);
  c.powerLawLambda.useTimeFactor = false;
  c.powerLawLambda.useRegionFactor = true;
  expect("power_law_lambda region_factor", c, false);
  c = PD::ModelConfiguration();
  c.model = PD::ModelId::BrokenRigidityKappa;
  expect("broken_rigidity_kappa", c, true);
  c.brokenRigidityKappa.useRadialFactor = true;
  expect("broken_rigidity_kappa radial_factor", c, false);
  c.brokenRigidityKappa.useRadialFactor = false;
  c.brokenRigidityKappa.useRegionFactor = true;
  expect("broken_rigidity_kappa region_factor", c, false);
  c = PD::ModelConfiguration();
  c.model = PD::ModelId::Bohm;
  expect("bohm mean_field", c, true);
  c.bohm.fieldDefinition = PD::BohmFieldDefinition::EffectiveField;
  expect("bohm effective_field", c, false);
  c = PD::ModelConfiguration();
  c.model = PD::ModelId::TabulatedParallel;
  c.table.axes = {PD::TableAxis::Rigidity, PD::TableAxis::HeliocentricRadius};
  expect("tabulated_parallel rigidity/radius", c, true);
  c.table.axes = {PD::TableAxis::EnergyPerNucleon};
  expect("tabulated_parallel energy_per_nucleon", c, false);

  Add(tests, "PD13-HOST-GATE", "registry", pass,
      "host-declared optional-input availability",
      pass ? "models needing undeclared inputs are rejected; others accepted"
           : std::string("host gate misclassified ") + failed);
  return pass;
}

bool TestKStarConverter(std::vector<TestRecord>* tests) {
  // Verify the convention through K_star*beta=kappa_reference rather than by
  // copying the converter's division.  The particle supplies beta because it
  // is not species independent at fixed rigidity.
  const PD::ParticleState particle = TenMeVProton();
  PD::ParticleKinematics kinematics;
  PD::ComputeParticleKinematics(particle, &kinematics);
  double kStar = 0.0;
  const double referenceKappa = 8.0e17;
  const PD::Status status =
      PD::KStarFromReferenceKappa(referenceKappa, particle, &kStar);
  const bool pass = status.ok() &&
      Relative(kStar * kinematics.beta, referenceKappa, 2.0e-15);
  Add(tests, "PD02-KSTAR-CONVERTER", "broken_rigidity_kappa", pass,
      "Section 6.1 explicit actual-kappa to K_star conversion",
      pass ? "K_star=reference_kappa/beta"
           : "K_star convention conversion failed");
  return pass;
}

}  // namespace

int main(int argc, char** argv) {
  std::string jsonPath;
  if (argc == 3 && std::string(argv[1]) == "--json") jsonPath = argv[2];
  else if (argc != 1) {
    std::cerr << "usage: test_parallel_diffusion [--json FILE]\n";
    return 2;
  }

  std::vector<TestRecord> tests;
  bool passed = true;
  // Must precede every test that installs an active model.
  passed = TestUnconfiguredBoundDispatch(&tests) && passed;
  passed = TestSectionParser(&tests) && passed;
  passed = TestKinematics(&tests) && passed;
  passed = TestConstants(&tests) && passed;
  passed = TestPowerLaw(&tests) && passed;
  passed = TestBrokenRigidity(&tests) && passed;
  passed = TestBohm(&tests) && passed;
  passed = TestParserAndDispatch(&tests) && passed;
  passed = TestBoundDispatch(&tests) && passed;
  passed = TestHostInputAvailability(&tests) && passed;
  passed = TestKStarConverter(&tests) && passed;

  for (const TestRecord& test : tests) {
    std::cout << (test.passed ? "PASS " : "FAIL ") << test.id << " ["
              << test.model << "]: " << test.detail << '\n';
  }
  if (!jsonPath.empty()) {
    std::ofstream output(jsonPath);
    if (!output) {
      std::cerr << "cannot write JSON report: " << jsonPath << '\n';
      return 2;
    }
    output << "{\n  \"schema\": \"parallel-diffusion-test-report-v1\",\n"
           << "  \"passed\": " << (passed ? "true" : "false") << ",\n"
           << "  \"tests\": [\n";
    for (std::size_t i = 0; i < tests.size(); ++i) {
      const TestRecord& test = tests[i];
      output << "    {\"id\": \"" << JsonEscape(test.id)
             << "\", \"model\": \"" << JsonEscape(test.model)
             << "\", \"status\": \"" << (test.passed ? "PASS" : "FAIL")
             << "\", \"fixture\": \"" << JsonEscape(test.fixture)
             << "\", \"detail\": \"" << JsonEscape(test.detail) << "\"}"
             << (i + 1 == tests.size() ? "\n" : ",\n");
    }
    output << "  ]\n}\n";
  }
  std::cout << (passed ? "parallel_diffusion: PASS" : "parallel_diffusion: FAIL")
            << " (" << tests.size() << " checks)\n";
  return passed ? 0 : 1;
}
