#include "../ForwardModel.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <limits>
#include <string>
#include <vector>

namespace FM = Enceladus::ForwardModel;

namespace {

int Fail(const std::string& reason) {
  std::cerr << "FAIL: " << reason << '\n';
  return EXIT_FAILURE;
}

void Value(const std::string& key, double value) {
  std::cout << key << '=' << std::scientific << std::setprecision(17)
            << value << '\n';
}

FM::Configuration TestConfiguration(const std::string& sourceFile) {
  FM::Configuration c;
  c.schemaVersion = 2;
  c.sourceDistributionFile = sourceFile;
  c.simulatedSpecies = {"H2O", "CO2"};
  c.injectedSpecies = {"H2O", "CO2"};
  c.chemistrySpecies = {"H2O"};
  c.selectedReactions = {"burger2007.photo.h2o_to_oh_h"};
  c.totalSourceRatePerS = {{"H2O", 1.0e28}, {"CO2", 2.0e26}};
  c.characteristicSpeedMPerS = {{"H2O", 2000.0}, {"CO2", 1000.0}};
  c.courantNumber = {{"H2O", 0.5}, {"CO2", 0.5}};
  c.primaryMacroparticlesPerStep = {{"H2O", 100}, {"CO2", 50}};
  c.sampleLengthSteps = 200;
  c.totalIterations = 800;
  c.randomSeed = 19;
  c.populationControl = FM::PopulationControl::Off;
  c.solarActivity = FM::SolarActivity::Quiet;
  c.plasmaBackgroundMode = FM::PlasmaBackgroundMode::CacheOnly;
  c.plasmaCacheDirectory = ".";
  return c;
}

int SourceCase(const std::string& sourceFile) {
  FM::SourceDistribution distribution;
  FM::Status status = FM::ParseSourceDistribution(sourceFile, &distribution);
  if (!status.success) return Fail(status.message);
  if (distribution.sources.size() != 15) return Fail("source count is not 15");
  if (distribution.stripeSurfaceMapping !=
          FM::StripeSurfaceMapping::LocalVerticalToReferenceSphere ||
      distribution.stripeReferenceRadiusM != 252100.0)
    return Fail("stripe surface projection receipt is not explicit");
  const FM::Source* baghdad = nullptr;
  const FM::Source* alexandria = nullptr;
  for (const FM::Source& source : distribution.sources) {
    if (source.id == "I_baghdad") baghdad = &source;
    if (source.id == "alexandria") alexandria = &source;
  }
  if (baghdad == nullptr || alexandria == nullptr)
    return Fail("canonical sources were not parsed");
  const FM::Vec3 normal = FM::PointJetNormal(*baghdad);
  const FM::Vec3 axis = FM::PointJetAxis(*baghdad);
  const double normalNorm = std::sqrt(normal.x*normal.x + normal.y*normal.y +
                                      normal.z*normal.z);
  const double axisNorm = std::sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
  std::vector<double> stripeWeights;
  status = FM::StripeNodeWeights(*alexandria, &stripeWeights);
  if (!status.success) return Fail(status.message);
  double stripeSum = 0.0;
  for (double weight : stripeWeights) stripeSum += weight;
  double stripeMinimumRadius = std::numeric_limits<double>::infinity();
  double stripeMaximumRadius = 0.0;
  for (const FM::Source& source : distribution.sources)
    if (source.geometry == FM::SourceGeometry::PolylineStripe)
      for (const FM::Vec3& vertex : source.verticesM) {
        const double radius = std::sqrt(vertex.x*vertex.x + vertex.y*vertex.y +
                                        vertex.z*vertex.z);
        stripeMinimumRadius = std::min(stripeMinimumRadius, radius);
        stripeMaximumRadius = std::max(stripeMaximumRadius, radius);
      }
  Value("source_count", static_cast<double>(distribution.sources.size()));
  Value("stripe_weight_sum", stripeSum);
  Value("point_normal_norm", normalNorm);
  Value("point_axis_norm", axisNorm);
  Value("stripe_minimum_radius_m", stripeMinimumRadius);
  Value("stripe_maximum_radius_m", stripeMaximumRadius);
  Value("stripe_reference_radius_m", distribution.stripeReferenceRadiusM);
  return EXIT_SUCCESS;
}

int CatalogCase(const std::string& file, std::size_t expected,
                bool requireSubnormal = false) {
  FM::SourceDistribution distribution;
  const auto status=FM::ParseSourceDistribution(file,&distribution);
  if (!status.success) return Fail(status.message);
  if (distribution.sources.size()!=expected) return Fail("catalog source count differs");
  // This decimal is supplied independently by the Python wire fixture.
  // Require its actual binary64 value, rather than merely accepting a file
  // after silently rounding the tiny strength down to zero.
  if (requireSubnormal &&
      distribution.sources.front().strength.at("H2O") != 3.7637234148938442e-317)
    return Fail("subnormal source strength was not preserved");
  std::map<std::string,double> totals;
  for (const auto& source:distribution.sources) {
    if (source.geometry!=FM::SourceGeometry::PointJet) return Fail("catalog contains a stripe");
    for (const auto& item:source.strength) totals[item.first]+=item.second;
  }
  for (const auto& item:totals) {
    if (!(item.second>0) || !std::isfinite(item.second)) return Fail("invalid catalog total");
    if (distribution.mode==FM::SourceMode::Shape && std::abs(item.second-1)>1e-12)
      return Fail("catalog shape is not normalized");
  }
  Value("catalog_source_count",static_cast<double>(expected));
  return EXIT_SUCCESS;
}

int ChemistryCase() {
  const FM::Reaction* photo =
      FM::FindReaction("burger2007.photo.h2o_to_oh_h");
  const FM::Reaction* cx =
      FM::FindReaction("burger2007.cx.hplus_h2o_to_h_h2op");
  if (photo == nullptr || cx == nullptr) return Fail("reaction lookup failed");
  const FM::Reaction* resonant =
      FM::FindReaction("burger2007.cx.h2op_h2o_resonant");
  if (resonant == nullptr || resonant->products.size() != 2 ||
      resonant->products[0].species != "H2O" ||
      resonant->products[1].species != "H2O_PLUS")
    return Fail("resonant Burger channel must retain H2O and produce H2O_PLUS");
  double quiet = 0.0, active = 0.0, rate = 0.0;
  FM::Status status = FM::EvaluateReactionRate(
      *photo, FM::SolarActivity::Quiet, 0, 0, 0, &quiet);
  if (!status.success) return Fail(status.message);
  status = FM::EvaluateReactionRate(
      *photo, FM::SolarActivity::Active, 0, 0, 0, &active);
  if (!status.success) return Fail(status.message);
  status = FM::EvaluateReactionRate(
      *cx, FM::SolarActivity::Quiet, 7.0e7, 0.15, 2.0e4, &rate);
  if (!status.success) return Fail(status.message);
  double rejected = 0.0;
  status = FM::EvaluateReactionRate(
      *cx, FM::SolarActivity::Quiet, 7.0e7, 0.15, 2.7e4, &rejected);
  if (status.success) return Fail("Burger speed extrapolation was not rejected");
  double zeroRate = -1.0;
  status = FM::EvaluateReactionRate(
      *cx, FM::SolarActivity::Quiet, 0.0, 0.15, 0.0, &zeroRate);
  if (!status.success || zeroRate != 0.0)
    return Fail("zero reactant density did not produce an exact zero rate");
  const double alpha = (20000.0 - 15000.0) / 11000.0;
  const double sigma = cx->crossSectionAt15KmS_M2 + alpha *
      (cx->crossSectionAt26KmS_M2 - cx->crossSectionAt15KmS_M2);
  Value("photo_h2o_oh_h_quiet_s-1", quiet);
  Value("photo_h2o_oh_h_active_s-1", active);
  Value("cx_hplus_h2o_sigma_20kms_m2", sigma);
  Value("cx_hplus_h2o_rate_20kms_s-1", rate);
  Value("cx_zero_reactant_rate_s-1", zeroRate);
  return EXIT_SUCCESS;
}

long double DirectFluxIntegral(double cosine, double speedRatio,
                               double upper) {
  // Independent composite Simpson quadrature of the defining integral. It
  // intentionally does not reuse ForwardModel's error-function recurrence.
  constexpr int intervals = 200000;
  const long double h = static_cast<long double>(upper)/intervals;
  auto integrand = [=](long double t) {
    return t*t*t*std::exp(-t*t + 2.0L*speedRatio*cosine*t -
                          speedRatio*speedRatio);
  };
  long double sum = integrand(0.0L)+integrand(upper);
  for (int i = 1; i < intervals; ++i)
    sum += (i % 2 == 0 ? 2.0L : 4.0L)*integrand(i*h);
  return sum*h/3.0L;
}

int SourceLawCase() {
  // The first pair is close to the most adverse current hybrid-8 CO2 jet
  // direction (s*c=-1.68), and therefore exercises the cancellation-sensitive
  // side of the qualified recurrence.
  const double cosine = -0.35;
  const double speedRatio = 4.8;
  double moment = 0.0;
  FM::Status status = FM::Tenishev2010FluxMoment3(
      cosine, speedRatio, &moment);
  if (!status.success) return Fail(status.message);
  const long double direct = DirectFluxIntegral(cosine, speedRatio, 16.0);
  if (std::abs(moment-static_cast<double>(direct))/static_cast<double>(direct) >
      2.0e-10)
    return Fail("Tenishev K3 recurrence disagrees with direct quadrature");

  constexpr double probability = 0.73;
  double quantile = 0.0;
  status = FM::Tenishev2010SpeedQuantile(
      cosine, speedRatio, probability, &quantile);
  if (!status.success) return Fail(status.message);
  const long double partial = DirectFluxIntegral(cosine, speedRatio, quantile);
  if (std::abs(static_cast<double>(partial/direct)-probability) > 2.0e-9)
    return Fail("Tenishev inverse CDF disagrees with direct quadrature");

  // Tilted Porco sources reach s*c<-2. Keep the same independent defining
  // integral oracle; extending support must reproduce the physical tail,
  // rather than weakening the old guard or clipping small source fluxes.
  for (const auto& pair:std::vector<std::pair<double,double>>{
      {-0.5,5.0},{-0.95,6.0},{-1.0,8.0}}) {
    double tail=0,tailQuantile=0;
    status=FM::Tenishev2010FluxMoment3(pair.first,pair.second,&tail);
    if (!status.success) return Fail(status.message);
    const long double independent=DirectFluxIntegral(pair.first,pair.second,16);
    if (std::abs(tail/static_cast<double>(independent)-1)>2.0e-9)
      return Fail("negative-tail K3 differs from independent quadrature");
    for (double probability:{0.0001,0.1,0.73,0.9999}) {
      status=FM::Tenishev2010SpeedQuantile(pair.first,pair.second,probability,&tailQuantile);
      if (!status.success) return Fail(status.message);
      const long double cdf=DirectFluxIntegral(pair.first,pair.second,tailQuantile)/independent;
      if (std::abs(static_cast<double>(cdf)-probability)>2.0e-9)
        return Fail("negative-tail inverse CDF differs from independent quadrature");
    }
  }
  Value("tenishev_k3_c-0.35_s4.8", moment);
  Value("tenishev_q73_c-0.35_s4.8", quantile);
  return EXIT_SUCCESS;
}

int NumericsCase(const std::string& sourceFile) {
  FM::SourceDistribution distribution;
  FM::Status status = FM::ParseSourceDistribution(sourceFile, &distribution);
  if (!status.success) return Fail(status.message);
  const FM::Configuration configuration = TestConfiguration(sourceFile);
  FM::NumericalPlan plan;
  status = FM::BuildNumericalPlan(configuration, distribution, 100.0, &plan);
  if (!status.success) return Fail(status.message);
  if (plan.limitingSpecies != "H2O") return Fail("wrong Courant-limiting species");
  for (const std::string& species : configuration.injectedSpecies)
    if (plan.expectedMacroparticlesPerStep.at(species) !=
        static_cast<double>(configuration.primaryMacroparticlesPerStep.at(species)))
      return Fail("primary macro count is not the requested exact integer");
  Value("minimum_cell_size_m", 100.0);
  Value("global_time_step_s", plan.globalTimeStepS);
  Value("primary_weight_H2O", plan.statisticalWeight.at("H2O"));
  Value("primary_weight_CO2", plan.statisticalWeight.at("CO2"));
  if (plan.statisticalWeight.count("OH") != 0 ||
      plan.statisticalWeight.count("H") != 0)
    return Fail("loss-only photolysis planned daughter weights");
  Value("planned_daughter_species", 0.0);
  std::cout << "limiting_species=" << plan.limitingSpecies << '\n';

  FM::SourceDistribution absolute = distribution;
  absolute.mode = FM::SourceMode::Absolute;
  for (FM::Source& source : absolute.sources) {
    source.strength["H2O"] *= 1.0e28;
    source.strength["CO2"] *= 2.0e26;
  }
  FM::Configuration absoluteConfiguration = configuration;
  absoluteConfiguration.totalSourceRatePerS.clear();
  FM::NumericalPlan absolutePlan;
  status = FM::BuildNumericalPlan(
      absoluteConfiguration, absolute, 100.0, &absolutePlan);
  if (!status.success) return Fail(status.message);
  Value("absolute_rate_H2O_s-1",
        absolutePlan.physicalSourceRatePerS.at("H2O"));
  Value("absolute_rate_CO2_s-1",
        absolutePlan.physicalSourceRatePerS.at("CO2"));
  Value("absolute_primary_weight_H2O",
        absolutePlan.statisticalWeight.at("H2O"));
  return EXIT_SUCCESS;
}

int InputCase(const std::string& inputFile, const std::string& sourceFile) {
  FM::Configuration parsed;
  FM::Status status = FM::ParseConfiguration(inputFile, &parsed);
  if (!status.success) return Fail(status.message);
  FM::SourceDistribution distribution;
  status = FM::ParseSourceDistribution(sourceFile, &distribution);
  if (!status.success) return Fail(status.message);
  status = FM::ValidateConfiguration(parsed, distribution);
  if (!status.success) return Fail(status.message);
  if (parsed.sourceDistributionFile != sourceFile)
    return Fail("relative source path was not resolved against runtime input");
  if (parsed.flybyDensityOutput) {
    const std::filesystem::path expectedTrajectory =
        std::filesystem::path(inputFile).parent_path()/"trajectory.csv";
    if (parsed.flybyId != "E3" ||
        parsed.flybyTrajectoryFile !=
            std::filesystem::absolute(expectedTrajectory).lexically_normal().string() ||
        !parsed.flybyMeshRefinement ||
        parsed.flybyRefinementTargetCellSizeM != 20000.0 ||
        parsed.flybyRefinementRadiusM != 150000.0 ||
        parsed.flybyRefinementTimeMinS != -120.0 ||
        parsed.flybyRefinementTimeMaxS != 100.0 || parsed.volumeFieldOutput)
      return Fail("flyby output controls or relative trajectory path were not parsed");
    Value("parsed_flyby_refinement_target_cell_size_m",
          parsed.flybyRefinementTargetCellSizeM);
    Value("parsed_flyby_refinement_radius_m",
          parsed.flybyRefinementRadiusM);
  } else if (!parsed.flybyId.empty() || !parsed.flybyTrajectoryFile.empty()) {
    return Fail("disabled flyby output retained dependent controls");
  }
  if (!parsed.observations.empty()) {
    const std::filesystem::path expectedContext =
        std::filesystem::path(inputFile).parent_path()/"context.csv";
    if (parsed.observations.size() != 3 ||
        parsed.observations[0].id != "inms_test" ||
        parsed.observations[0].instrument != FM::ObservationInstrument::Inms ||
        parsed.observations[0].comparisonSpecies != "H2O" ||
        parsed.observations[1].id != "inms_context_test" ||
        parsed.observations[1].contextObservable !=
            "INMS_OSNB_H2O_counts_per_integration" ||
        parsed.observations[1].dataFile !=
            std::filesystem::absolute(expectedContext).lexically_normal().string() ||
        parsed.observations[2].id != "uvis_test" ||
        parsed.observations[2].instrument != FM::ObservationInstrument::Uvis ||
        parsed.observations[2].uvisRayGeometry !=
            FM::UvisRayGeometry::FiniteObserverToSun ||
        parsed.observations[2].publishedH2OColumnM2 != 1.4e20)
      return Fail("repeatable INMS/UVIS observation blocks were not parsed");
    Value("parsed_observation_requests",
          static_cast<double>(parsed.observations.size()));
  }
  Value("parsed_sample_length_steps", parsed.sampleLengthSteps);
  Value("parsed_total_iterations", parsed.totalIterations);
  if (parsed.populationControl == FM::PopulationControl::VelocityShift) {
    Value("parsed_split_position_cell_diagonal_fraction",
          parsed.splitPositionShiftCellDiagonalFraction);
    Value("parsed_split_velocity_thermal_speed_fraction",
          parsed.splitVelocityShiftThermalSpeedFraction);
    Value("parsed_nonuniform_split_position_shift",
          parsed.nonuniformSplitPositionShift ? 1.0 : 0.0);
  }
  return EXIT_SUCCESS;
}

int ValidationCase(const std::string& sourceFile) {
  FM::SourceDistribution distribution;
  FM::Status status = FM::ParseSourceDistribution(sourceFile, &distribution);
  if (!status.success) return Fail(status.message);

  FM::Configuration configuration = TestConfiguration(sourceFile);
  configuration.simulatedSpecies.push_back("H3O_PLUS");
  configuration.characteristicSpeedMPerS["H3O_PLUS"] = 1000.0;
  configuration.courantNumber["H3O_PLUS"] = 0.5;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("H3O_PLUS") == std::string::npos)
    return Fail("H3O_PLUS was not explicitly rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.selectedReactions = {"burger2007.photo.co2_invented"};
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("unknown") == std::string::npos)
    return Fail("unknown CO2 chemistry was not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.totalSourceRatePerS.erase("CO2");
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("shape mode") == std::string::npos)
    return Fail("shape mode without a CO2 total was not rejected");

  configuration = TestConfiguration(sourceFile);
  distribution.mode = FM::SourceMode::Absolute;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("forbids") == std::string::npos)
    return Fail("absolute mode with a competing total was not rejected");

  distribution.mode = FM::SourceMode::Shape;
  configuration = TestConfiguration(sourceFile);
  configuration.injectedSpecies = {"H2O"};
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("not injected") == std::string::npos)
    return Fail("inactive-species run controls were not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.plasmaBackgroundMode = FM::PlasmaBackgroundMode::CacheOrIces;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success ||
      status.message.find("ices_installation_directory") == std::string::npos)
    return Fail("cache_or_ices without its explicit ICES authority was not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.plasmaBackgroundMode = FM::PlasmaBackgroundMode::None;
  configuration.plasmaCacheDirectory.clear();
  status = FM::ValidateConfiguration(configuration, distribution);
  if (!status.success)
    return Fail("plasma_background_mode=none with no plasma settings was rejected");

  configuration.plasmaCacheDirectory = ".";
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success ||
      status.message.find("none forbids") == std::string::npos)
    return Fail("plasma_background_mode=none accepted a cache directory");

  configuration = TestConfiguration(sourceFile);
  configuration.selectedReactions.clear();
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("chemistry_species must be none") ==
                            std::string::npos)
    return Fail("chemistry species without reactions were not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.injectedSpecies = {"CO2"};
  configuration.primaryMacroparticlesPerStep.erase("H2O");
  configuration.totalSourceRatePerS.erase("H2O");
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("not injected") == std::string::npos)
    return Fail("a non-primary simulated reaction parent was not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.plasmaIonFraction["H_PLUS"] = 0.15;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success ||
      status.message.find("absent from selected_reactions") == std::string::npos)
    return Fail("an unused plasma ion fraction was not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.chemistryWeightReferenceNumberDensityM3 = 7.0e7;
  configuration.chemistryWeightReferenceRelativeSpeedMPerS = 2.0e4;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("forbidden") == std::string::npos)
    return Fail("charge-exchange references without charge exchange were not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.splitPositionShiftCellDiagonalFraction = 0.2;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success ||
      status.message.find("splitting_merging=off") == std::string::npos)
    return Fail("split controls were accepted while population control was off");

  configuration = TestConfiguration(sourceFile);
  configuration.populationControl = FM::PopulationControl::Scatter;
  configuration.minParticlesPerCell = 10;
  configuration.maxParticlesPerCell = 20;
  configuration.splitVelocityShiftThermalSpeedFraction = 0.01;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("scatter forbids") == std::string::npos)
    return Fail("velocity-shift controls were accepted in scatter mode");

  configuration = TestConfiguration(sourceFile);
  configuration.plasmaBackgroundMode = FM::PlasmaBackgroundMode::CacheOrIces;
  configuration.icesInstallationDirectory = "/safe/ICES/Models";
  configuration.icesCase = "case;touch_bad";
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success || status.message.find("ices_case") == std::string::npos)
    return Fail("a shell-unsafe ICES case was not rejected");

  configuration = TestConfiguration(sourceFile);
  configuration.selectedReactions = {
      "burger2007.cx.hplus_h2o_to_h_h2op"};
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success ||
      status.message.find("only photolytic reactions") == std::string::npos)
    return Fail("charge exchange was not rejected by loss-only schema 2");

  configuration = TestConfiguration(sourceFile);
  configuration.flybyMeshRefinement = true;
  configuration.flybyRefinementTargetCellSizeM = 20000.0;
  configuration.flybyRefinementRadiusM = 150000.0;
  configuration.flybyRefinementTimeMinS = -120.0;
  configuration.flybyRefinementTimeMaxS = 100.0;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success ||
      status.message.find("requires flyby_density_output") == std::string::npos)
    return Fail("flyby mesh refinement was accepted without flyby output");

  configuration.flybyDensityOutput = true;
  configuration.flybyId = "E3";
  configuration.flybyTrajectoryFile = "trajectory.csv";
  configuration.flybyRefinementTimeMaxS = -120.0;
  status = FM::ValidateConfiguration(configuration, distribution);
  if (status.success ||
      status.message.find("time_max_s greater") == std::string::npos)
    return Fail("an empty flyby refinement interval was accepted");

  Value("negative_validation_checks", 17.0);
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) return Fail("usage: forward_model_tests CASE [arguments]");
  const std::string test = argv[1];
  if (test == "source" && argc == 3) return SourceCase(argv[2]);
  if (test == "catalog" && argc == 4) return CatalogCase(argv[2],std::stoul(argv[3]));
  if (test == "catalog_subnormal" && argc == 4)
    return CatalogCase(argv[2],std::stoul(argv[3]),true);
  if (test == "chemistry" && argc == 2) return ChemistryCase();
  if (test == "source_law" && argc == 2) return SourceLawCase();
  if (test == "numerics" && argc == 3) return NumericsCase(argv[2]);
  if (test == "input" && argc == 4) return InputCase(argv[2], argv[3]);
  if (test == "validation" && argc == 3) return ValidationCase(argv[2]);
  return Fail("invalid case or argument count");
}
