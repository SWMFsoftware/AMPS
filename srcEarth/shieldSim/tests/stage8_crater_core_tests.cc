/* Deterministic scientific-boundary tests for the axial CRaTER benchmark. */

#include "CraterBenchmark.hh"
#include "StrictJson.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

namespace {
int failures = 0;
void Check(bool condition, const std::string& description) {
  if (!condition) { ++failures; std::cerr << "FAIL: " << description << '\n'; }
}
bool Near(double actual, double expected, double relative = 1.0e-11) {
  return std::abs(actual-expected) <= relative *
      std::max({1.0, std::abs(actual), std::abs(expected)});
}
}

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: stage8_crater_core_tests SOURCE_ROOT EVIDENCE_DIR "
                 "CRATER_REFERENCE_DIR\n";
    return 2;
  }
  namespace fs = std::filesystem;
  const fs::path root = fs::absolute(argv[1]);
  const fs::path evidence = fs::absolute(argv[2]);
  const fs::path reference = fs::absolute(argv[3]);
  fs::create_directories(evidence);
  shield_sim::Json::Array comparisons;

  const auto& stack = shield_sim::PublishedCraterAxialStack();
  const double total = std::accumulate(stack.begin(), stack.end(), 0.0,
      [](double sum, const auto& element) { return sum + element.thicknessMm; });
  Check(stack.size() == 19, "all 19 published solid/gap axial entries retained");
  Check(Near(total, 132.02606), "published axial thickness sum");
  Check(std::count_if(stack.begin(), stack.end(),
      [](const auto& element){ return element.sensitive; }) == 6,
      "exactly six sensitive silicon detectors");
  comparisons.emplace_back(shield_sim::Json::Object{
      {"id", "published_axial_stack_total"}, {"unit", "mm"},
      {"expected", 132.02606}, {"actual", total}});

  const auto& ions = shield_sim::PublishedCraterPrimaryIons();
  Check(ions.size() == 28, "BON2020 maps exactly one primary isotope per Z");
  for (std::size_t i = 0; i < ions.size(); ++i)
    Check(ions[i].z == static_cast<int>(i+1), "ion mapping is ordered Z=1..28");
  Check(ions.front().a == 1 && ions[1].a == 4 && ions[25].a == 56 &&
        ions.back().a == 58, "archive-derived H, He, Fe, and Ni masses");

  shield_sim::CraterBonSpectrum analytic;
  analytic.energyMeVPerNucleon = {10.0, 100.0};
  analytic.fluxByEnergy = {std::vector<double>(28, 10.0),
                           std::vector<double>(28, 1.0)};
  const double interpolated = shield_sim::InterpolateCraterFlux(
      analytic, 1, std::sqrt(1000.0));
  Check(Near(interpolated, std::sqrt(10.0)),
        "log-log source interpolation analytic power law");
  comparisons.emplace_back(shield_sim::Json::Object{
      {"id", "log_log_flux_interpolation"}, {"unit", "arbitrary"},
      {"expected", std::sqrt(10.0)}, {"actual", interpolated}});

  const double ratio = std::pow(10.0, 1.0/80.0);
  const auto edges = shield_sim::CraterObservationBinEdges({1.0, ratio, ratio*ratio});
  Check(edges.size() == 4 && Near(edges.front(), 1.0/std::sqrt(ratio)) &&
        Near(edges.back(), ratio*ratio*std::sqrt(ratio)),
        "published 80-bin/decade centers reconstruct exact edges");

  shield_sim::CraterDeposits low;
  low.kevPerUm = {0.0, 0.30, 0.0, 0.30, 0.30, 0.30};
  Check(shield_sim::CraterCoincidence(low), "low-LET coincidence");
  Check(Near(shield_sim::CraterCombinedObservable(low), 0.30),
        "D6 supplies low normalized-deposit branch");
  Check(shield_sim::CraterCombinedPublishedCut(low),
        "published D6 cut accepts a diagonal MIP-like event");

  shield_sim::CraterDeposits high;
  high.kevPerUm = {0.0, 5.0, 0.0, 5.0, 50.0, 100.0};
  Check(Near(shield_sim::CraterCombinedObservable(high), 50.0),
        "D5 supplies branch when D6 exceeds 20 keV/um");
  Check(shield_sim::CraterCombinedPublishedCut(high),
        "published D5 logic accepts a foreground-like high-deposit event");
  shield_sim::CraterDeposits rejected = high;
  rejected.kevPerUm[1] = 0.5;
  Check(!shield_sim::CraterCombinedPublishedCut(rejected),
        "published D5 logic rejects low-D2 background-like event");

  const auto definition = shield_sim::LoadCraterBenchmarkCase(
      (root/"examples/stage8_crater_foreground_axial.json").string());
  Check(definition.targetRelativeStatisticalUncertainty == 0.10,
        "Stage-8 example explicitly records the approved 10% target");
  const auto source = shield_sim::LoadCraterBonSpectrum(
      (reference/"glace/solmin_gcr.json").string(), definition.sourceSha256);
  Check(source.energyMeVPerNucleon.size() == 100 &&
        source.fluxByEnergy.size() == 100 &&
        source.fluxByEnergy.front().size() == 28,
        "checksum-pinned BON2020 source dimensions");
  const auto observation = shield_sim::LoadCraterObservation(
      definition.observationPath, definition.observationSha256,
      definition.observationColumn);
  Check(observation.centersKevPerUm.size() == 291,
        "checksum-pinned CRaTER observation has 291 bins");

  shield_sim::Json document(shield_sim::Json::Object{
      {"schema_version", "shieldSim.validation_evidence.v1"},
      {"test_case", "stage8_crater_foreground_axial_core"},
      {"benchmark_scope", shield_sim::kCraterBenchmarkScope},
      {"comparisons", shield_sim::Json(std::move(comparisons))},
      {"scientific_status", "deterministic_boundary_pass_not_observational_validation"},
  });
  std::ofstream output(evidence/"stage8_crater_core.json");
  output << document.DumpCanonical() << '\n';
  if (failures) {
    std::cerr << failures << " CRaTER core test(s) failed\n";
    return 1;
  }
  std::cout << "All foreground-only axial CRaTER core tests passed\n";
  return 0;
}
