#include "PublishedSourceAdapters.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

int failures = 0;

void Check(bool condition, const std::string& label) {
  if (!condition) {
    std::cerr << "FAIL " << label << '\n';
    ++failures;
  } else {
    std::cout << "PASS " << label << '\n';
  }
}

bool Near(double actual, double expected, double relative = 1.0e-10) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  return std::abs(actual - expected) <= relative * scale;
}

void CheckThrows(const std::function<void()>& operation, const std::string& label) {
  try {
    operation();
    Check(false, label + " (did not throw)");
  } catch (const std::exception&) {
    Check(true, label);
  }
}

shield_sim::GoesSgpsRecord Record(std::int64_t start, double flux) {
  shield_sim::GoesSgpsRecord record;
  record.startUnixSeconds = start;
  record.midpointUnixSeconds = start + 150;
  record.endUnixSeconds = start + 300;
  record.sensorUnit = 0;
  record.spacecraftAxis = "-X";
  record.nominalLook = "westward";
  record.channel = "PTEST";
  record.energyLowMeV = 10.0;
  record.energyHighMeV = 20.0;
  record.l1bRecordsInAverage = 300;
  record.validOneSecondSamples = 300;
  record.hasFlux = true;
  record.fluxPerCm2SrSMeV = flux;
  record.hasUncertainty = true;
  record.uncertaintyPerCm2SrSMeV = 0.1;
  return record;
}

shield_sim::GoesEventWindow Window(std::int64_t start, std::int64_t end) {
  shield_sim::GoesEventWindow window;
  window.id = "analytic";
  window.role = "primary";
  window.startUnixSeconds = start;
  window.endUnixSeconds = end;
  window.startUtc = "analytic-start";
  window.endUtc = "analytic-end";
  window.definition = "analytic test";
  window.sourceUrl = "local analytic fixture";
  return window;
}

} // namespace

int main(int argc, char** argv) {
  if (argc > 2) {
    std::cerr << "usage: stage2_source_adapter_tests [validation-evidence-directory]\n";
    return 2;
  }
  const std::string root = SHIELDSIM_SOURCE_DIR;
  shield_sim::Json::Array evidenceComparisons;
  auto Evidence = [&](const std::string& id, double expected, double actual,
                      const std::string& unit, const std::string& basis) {
    evidenceComparisons.emplace_back(shield_sim::Json::Object{
        {"actual", actual}, {"basis", basis}, {"expected", expected},
        {"id", id}, {"unit", unit}});
  };

  const auto ams = shield_sim::AMS02HapiAdapter::Load(
      root + "/tests/reference_spectra/gcr/ams02_20131118_hapi_info.json",
      root + "/tests/reference_spectra/gcr/ams02_20131118_hapi_data.csv",
      root + "/tests/reference_spectra/gcr/ams02_20131118_interpretation.v1.json");
  Check(ams.components.size() == 2, "AMS adapter retains separate proton and helium grids");
  Check(ams.components[0].binLower.size() == 45 && ams.components[1].binLower.size() == 40,
        "AMS adapter retains published grid sizes");
  Check(ams.components[1].speciesId == "helium_mixed_isotopes" &&
        ams.components[1].massNumber == 0 && ams.components[1].geant4Particle.empty(),
        "AMS inclusive helium is not silently relabeled as alpha particles");
  Check(Near(ams.components[0].binLower.front(), 432.9881050056602, 1.0e-12) &&
        Near(ams.components[1].binLower.front(), 404.0969061796101, 1.0e-12),
        "AMS rigidity edges convert to kinetic energy per nucleon analytically");
  Evidence("ams_first_proton_lower_edge", 432.9881050056602,
           ams.components[0].binLower.front(), "MeV/n", "analytic");
  Evidence("ams_first_helium_lower_edge", 404.0969061796101,
           ams.components[1].binLower.front(), "MeV/n", "analytic");

  // Primary Table SM XXXV gives 395/(m2 sr s GV) in the first proton bin
  // [1.00,1.16] GV and 46.58 in the first helium bin [1.92,2.15] GV. A
  // bin-average Jacobian J_T = J_R*dR/dT reproduces the HAPI 524 and 127
  // values. This independent spot check is what resolves the HAPI coordinate
  // array's misleading generic JSON key "centers".
  const double pExpected = 395.0 * 0.16 /
      ((ams.components[0].binUpper[0] - ams.components[0].binLower[0]) / 1000.0);
  const double heExpected = 46.58 * 0.23 /
      ((ams.components[1].binUpper[0] - ams.components[1].binLower[0]) / 1000.0);
  Check(Near(ams.components[0].value[0] / 1.0e-7, pExpected, 1.0e-3),
        "AMS first proton value matches primary-table Jacobian");
  Check(Near(ams.components[1].value[0] / 1.0e-7, heExpected, 2.0e-3),
        "AMS first helium value matches primary-table Jacobian");
  Evidence("ams_first_proton_jacobian", pExpected,
           ams.components[0].value[0] / 1.0e-7,
           "1/(m2 sr s (GeV/n))", "analytic_primary_table");
  Evidence("ams_first_helium_jacobian", heExpected,
           ams.components[1].value[0] / 1.0e-7,
           "1/(m2 sr s (GeV/n))", "analytic_primary_table");
  Check(shield_sim::Json::Parse(ams.canonicalMetadata.DumpCanonical()).DumpCanonical() ==
        ams.canonicalMetadata.DumpCanonical(), "AMS canonical dataset round-trips exactly");

  const std::string amsCycleInterpretation =
      root + "/tests/reference_spectra/gcr/ams02_2025/interpretation.v1.json";
  const auto amsCycle =
      shield_sim::AMS02SolarCycleAdapter::LoadRotation(amsCycleInterpretation, 2581);
  Check(amsCycle.components.size() == 7,
        "AMS solar-cycle adapter retains all published elements He through O");
  Check(std::all_of(amsCycle.components.begin(), amsCycle.components.end(),
                    [](const auto& component) {
                      return component.binLower.size() == 40 &&
                             component.value.size() == 40 &&
                             component.statisticalUncertainty.size() == 40 &&
                             component.timeDependentSystematicUncertainty.size() == 40 &&
                             component.timeIndependentSystematicUncertainty.size() == 40 &&
                             component.massNumber == 0 && component.geant4Particle.empty();
                    }),
        "AMS solar-cycle grids and three uncertainty components remain lossless");
  Check(amsCycle.energyVariable == "rigidity" && amsCycle.energyUnit == "GV" &&
        amsCycle.components.front().binLower.front() == 1.92 &&
        amsCycle.components.front().binUpper.back() == 60.3,
        "AMS inclusive elements remain on their published rigidity grid");
  Check(Near(amsCycle.components.front().value.front(), 77.56, 1.0e-14) &&
        Near(amsCycle.components.back().value.front(), 1.631, 1.0e-14),
        "AMS held-out rotation matches independent helium and oxygen row values");
  Evidence("ams2025_br2581_helium_first_bin", 77.56,
           amsCycle.components.front().value.front(), "1/(m2 sr s GV)",
           "independent_observation_oracle");
  Evidence("ams2025_br2581_oxygen_first_bin", 1.631,
           amsCycle.components.back().value.front(), "1/(m2 sr s GV)",
           "independent_observation_oracle");
  CheckThrows([&] {
    shield_sim::AMS02SolarCycleAdapter::LoadRotation(amsCycleInterpretation, 2472);
  }, "AMS missing published rotation is rejected rather than interpolated");
  Check(shield_sim::Json::Parse(amsCycle.canonicalMetadata.DumpCanonical()).DumpCanonical() ==
        amsCycle.canonicalMetadata.DumpCanonical(),
        "AMS solar-cycle canonical dataset round-trips exactly");

  const auto cris = shield_sim::ACECRISHapiAdapter::LoadRecord(
      root + "/tests/reference_spectra/gcr/ace_cris/ac_h3_cris_hapi_info.json",
      root + "/tests/reference_spectra/gcr/ace_cris/ac_h3_cris_20191201_all.csv",
      root + "/tests/reference_spectra/gcr/ace_cris/interpretation.v1.json");
  Check(cris.components.size() == 24 && cris.components.front().atomicNumber == 5 &&
        cris.components.back().atomicNumber == 28,
        "ACE/CRIS adapter retains every published element B through Ni");
  Check(std::all_of(cris.components.begin(), cris.components.end(),
                    [](const auto& component) {
                      return component.binLower.size() == 7 && component.binUpper.size() == 7 &&
                             component.value.size() == 7 && component.counts.size() == 7 &&
                             component.statisticalUncertainty.size() == 7 &&
                             component.valid.size() == 7 && component.massNumber == 0 &&
                             component.geant4Particle.empty();
                    }),
        "ACE/CRIS bin bounds, counts, masks, and statistical errors remain lossless");
  const auto& oxygen = cris.components[3];
  Check(oxygen.binLower.front() == 66.3 && oxygen.binUpper.back() == 235.8,
        "ACE/CRIS oxygen bounds come from the pinned V04 CDF labels");
  Check(std::abs(oxygen.value.front() - 7.375e-7) < 1.0e-18 &&
        oxygen.counts.front() == 139.0 &&
        std::abs(oxygen.statisticalUncertainty.front() - 6.255393356389787e-8) < 1.0e-18,
        "ACE/CRIS oxygen value/count pair gives the documented Poisson error");
  Evidence("ace_cris_oxygen_first_bin_flux", 7.375e-7,
           oxygen.value.front(), "1/(cm2 sr s (MeV/n))", "independent_observation_oracle");
  Evidence("ace_cris_oxygen_first_bin_statistical_uncertainty", 6.255393356389787e-8,
           oxygen.statisticalUncertainty.front(), "1/(cm2 sr s (MeV/n))", "analytic_counts");
  std::size_t invalidCrisBins = 0;
  for (const auto& component : cris.components)
    invalidCrisBins += static_cast<std::size_t>(
        std::count(component.valid.begin(), component.valid.end(), false));
  Check(invalidCrisBins == 19,
        "ACE/CRIS zero-count bins are explicitly masked instead of assigned zero error");
  Check(shield_sim::Json::Parse(cris.canonicalMetadata.DumpCanonical()).DumpCanonical() ==
        cris.canonicalMetadata.DumpCanonical(),
        "ACE/CRIS canonical dataset round-trips exactly");

  const auto goes = shield_sim::GOESRSGPSAdapter::Load(
      root + "/tests/reference_spectra/sep/goes16_sgps_gle72_timeseries.csv",
      root + "/tests/reference_spectra/sep/goes16_gle72_event_definition.v1.json");
  Check(goes.Records().size() == 1181u * 2u * 13u,
        "GOES adapter retains every time/sensor/channel record");
  Check(goes.Records().front().midpointUnixSeconds == goes.Records().front().startUnixSeconds + 150 &&
        goes.Records().front().endUnixSeconds == goes.Records().front().startUnixSeconds + 300,
        "GOES source start timestamp maps to the requested midpoint exactly");
  Check(goes.Windows().size() == 2 &&
        goes.Windows()[0].endUnixSeconds - goes.Windows()[0].startUnixSeconds == 195300 &&
        goes.Windows()[1].endUnixSeconds - goes.Windows()[1].startUnixSeconds == 353100,
        "GOES primary and sensitivity windows retain published boundaries");
  const auto fluence = goes.Integrate();
  Check(fluence.channels.size() == 52, "GOES integration keeps both sensors and all channels separate");
  std::set<std::string> directions;
  bool foundIncomplete = false;
  bool foundRejected = false;
  bool checkedPublishedSpot = false;
  double publishedP7Fluence = 0.0;
  double publishedP7ValidSeconds = 0.0;
  double publishedP7RejectedSeconds = 0.0;
  for (const auto& channel : fluence.channels) {
    directions.insert(channel.spacecraftAxis);
    foundIncomplete = foundIncomplete || !channel.complete;
    foundRejected = foundRejected || channel.rejectedSeconds > 0.0;
    Check(channel.validSeconds + channel.rejectedSeconds + channel.missingSeconds ==
          channel.expectedSeconds, "GOES time ledger closes for " + channel.windowId + '/' +
          channel.sensorId + '/' + channel.channel);
    if (channel.windowId == "gle72_gt100MeV_primary" &&
        channel.sensorId == "sgps_0" && channel.channel == "P7") {
      // Independently summed from the pinned CSV by the dependency-free
      // Python reference test. This catches production-parser/integrator
      // agreement errors that a purely synthetic test could miss.
      Check(Near(channel.partialFluencePerCm2SrMeV, 402238.9858300105, 1.0e-11) &&
            Near(channel.validSeconds, 193800.0) &&
            Near(channel.rejectedSeconds, 1500.0) && Near(channel.missingSeconds, 0.0),
            "GOES published-row primary-window spot check");
      publishedP7Fluence = channel.partialFluencePerCm2SrMeV;
      publishedP7ValidSeconds = channel.validSeconds;
      publishedP7RejectedSeconds = channel.rejectedSeconds;
      checkedPublishedSpot = true;
    }
  }
  Check(directions == std::set<std::string>({"+X","-X"}),
        "GOES directions survive event integration");
  Check(foundIncomplete && foundRejected,
        "GOES flagged samples produce explicit incomplete partial fluences");
  Check(checkedPublishedSpot, "GOES published spot-check channel is present");
  Evidence("goes_primary_sgps0_p7_fluence", 402238.9858300105,
           publishedP7Fluence, "1/(cm2 sr MeV)", "independent_observation_oracle");
  Evidence("goes_primary_sgps0_p7_valid_time", 193800.0,
           publishedP7ValidSeconds, "s", "independent_observation_oracle");
  Evidence("goes_primary_sgps0_p7_rejected_time", 1500.0,
           publishedP7RejectedSeconds, "s", "independent_observation_oracle");
  Check(shield_sim::Json::Parse(fluence.canonicalJson.DumpCanonical()).DumpCanonical() ==
        fluence.canonicalJson.DumpCanonical(), "GOES fluence result round-trips exactly");

  // Exact analytic integration checks exercise the production kernel without
  // relying on the published event values. A five-minute average represents
  // its entire support, so fluence is value times the overlap duration.
  const auto constant = shield_sim::GOESRSGPSAdapter::IntegrateRecords(
      {Record(1000, 2.0)}, {Window(1100, 1200)});
  Check(constant.size() == 1 && Near(constant[0].partialFluencePerCm2SrMeV, 200.0) &&
        constant[0].complete,
        "analytic partial boundary bin uses exact overlap duration");
  Evidence("goes_partial_boundary_fluence", 200.0,
           constant[0].partialFluencePerCm2SrMeV, "1/(cm2 sr MeV)", "analytic");

  const auto linear = shield_sim::GOESRSGPSAdapter::IntegrateRecords(
      {Record(0, 1.0), Record(300, 3.0)}, {Window(0, 600)});
  Check(linear.size() == 1 && Near(linear[0].partialFluencePerCm2SrMeV, 1200.0) &&
        linear[0].complete,
        "analytic midpoint samples integrate a linear trend exactly");
  Evidence("goes_linear_midpoint_fluence", 1200.0,
           linear[0].partialFluencePerCm2SrMeV, "1/(cm2 sr MeV)", "analytic");

  const auto missing = shield_sim::GOESRSGPSAdapter::IntegrateRecords(
      {Record(0, 1.0)}, {Window(0, 600)});
  Check(Near(missing[0].validSeconds, 300.0) && Near(missing[0].missingSeconds, 300.0) &&
        !missing[0].complete,
        "missing interval is reported and never interpolated");
  Evidence("goes_missing_interval_ledger", 300.0, missing[0].missingSeconds,
           "s", "analytic");

  auto rejectedRecord = Record(300, 3.0);
  rejectedRecord.dynamicErrorFlaggedSamples = 1;
  const auto rejected = shield_sim::GOESRSGPSAdapter::IntegrateRecords(
      {Record(0, 1.0), rejectedRecord}, {Window(0, 600)});
  Check(Near(rejected[0].validSeconds, 300.0) && Near(rejected[0].rejectedSeconds, 300.0) &&
        Near(rejected[0].missingSeconds, 0.0) && !rejected[0].complete,
        "flagged interval is rejected without being mislabeled as absent");

  // Every individual quality requirement is independently capable of
  // rejecting a record. This guards against accidentally weakening a
  // conjunction while refactoring the published-data policy.
  std::vector<shield_sim::GoesSgpsRecord> qualityMutations;
  auto shortAverage = Record(0, 1.0); shortAverage.l1bRecordsInAverage = 299;
  auto shortValid = Record(0, 1.0); shortValid.validOneSecondSamples = 299;
  auto deadTime = Record(0, 1.0); deadTime.deadTimeFlaggedSamples = 1;
  auto outOfBand = Record(0, 1.0); outOfBand.outOfBandFlaggedSamples = 1;
  auto dynamicError = Record(0, 1.0); dynamicError.dynamicErrorFlaggedSamples = 1;
  auto noFlux = Record(0, 1.0); noFlux.hasFlux = false;
  auto noUncertainty = Record(0, 1.0); noUncertainty.hasUncertainty = false;
  qualityMutations = {shortAverage, shortValid, deadTime, outOfBand,
                      dynamicError, noFlux, noUncertainty};
  Check(std::all_of(qualityMutations.begin(), qualityMutations.end(),
                    [](const auto& record) { return !record.IsUsable(); }),
        "every GOES quality requirement independently rejects a sample");

  CheckThrows([&] {
    const auto duplicate = Record(0, 2.0);
    shield_sim::GOESRSGPSAdapter::IntegrateRecords({Record(0, 1.0), duplicate}, {Window(0, 300)});
  }, "duplicate time/sensor/channel bin is rejected");

  if (failures) {
    std::cerr << failures << " Stage-2 source-adapter test(s) failed\n";
    return 1;
  }
  if (argc == 2) {
    namespace fs = std::filesystem;
    fs::create_directories(argv[1]);
    shield_sim::Json evidence(shield_sim::Json::Object{
        {"comparisons", shield_sim::Json(std::move(evidenceComparisons))},
        {"schema_version", "shieldSim.validation_evidence.v1"},
        {"test_case", "stage2_source_adapters"}});
    std::ofstream output(fs::path(argv[1]) / "stage2_source_adapters.json",
                         std::ios::binary | std::ios::trunc);
    if (!output) {
      std::cerr << "cannot write Stage-2 validation evidence\n";
      return 2;
    }
    output << evidence.DumpCanonical() << '\n';
    if (!output) {
      std::cerr << "failed writing Stage-2 validation evidence\n";
      return 2;
    }
  }
  std::cout << "All Stage-2 source-adapter tests passed\n";
  return 0;
}
