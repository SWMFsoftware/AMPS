#include "CaseDefinition.hh"
#include "SpectrumDataset.hh"
#include "StrictJson.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

int failures = 0;

void Check(bool condition, const std::string& message) {
  if (condition) return;
  ++failures;
  std::cerr << "FAIL: " << message << '\n';
}

bool Near(double actual, double expected, double relative = 1.0e-12) {
  return std::abs(actual - expected) <=
      relative * std::max({1.0, std::abs(actual), std::abs(expected)});
}

void WriteJson(const fs::path& path, const shield_sim::Json& value) {
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot write " + path.string());
  output << value.DumpCanonical() << '\n';
}

template<class Function>
void ExpectFailure(Function&& function, const std::string& expected,
                   const std::string& message) {
  try {
    function();
    Check(false, message + " (no exception)");
  } catch (const std::exception& error) {
    Check(std::string(error.what()).find(expected) != std::string::npos,
          message + " (unexpected message: " + error.what() + ')');
  }
}

shield_sim::Json ExchangeCaseWith(const fs::path& sourceDirectory,
                                  const fs::path& spectrumPath,
                                  const std::string& componentId) {
  shield_sim::Json value = shield_sim::Json::ParseFile(
      (sourceDirectory / "examples/stage11_spectrum_exchange_case.json").string());
  auto& root = value.AsObject("case");
  auto& source = root.at("source").AsObject("case.source");
  auto& spectrum = source.at("spectrum").AsObject("case.source.spectrum");
  spectrum["path"] = spectrumPath.string();
  spectrum["component_id"] = componentId;
  return value;
}

} // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: stage11_spectrum_exchange_tests SOURCE_DIR EVIDENCE_DIR\n";
    return 2;
  }
  const fs::path sourceDirectory = fs::absolute(argv[1]);
  const fs::path evidenceDirectory = fs::absolute(argv[2]);
  fs::create_directories(evidenceDirectory);
  const fs::path fixture = sourceDirectory /
      "tests/data/stage11_spectrum_exchange.v1.json";

  const auto dataset = shield_sim::SpectrumDataset::LoadSpectrumExchangeV1(
      fixture.string());
  Check(dataset.id == "stage11_two_component_verification_fixture",
        "dataset identity is retained");
  Check(dataset.components.size() == 2, "both source components are retained");
  Check(dataset.qualityStatus == "accepted" &&
        Near(dataset.maximumUnresolvedFraction, 0.01),
        "aggregate AMPS trajectory quality is retained");
  const auto& proton = dataset.components.at(0);
  Check(proton.representative.size() == 3 && proton.binLower.empty(),
        "point samples are not reinterpreted as bins during ingestion");
  Check(proton.lowerBound.size() == 3 && proton.upperBound.size() == 3 &&
        proton.statisticalUncertainty.empty(),
        "conservative envelopes remain distinct from statistical uncertainty");
  Check(proton.sampledDirections.at(1) == 100.0 &&
        proton.resolvedDirections.at(1) == 99.0 &&
        Near(proton.unresolvedFraction.at(1), 0.01),
        "per-node trajectory resolution diagnostics are retained");

  const fs::path protonCasePath = evidenceDirectory / "stage11_proton_case.json";
  WriteJson(protonCasePath, ExchangeCaseWith(sourceDirectory, fixture,
                                             "local_proton_points"));
  const auto protonCase = shield_sim::CaseDefinition::Load(protonCasePath.string());
  Check(protonCase->sourceFormat == "spectrum_exchange_v1" &&
        protonCase->sourceComponentId == "local_proton_points",
        "case parser binds the requested exchange component");
  Check(protonCase->spectrum.transportInterpolation ==
        "log_log_positive_linear_zero",
        "point interpolation survives transport admission");
  Check(protonCase->spectrum.transportBinLowerMeV.size() == 2 &&
        Near(protonCase->spectrum.transportBinLowerMeV.at(0), 10.0) &&
        Near(protonCase->spectrum.transportBinUpperMeV.at(1), 1000.0),
        "point support becomes exactly two declared transport segments");
  Check(Near(protonCase->spectrum.transportValueLowerPerMeV.at(0), 0.1) &&
        Near(protonCase->spectrum.transportValueUpperPerMeV.at(1), 0.0),
        "transport retains the exact point intensities, including physical zero");

  // The alpha component is stated per nucleon.  Admission must map [10,20]
  // MeV/n and 40 m^-2/(MeV/n) to [40,80] total MeV and
  // 0.001 cm^-2/MeV.  The physical integral is 400 m^-2 = 0.04 cm^-2.
  const fs::path alphaCasePath = evidenceDirectory / "stage11_alpha_case.json";
  WriteJson(alphaCasePath, ExchangeCaseWith(sourceDirectory, fixture,
                                            "local_alpha_bins"));
  const auto alphaCase = shield_sim::CaseDefinition::Load(alphaCasePath.string());
  Check(alphaCase->spectrum.transportInterpolation == "piecewise_constant" &&
        alphaCase->spectrum.transportBinLowerMeV.size() == 1,
        "bin-average alpha component uses piecewise-constant transport");
  Check(Near(alphaCase->spectrum.transportBinLowerMeV.at(0), 40.0) &&
        Near(alphaCase->spectrum.transportBinUpperMeV.at(0), 80.0) &&
        Near(alphaCase->spectrum.alphaValue.at(0), 0.001) &&
        Near(alphaCase->spectrum.alphaValue.at(0) *
             (alphaCase->spectrum.transportBinUpperMeV.at(0) -
              alphaCase->spectrum.transportBinLowerMeV.at(0)), 0.04),
        "per-nucleon and SI-area Jacobians conserve the alpha intensity integral");

  shield_sim::Json wrongVersion = shield_sim::Json::ParseFile(fixture.string());
  wrongVersion.AsObject("exchange")["schema_version"] =
      "shieldSim.spectrum_exchange.v2";
  const fs::path wrongVersionPath = evidenceDirectory / "stage11_wrong_version.json";
  WriteJson(wrongVersionPath, wrongVersion);
  ExpectFailure([&] {
    shield_sim::SpectrumDataset::LoadSpectrumExchangeV1(wrongVersionPath.string());
  }, "schema_version", "unknown exchange versions fail closed");

  shield_sim::Json ambiguous = shield_sim::Json::ParseFile(fixture.string());
  auto& ambiguousEnergy = ambiguous.AsObject("exchange").at("components")
      .AsArray("components").at(0).AsObject("component").at("energy")
      .AsObject("energy");
  ambiguousEnergy.erase("interpolation");
  const fs::path ambiguousPath = evidenceDirectory / "stage11_missing_interpolation.json";
  WriteJson(ambiguousPath, ambiguous);
  ExpectFailure([&] {
    shield_sim::SpectrumDataset::LoadSpectrumExchangeV1(ambiguousPath.string());
  }, "interpolation", "missing interpolation is never inferred");

  shield_sim::Json invalidUtc = shield_sim::Json::ParseFile(fixture.string());
  invalidUtc.AsObject("exchange").at("context").AsObject("context").at("time")
      .AsObject("time")["reference_utc"] = "2026-02-30T00:00:00Z";
  const fs::path invalidUtcPath = evidenceDirectory / "stage11_invalid_utc.json";
  WriteJson(invalidUtcPath, invalidUtc);
  ExpectFailure([&] {
    shield_sim::SpectrumDataset::LoadSpectrumExchangeV1(invalidUtcPath.string());
  }, "invalid UTC calendar", "calendar-invalid UTC is rejected");

  shield_sim::Json invalid = shield_sim::Json::ParseFile(fixture.string());
  invalid.AsObject("exchange").at("components").AsArray("components").at(0)
      .AsObject("component").at("data").AsArray("data").at(1)
      .AsObject("row")["valid"] = false;
  const fs::path invalidPath = evidenceDirectory / "stage11_invalid_node.json";
  WriteJson(invalidPath, invalid);
  const fs::path invalidCasePath = evidenceDirectory / "stage11_invalid_case.json";
  WriteJson(invalidCasePath, ExchangeCaseWith(sourceDirectory, invalidPath,
                                              "local_proton_points"));
  ExpectFailure([&] { shield_sim::CaseDefinition::Load(invalidCasePath.string()); },
                "invalid sample", "transport never drops or bridges invalid nodes");

  shield_sim::Json::Object evidence{
    {"schema_version", "shieldSim.validation_evidence.stage11.v1"},
    {"status", failures == 0 ? "PASS" : "FAIL"},
    {"checks", 15},
    {"fixture", fixture.string()},
    {"fixture_content_sha256", dataset.contentHash},
    {"proton_segment_count", static_cast<double>(
        protonCase->spectrum.transportBinLowerMeV.size())},
    {"alpha_integral_input_m^-2_s^-1_sr^-1", 400.0},
    {"alpha_integral_internal_cm^-2_s^-1_sr^-1", 0.04}
  };
  WriteJson(evidenceDirectory / "stage11_spectrum_exchange_core.json",
            shield_sim::Json(std::move(evidence)));
  if (failures != 0) return 1;
  std::cout << "Stage-11 spectrum exchange tests passed\n";
  return 0;
}
