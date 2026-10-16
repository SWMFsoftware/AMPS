#include "SpectrumDataset.hh"

#include "ContentHash.hh"

#include <cmath>
#include <fstream>
#include <sstream>

namespace shield_sim {
namespace {

std::string RequiredString(const Json::Object& object, const std::string& field,
                           const std::string& context) {
  const auto found = object.find(field);
  if (found == object.end()) throw JsonError(context + " is missing required field '" + field + "'");
  const std::string result = found->second.AsString(context + '.' + field);
  if (result.empty()) throw JsonError(context + '.' + field + " must not be empty");
  return result;
}

double PositiveNumber(const Json& value, const std::string& context) {
  const double result = value.AsNumber(context);
  if (!std::isfinite(result) || result <= 0.0)
    throw JsonError(context + " must be finite and > 0");
  return result;
}

} // namespace

SpectrumDataset SpectrumDataset::LoadLegacyThreeColumn(const std::string& spectrumPath,
                                                       const std::string& metadataPath) {
  if (metadataPath.empty())
    throw JsonError("legacy three-column spectrum '" + spectrumPath + "' requires an explicit metadata sidecar");
  SpectrumDataset dataset;
  dataset.spectrumPath = spectrumPath;
  dataset.metadataPath = metadataPath;
  dataset.canonicalMetadata = Json::ParseFile(metadataPath);
  const auto& root = dataset.canonicalMetadata.AsObject("spectrum metadata");
  RejectUnknownFields(root,
    {"schema_version","dataset_id","format","columns","energy_variable","energy_unit",
     "value_type","value_unit","angular_model","bin_representation","edge_rule",
     "single_bin_width","science_status","provenance"}, "spectrum metadata");
  if (RequiredString(root, "schema_version", "spectrum metadata") != "shieldSim.spectrum_metadata.v1")
    throw JsonError("spectrum metadata.schema_version must be 'shieldSim.spectrum_metadata.v1'");
  if (RequiredString(root, "format", "spectrum metadata") != "legacy_three_column")
    throw JsonError("spectrum metadata.format must be 'legacy_three_column'");
  dataset.id = RequiredString(root, "dataset_id", "spectrum metadata");
  dataset.energyVariable = RequiredString(root, "energy_variable", "spectrum metadata");
  if (dataset.energyVariable != "total_kinetic_energy")
    throw JsonError("legacy Stage-1 adapter requires energy_variable='total_kinetic_energy'");
  dataset.energyUnit = RequiredString(root, "energy_unit", "spectrum metadata");
  if (dataset.energyUnit != "MeV") throw JsonError("legacy Stage-1 adapter requires energy_unit='MeV'");
  dataset.valueType = RequiredString(root, "value_type", "spectrum metadata");
  if (dataset.valueType != "differential_rate" &&
      dataset.valueType != "differential_intensity" &&
      dataset.valueType != "arbitrary_test_weight")
    throw JsonError("spectrum metadata.value_type is unsupported");
  dataset.valueUnit = RequiredString(root, "value_unit", "spectrum metadata");
  dataset.angularModel = RequiredString(root, "angular_model", "spectrum metadata");
  if (dataset.angularModel != "beam" && dataset.angularModel != "isotropic" &&
      dataset.angularModel != "declared_by_case")
    throw JsonError("spectrum metadata.angular_model must be beam, isotropic, or declared_by_case");
  dataset.scienceStatus = RequiredString(root, "science_status", "spectrum metadata");
  if (dataset.scienceStatus != "synthetic_test_only" &&
      dataset.scienceStatus != "published_adapter_pending" &&
      dataset.scienceStatus != "user_supplied_unvalidated")
    throw JsonError("spectrum metadata.science_status is unsupported");

  const auto& columns = dataset.canonicalMetadata.At("columns", "spectrum metadata")
                            .AsArray("spectrum metadata.columns");
  const std::vector<std::string> expected = {"energy", "proton", "alpha"};
  if (columns.size() != expected.size())
    throw JsonError("spectrum metadata.columns must be exactly [energy, proton, alpha]");
  for (std::size_t i = 0; i < expected.size(); ++i)
    if (columns[i].AsString("spectrum metadata.columns[]") != expected[i])
      throw JsonError("spectrum metadata.columns must be exactly [energy, proton, alpha]");
  if (RequiredString(root, "bin_representation", "spectrum metadata") != "centers")
    throw JsonError("legacy Stage-1 adapter supports only bin_representation='centers'");
  if (RequiredString(root, "edge_rule", "spectrum metadata") != "geometric_midpoint_edges")
    throw JsonError("legacy Stage-1 adapter requires edge_rule='geometric_midpoint_edges'");
  if (!root.count("provenance") || !root.at("provenance").IsObject())
    throw JsonError("spectrum metadata.provenance must be an object");

  std::ifstream input(spectrumPath);
  if (!input) throw JsonError("cannot open spectrum file '" + spectrumPath + "'");
  std::string line;
  std::size_t lineNumber = 0;
  bool anyPositive = false;
  while (std::getline(input, line)) {
    ++lineNumber;
    const std::size_t first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') continue;
    std::istringstream row(line);
    double energy = 0.0, proton = 0.0, alpha = 0.0;
    std::string extra;
    if (!(row >> energy >> proton >> alpha) || (row >> extra))
      throw JsonError(spectrumPath + ':' + std::to_string(lineNumber) +
                      ": expected exactly three numeric columns");
    if (!std::isfinite(energy) || energy <= 0.0)
      throw JsonError(spectrumPath + ':' + std::to_string(lineNumber) + ": energy must be finite and > 0");
    if (!std::isfinite(proton) || proton < 0.0 || !std::isfinite(alpha) || alpha < 0.0)
      throw JsonError(spectrumPath + ':' + std::to_string(lineNumber) +
                      ": spectrum values must be finite and >= 0; fill values are not accepted");
    if (!dataset.energyMeV.empty() && energy <= dataset.energyMeV.back())
      throw JsonError(spectrumPath + ':' + std::to_string(lineNumber) +
                      ": energy centers must be strictly increasing; rows are never silently reordered");
    dataset.energyMeV.push_back(energy);
    dataset.protonValue.push_back(proton);
    dataset.alphaValue.push_back(alpha);
    anyPositive = anyPositive || proton > 0.0 || alpha > 0.0;
  }
  if (dataset.energyMeV.empty()) throw JsonError("spectrum file '" + spectrumPath + "' has no data rows");
  if (!anyPositive) throw JsonError("spectrum file '" + spectrumPath + "' has zero weight for every species");

  dataset.binWidthMeV.resize(dataset.energyMeV.size());
  std::vector<double> canonicalEdges(dataset.energyMeV.size() + 1, 0.0);
  if (dataset.energyMeV.size() == 1) {
    if (!root.count("single_bin_width"))
      throw JsonError("single-row spectrum requires spectrum metadata.single_bin_width");
    const auto& width = root.at("single_bin_width").AsObject("spectrum metadata.single_bin_width");
    RejectUnknownFields(width, {"value","unit"}, "spectrum metadata.single_bin_width");
    if (RequiredString(width, "unit", "spectrum metadata.single_bin_width") != "MeV")
      throw JsonError("spectrum metadata.single_bin_width.unit must be 'MeV'");
    dataset.binWidthMeV[0] = PositiveNumber(width.at("value"), "spectrum metadata.single_bin_width.value");
    canonicalEdges[0] = dataset.energyMeV[0] - dataset.binWidthMeV[0] / 2.0;
    canonicalEdges[1] = dataset.energyMeV[0] + dataset.binWidthMeV[0] / 2.0;
    if (canonicalEdges[0] < 0.0)
      throw JsonError("single-row spectrum width extends below zero energy");
  } else {
    for (std::size_t i = 1; i < dataset.energyMeV.size(); ++i)
      canonicalEdges[i] = std::sqrt(dataset.energyMeV[i - 1] * dataset.energyMeV[i]);
    canonicalEdges.front() = dataset.energyMeV.front() * dataset.energyMeV.front() / canonicalEdges[1];
    canonicalEdges.back() = dataset.energyMeV.back() * dataset.energyMeV.back() /
                            canonicalEdges[canonicalEdges.size() - 2];
    for (std::size_t i = 0; i < dataset.energyMeV.size(); ++i) {
      dataset.binWidthMeV[i] = canonicalEdges[i + 1] - canonicalEdges[i];
      if (!std::isfinite(dataset.binWidthMeV[i]) || dataset.binWidthMeV[i] <= 0.0)
        throw JsonError("derived spectrum bin width is invalid at row " + std::to_string(i));
    }
  }

  // Populate the general component representation as well as the historical
  // shared-grid vectors.  Later adapters use only components, which supports
  // independent species grids and directional components without information
  // loss.  Keeping this conversion here gives every source reader the same
  // canonical boundary while preserving the Stage-1 generator API.
  auto makeComponent = [&](const std::string& componentId,
                           const std::string& speciesId,
                           const std::string& particle,
                           int z, int a,
                           const std::vector<double>& values) {
    SpectrumComponent component;
    component.componentId = componentId;
    component.speciesId = speciesId;
    component.geant4Particle = particle;
    component.atomicNumber = z;
    component.massNumber = a;
    component.energyVariable = dataset.energyVariable;
    component.energyUnit = dataset.energyUnit;
    component.valueType = dataset.valueType;
    component.valueUnit = dataset.valueUnit;
    component.angularRepresentation = dataset.angularModel;
    component.directionId = dataset.angularModel;
    component.uncertaintySemantics = "not_supplied";
    component.representative = dataset.energyMeV;
    component.value = values;
    component.binLower.assign(canonicalEdges.begin(), canonicalEdges.end() - 1);
    component.binUpper.assign(canonicalEdges.begin() + 1, canonicalEdges.end());
    return component;
  };
  dataset.components.push_back(makeComponent("proton", "proton", "proton", 1, 1,
                                             dataset.protonValue));
  dataset.components.push_back(makeComponent("alpha", "alpha", "alpha", 2, 4,
                                             dataset.alphaValue));

  std::ostringstream numeric;
  numeric.precision(17);
  for (std::size_t i = 0; i < dataset.energyMeV.size(); ++i)
    numeric << dataset.energyMeV[i] << ',' << dataset.protonValue[i] << ','
            << dataset.alphaValue[i] << ',' << dataset.binWidthMeV[i] << '\n';
  dataset.contentHash = Sha256(dataset.canonicalMetadata.DumpCanonical() + "\n" + numeric.str());
  return dataset;
}

} // namespace shield_sim
