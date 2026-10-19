#include "SpectrumDataset.hh"

#include "ContentHash.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <regex>
#include <set>
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

double NonNegativeNumber(const Json& value, const std::string& context) {
  const double result = value.AsNumber(context);
  if (!std::isfinite(result) || result < 0.0)
    throw JsonError(context + " must be finite and >= 0");
  return result;
}

long NonNegativeInteger(const Json& value, const std::string& context) {
  const double result = NonNegativeNumber(value, context);
  if (std::floor(result) != result ||
      result > static_cast<double>(std::numeric_limits<long>::max()))
    throw JsonError(context + " must be a non-negative integer in range");
  return static_cast<long>(result);
}

long double UtcOrdinalSeconds(const std::string& value,
                              const std::string& context) {
  // The exchange boundary deliberately accepts only UTC, not a local offset
  // whose interpretation could change when the file moves to another host.
  // Validate calendar fields as well as spelling: a regex-only check would
  // accidentally accept dates such as month 19 or February 31.
  static const std::regex utc(
      R"(^([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2})(\.[0-9]+)?Z$)");
  std::smatch match;
  if (!std::regex_match(value, match, utc))
    throw JsonError(context + " must be an ISO-8601 UTC timestamp ending in Z");
  const int year = std::stoi(match[1].str());
  const int month = std::stoi(match[2].str());
  const int day = std::stoi(match[3].str());
  const int hour = std::stoi(match[4].str());
  const int minute = std::stoi(match[5].str());
  const int second = std::stoi(match[6].str());
  const bool leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
  const int daysByMonth[12] = {31, leap ? 29 : 28, 31, 30, 31, 30,
                               31, 31, 30, 31, 30, 31};
  if (year == 0 || month < 1 || month > 12 || day < 1 ||
      day > daysByMonth[month - 1] || hour > 23 || minute > 59 || second > 59)
    throw JsonError(context + " contains an invalid UTC calendar/time field");

  // Gregorian ordinal sufficient for ordering the explicitly UTC instants.
  long long days = 0;
  for (int y = 1; y < year; ++y)
    days += (y % 4 == 0 && (y % 100 != 0 || y % 400 == 0)) ? 366 : 365;
  for (int m = 1; m < month; ++m) days += daysByMonth[m - 1];
  days += day - 1;
  long double fraction = 0.0L;
  if (match[7].matched) fraction = std::stold("0" + match[7].str());
  return static_cast<long double>(days) * 86400.0L + hour * 3600.0L +
         minute * 60.0L + second + fraction;
}

std::string RequiredUtc(const Json::Object& object, const std::string& field,
                        const std::string& context) {
  const std::string value = RequiredString(object, field, context);
  UtcOrdinalSeconds(value, context + '.' + field);
  return value;
}

bool IsLowerHexSha256(const std::string& value) {
  return value.size() == 64 &&
      std::all_of(value.begin(), value.end(), [](unsigned char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
      });
}

void CheckOptionalFingerprint(const Json::Object& object, const std::string& field,
                              const std::string& context) {
  const auto found = object.find(field);
  if (found == object.end()) return;
  const std::string value = found->second.AsString(context + '.' + field);
  if (!IsLowerHexSha256(value))
    throw JsonError(context + '.' + field +
                    " must be a 64-character lowercase SHA-256 hex string");
}

void CheckParticleIdentity(const std::string& name, int z, int a,
                           const std::string& context) {
  if (a < z) throw JsonError(context + " requires A >= Z");
  if (name == "proton" && (z != 1 || a != 1))
    throw JsonError(context + " proton identity requires Z=1 and A=1");
  if (name == "alpha" && (z != 2 || a != 4))
    throw JsonError(context + " alpha identity requires Z=2 and A=4");
  if (name == "neutron" && (z != 0 || a != 1))
    throw JsonError(context + " neutron identity requires Z=0 and A=1");
  if (name == "gamma" && (z != 0 || a != 0))
    throw JsonError(context + " gamma identity requires Z=0 and A=0");
  if (name == "ion" && (z <= 0 || a <= 0))
    throw JsonError(context + " ion identity requires explicit positive Z and A");
  if (name != "proton" && name != "alpha" && name != "neutron" &&
      name != "gamma" && name != "ion")
    throw JsonError(context +
                    " name must be proton, alpha, neutron, gamma, or ion");
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

SpectrumDataset SpectrumDataset::LoadSpectrumExchangeV1(
    const std::string& spectrumPath) {
  SpectrumDataset dataset;
  dataset.spectrumPath = spectrumPath;
  dataset.canonicalMetadata = Json::ParseFile(spectrumPath);
  const auto& root = dataset.canonicalMetadata.AsObject("spectrum exchange");
  RejectUnknownFields(root,
      {"schema_version","dataset_id","science_status","producer","context",
       "quality","components"}, "spectrum exchange");
  if (RequiredString(root, "schema_version", "spectrum exchange") !=
      "shieldSim.spectrum_exchange.v1")
    throw JsonError("spectrum exchange.schema_version must be "
                    "'shieldSim.spectrum_exchange.v1'");
  dataset.id = RequiredString(root, "dataset_id", "spectrum exchange");
  dataset.scienceStatus = RequiredString(root, "science_status", "spectrum exchange");
  if (dataset.scienceStatus != "calculated_unvalidated" &&
      dataset.scienceStatus != "calculated_validation_qualified" &&
      dataset.scienceStatus != "verification_fixture")
    throw JsonError("spectrum exchange.science_status is unsupported");

  const auto& producer = root.at("producer").AsObject("spectrum exchange.producer");
  RejectUnknownFields(producer,
      {"software","software_version","code_revision","run_id","created_utc",
       "calculation_mode","source_product","input_fingerprint",
       "field_snapshot_id","field_snapshot_fingerprint",
       "boundary_spectrum_fingerprint","product_control_fingerprint","description"},
      "spectrum exchange.producer");
  RequiredString(producer, "software", "spectrum exchange.producer");
  RequiredString(producer, "software_version", "spectrum exchange.producer");
  RequiredString(producer, "code_revision", "spectrum exchange.producer");
  RequiredString(producer, "run_id", "spectrum exchange.producer");
  RequiredUtc(producer, "created_utc", "spectrum exchange.producer");
  const std::string calculationMode = RequiredString(
      producer, "calculation_mode", "spectrum exchange.producer");
  if (calculationMode != "gridless" && calculationMode != "mode3d" &&
      calculationMode != "verification_fixture")
    throw JsonError("spectrum exchange.producer.calculation_mode must be gridless, mode3d, or verification_fixture");
  RequiredString(producer, "source_product", "spectrum exchange.producer");
  CheckOptionalFingerprint(producer, "input_fingerprint", "spectrum exchange.producer");
  CheckOptionalFingerprint(producer, "field_snapshot_fingerprint", "spectrum exchange.producer");
  CheckOptionalFingerprint(producer, "boundary_spectrum_fingerprint", "spectrum exchange.producer");
  CheckOptionalFingerprint(producer, "product_control_fingerprint", "spectrum exchange.producer");
  if (calculationMode != "verification_fixture") {
    if (!producer.count("input_fingerprint") ||
        !producer.count("boundary_spectrum_fingerprint"))
      throw JsonError("gridless/mode3d spectrum producers require input_fingerprint and boundary_spectrum_fingerprint");
  }
  if (calculationMode == "mode3d") {
    RequiredString(producer, "field_snapshot_id", "spectrum exchange.producer");
    if (!producer.count("field_snapshot_fingerprint") ||
        !producer.count("product_control_fingerprint"))
      throw JsonError("mode3d spectrum producers require field_snapshot_id, field_snapshot_fingerprint, and product_control_fingerprint");
  } else if (producer.count("field_snapshot_id") ||
             producer.count("field_snapshot_fingerprint")) {
    throw JsonError("field snapshot identity/fingerprint are permitted only for calculation_mode='mode3d'");
  }
  if (producer.count("description"))
    RequiredString(producer, "description", "spectrum exchange.producer");

  const auto& context = root.at("context").AsObject("spectrum exchange.context");
  RejectUnknownFields(context, {"time","location"}, "spectrum exchange.context");
  const auto& time = context.at("time").AsObject("spectrum exchange.context.time");
  RejectUnknownFields(time, {"interpretation","reference_utc","start_utc","end_utc"},
                      "spectrum exchange.context.time");
  const std::string timeInterpretation = RequiredString(
      time, "interpretation", "spectrum exchange.context.time");
  const std::string referenceUtc = RequiredUtc(
      time, "reference_utc", "spectrum exchange.context.time");
  if (timeInterpretation == "instantaneous_quasi_static") {
    if (time.count("start_utc") || time.count("end_utc"))
      throw JsonError("instantaneous_quasi_static spectrum time must not declare an interval");
  } else if (timeInterpretation == "interval_average") {
    const std::string startUtc = RequiredUtc(
        time, "start_utc", "spectrum exchange.context.time");
    const std::string endUtc = RequiredUtc(
        time, "end_utc", "spectrum exchange.context.time");
    const long double start = UtcOrdinalSeconds(
        startUtc, "spectrum exchange.context.time.start_utc");
    const long double end = UtcOrdinalSeconds(
        endUtc, "spectrum exchange.context.time.end_utc");
    const long double reference = UtcOrdinalSeconds(
        referenceUtc, "spectrum exchange.context.time.reference_utc");
    if (!(start < end) || reference < start || reference > end)
      throw JsonError("interval_average requires start_utc < end_utc and reference_utc inside the closed interval");
  } else {
    throw JsonError("spectrum exchange.context.time.interpretation must be "
                    "instantaneous_quasi_static or interval_average");
  }
  const auto& location = context.at("location").AsObject(
      "spectrum exchange.context.location");
  RejectUnknownFields(location, {"location_id","frame","unit","position"},
                      "spectrum exchange.context.location");
  RequiredString(location, "location_id", "spectrum exchange.context.location");
  RequiredString(location, "frame", "spectrum exchange.context.location");
  if (RequiredString(location, "unit", "spectrum exchange.context.location") != "km")
    throw JsonError("spectrum exchange.context.location.unit must be 'km'");
  const auto& position = location.at("position").AsArray(
      "spectrum exchange.context.location.position");
  if (position.size() != 3)
    throw JsonError("spectrum exchange.context.location.position must contain x, y, z");
  for (std::size_t i = 0; i < position.size(); ++i) {
    const double coordinate = position[i].AsNumber(
        "spectrum exchange.context.location.position[]");
    if (!std::isfinite(coordinate))
      throw JsonError("spectrum exchange.context.location.position values must be finite");
  }

  const auto& quality = root.at("quality").AsObject("spectrum exchange.quality");
  RejectUnknownFields(quality,
      {"status","unresolved_semantics","maximum_unresolved_fraction",
       "acceptance_limit"}, "spectrum exchange.quality");
  dataset.qualityStatus = RequiredString(quality, "status", "spectrum exchange.quality");
  if (dataset.qualityStatus != "accepted" && dataset.qualityStatus != "rejected")
    throw JsonError("spectrum exchange.quality.status must be accepted or rejected");
  if (RequiredString(quality, "unresolved_semantics", "spectrum exchange.quality") !=
      "lower_forbidden_upper_allowed")
    throw JsonError("spectrum exchange.quality.unresolved_semantics must be "
                    "'lower_forbidden_upper_allowed'");
  dataset.maximumUnresolvedFraction = NonNegativeNumber(
      quality.at("maximum_unresolved_fraction"),
      "spectrum exchange.quality.maximum_unresolved_fraction");
  dataset.unresolvedFractionAcceptanceLimit = NonNegativeNumber(
      quality.at("acceptance_limit"), "spectrum exchange.quality.acceptance_limit");
  if (dataset.maximumUnresolvedFraction > 1.0 ||
      dataset.unresolvedFractionAcceptanceLimit > 1.0)
    throw JsonError("spectrum exchange unresolved fractions must be <= 1");
  if ((dataset.qualityStatus == "accepted") !=
      (dataset.maximumUnresolvedFraction <= dataset.unresolvedFractionAcceptanceLimit))
    throw JsonError("spectrum exchange quality status conflicts with its unresolved limit");

  const auto& components = root.at("components").AsArray("spectrum exchange.components");
  if (components.empty())
    throw JsonError("spectrum exchange.components must not be empty");
  std::set<std::string> componentIds;
  double observedMaximumUnresolved = 0.0;
  for (std::size_t componentIndex = 0; componentIndex < components.size(); ++componentIndex) {
    const std::string componentContext = "spectrum exchange.components[" +
        std::to_string(componentIndex) + ']';
    const auto& object = components[componentIndex].AsObject(componentContext);
    RejectUnknownFields(object,
        {"component_id","particle","energy","quantity","uncertainty","data"},
        componentContext);
    SpectrumComponent component;
    component.componentId = RequiredString(object, "component_id", componentContext);
    if (!componentIds.insert(component.componentId).second)
      throw JsonError(componentContext + ".component_id is duplicated");

    const auto& particle = object.at("particle").AsObject(componentContext + ".particle");
    RejectUnknownFields(particle, {"name","Z","A"}, componentContext + ".particle");
    component.speciesId = RequiredString(particle, "name", componentContext + ".particle");
    component.geant4Particle = component.speciesId;
    component.atomicNumber = static_cast<int>(NonNegativeInteger(
        particle.at("Z"), componentContext + ".particle.Z"));
    component.massNumber = static_cast<int>(NonNegativeInteger(
        particle.at("A"), componentContext + ".particle.A"));
    CheckParticleIdentity(component.speciesId, component.atomicNumber,
                          component.massNumber, componentContext + ".particle");

    const auto& energy = object.at("energy").AsObject(componentContext + ".energy");
    RejectUnknownFields(energy,
        {"coordinate","unit","representation","interpolation","support"},
        componentContext + ".energy");
    component.energyVariable = RequiredString(energy, "coordinate", componentContext + ".energy");
    component.energyUnit = RequiredString(energy, "unit", componentContext + ".energy");
    const std::string representation = RequiredString(
        energy, "representation", componentContext + ".energy");
    const std::string interpolation = RequiredString(
        energy, "interpolation", componentContext + ".energy");
    if (RequiredString(energy, "support", componentContext + ".energy") !=
        "closed_interval_no_extrapolation")
      throw JsonError(componentContext +
                      ".energy.support must be 'closed_interval_no_extrapolation'");
    if (component.energyVariable == "total_kinetic_energy") {
      if (component.energyUnit != "MeV")
        throw JsonError(componentContext + " total kinetic energy requires unit MeV");
    } else if (component.energyVariable == "kinetic_energy_per_nucleon") {
      if (component.energyUnit != "MeV/nucleon" || component.massNumber <= 0 ||
          component.speciesId == "gamma")
        throw JsonError(componentContext +
                        " per-nucleon energy requires a massive explicit-A particle and MeV/nucleon");
    } else {
      throw JsonError(componentContext + ".energy.coordinate is unsupported");
    }
    if ((representation == "point_samples" &&
         interpolation != "log_log_positive_linear_zero") ||
        (representation == "bin_averages" && interpolation != "piecewise_constant") ||
        (representation != "point_samples" && representation != "bin_averages"))
      throw JsonError(componentContext +
                      " has an inconsistent energy representation/interpolation pair");

    const auto& quantity = object.at("quantity").AsObject(componentContext + ".quantity");
    RejectUnknownFields(quantity, {"type","unit","angular_representation"},
                        componentContext + ".quantity");
    component.valueType = RequiredString(quantity, "type", componentContext + ".quantity");
    component.valueUnit = RequiredString(quantity, "unit", componentContext + ".quantity");
    component.angularRepresentation = RequiredString(
        quantity, "angular_representation", componentContext + ".quantity");
    component.directionId = component.angularRepresentation;
    if (component.valueType != "differential_directional_intensity")
      throw JsonError(componentContext +
                      ".quantity.type must be differential_directional_intensity");
    const std::string expectedUnit = component.energyVariable == "total_kinetic_energy"
        ? "m^-2 s^-1 sr^-1 MeV^-1"
        : "m^-2 s^-1 sr^-1 (MeV/nucleon)^-1";
    if (component.valueUnit != expectedUnit)
      throw JsonError(componentContext + ".quantity.unit must be '" + expectedUnit + "'");
    if (component.angularRepresentation != "isotropic")
      throw JsonError(componentContext +
                      ".quantity.angular_representation must be isotropic in schema v1");

    const auto& uncertainty = object.at("uncertainty").AsObject(
        componentContext + ".uncertainty");
    RejectUnknownFields(uncertainty, {"semantics","source"},
                        componentContext + ".uncertainty");
    component.uncertaintySemantics = RequiredString(
        uncertainty, "semantics", componentContext + ".uncertainty");
    if (component.uncertaintySemantics != "conservative_bounds" &&
        component.uncertaintySemantics != "not_supplied")
      throw JsonError(componentContext + ".uncertainty.semantics is unsupported");
    RequiredString(uncertainty, "source", componentContext + ".uncertainty");

    const auto& data = object.at("data").AsArray(componentContext + ".data");
    if ((representation == "point_samples" && data.size() < 2) || data.empty())
      throw JsonError(componentContext +
                      ".data requires at least two point samples or one bin");
    double previousCoordinate = -1.0;
    double previousUpper = -1.0;
    for (std::size_t rowIndex = 0; rowIndex < data.size(); ++rowIndex) {
      const std::string rowContext = componentContext + ".data[" +
          std::to_string(rowIndex) + ']';
      const auto& row = data[rowIndex].AsObject(rowContext);
      const std::vector<std::string> coordinateFields = representation == "point_samples"
          ? std::vector<std::string>{"energy","value","lower","upper","valid",
              "unresolved_fraction","n_sampled","n_resolved","n_allowed"}
          : std::vector<std::string>{"energy_lower","energy_upper","value","lower","upper",
              "valid","unresolved_fraction","n_sampled","n_resolved","n_allowed"};
      RejectUnknownFields(row, coordinateFields, rowContext);
      double lowerCoordinate = 0.0, upperCoordinate = 0.0;
      if (representation == "point_samples") {
        const double coordinate = PositiveNumber(row.at("energy"), rowContext + ".energy");
        if (coordinate <= previousCoordinate)
          throw JsonError(rowContext + ".energy must be strictly increasing");
        previousCoordinate = coordinate;
        component.representative.push_back(coordinate);
      } else {
        lowerCoordinate = PositiveNumber(row.at("energy_lower"), rowContext + ".energy_lower");
        upperCoordinate = PositiveNumber(row.at("energy_upper"), rowContext + ".energy_upper");
        if (upperCoordinate <= lowerCoordinate)
          throw JsonError(rowContext + " requires energy_upper > energy_lower");
        if (previousUpper > lowerCoordinate)
          throw JsonError(rowContext + " overlaps the preceding energy bin");
        previousUpper = upperCoordinate;
        component.binLower.push_back(lowerCoordinate);
        component.binUpper.push_back(upperCoordinate);
      }
      const double value = NonNegativeNumber(row.at("value"), rowContext + ".value");
      component.value.push_back(value);
      const bool valid = row.at("valid").AsBool(rowContext + ".valid");
      component.valid.push_back(valid);
      if (component.uncertaintySemantics == "conservative_bounds") {
        const double lower = NonNegativeNumber(row.at("lower"), rowContext + ".lower");
        const double upper = NonNegativeNumber(row.at("upper"), rowContext + ".upper");
        if (lower > value || value > upper)
          throw JsonError(rowContext + " requires lower <= value <= upper");
        component.lowerBound.push_back(lower);
        component.upperBound.push_back(upper);
      } else if (row.count("lower") || row.count("upper")) {
        throw JsonError(rowContext +
                        " lower/upper require uncertainty.semantics=conservative_bounds");
      }

      const double unresolved = NonNegativeNumber(
          row.at("unresolved_fraction"), rowContext + ".unresolved_fraction");
      const long sampled = NonNegativeInteger(row.at("n_sampled"), rowContext + ".n_sampled");
      const long resolved = NonNegativeInteger(row.at("n_resolved"), rowContext + ".n_resolved");
      if (sampled <= 0 || resolved > sampled)
        throw JsonError(rowContext + " requires 0 <= n_resolved <= n_sampled and n_sampled > 0");
      const double expectedUnresolved = static_cast<double>(sampled - resolved) /
                                        static_cast<double>(sampled);
      if (unresolved > 1.0 || std::abs(unresolved - expectedUnresolved) > 1.0e-12)
        throw JsonError(rowContext +
                        ".unresolved_fraction conflicts with n_sampled/n_resolved");
      component.unresolvedFraction.push_back(unresolved);
      component.sampledDirections.push_back(static_cast<double>(sampled));
      component.resolvedDirections.push_back(static_cast<double>(resolved));
      if (row.count("n_allowed")) {
        const long allowed = NonNegativeInteger(row.at("n_allowed"), rowContext + ".n_allowed");
        if (allowed > resolved)
          throw JsonError(rowContext + ".n_allowed must not exceed n_resolved");
        component.allowedDirections.push_back(static_cast<double>(allowed));
      } else if (!component.allowedDirections.empty()) {
        throw JsonError(componentContext + ".data must provide n_allowed for every row or none");
      }
      observedMaximumUnresolved = std::max(observedMaximumUnresolved, unresolved);
    }
    if (!component.allowedDirections.empty() &&
        component.allowedDirections.size() != data.size())
      throw JsonError(componentContext + ".data must provide n_allowed for every row or none");
    dataset.components.push_back(std::move(component));
  }

  if (std::abs(observedMaximumUnresolved - dataset.maximumUnresolvedFraction) > 1.0e-12)
    throw JsonError("spectrum exchange.quality.maximum_unresolved_fraction does not match data");

  dataset.energyVariable = "component_specific";
  dataset.energyUnit = "component_specific";
  dataset.valueType = "differential_directional_intensity";
  dataset.valueUnit = "component_specific";
  dataset.angularModel = "isotropic";
  dataset.contentHash = Sha256(dataset.canonicalMetadata.DumpCanonical());
  return dataset;
}

} // namespace shield_sim
