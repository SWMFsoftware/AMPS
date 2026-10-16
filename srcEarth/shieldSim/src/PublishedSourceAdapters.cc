#include "PublishedSourceAdapters.hh"

#include "ContentHash.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <tuple>

namespace shield_sim {
namespace {

namespace fs = std::filesystem;

std::string ReadBytes(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw JsonError("cannot open published-source file '" + path + "'");
  std::ostringstream bytes;
  bytes << stream.rdbuf();
  if (!stream.good() && !stream.eof())
    throw JsonError("failed while reading published-source file '" + path + "'");
  return bytes.str();
}

std::string RequiredString(const Json::Object& object, const std::string& field,
                           const std::string& context) {
  const auto found = object.find(field);
  if (found == object.end()) throw JsonError(context + " is missing required field '" + field + "'");
  const std::string value = found->second.AsString(context + '.' + field);
  if (value.empty()) throw JsonError(context + '.' + field + " must not be empty");
  return value;
}

double Finite(const Json& value, const std::string& context) {
  const double result = value.AsNumber(context);
  if (!std::isfinite(result)) throw JsonError(context + " must be finite");
  return result;
}

int ExactInt(const Json& value, const std::string& context) {
  const double number = Finite(value, context);
  if (std::floor(number) != number ||
      number < static_cast<double>(std::numeric_limits<int>::min()) ||
      number > static_cast<double>(std::numeric_limits<int>::max()))
    throw JsonError(context + " must be an integer in range");
  return static_cast<int>(number);
}

std::vector<double> NumberArray(const Json& value, const std::string& context,
                                std::size_t minimumSize = 1) {
  std::vector<double> result;
  for (const auto& item : value.AsArray(context)) result.push_back(Finite(item, context + "[]"));
  if (result.size() < minimumSize)
    throw JsonError(context + " must contain at least " + std::to_string(minimumSize) + " values");
  for (std::size_t i = 0; i < result.size(); ++i) {
    if (result[i] <= 0.0) throw JsonError(context + " values must be > 0");
    if (i && result[i] <= result[i - 1])
      throw JsonError(context + " must be strictly increasing");
  }
  return result;
}

std::vector<std::string> ParseCsvLine(const std::string& line,
                                      const std::string& context) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (quoted) {
      if (c == '"') {
        if (i + 1 < line.size() && line[i + 1] == '"') {
          field.push_back('"');
          ++i;
        } else {
          quoted = false;
        }
      } else {
        field.push_back(c);
      }
    } else if (c == ',') {
      fields.push_back(field);
      field.clear();
    } else if (c == '"') {
      if (!field.empty()) throw JsonError(context + ": quote begins inside an unquoted CSV field");
      quoted = true;
    } else if (c != '\r') {
      field.push_back(c);
    }
  }
  if (quoted) throw JsonError(context + ": unterminated quoted CSV field");
  fields.push_back(field);
  return fields;
}

double CsvDouble(const std::string& text, const std::string& context) {
  std::size_t consumed = 0;
  double result = 0.0;
  try { result = std::stod(text, &consumed); }
  catch (const std::exception&) { throw JsonError(context + " must be a number"); }
  if (consumed != text.size() || !std::isfinite(result))
    throw JsonError(context + " must be a finite number");
  return result;
}

int CsvInt(const std::string& text, const std::string& context) {
  std::size_t consumed = 0;
  long long result = 0;
  try { result = std::stoll(text, &consumed); }
  catch (const std::exception&) { throw JsonError(context + " must be an integer"); }
  if (consumed != text.size() || result < std::numeric_limits<int>::min() ||
      result > std::numeric_limits<int>::max())
    throw JsonError(context + " must be an integer in range");
  return static_cast<int>(result);
}

// Days since 1970-01-01. This is the public-domain civil-calendar algorithm
// by Howard Hinnant, written out here so timestamp validation does not depend
// on the process timezone or a non-standard strptime implementation.
std::int64_t DaysFromCivil(int year, unsigned month, unsigned day) {
  year -= month <= 2;
  const int era = (year >= 0 ? year : year - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(year - era * 400);
  const unsigned doy = (153 * (month + (month > 2 ? -3 : 9)) + 2) / 5 + day - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return static_cast<std::int64_t>(era) * 146097 + static_cast<int>(doe) - 719468;
}

bool LeapYear(int year) {
  return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

std::int64_t ParseUtcSeconds(const std::string& value, const std::string& context) {
  if (value.size() != 20 || value[4] != '-' || value[7] != '-' || value[10] != 'T' ||
      value[13] != ':' || value[16] != ':' || value[19] != 'Z')
    throw JsonError(context + " must use canonical YYYY-MM-DDTHH:MM:SSZ UTC format");
  auto digits = [&](std::size_t offset, std::size_t count) {
    int number = 0;
    for (std::size_t i = 0; i < count; ++i) {
      const char c = value[offset + i];
      if (c < '0' || c > '9') throw JsonError(context + " contains a non-digit date field");
      number = number * 10 + c - '0';
    }
    return number;
  };
  const int year = digits(0, 4);
  const int month = digits(5, 2);
  const int day = digits(8, 2);
  const int hour = digits(11, 2);
  const int minute = digits(14, 2);
  const int second = digits(17, 2);
  static const int monthDays[] = {31,28,31,30,31,30,31,31,30,31,30,31};
  if (month < 1 || month > 12 || hour > 23 || minute > 59 || second > 59)
    throw JsonError(context + " is outside the supported civil-time range");
  const int maximumDay = monthDays[month - 1] + (month == 2 && LeapYear(year) ? 1 : 0);
  if (day < 1 || day > maximumDay) throw JsonError(context + " has an invalid calendar day");
  return DaysFromCivil(year, static_cast<unsigned>(month), static_cast<unsigned>(day)) * 86400 +
         hour * 3600 + minute * 60 + second;
}

Json DoubleArrayJson(const std::vector<double>& values) {
  Json::Array result;
  for (double value : values) result.emplace_back(value);
  return Json(std::move(result));
}

Json ComponentJson(const SpectrumComponent& component) {
  Json::Object object;
  object["angular_representation"] = component.angularRepresentation;
  object["atomic_number"] = component.atomicNumber;
  object["bin_lower"] = DoubleArrayJson(component.binLower);
  object["bin_upper"] = DoubleArrayJson(component.binUpper);
  object["component_id"] = component.componentId;
  object["direction_id"] = component.directionId;
  object["energy_unit"] = component.energyUnit;
  object["energy_variable"] = component.energyVariable;
  object["geant4_particle"] = component.geant4Particle;
  object["mass_number"] = component.massNumber;
  object["species_id"] = component.speciesId;
  object["total_uncertainty"] = DoubleArrayJson(component.totalUncertainty);
  object["uncertainty_semantics"] = component.uncertaintySemantics;
  object["value"] = DoubleArrayJson(component.value);
  object["value_type"] = component.valueType;
  object["value_unit"] = component.valueUnit;
  return Json(std::move(object));
}

const Json::Object& FindParameter(const Json::Array& parameters,
                                  const std::string& name) {
  const Json::Object* found = nullptr;
  for (const auto& parameterValue : parameters) {
    const auto& parameter = parameterValue.AsObject("AMS HAPI parameter");
    const auto it = parameter.find("name");
    if (it != parameter.end() && it->second.AsString("AMS HAPI parameter.name") == name) {
      if (found) throw JsonError("AMS HAPI metadata repeats parameter '" + name + "'");
      found = &parameter;
    }
  }
  if (!found) throw JsonError("AMS HAPI metadata lacks parameter '" + name + "'");
  return *found;
}

std::vector<double> HapiCoordinates(const Json::Object& parameter,
                                    const std::string& expectedBinName,
                                    std::size_t expectedSize,
                                    const std::string& context) {
  const auto& size = parameter.at("size").AsArray(context + ".size");
  if (size.size() != 1 || ExactInt(size[0], context + ".size[0]") != static_cast<int>(expectedSize))
    throw JsonError(context + ".size does not match the primary-table bin count");
  const auto& bins = parameter.at("bins").AsArray(context + ".bins");
  if (bins.size() != 1) throw JsonError(context + ".bins must contain exactly one coordinate");
  const auto& coordinate = bins[0].AsObject(context + ".bins[0]");
  if (RequiredString(coordinate, "name", context + ".bins[0]") != expectedBinName)
    throw JsonError(context + " has an unexpected bin-coordinate name");
  if (RequiredString(coordinate, "units", context + ".bins[0]") != "GeV/n")
    throw JsonError(context + " coordinate unit must be GeV/n");
  const auto result = NumberArray(coordinate.at("centers"), context + ".bins[0].centers");
  if (result.size() != expectedSize) throw JsonError(context + " coordinate count changed");
  return result;
}

double KineticEnergyPerNucleonGeV(double rigidityGV, double chargeOverMass,
                                  double nucleonMassGeV) {
  const double momentumPerNucleon = rigidityGV * chargeOverMass;
  return std::sqrt(momentumPerNucleon * momentumPerNucleon +
                   nucleonMassGeV * nucleonMassGeV) - nucleonMassGeV;
}

std::vector<double> ParseCsvSingleRow(const std::vector<std::string>& fields,
                                      std::size_t offset, std::size_t count,
                                      const std::string& context) {
  std::vector<double> result;
  for (std::size_t i = 0; i < count; ++i) {
    const double value = CsvDouble(fields.at(offset + i), context + '[' + std::to_string(i) + ']');
    if (value < 0.0) throw JsonError(context + " contains a fill/negative value");
    result.push_back(value);
  }
  return result;
}

} // namespace

SpectrumDataset AMS02HapiAdapter::Load(const std::string& infoPath,
                                       const std::string& dataPath,
                                       const std::string& interpretationPath) {
  const std::string infoBytes = ReadBytes(infoPath);
  const std::string dataBytes = ReadBytes(dataPath);
  const Json interpretation = Json::ParseFile(interpretationPath);
  const auto& contract = interpretation.AsObject("AMS adapter interpretation");
  RejectUnknownFields(contract,
      {"schema_version","dataset_id","adapter_version","info_sha256","data_sha256",
       "supplement","record_start_utc","record_stop_utc","bartels_rotation",
       "nucleon_rest_mass_GeV",
       "proton_rigidity_edges_GV","helium_rigidity_edges_GV","science_status"},
      "AMS adapter interpretation");
  if (RequiredString(contract, "schema_version", "AMS adapter interpretation") !=
      "shieldSim.ams02_hapi_interpretation.v1")
    throw JsonError("unsupported AMS adapter interpretation schema");
  if (RequiredString(contract, "adapter_version", "AMS adapter interpretation") !=
      "AMS02HapiAdapter.v1")
    throw JsonError("unsupported AMS adapter version");
  if (RequiredString(contract, "info_sha256", "AMS adapter interpretation") != Sha256(infoBytes))
    throw JsonError("AMS HAPI info SHA-256 does not match its interpretation record");
  if (RequiredString(contract, "data_sha256", "AMS adapter interpretation") != Sha256(dataBytes))
    throw JsonError("AMS HAPI data SHA-256 does not match its interpretation record");

  const fs::path interpretationDirectory = fs::absolute(interpretationPath).parent_path();
  const auto& supplement = contract.at("supplement").AsObject("AMS adapter interpretation.supplement");
  RejectUnknownFields(supplement, {"path","sha256","url","table"},
                      "AMS adapter interpretation.supplement");
  const fs::path supplementPath =
      (interpretationDirectory / RequiredString(supplement, "path",
                                                "AMS adapter interpretation.supplement")).lexically_normal();
  if (Sha256(ReadBytes(supplementPath.string())) !=
      RequiredString(supplement, "sha256", "AMS adapter interpretation.supplement"))
    throw JsonError("AMS supplemental-table SHA-256 does not match its interpretation record");

  const std::vector<double> protonRigidity = NumberArray(
      contract.at("proton_rigidity_edges_GV"), "AMS proton rigidity edges", 2);
  const std::vector<double> heliumRigidity = NumberArray(
      contract.at("helium_rigidity_edges_GV"), "AMS helium rigidity edges", 2);
  if (protonRigidity.size() != 46 || heliumRigidity.size() != 41)
    throw JsonError("AMS primary-table grids must contain 46 proton and 41 helium edges");
  const double nucleonMassGeV = Finite(contract.at("nucleon_rest_mass_GeV"),
                                       "AMS nucleon_rest_mass_GeV");
  if (nucleonMassGeV <= 0.0) throw JsonError("AMS nucleon_rest_mass_GeV must be > 0");

  const Json info = Json::Parse(infoBytes, infoPath);
  const auto& infoObject = info.AsObject("AMS HAPI info");
  if (RequiredString(infoObject, "HAPI", "AMS HAPI info") != "2.0")
    throw JsonError("AMS fixture must use HAPI 2.0 metadata");
  const auto& status = infoObject.at("status").AsObject("AMS HAPI status");
  if (ExactInt(status.at("code"), "AMS HAPI status.code") != 1200)
    throw JsonError("AMS HAPI metadata status is not OK");
  const auto& parameters = infoObject.at("parameters").AsArray("AMS HAPI parameters");
  const auto& proton = FindParameter(parameters, "ProFlux");
  const auto& protonError = FindParameter(parameters, "ProFlux_TotErr");
  const auto& helium = FindParameter(parameters, "HeFlux");
  const auto& heliumError = FindParameter(parameters, "HeFlux_TotErr");
  const std::string hapiUnit = "1/[m**2-sr-s-GeV/n]";
  for (const auto* parameter : {&proton, &protonError, &helium, &heliumError}) {
    if (RequiredString(*parameter, "units", "AMS HAPI parameter") != hapiUnit)
      throw JsonError("AMS HAPI flux/error units changed");
  }
  const std::vector<double> hapiProtonLow =
      HapiCoordinates(proton, "BinProLow", 45, "AMS ProFlux");
  const std::vector<double> hapiHeliumLow =
      HapiCoordinates(helium, "BinHeLow", 40, "AMS HeFlux");
  const std::vector<double> hapiProtonErrorLow =
      HapiCoordinates(protonError, "BinProLow", 45, "AMS ProFlux_TotErr");
  const std::vector<double> hapiHeliumErrorLow =
      HapiCoordinates(heliumError, "BinHeLow", 40, "AMS HeFlux_TotErr");
  if (hapiProtonErrorLow != hapiProtonLow || hapiHeliumErrorLow != hapiHeliumLow)
    throw JsonError("AMS flux and total-uncertainty coordinate grids differ");

  auto transformedEdges = [&](const std::vector<double>& rigidity, double zOverA) {
    std::vector<double> result;
    for (double edge : rigidity)
      result.push_back(1000.0 * KineticEnergyPerNucleonGeV(edge, zOverA, nucleonMassGeV));
    return result;
  };
  const std::vector<double> protonEdges = transformedEdges(protonRigidity, 1.0);
  const std::vector<double> heliumEdges = transformedEdges(heliumRigidity, 0.5);
  // The HAPI coordinate is rounded to two or three significant figures.  It is
  // named Bin*Low and agrees with the lower edges transformed from the primary
  // rigidity tables.  This explicit check resolves, rather than ignores, the
  // contradictory generic HAPI JSON key called "centers".
  for (std::size_t i = 0; i < hapiProtonLow.size(); ++i)
    if (std::abs(hapiProtonLow[i] * 1000.0 - protonEdges[i]) > 60.0)
      throw JsonError("AMS proton HAPI coordinate is inconsistent with transformed rigidity lower edges");
  for (std::size_t i = 0; i < hapiHeliumLow.size(); ++i)
    if (std::abs(hapiHeliumLow[i] * 1000.0 - heliumEdges[i]) > 60.0)
      throw JsonError("AMS helium HAPI coordinate is inconsistent with transformed rigidity lower edges");

  std::vector<std::string> dataRows;
  std::istringstream dataStream(dataBytes);
  std::string line;
  while (std::getline(dataStream, line)) {
    if (!line.empty() && line[0] != '#') dataRows.push_back(line);
  }
  if (dataRows.size() != 1) throw JsonError("AMS adapter fixture must contain exactly one HAPI record");
  const auto fields = ParseCsvLine(dataRows[0], "AMS HAPI CSV row");
  if (fields.size() != 172) throw JsonError("AMS HAPI selected-field row must contain 172 columns");
  const std::string expectedStart = RequiredString(contract, "record_start_utc", "AMS adapter interpretation");
  const std::string expectedStop = RequiredString(contract, "record_stop_utc", "AMS adapter interpretation");
  const std::int64_t startSeconds = ParseUtcSeconds(expectedStart, "AMS record_start_utc");
  const std::int64_t stopSeconds = ParseUtcSeconds(expectedStop, "AMS record_stop_utc");
  if (startSeconds >= stopSeconds) throw JsonError("AMS record interval must have start < stop");
  if (fields[0] != expectedStart.substr(0, 19) + ".000000000Z")
    throw JsonError("AMS HAPI record timestamp changed");
  const int expectedBartelsRotation = ExactInt(contract.at("bartels_rotation"),
                                                "AMS bartels_rotation");
  // This HAPI response leaves the requested BR scalar empty. The pinned AMS
  // supplemental table supplies BR 2460, so retain that provenance in the
  // interpretation record. If CDAWeb later populates the scalar it must agree;
  // an empty field is never invented from the numeric payload itself.
  if (!fields[1].empty() &&
      CsvInt(fields[1], "AMS HAPI Bartels rotation") != expectedBartelsRotation)
    throw JsonError("AMS HAPI Bartels rotation changed");

  const std::vector<double> protonFlux = ParseCsvSingleRow(fields, 2, 45, "AMS proton flux");
  const std::vector<double> protonUncertainty = ParseCsvSingleRow(fields, 47, 45, "AMS proton uncertainty");
  const std::vector<double> heliumFlux = ParseCsvSingleRow(fields, 92, 40, "AMS helium flux");
  const std::vector<double> heliumUncertainty = ParseCsvSingleRow(fields, 132, 40, "AMS helium uncertainty");
  constexpr double kPerM2PerGeVToPerCm2PerMeV = 1.0e-7;

  auto makeComponent = [&](const std::string& id, const std::string& species,
                           const std::string& particle, int z, int a,
                           const std::vector<double>& edges,
                           const std::vector<double>& flux,
                           const std::vector<double>& uncertainty) {
    SpectrumComponent component;
    component.componentId = id;
    component.speciesId = species;
    component.geant4Particle = particle;
    component.atomicNumber = z;
    component.massNumber = a;
    component.energyVariable = "kinetic_energy_per_nucleon";
    component.energyUnit = "MeV/n";
    component.valueType = "differential_intensity";
    component.valueUnit = "1/(cm2 sr s (MeV/n))";
    component.angularRepresentation = "isotropic_intensity_at_instrument";
    component.directionId = "AMS02_acceptance";
    component.uncertaintySemantics = "published_total_statistical_and_systematic_quadrature";
    component.binLower.assign(edges.begin(), edges.end() - 1);
    component.binUpper.assign(edges.begin() + 1, edges.end());
    for (double value : flux) component.value.push_back(value * kPerM2PerGeVToPerCm2PerMeV);
    for (double value : uncertainty)
      component.totalUncertainty.push_back(value * kPerM2PerGeVToPerCm2PerMeV);
    return component;
  };

  SpectrumDataset result;
  result.id = RequiredString(contract, "dataset_id", "AMS adapter interpretation");
  result.energyVariable = "kinetic_energy_per_nucleon";
  result.energyUnit = "MeV/n";
  result.valueType = "differential_intensity";
  result.valueUnit = "1/(cm2 sr s (MeV/n))";
  result.angularModel = "isotropic_intensity_at_instrument";
  result.scienceStatus = RequiredString(contract, "science_status", "AMS adapter interpretation");
  result.spectrumPath = dataPath;
  result.metadataPath = interpretationPath;
  result.components.push_back(makeComponent("proton", "proton", "proton", 1, 1,
                                            protonEdges, protonFlux, protonUncertainty));
  // AMS publishes inclusive Z=2 helium.  A=0 and an empty Geant4 particle are
  // intentional: treating this population as pure alpha particles would be an
  // unsupported isotope assumption.
  result.components.push_back(makeComponent("helium", "helium_mixed_isotopes", "", 2, 0,
                                            heliumEdges, heliumFlux, heliumUncertainty));

  Json::Array componentJson;
  for (const auto& component : result.components) componentJson.push_back(ComponentJson(component));
  Json::Object canonical;
  canonical["adapter"] = "AMS02HapiAdapter.v1";
  canonical["components"] = Json(std::move(componentJson));
  canonical["dataset_id"] = result.id;
  canonical["info_sha256"] = Sha256(infoBytes);
  canonical["observation_start_utc"] = expectedStart;
  canonical["observation_stop_utc"] = expectedStop;
  canonical["science_status"] = result.scienceStatus;
  canonical["source_data_sha256"] = Sha256(dataBytes);
  canonical["source_interpretation"] = interpretation;
  canonical["schema_version"] = "shieldSim.spectrum_dataset.v2";
  result.canonicalMetadata = Json(std::move(canonical));
  result.contentHash = Sha256(result.canonicalMetadata.DumpCanonical());
  return result;
}

bool GoesSgpsRecord::IsUsable() const {
  return l1bRecordsInAverage == 300 && validOneSecondSamples == 300 &&
         deadTimeFlaggedSamples == 0 && outOfBandFlaggedSamples == 0 &&
         dynamicErrorFlaggedSamples == 0 && hasFlux && hasUncertainty;
}

GOESRSGPSAdapter GOESRSGPSAdapter::Load(const std::string& csvPath,
                                        const std::string& eventDefinitionPath) {
  GOESRSGPSAdapter result;
  const std::string csvBytes = ReadBytes(csvPath);
  result.canonicalMetadata_ = Json::ParseFile(eventDefinitionPath);
  const auto& definition = result.canonicalMetadata_.AsObject("GOES event definition");
  RejectUnknownFields(definition,
      {"schema_version","dataset_id","adapter_version","source_file","timestamp_semantics",
       "quality_policy","windows","markers","science_status"}, "GOES event definition");
  if (RequiredString(definition, "schema_version", "GOES event definition") !=
      "shieldSim.goes_sgps_event_definition.v1")
    throw JsonError("unsupported GOES event-definition schema");
  if (RequiredString(definition, "adapter_version", "GOES event definition") !=
      "GOESRSGPSAdapter.v1")
    throw JsonError("unsupported GOES adapter version");
  RequiredString(definition, "dataset_id", "GOES event definition");
  if (RequiredString(definition, "science_status", "GOES event definition") !=
      "published_directional_geosynchronous_event_ingestion_not_transport_boundary")
    throw JsonError("GOES science status must retain the transport-boundary prohibition");

  const auto& source = definition.at("source_file").AsObject("GOES source_file");
  RejectUnknownFields(source, {"url","sha256","extracted_csv_sha256","license","product_id"},
                      "GOES source_file");
  if (RequiredString(source, "extracted_csv_sha256", "GOES source_file") != Sha256(csvBytes))
    throw JsonError("GOES extracted CSV SHA-256 does not match event definition");
  const auto& timing = definition.at("timestamp_semantics").AsObject("GOES timestamp_semantics");
  RejectUnknownFields(timing, {"source_coordinate","canonical_coordinate","cadence_seconds",
                               "sample_support","boundary_policy"}, "GOES timestamp_semantics");
  if (RequiredString(timing, "source_coordinate", "GOES timestamp_semantics") != "interval_start" ||
      RequiredString(timing, "canonical_coordinate", "GOES timestamp_semantics") != "midpoint" ||
      ExactInt(timing.at("cadence_seconds"), "GOES timestamp_semantics.cadence_seconds") != 300 ||
      RequiredString(timing, "sample_support", "GOES timestamp_semantics") !=
          "[midpoint-150s,midpoint+150s)" ||
      RequiredString(timing, "boundary_policy", "GOES timestamp_semantics") !=
          "duration_weighted_overlap")
    throw JsonError("GOES timestamp semantics differ from the approved Stage-2 contract");
  const auto& quality = definition.at("quality_policy").AsObject("GOES quality_policy");
  RejectUnknownFields(quality, {"required_l1b_records","required_valid_samples",
                                "require_zero_dead_time_flags","require_zero_out_of_band_flags",
                                "require_zero_dynamic_error_flags","missing_policy",
                                "uncertainty_policy"}, "GOES quality_policy");
  if (ExactInt(quality.at("required_l1b_records"), "GOES quality_policy.required_l1b_records") != 300 ||
      ExactInt(quality.at("required_valid_samples"), "GOES quality_policy.required_valid_samples") != 300 ||
      !quality.at("require_zero_dead_time_flags").AsBool("GOES quality_policy") ||
      !quality.at("require_zero_out_of_band_flags").AsBool("GOES quality_policy") ||
      !quality.at("require_zero_dynamic_error_flags").AsBool("GOES quality_policy") ||
      RequiredString(quality, "missing_policy", "GOES quality_policy") != "no_interpolation" ||
      RequiredString(quality, "uncertainty_policy", "GOES quality_policy") !=
          "retain_total_do_not_combine_without_covariance")
    throw JsonError("GOES quality policy differs from the approved Stage-2 contract");

  std::set<std::string> windowIds;
  std::set<std::string> windowRoles;
  for (const auto& windowValue : definition.at("windows").AsArray("GOES windows")) {
    const auto& windowObject = windowValue.AsObject("GOES window");
    RejectUnknownFields(windowObject, {"id","role","start_utc","end_utc","definition","source_url"},
                        "GOES window");
    GoesEventWindow window;
    window.id = RequiredString(windowObject, "id", "GOES window");
    if (!windowIds.insert(window.id).second) throw JsonError("duplicate GOES event-window id");
    window.role = RequiredString(windowObject, "role", "GOES window");
    if (window.role != "primary" && window.role != "sensitivity")
      throw JsonError("GOES window.role must be primary or sensitivity");
    if (!windowRoles.insert(window.role).second)
      throw JsonError("GOES Stage-2 definition must contain one primary and one sensitivity window");
    window.startUtc = RequiredString(windowObject, "start_utc", "GOES window");
    window.endUtc = RequiredString(windowObject, "end_utc", "GOES window");
    window.startUnixSeconds = ParseUtcSeconds(window.startUtc, "GOES window.start_utc");
    window.endUnixSeconds = ParseUtcSeconds(window.endUtc, "GOES window.end_utc");
    if (window.startUnixSeconds >= window.endUnixSeconds)
      throw JsonError("GOES event window must have start < end");
    window.definition = RequiredString(windowObject, "definition", "GOES window");
    window.sourceUrl = RequiredString(windowObject, "source_url", "GOES window");
    result.windows_.push_back(std::move(window));
  }
  if (result.windows_.size() != 2)
    throw JsonError("GOES Stage-2 definition requires primary and sensitivity windows");

  // Markers do not control the integration bounds, but they are part of the
  // auditable event definition. Validate them instead of treating arbitrary
  // nested JSON as trustworthy metadata.
  const auto& markers = definition.at("markers").AsArray("GOES markers");
  if (markers.empty()) throw JsonError("GOES event definition requires an independent onset marker");
  std::set<std::string> markerIds;
  for (const auto& markerValue : markers) {
    const auto& marker = markerValue.AsObject("GOES marker");
    RejectUnknownFields(marker, {"id","time_utc","role","source_url"}, "GOES marker");
    const std::string markerId = RequiredString(marker, "id", "GOES marker");
    if (!markerIds.insert(markerId).second) throw JsonError("duplicate GOES marker id");
    ParseUtcSeconds(RequiredString(marker, "time_utc", "GOES marker"), "GOES marker.time_utc");
    if (RequiredString(marker, "role", "GOES marker") !=
        "independent_onset_marker_not_an_integration_boundary")
      throw JsonError("GOES marker must not be used as an integration boundary");
    RequiredString(marker, "source_url", "GOES marker");
  }

  std::istringstream stream(csvBytes);
  std::string line;
  std::vector<std::string> header;
  std::size_t lineNumber = 0;
  std::set<std::tuple<std::int64_t,int,std::string>> uniqueRows;
  const std::vector<std::string> expectedHeader = {
      "sample_start_utc","sample_midpoint_utc","sample_end_utc","l1b_records_in_average",
      "sensor_unit","spacecraft_axis","upright_nominal_look","yaw_flip_flag","channel",
      "energy_low_MeV","energy_high_MeV","proton_flux_per_cm2_sr_s_MeV",
      "flux_uncertainty_per_cm2_sr_s_MeV","valid_1s_samples","dead_time_flagged_samples",
      "out_of_band_flagged_samples","dynamic_error_flagged_samples"};
  while (std::getline(stream, line)) {
    ++lineNumber;
    if (line.empty() || line[0] == '#') continue;
    const auto fields = ParseCsvLine(line, csvPath + ':' + std::to_string(lineNumber));
    if (header.empty()) {
      header = fields;
      if (header != expectedHeader) throw JsonError("GOES extracted CSV header/schema changed");
      continue;
    }
    if (fields.size() != expectedHeader.size())
      throw JsonError(csvPath + ':' + std::to_string(lineNumber) + ": wrong field count");
    auto context = [&](std::size_t index) {
      return csvPath + ':' + std::to_string(lineNumber) + ':' + expectedHeader[index];
    };
    GoesSgpsRecord record;
    record.startUtc = fields[0];
    record.midpointUtc = fields[1];
    record.endUtc = fields[2];
    record.startUnixSeconds = ParseUtcSeconds(record.startUtc, context(0));
    record.midpointUnixSeconds = ParseUtcSeconds(record.midpointUtc, context(1));
    record.endUnixSeconds = ParseUtcSeconds(record.endUtc, context(2));
    if (record.midpointUnixSeconds != record.startUnixSeconds + 150 ||
        record.endUnixSeconds != record.startUnixSeconds + 300)
      throw JsonError(context(1) + ": midpoint/end do not describe the source five-minute interval");
    record.l1bRecordsInAverage = CsvInt(fields[3], context(3));
    record.sensorUnit = CsvInt(fields[4], context(4));
    if (record.sensorUnit != 0 && record.sensorUnit != 1) throw JsonError(context(4) + " must be 0 or 1");
    record.spacecraftAxis = fields[5];
    record.nominalLook = fields[6];
    if ((record.sensorUnit == 0 && (record.spacecraftAxis != "-X" || record.nominalLook != "westward")) ||
        (record.sensorUnit == 1 && (record.spacecraftAxis != "+X" || record.nominalLook != "eastward")))
      throw JsonError(context(5) + ": sensor direction mapping changed");
    record.yawFlipFlag = CsvInt(fields[7], context(7));
    record.channel = fields[8];
    if (record.channel.empty()) throw JsonError(context(8) + " must not be empty");
    record.energyLowMeV = CsvDouble(fields[9], context(9));
    record.energyHighMeV = CsvDouble(fields[10], context(10));
    if (!(record.energyLowMeV > 0.0 && record.energyHighMeV > record.energyLowMeV))
      throw JsonError(context(9) + ": invalid energy bounds");
    if (fields[11] != "NA") {
      record.fluxPerCm2SrSMeV = CsvDouble(fields[11], context(11));
      if (record.fluxPerCm2SrSMeV < 0.0) throw JsonError(context(11) + " must be >= 0 or NA");
      record.hasFlux = true;
    }
    if (fields[12] != "NA") {
      record.uncertaintyPerCm2SrSMeV = CsvDouble(fields[12], context(12));
      if (record.uncertaintyPerCm2SrSMeV < 0.0) throw JsonError(context(12) + " must be >= 0 or NA");
      record.hasUncertainty = true;
    }
    record.validOneSecondSamples = CsvInt(fields[13], context(13));
    record.deadTimeFlaggedSamples = CsvInt(fields[14], context(14));
    record.outOfBandFlaggedSamples = CsvInt(fields[15], context(15));
    record.dynamicErrorFlaggedSamples = CsvInt(fields[16], context(16));
    for (int count : {record.l1bRecordsInAverage, record.yawFlipFlag,
                      record.validOneSecondSamples,
                      record.deadTimeFlaggedSamples, record.outOfBandFlaggedSamples,
                      record.dynamicErrorFlaggedSamples})
      if (count < 0 || count > 301) throw JsonError("GOES sample/flag count is outside [0,301]");
    if (record.validOneSecondSamples > record.l1bRecordsInAverage)
      throw JsonError(context(13) + ": valid samples exceed records in average");
    for (int flagged : {record.deadTimeFlaggedSamples, record.outOfBandFlaggedSamples,
                        record.dynamicErrorFlaggedSamples})
      if (flagged > record.validOneSecondSamples)
        throw JsonError(context(13) + ": flagged samples exceed valid samples");
    const auto key = std::make_tuple(record.startUnixSeconds, record.sensorUnit, record.channel);
    if (!uniqueRows.insert(key).second) throw JsonError("duplicate GOES time/sensor/channel record");
    result.records_.push_back(std::move(record));
  }
  if (header.empty() || result.records_.empty()) throw JsonError("GOES extracted CSV contains no records");
  std::sort(result.records_.begin(), result.records_.end(), [](const auto& a, const auto& b) {
    return std::tie(a.startUnixSeconds,a.sensorUnit,a.channel) <
           std::tie(b.startUnixSeconds,b.sensorUnit,b.channel);
  });
  result.contentHash_ = Sha256(result.canonicalMetadata_.DumpCanonical() + "\n" + Sha256(csvBytes));
  return result;
}

std::vector<GoesFluenceChannel> GOESRSGPSAdapter::IntegrateRecords(
    const std::vector<GoesSgpsRecord>& records,
    const std::vector<GoesEventWindow>& windows) {
  using Key = std::tuple<int,std::string>;
  std::map<Key, const GoesSgpsRecord*> definitions;
  std::map<Key, std::set<std::int64_t>> starts;
  for (const auto& record : records) {
    const Key key(record.sensorUnit, record.channel);
    const auto inserted = definitions.emplace(key, &record);
    if (!inserted.second) {
      const auto& reference = *inserted.first->second;
      if (record.energyLowMeV != reference.energyLowMeV ||
          record.energyHighMeV != reference.energyHighMeV ||
          record.spacecraftAxis != reference.spacecraftAxis ||
          record.nominalLook != reference.nominalLook)
        throw JsonError("GOES channel definition changes within the event time series");
    }
    if (!starts[key].insert(record.startUnixSeconds).second)
      throw JsonError("GOES integration encountered a duplicate time bin");
  }
  if (definitions.empty()) throw JsonError("event integration requires at least one channel");
  std::vector<GoesFluenceChannel> outputs;
  for (const auto& window : windows) {
    const double expected = static_cast<double>(window.endUnixSeconds - window.startUnixSeconds);
    for (const auto& definition : definitions) {
      const int sensor = std::get<0>(definition.first);
      const std::string& channel = std::get<1>(definition.first);
      const auto& reference = *definition.second;
      GoesFluenceChannel output;
      output.windowId = window.id;
      output.sensorId = "sgps_" + std::to_string(sensor);
      output.spacecraftAxis = reference.spacecraftAxis;
      output.nominalLook = reference.nominalLook;
      output.channel = channel;
      output.energyLowMeV = reference.energyLowMeV;
      output.energyHighMeV = reference.energyHighMeV;
      output.expectedSeconds = expected;
      output.uncertaintyStatus = "not_assessable_temporal_covariance_unpublished";
      for (const auto& record : records) {
        if (record.sensorUnit != sensor || record.channel != channel) continue;
        const std::int64_t overlapStart = std::max(record.startUnixSeconds, window.startUnixSeconds);
        const std::int64_t overlapEnd = std::min(record.endUnixSeconds, window.endUnixSeconds);
        if (overlapStart >= overlapEnd) continue;
        const double seconds = static_cast<double>(overlapEnd - overlapStart);
        if (record.IsUsable()) {
          output.validSeconds += seconds;
          output.partialFluencePerCm2SrMeV += record.fluxPerCm2SrSMeV * seconds;
        } else {
          output.rejectedSeconds += seconds;
        }
      }
      output.missingSeconds = output.expectedSeconds - output.validSeconds - output.rejectedSeconds;
      if (output.missingSeconds < 0.0 && output.missingSeconds > -1e-9) output.missingSeconds = 0.0;
      if (output.missingSeconds < 0.0)
        throw JsonError("GOES time intervals overlap within a sensor/channel");
      output.coverageFraction = output.validSeconds / output.expectedSeconds;
      output.complete = output.validSeconds == output.expectedSeconds;
      outputs.push_back(output);
    }
  }
  return outputs;
}

GoesFluenceResult GOESRSGPSAdapter::Integrate() const {
  GoesFluenceResult result;
  result.windows = windows_;
  result.channels = IntegrateRecords(records_, windows_);
  if (result.channels.size() != windows_.size() * 26)
    throw JsonError("GOES event series must retain 2 sensors x 13 channels");

  Json::Array windowResults;
  for (const auto& window : windows_) {
    Json::Array channelResults;
    for (const auto& output : result.channels) {
      if (output.windowId != window.id) continue;
      Json::Object json;
      json["channel"] = output.channel;
      json["complete"] = output.complete;
      json["coverage_fraction"] = output.coverageFraction;
      json["energy_high_MeV"] = output.energyHighMeV;
      json["energy_low_MeV"] = output.energyLowMeV;
      json["expected_seconds"] = output.expectedSeconds;
      json["missing_seconds"] = output.missingSeconds;
      json["nominal_look"] = output.nominalLook;
      json["partial_fluence_per_cm2_sr_MeV"] = output.partialFluencePerCm2SrMeV;
      json["rejected_seconds"] = output.rejectedSeconds;
      json["sensor_id"] = output.sensorId;
      json["spacecraft_axis"] = output.spacecraftAxis;
      json["uncertainty_status"] = output.uncertaintyStatus;
      json["valid_seconds"] = output.validSeconds;
      channelResults.emplace_back(std::move(json));
    }
    Json::Object windowJson;
    windowJson["channels"] = Json(std::move(channelResults));
    windowJson["definition"] = window.definition;
    windowJson["end_utc_exclusive"] = window.endUtc;
    windowJson["id"] = window.id;
    windowJson["role"] = window.role;
    windowJson["source_url"] = window.sourceUrl;
    windowJson["start_utc_inclusive"] = window.startUtc;
    windowResults.emplace_back(std::move(windowJson));
  }
  Json::Object canonical;
  canonical["adapter"] = "GOESRSGPSAdapter.v1";
  canonical["event_definition"] = canonicalMetadata_;
  canonical["schema_version"] = "shieldSim.goes_event_fluence.v1";
  canonical["source_content_hash"] = contentHash_;
  canonical["windows"] = Json(std::move(windowResults));
  result.canonicalJson = Json(std::move(canonical));
  result.contentHash = Sha256(result.canonicalJson.DumpCanonical());
  return result;
}

} // namespace shield_sim
