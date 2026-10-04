#include "sep_coronal_cme/configuration_parser.h"

#include "sha256.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

namespace SEP {
namespace CoronalCME {
namespace {

struct GrammarEntry {
  const char* sectionPattern;
  const char* key;
  const char* expression;
};

constexpr GrammarEntry kSchema5Grammar[] = {
#include "schema5_registry.inc"
};

struct Assignment {
  std::string section;
  std::string key;
  std::string value;
  std::size_t line = 0;
};

std::string Trim(const std::string& input) {
  const std::size_t first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  const std::size_t last = input.find_last_not_of(" \t\r\n");
  return input.substr(first, last - first + 1);
}

std::vector<std::string> Split(const std::string& input,
                               const std::string& delimiter) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (true) {
    const std::size_t position = input.find(delimiter, start);
    fields.push_back(Trim(input.substr(
        start, position == std::string::npos ? position : position - start)));
    if (position == std::string::npos) break;
    start = position + delimiter.size();
  }
  return fields;
}

bool ParseFiniteDouble(const std::string& input, double* value) {
  if (value == nullptr || input.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(input.c_str(), &end);
  if (errno == ERANGE || end == input.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ParseUnsigned64(const std::string& input, std::uint64_t* value) {
  if (value == nullptr || input.empty() || input.front() == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(input.c_str(), &end, 10);
  if (errno == ERANGE || end == input.c_str() || *end != '\0') return false;
  *value = static_cast<std::uint64_t>(parsed);
  return true;
}

bool SectionMatches(const std::string& pattern, const std::string& section) {
  const std::string marker = ".ID";
  const std::size_t wildcard = pattern.find(marker);
  if (wildcard == std::string::npos) return pattern == section;
  const std::string prefix = pattern.substr(0, wildcard + 1u);
  if (section.size() <= prefix.size() || section.compare(0, prefix.size(), prefix) != 0)
    return false;
  return section.find('.', prefix.size()) == std::string::npos;
}

const GrammarEntry* FindGrammar(const std::string& section,
                                const std::string& key) {
  for (const GrammarEntry& entry : kSchema5Grammar) {
    if (entry.key == key && SectionMatches(entry.sectionPattern, section)) {
      return &entry;
    }
  }
  return nullptr;
}

bool IsNumericKey(const std::string& key, const std::string& expression) {
  if (expression == "REQUIRED_OR_ZERO" || expression == "0" || expression == "5")
    return true;
  const std::vector<std::string> fragments = {
      "_m", "_s", "_kg", "_c", "_k", "_tesla", "_pa", "_j",
      "_m3", "_rad", "_fraction", "_ratio", "_index", "_degree",
      "_points", "_level", "_steps", "_u64", "_slot", "_samples",
      "_cadence", "_order", "_power", "_tolerance", "_quantile",
      "_rate", "_count", "_number", "_weight", "_efficiency"};
  for (const std::string& fragment : fragments) {
    if (key.size() >= fragment.size() &&
        key.compare(key.size() - fragment.size(), fragment.size(), fragment) == 0)
      return true;
  }
  return false;
}

bool IsNumericListKey(const std::string& key) {
  return key.find("radii_m") != std::string::npos ||
         key.find("ensemble_m") != std::string::npos ||
         key.find("coefficients") != std::string::npos;
}

Core::Status ValidateNumber(const std::string& dottedKey,
                            const std::string& value, bool allowList) {
  const std::vector<std::string> fields = allowList ? Split(value, ",")
                                                     : std::vector<std::string>{value};
  if (fields.empty()) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 dottedKey + " requires a bare SI scalar");
  }
  for (const std::string& field : fields) {
    double parsed = 0.0;
    if (!ParseFiniteDouble(field, &parsed)) {
      return Core::Status::Failure(
          Core::StatusCode::InvalidConfiguration,
          dottedKey + " must contain only finite bare SI scalar(s); got '" +
              value + "'");
    }
  }
  return Core::Status::Success();
}

Core::Status ValidateValue(const GrammarEntry& grammar,
                           const Assignment& assignment) {
  const std::string dotted = assignment.section + "." + assignment.key;
  const std::string expression = grammar.expression;
  const std::vector<std::string> alternatives = Split(expression, "|");
  bool permitsRequired = false;
  for (const std::string& alternative : alternatives) {
    if (alternative == assignment.value) return Core::Status::Success();
    if (alternative == "REQUIRED" || alternative == "REQUIRED_OR_ZERO")
      permitsRequired = true;
  }
  if (!permitsRequired) {
    return Core::Status::Failure(
        Core::StatusCode::InvalidConfiguration,
        "line " + std::to_string(assignment.line) + ": invalid value '" +
            assignment.value + "' for " + dotted + "; expected " + expression);
  }
  if (assignment.value.empty() || assignment.value == "REQUIRED" ||
      assignment.value == "REQUIRED_OR_ZERO") {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 dotted + " is incomplete");
  }
  if (IsNumericKey(assignment.key, expression)) {
    return ValidateNumber(dotted, assignment.value,
                          IsNumericListKey(assignment.key));
  }
  // Non-numeric REQUIRED fields are identifiers, paths, frames, epochs, or
  // checksums.  Whitespace is forbidden so a unit-bearing number such as
  // "500 km/s" cannot pass by being misclassified as an opaque string.
  if (assignment.value.find_first_of(" \t\r\n") != std::string::npos) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 dotted + " must be one token");
  }
  return Core::Status::Success();
}

Core::Result<std::vector<Assignment>> Lex(const std::string& input) {
  std::vector<Assignment> assignments;
  std::istringstream stream(input);
  std::string section;
  std::string line;
  std::size_t lineNumber = 0;
  while (std::getline(stream, line)) {
    ++lineNumber;
    const std::string trimmed = Trim(line);
    if (trimmed.empty() || trimmed.front() == '#' || trimmed.front() == '!')
      continue;
    if (trimmed.front() == '[' && trimmed.back() == ']') {
      section = Trim(trimmed.substr(1u, trimmed.size() - 2u));
      if (section.empty()) {
        return Core::Result<std::vector<Assignment>>::Failure(
            Core::StatusCode::InvalidConfiguration,
            "line " + std::to_string(lineNumber) + ": empty section");
      }
      continue;
    }
    const std::size_t equals = trimmed.find('=');
    if (section.empty() || equals == std::string::npos ||
        trimmed.find('=', equals + 1u) != std::string::npos) {
      return Core::Result<std::vector<Assignment>>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "line " + std::to_string(lineNumber) +
              ": expected one key = value assignment inside a section");
    }
    Assignment assignment;
    assignment.section = section;
    assignment.key = Trim(trimmed.substr(0u, equals));
    assignment.value = Trim(trimmed.substr(equals + 1u));
    assignment.line = lineNumber;
    if (assignment.key.empty() || assignment.value.empty()) {
      return Core::Result<std::vector<Assignment>>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "line " + std::to_string(lineNumber) + ": empty key or value");
    }
    assignments.push_back(std::move(assignment));
  }
  return Core::Result<std::vector<Assignment>>::Success(std::move(assignments));
}

template <typename Enum>
Core::Status ParseEnum(const std::map<std::string, std::string>& assignments,
                       const std::string& key,
                       const std::vector<std::pair<std::string, Enum>>& choices,
                       Enum* result) {
  const auto found = assignments.find(key);
  if (found == assignments.end()) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "missing required key " + key);
  }
  for (const auto& choice : choices) {
    if (found->second == choice.first) {
      *result = choice.second;
      return Core::Status::Success();
    }
  }
  return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                               "invalid enum value for " + key);
}

double RequiredDouble(const std::map<std::string, std::string>& assignments,
                      const std::string& key, Core::Status* status) {
  const auto found = assignments.find(key);
  double value = 0.0;
  if (found == assignments.end() || !ParseFiniteDouble(found->second, &value)) {
    *status = Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                    "missing or invalid numeric key " + key);
  }
  return value;
}

Core::Status BuildTyped(const std::map<std::string, std::string>& assignments,
                        ModelConfiguration* configuration) {
  if (configuration == nullptr) {
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                 "null configuration output");
  }
  configuration->schemaVersion = 5;
  configuration->assignments = assignments;
  Core::Status status;
  if (!(status = ParseEnum(assignments, "run.intent", {
          {"production-shock-injection", RunIntent::ProductionShockInjection},
          {"analytic-verification", RunIntent::AnalyticVerification}},
          &configuration->intent)).ok()) return status;
  if (!(status = ParseEnum(assignments, "run.transport", {
          {"ballistic-verification", TransportModel::BallisticVerification},
          {"parker", TransportModel::Parker},
          {"focused-pitch-angle-diffusion", TransportModel::FocusedPitchAngleDiffusion},
          {"focused-discrete-scattering", TransportModel::FocusedDiscreteScattering}},
          &configuration->transport)).ok()) return status;
  if (!(status = ParseEnum(assignments, "run.transport_frame", {
          {"inertial", TransportFrame::Inertial},
          {"rigid-corotating", TransportFrame::RigidCorotating}},
          &configuration->transportFrame)).ok()) return status;
  if (!(status = ParseEnum(assignments, "solar_rotation.model", {
          {"rigid", SolarRotationModel::Rigid},
          {"latitude-dependent-verification",
           SolarRotationModel::LatitudeDependentVerification}},
          &configuration->rotationModel)).ok()) return status;
  if (!(status = ParseEnum(assignments, "solar_wind.model", {
          {"flux-tube-polytropic", WindModel::FluxTubePolytropic},
          {"empirical-kinematic", WindModel::EmpiricalKinematic}},
          &configuration->windModel)).ok()) return status;
  if (!(status = ParseEnum(assignments, "solar_wind.energy_closure", {
          {"isothermal", WindEnergyClosure::Isothermal},
          {"polytropic", WindEnergyClosure::Polytropic},
          {"empirical-profile", WindEnergyClosure::EmpiricalProfile}},
          &configuration->windEnergyClosure)).ok()) return status;
  if (!(status = ParseEnum(assignments, "closed_field_plasma.model", {
          {"isothermal-hydrostatic", ClosedFieldModel::IsothermalHydrostatic},
          {"polytropic-hydrostatic", ClosedFieldModel::PolytropicHydrostatic}},
          &configuration->closedFieldModel)).ok()) return status;
  if (!(status = ParseEnum(assignments, "open_closed_interface.representation", {
          {"sharp-one-sided", InterfaceRepresentation::SharpOneSided},
          {"finite-width-volume", InterfaceRepresentation::FiniteWidthVolume}},
          &configuration->interfaceRepresentation)).ok()) return status;
  if (!(status = ParseEnum(assignments, "open_closed_interface.policy", {
          {"diagnostic-kinematic", InterfacePolicy::DiagnosticKinematic},
          {"bounded-approximation", InterfacePolicy::BoundedApproximation},
          {"stationary-td-equilibrium",
           InterfacePolicy::StationaryTangentialDiscontinuity}},
          &configuration->interfacePolicy)).ok()) return status;

  configuration->startTimeS = RequiredDouble(assignments, "run.start_time_s", &status);
  if (!status.ok()) return status;
  configuration->endTimeS = RequiredDouble(assignments, "run.end_time_s", &status);
  if (!status.ok()) return status;
  const auto seed = assignments.find("run.campaign_seed_u64");
  if (seed == assignments.end() || !ParseUnsigned64(seed->second,
                                                     &configuration->campaignSeed)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "run.campaign_seed_u64 must be an unsigned integer");
  }
  configuration->solarRadiusM = RequiredDouble(assignments, "domain.solar_radius_m", &status);
  if (!status.ok()) return status;
  configuration->outerRadiusM = RequiredDouble(assignments, "domain.outer_radius_m", &status);
  if (!status.ok()) return status;
  configuration->pfssOuterBoundaryRadiusM =
      RequiredDouble(assignments, "pfss.outer_boundary_radius_m", &status);
  if (!status.ok()) return status;
  configuration->currentSheetInterfaceRadiusM =
      RequiredDouble(assignments, "current_sheet.interface_radius_m", &status);
  if (!status.ok()) return status;
  configuration->currentSheetOuterRadiusM =
      RequiredDouble(assignments, "current_sheet.outer_radial_radius_m", &status);
  if (!status.ok()) return status;
  configuration->siderealRotationRateRadPerS = RequiredDouble(
      assignments, "solar_rotation.rigid_input_rotation_rate_rad_per_s", &status);
  if (!status.ok()) return status;
  configuration->gammaWind =
      RequiredDouble(assignments, "solar_wind.polytropic_index", &status);
  if (!status.ok()) return status;
  configuration->gammaAdiabatic =
      RequiredDouble(assignments, "plasma_eos.adiabatic_index", &status);
  if (!status.ok()) return status;
  configuration->gammaClosed = RequiredDouble(
      assignments, "closed_field_plasma.closed_polytropic_index", &status);
  if (!status.ok()) return status;
  const auto electronMass = assignments.find("plasma_eos.electron_mass_in_density");
  configuration->includeElectronMassInDensity =
      electronMass != assignments.end() && electronMass->second == "include";
  configuration->physicsFingerprint = ComputePhysicsFingerprint(assignments);
  return ValidateConfiguration(*configuration);
}

}  // namespace

Core::Result<VersionedConfiguration> ParseConfiguration(
    const std::string& inputBytes) {
  const Core::Result<std::vector<Assignment>> lexical = Lex(inputBytes);
  if (!lexical.ok()) {
    return Core::Result<VersionedConfiguration>::Failure(
        lexical.status.code, lexical.status.message);
  }

  // Pass one is intentionally ignorant of schema-5 keys.  It finds exactly
  // one version selector and cannot accidentally validate a new deck with an
  // old grammar.
  int schemaVersion = 0;
  std::size_t versionCount = 0;
  for (const Assignment& assignment : lexical.value) {
    if (assignment.section == "run" && assignment.key == "schema_version") {
      ++versionCount;
      double parsed = 0.0;
      if (!ParseFiniteDouble(assignment.value, &parsed) ||
          std::floor(parsed) != parsed || parsed < 1.0 || parsed > 5.0) {
        return Core::Result<VersionedConfiguration>::Failure(
            Core::StatusCode::InvalidConfiguration,
            "line " + std::to_string(assignment.line) +
                ": run.schema_version must be an integer in [1,5]");
      }
      schemaVersion = static_cast<int>(parsed);
    }
  }
  if (versionCount != 1u) {
    return Core::Result<VersionedConfiguration>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "configuration requires exactly one run.schema_version assignment");
  }

  VersionedConfiguration resolved;
  if (schemaVersion < 5) {
    resolved.disposition = ParseDisposition::LegacyPassThrough;
    resolved.legacySchemaVersion = schemaVersion;
    resolved.legacyBytes = inputBytes;
    resolved.fingerprint = Internal::Sha256Hex(inputBytes);
    return Core::Result<VersionedConfiguration>::Success(std::move(resolved));
  }

  std::map<std::string, std::string> normalized;
  std::set<std::string> presentSections;
  for (const Assignment& assignment : lexical.value) {
    const GrammarEntry* grammar = FindGrammar(assignment.section, assignment.key);
    if (grammar == nullptr) {
      return Core::Result<VersionedConfiguration>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "line " + std::to_string(assignment.line) + ": unknown schema-5 key " +
              assignment.section + "." + assignment.key);
    }
    const std::string dotted = assignment.section + "." + assignment.key;
    if (normalized.count(dotted) != 0u) {
      return Core::Result<VersionedConfiguration>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "line " + std::to_string(assignment.line) + ": duplicate key " + dotted);
    }
    const Core::Status valueStatus = ValidateValue(*grammar, assignment);
    if (!valueStatus.ok()) {
      return Core::Result<VersionedConfiguration>::Failure(
          valueStatus.code, valueStatus.message);
    }
    normalized.emplace(dotted, assignment.value);
    presentSections.insert(assignment.section);
  }

  // Every key in a present schema-5 section is mandatory.  This is how the
  // grammar distinguishes an explicit inactive zero/none from an omitted
  // physical choice.  Repeated ID sections are validated per concrete ID.
  for (const std::string& concreteSection : presentSections) {
    for (const GrammarEntry& grammar : kSchema5Grammar) {
      if (SectionMatches(grammar.sectionPattern, concreteSection)) {
        const std::string dotted = concreteSection + "." + grammar.key;
        if (normalized.count(dotted) == 0u) {
          return Core::Result<VersionedConfiguration>::Failure(
              Core::StatusCode::InvalidConfiguration,
              "incomplete section [" + concreteSection + "]: missing " + grammar.key);
        }
      }
    }
  }

  // Core sections are mandatory even when a later provider is unavailable.
  const std::vector<std::string> mandatorySections = {
      "run", "domain", "solar_rotation", "pfss", "current_sheet",
      "closed_field_plasma", "open_closed_interface", "solar_wind",
      "plasma_eos"};
  for (const std::string& section : mandatorySections) {
    if (presentSections.count(section) == 0u) {
      return Core::Result<VersionedConfiguration>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "missing mandatory schema-5 section [" + section + "]");
    }
  }

  const Core::Status typedStatus = BuildTyped(normalized, &resolved.schema5);
  if (!typedStatus.ok()) {
    return Core::Result<VersionedConfiguration>::Failure(
        typedStatus.code, typedStatus.message);
  }
  resolved.disposition = ParseDisposition::ResolvedSchema5;
  resolved.fingerprint = resolved.schema5.physicsFingerprint;
  return Core::Result<VersionedConfiguration>::Success(std::move(resolved));
}

Core::Status ValidateConfiguration(const ModelConfiguration& configuration) {
  if (configuration.schemaVersion != 5) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "typed SCCM configuration must use schema 5");
  }
  if (!(configuration.startTimeS < configuration.endTimeS)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "run.start_time_s must be less than run.end_time_s");
  }
  if (!(configuration.solarRadiusM > 0.0 &&
        configuration.pfssOuterBoundaryRadiusM > configuration.solarRadiusM &&
        configuration.outerRadiusM > configuration.pfssOuterBoundaryRadiusM)) {
    return Core::Status::Failure(
        Core::StatusCode::InvalidConfiguration,
        "radii must satisfy 0 < R_sun < R_b < domain.outer_radius_m");
  }
  const auto sheetModel = configuration.assignments.find("current_sheet.model");
  if (sheetModel != configuration.assignments.end() &&
      sheetModel->second == "finite-shell-schatten") {
    if (!(configuration.currentSheetInterfaceRadiusM >= configuration.solarRadiusM &&
          configuration.currentSheetInterfaceRadiusM <=
              configuration.pfssOuterBoundaryRadiusM &&
          configuration.currentSheetOuterRadiusM >
              configuration.currentSheetInterfaceRadiusM &&
          configuration.currentSheetOuterRadiusM < configuration.outerRadiusM)) {
      return Core::Status::Failure(
          Core::StatusCode::InvalidConfiguration,
          "finite SCS requires R_sun <= R_i <= R_b and R_i < R_scs < R_out");
    }
  } else if (sheetModel != configuration.assignments.end() &&
             sheetModel->second == "none" &&
             (configuration.currentSheetInterfaceRadiusM != 0.0 ||
              configuration.currentSheetOuterRadiusM != 0.0)) {
    return Core::Status::Failure(
        Core::StatusCode::InvalidConfiguration,
        "inactive current-sheet radii must be exactly zero");
  }
  if (!(configuration.gammaAdiabatic > 1.0)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "plasma_eos.adiabatic_index must exceed one");
  }
  if (configuration.windEnergyClosure == WindEnergyClosure::Polytropic &&
      !(configuration.gammaWind > 1.0)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "polytropic wind requires gamma_w > 1");
  }
  if (configuration.windEnergyClosure == WindEnergyClosure::Isothermal &&
      configuration.gammaWind != 0.0) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "isothermal wind requires inactive gamma_w=0");
  }
  if (configuration.closedFieldModel == ClosedFieldModel::PolytropicHydrostatic &&
      !(configuration.gammaClosed > 1.0)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "polytropic closed plasma requires gamma_c > 1");
  }
  if (configuration.closedFieldModel == ClosedFieldModel::IsothermalHydrostatic &&
      configuration.gammaClosed != 0.0) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "isothermal closed plasma requires inactive gamma_c=0");
  }

  const auto rotationConvention =
      configuration.assignments.find("solar_rotation.input_rate_convention");
  const auto ephemeris = configuration.assignments.find(
      "solar_rotation.synodic_conversion_ephemeris_file");
  const auto differential = configuration.assignments.find(
      "solar_rotation.differential_rotation_coefficients_file");
  if (configuration.rotationModel == SolarRotationModel::Rigid) {
    if (!(configuration.siderealRotationRateRadPerS > 0.0)) {
      return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                   "rigid rotation requires a positive rate");
    }
    if (differential != configuration.assignments.end() &&
        differential->second != "none") {
      return Core::Status::Failure(
          Core::StatusCode::InvalidConfiguration,
          "rigid rotation requires inactive differential coefficients");
    }
    if (rotationConvention != configuration.assignments.end() &&
        rotationConvention->second == "sidereal" &&
        ephemeris != configuration.assignments.end() && ephemeris->second != "none") {
      return Core::Status::Failure(
          Core::StatusCode::InvalidConfiguration,
          "sidereal rotation requires no synodic conversion asset");
    }
    if (rotationConvention != configuration.assignments.end() &&
        rotationConvention->second == "synodic" &&
        (ephemeris == configuration.assignments.end() || ephemeris->second == "none")) {
      return Core::Status::Failure(
          Core::StatusCode::InvalidConfiguration,
          "synodic rotation requires an ephemeris conversion asset");
    }
  } else {
    if (configuration.intent != RunIntent::AnalyticVerification) {
      return Core::Status::Failure(
          Core::StatusCode::InvalidConfiguration,
          "latitude-dependent rotation is verification-only in schema 5");
    }
    if (configuration.siderealRotationRateRadPerS != 0.0 ||
        differential == configuration.assignments.end() ||
        differential->second == "none") {
      return Core::Status::Failure(
          Core::StatusCode::InvalidConfiguration,
          "differential rotation requires zero rigid rate and a coefficient asset");
    }
  }

  // R1 radialization policy is validated even before the Stage-3 SCS provider
  // exists.  Stage 0 freezes the selector and inactive-value semantics; Stage
  // 3 will evaluate the actual zonal-power and latitude diagnostics.
  const auto radialization =
      configuration.assignments.find("current_sheet.radialization_gate");
  if (radialization != configuration.assignments.end()) {
    const auto requiredNumber = [&](const std::string& key, double* value) {
      const auto item = configuration.assignments.find(key);
      return item != configuration.assignments.end() &&
             ParseFiniteDouble(item->second, value);
    };
    double maximumPower = 0.0;
    double minimumCoverage = 0.0;
    double maximumRms = 0.0;
    double maximumRatio = 0.0;
    if (!requiredNumber(
            "current_sheet.maximum_outer_zonal_nonmonopole_power_fraction",
            &maximumPower) ||
        !requiredNumber(
            "current_sheet.latitude_minimum_unmasked_longitude_fraction",
            &minimumCoverage) ||
        !requiredNumber("current_sheet.maximum_unsigned_radial_flux_rms_fraction",
                        &maximumRms) ||
        !requiredNumber(
            "current_sheet.maximum_unsigned_radial_flux_p95_to_p05_ratio",
            &maximumRatio)) {
      return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                   "SCS radialization thresholds are incomplete");
    }
    if (radialization->second == "diagnostic-only") {
      if (configuration.intent != RunIntent::AnalyticVerification ||
          maximumPower != 0.0 || minimumCoverage != 0.0 || maximumRms != 0.0 ||
          maximumRatio != 0.0) {
        return Core::Status::Failure(
            Core::StatusCode::InvalidConfiguration,
            "diagnostic-only SCS requires analytic intent and zero acceptance bounds");
      }
    } else if (radialization->second ==
               "outer-zonal-power-and-latitude-flatness") {
      if (!(maximumPower > 0.0 && maximumPower < 1.0 &&
            minimumCoverage > 0.0 && minimumCoverage <= 1.0 &&
            maximumRms > 0.0 && maximumRatio > 1.0)) {
        return Core::Status::Failure(
            Core::StatusCode::InvalidConfiguration,
            "production radialization bounds have invalid ranges");
      }
    } else if (radialization->second == "not-applicable") {
      if (sheetModel != configuration.assignments.end() &&
          sheetModel->second != "none") {
        return Core::Status::Failure(
            Core::StatusCode::InvalidConfiguration,
            "not-applicable radialization requires current_sheet.model=none");
      }
    }
  }
  return Core::Status::Success();
}

std::string ComputePhysicsFingerprint(
    const std::map<std::string, std::string>& normalizedAssignments) {
  std::string bytes = "sep-coronal-cme-schema5-fingerprint-v1\n";
  for (const auto& entry : normalizedAssignments) {
    const std::string& key = entry.first;
    // File-system spelling is not physics.  Asset ingest appends a derived
    // ``*.content_checksum`` entry, which *is* retained by this serializer.
    if (key.size() >= 5u && key.compare(key.size() - 5u, 5u, "_file") == 0)
      continue;
    if (key.find("comment") != std::string::npos) continue;
    bytes += std::to_string(key.size()) + ":" + key;
    bytes += std::to_string(entry.second.size()) + ":" + entry.second + "\n";
  }
  return Internal::Sha256Hex(bytes);
}

std::string ComputeContentChecksum(const std::string& bytes) {
  return Internal::Sha256Hex(bytes);
}

Core::Status CheckRestartIdentity(const std::string& expectedFingerprint,
                                  const std::string& actualFingerprint,
                                  const std::string& identityCategory) {
  if (expectedFingerprint.empty() || actualFingerprint.empty()) {
    return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                 identityCategory + " fingerprint is missing");
  }
  if (expectedFingerprint != actualFingerprint) {
    return Core::Status::Failure(
        Core::StatusCode::DataIntegrityFailure,
        identityCategory + " fingerprint differs from restart identity");
  }
  return Core::Status::Success();
}

Core::Status CheckCapabilityAvailability(
    const ModelConfiguration& configuration, int completedStage) {
  if (completedStage < 0) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                 "completed stage cannot be negative");
  }
  if (completedStage < 1) {
    return Core::Status::Failure(
        Core::StatusCode::NotImplemented,
        "PFSS provider is registered but unavailable before Stage 1");
  }
  if (completedStage < 2) {
    return Core::Status::Failure(
        Core::StatusCode::NotImplemented,
        "wind and closed-plasma providers are unavailable before Stage 2");
  }
  const auto sheet = configuration.assignments.find("current_sheet.model");
  if (sheet != configuration.assignments.end() &&
      sheet->second == "finite-shell-schatten" && completedStage < 3) {
    return Core::Status::Failure(
        Core::StatusCode::NotImplemented,
        "finite-shell SCS is a Stage-3 provider; no Parker/PFSS fallback is allowed");
  }
  return Core::Status::Success();
}

}  // namespace CoronalCME
}  // namespace SEP
