#include "ForwardModel.h"

#include "../src/general/input/amps_input.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace Enceladus {
namespace ForwardModel {
namespace {

namespace fs = std::filesystem;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kElectronVoltJ = 1.602176634e-19;

std::string Trim(const std::string& value) {
  const std::size_t first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return std::string();
  const std::size_t last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string Unquote(const std::string& value) {
  const std::string text = Trim(value);
  if (text.size() >= 2 &&
      ((text.front() == '"' && text.back() == '"') ||
       (text.front() == '\'' && text.back() == '\'')))
    return text.substr(1, text.size() - 2);
  return text;
}

bool ParseFinite(const std::string& text, double* output) {
  const std::string scalar = Unquote(text);
  errno = 0;
  char* end = nullptr;
  const double value = std::strtod(scalar.c_str(), &end);
  if (errno == ERANGE || end == scalar.c_str() || *end != '\0' ||
      !std::isfinite(value)) return false;
  *output = value;
  return true;
}

bool ParsePositive(const std::string& text, double* output) {
  return ParseFinite(text, output) && *output > 0.0;
}

bool ParseNonnegative(const std::string& text, double* output) {
  return ParseFinite(text, output) && *output >= 0.0;
}

bool ParsePositiveU64(const std::string& text, std::uint64_t* output) {
  const std::string scalar = Unquote(text);
  if (scalar.empty() || scalar.front() == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long value = std::strtoull(scalar.c_str(), &end, 10);
  if (errno == ERANGE || end == scalar.c_str() || *end != '\0' || value == 0)
    return false;
  *output = static_cast<std::uint64_t>(value);
  return true;
}

bool ParseNonnegativeU64(const std::string& text, std::uint64_t* output) {
  const std::string scalar = Unquote(text);
  if (scalar.empty() || scalar.front() == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long value = std::strtoull(scalar.c_str(), &end, 10);
  if (errno == ERANGE || end == scalar.c_str() || *end != '\0') return false;
  *output = static_cast<std::uint64_t>(value);
  return true;
}

bool ParsePositiveLong(const std::string& text, long* output) {
  const std::string scalar = Unquote(text);
  errno = 0;
  char* end = nullptr;
  const long value = std::strtol(scalar.c_str(), &end, 10);
  if (errno == ERANGE || end == scalar.c_str() || *end != '\0' || value <= 0)
    return false;
  *output = value;
  return true;
}

bool ParsePositiveInt(const std::string& text, int* output) {
  long value = 0;
  if (!ParsePositiveLong(text, &value) ||
      value > std::numeric_limits<int>::max()) return false;
  *output = static_cast<int>(value);
  return true;
}

bool ParseBoolean(const std::string& text, bool* output) {
  const std::string value = Unquote(text);
  if (value == "true") {
    *output = true;
    return true;
  }
  if (value == "false") {
    *output = false;
    return true;
  }
  return false;
}

std::vector<std::string> ParseList(const std::string& text) {
  std::vector<std::string> result;
  std::stringstream stream(Unquote(text));
  std::string field;
  while (std::getline(stream, field, ',')) {
    field = Trim(field);
    if (!field.empty()) result.push_back(field);
  }
  return result;
}

std::vector<std::string> ParseOptionalList(const std::string& text) {
  if (Unquote(text) == "none") return {};
  return ParseList(text);
}

Status InputError(const AMPS::Input::Assignment& item,
                  const std::string& reason) {
  return Status::Error(AMPS::Input::ErrorAt(item.where, reason).message);
}

Status BlockError(const AMPS::Input::Block& block,
                  const std::string& reason) {
  return Status::Error(AMPS::Input::ErrorAt(block.begin, reason).message);
}

bool Contains(const std::vector<std::string>& values,
              const std::string& value) {
  return std::find(values.begin(), values.end(), value) != values.end();
}

bool IsPortableName(const std::string& value) {
  if (value.empty()) return false;
  for (const unsigned char character : value)
    if (!std::isalnum(character) && character != '_' && character != '-' &&
        character != '.') return false;
  return true;
}

bool IsSpeciesSymbol(const std::string& value) {
  if (value.empty() || !std::isalpha(static_cast<unsigned char>(value.front())))
    return false;
  for (const unsigned char character : value)
    if (!std::isalnum(character) && character != '_') return false;
  return true;
}

// PIC::CPLR::DATAFILE::ICES::retriveSWMFdata composes its installation path
// and case into an unquoted POSIX shell command. Restrict the explicitly
// selected fallback authority to characters that are both shell-literal and
// valid in the path layout used by that adapter. cache_only never reaches it.
bool IsSafeIcesInstallationPath(const std::string& value) {
  if (value.empty()) return false;
  for (const unsigned char character : value)
    if (!std::isalnum(character) && character != '/' && character != '_' &&
        character != '-' && character != '.') return false;
  return true;
}

bool IsSafeIcesCase(const std::string& value) {
  return value != "." && value != ".." && IsPortableName(value);
}

Status RequireUniqueNonempty(const std::vector<std::string>& values,
                             const std::string& label) {
  if (values.empty()) return Status::Error(label + " must not be empty");
  std::set<std::string> unique;
  for (const std::string& value : values) {
    if (value.empty() || !unique.insert(value).second)
      return Status::Error(label + " contains an empty or duplicate value: " +
                           value);
  }
  return Status::Ok();
}

double Norm(const Vec3& value) {
  return std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
}

Vec3 Normalize(const Vec3& value) {
  const double norm = Norm(value);
  return Vec3{value.x / norm, value.y / norm, value.z / norm};
}

Status ResolveInputPath(const fs::path& inputFile, const std::string& value,
                        std::string* output) {
  fs::path path = Unquote(value);
  if (path.empty()) return Status::Error("path must not be empty");
  if (path.is_relative()) path = inputFile.parent_path() / path;
  std::error_code ec;
  path = fs::absolute(path, ec).lexically_normal();
  if (ec) return Status::Error("cannot resolve path: " + path.string());
  *output = path.string();
  return Status::Ok();
}

const std::vector<std::string>& HybridSourceIds() {
  static const std::vector<std::string> ids = {
      "I_baghdad", "II_damascus", "III_damascus", "IV_alexandria",
      "V_cairo", "VI_baghdad", "VII_baghdad", "VIII_cairo",
      "alexandria", "baghdad", "cairo", "damascus", "baghdad_extra",
      "cairo_extra", "damascus_extra"};
  return ids;
}

std::map<std::string, const AMPS::Input::Assignment*> AssignmentMap(
    const AMPS::Input::Block& block) {
  std::map<std::string, const AMPS::Input::Assignment*> result;
  for (const auto& item : block.assignments) result[item.key] = &item;
  return result;
}

Status RequireKeys(
    const AMPS::Input::Block& block,
    const std::map<std::string, const AMPS::Input::Assignment*>& settings,
    const std::vector<std::string>& keys) {
  for (const std::string& key : keys) {
    if (settings.find(key) == settings.end())
      return BlockError(block, "missing required key '" + key + "'");
  }
  return Status::Ok();
}

Status ParseVertex(const AMPS::Input::Block& block, Vec3* vertex) {
  if (!block.blocks.empty())
    return BlockError(block, "vertex blocks cannot contain child blocks");
  const auto settings = AssignmentMap(block);
  Status status = RequireKeys(block, settings, {"x_m", "y_m", "z_m"});
  if (!status.success) return status;
  for (const auto& item : settings) {
    if (item.first != "x_m" && item.first != "y_m" && item.first != "z_m")
      return InputError(*item.second, "unknown vertex key '" + item.first + "'");
  }
  if (!ParseFinite(settings.at("x_m")->value, &vertex->x) ||
      !ParseFinite(settings.at("y_m")->value, &vertex->y) ||
      !ParseFinite(settings.at("z_m")->value, &vertex->z))
    return BlockError(block, "vertex coordinates must be finite SI metres");
  if (!(Norm(*vertex) > 0.0))
    return BlockError(block, "vertex cannot be at the body centre");
  return Status::Ok();
}

Status ParseSourceBlock(const AMPS::Input::Block& block, Source* source) {
  if (block.name.empty()) return BlockError(block, "source id is required");
  source->id = block.name;
  const auto settings = AssignmentMap(block);
  Status status = RequireKeys(
      block, settings, {"geometry", "temperature_K", "bulk_speed_m_s"});
  if (!status.success) return status;
  const std::string geometry = Unquote(settings.at("geometry")->value);
  if (geometry == "point_jet") source->geometry = SourceGeometry::PointJet;
  else if (geometry == "polyline_stripe")
    source->geometry = SourceGeometry::PolylineStripe;
  else return InputError(*settings.at("geometry"),
                         "geometry must be point_jet or polyline_stripe");
  if (!ParsePositive(settings.at("temperature_K")->value,
                     &source->temperatureK))
    return InputError(*settings.at("temperature_K"),
                      "temperature_K must be finite and positive");
  if (!ParseNonnegative(settings.at("bulk_speed_m_s")->value,
                        &source->bulkSpeedMPerS))
    return InputError(*settings.at("bulk_speed_m_s"),
                      "bulk_speed_m_s must be finite and nonnegative");

  const std::set<std::string> common = {
      "geometry", "temperature_K", "bulk_speed_m_s"};
  const std::set<std::string> point = {
      "latitude_deg", "west_longitude_deg", "tilt_deg", "azimuth_deg"};
  for (const auto& item : settings) {
    if (common.count(item.first) != 0) continue;
    if (item.first.rfind("strength.", 0) == 0) {
      const std::string species = item.first.substr(9);
      double value = 0.0;
      if (!IsSpeciesSymbol(species) ||
          !ParseNonnegative(item.second->value, &value))
        return InputError(*item.second,
            "source strength species must match [A-Za-z][A-Za-z0-9_]* and "
            "its value must be finite and nonnegative");
      source->strength[species] = value;
      continue;
    }
    if (source->geometry == SourceGeometry::PointJet &&
        point.count(item.first) != 0) continue;
    return InputError(*item.second, "unknown source key '" + item.first + "'");
  }
  if (source->strength.empty())
    return BlockError(block, "at least one strength.SPECIES key is required");

  if (source->geometry == SourceGeometry::PointJet) {
    for (const std::string& key : point)
      if (settings.find(key) == settings.end())
        return BlockError(block, "point_jet is missing required key '" + key + "'");
    if (!block.blocks.empty())
      return BlockError(block, "point_jet cannot contain vertex blocks");
    if (!ParseFinite(settings.at("latitude_deg")->value, &source->latitudeDeg) ||
        source->latitudeDeg < -90.0 || source->latitudeDeg > 90.0)
      return InputError(*settings.at("latitude_deg"),
                        "latitude_deg must be in [-90,90]");
    if (!ParseFinite(settings.at("west_longitude_deg")->value,
                     &source->westLongitudeDeg) ||
        source->westLongitudeDeg < 0.0 || source->westLongitudeDeg >= 360.0)
      return InputError(*settings.at("west_longitude_deg"),
                        "west_longitude_deg must be in [0,360)");
    if (!ParseFinite(settings.at("tilt_deg")->value, &source->tiltDeg) ||
        source->tiltDeg < 0.0 || source->tiltDeg > 90.0)
      return InputError(*settings.at("tilt_deg"),
                        "tilt_deg must be in [0,90]");
    if (!ParseFinite(settings.at("azimuth_deg")->value, &source->azimuthDeg) ||
        source->azimuthDeg < 0.0 || source->azimuthDeg >= 360.0)
      return InputError(*settings.at("azimuth_deg"),
                        "azimuth_deg must be in [0,360)");
  } else {
    for (const std::string& key : point)
      if (settings.find(key) != settings.end())
        return InputError(*settings.at(key),
                          "polyline_stripe cannot use point-jet angles");
    const std::vector<const AMPS::Input::Block*> vertices =
        block.BlocksOfKind("vertex");
    if (vertices.size() < 2 || vertices.size() != block.blocks.size())
      return BlockError(block,
          "polyline_stripe requires at least two vertex child blocks and no other children");
    std::set<std::string> vertexIds;
    for (const AMPS::Input::Block* vertexBlock : vertices) {
      if (vertexBlock->name.empty() || !vertexIds.insert(vertexBlock->name).second)
        return BlockError(*vertexBlock, "vertex id must be nonempty and unique");
      Vec3 vertex;
      status = ParseVertex(*vertexBlock, &vertex);
      if (!status.success) return status;
      source->verticesM.push_back(vertex);
    }
    std::vector<double> weights;
    status = StripeNodeWeights(*source, &weights);
    if (!status.success) return BlockError(block, status.message);
  }
  return Status::Ok();
}

}  // namespace

const char* Name(SourceMode mode) {
  return mode == SourceMode::Shape ? "shape" : "absolute";
}
const char* Name(SourceGeometry geometry) {
  return geometry == SourceGeometry::PointJet ? "point_jet" : "polyline_stripe";
}
const char* Name(StripeSurfaceMapping mapping) {
  return mapping == StripeSurfaceMapping::LocalVerticalToReferenceSphere
      ? "local_vertical_to_reference_sphere" : "unknown";
}
const char* Name(SolarActivity activity) {
  return activity == SolarActivity::Quiet ? "quiet" : "active";
}
const char* Name(PopulationControl control) {
  switch (control) {
    case PopulationControl::Off: return "off";
    case PopulationControl::VelocityShift: return "velocity_shift";
    case PopulationControl::Scatter: return "scatter";
  }
  return "unknown";
}
const char* Name(ReactionProducts products) {
  return products == ReactionProducts::Discard ? "discard" : "unknown";
}
const char* Name(PlasmaBackgroundMode mode) {
  switch (mode) {
    case PlasmaBackgroundMode::None: return "none";
    case PlasmaBackgroundMode::CacheOnly: return "cache_only";
    case PlasmaBackgroundMode::CacheOrIces: return "cache_or_ices";
  }
  return "unknown";
}

Status ParseSourceDistribution(const std::string& file,
                               SourceDistribution* output) {
  if (output == nullptr) return Status::Error("source distribution output is null");
  AMPS::Input::Document document;
  const AMPS::Input::Status parsed = AMPS::Input::ParseFile(file, &document);
  if (!parsed.ok()) return Status::Error(parsed.message);
  AMPS::Input::Status generic = AMPS::Input::RequireNoAssignments(document.root);
  if (!generic.ok()) return Status::Error(generic.message);
  generic = AMPS::Input::RequireOnlyBlocks(
      document.root, {"source_distribution"}, "the Enceladus source reader");
  if (!generic.ok()) return Status::Error(generic.message);
  if (document.root.blocks.size() != 1)
    return Status::Error(file + ": exactly one source_distribution block is required");
  const AMPS::Input::Block& block = document.root.blocks.front();
  if (block.name.empty()) return BlockError(block, "distribution name is required");
  if (!IsPortableName(block.name))
    return BlockError(block,
        "distribution name may contain only letters, digits, '_', '-', and '.'");
  const auto settings = AssignmentMap(block);
  Status status = RequireKeys(block, settings,
      {"schema_version", "source_case", "frame", "mode", "source_law",
       "stripe_surface_mapping", "stripe_reference_radius_m"});
  if (!status.success) return status;
  for (const auto& item : settings) {
    const std::set<std::string> allowed = {
        "schema_version", "source_case", "frame", "mode", "source_law",
        "stripe_surface_mapping", "stripe_reference_radius_m"};
    if (allowed.count(item.first) == 0)
      return InputError(*item.second,
                        "unknown source_distribution key '" + item.first + "'");
  }

  SourceDistribution candidate;
  candidate.name = block.name;
  std::error_code ec;
  candidate.inputFile = fs::absolute(file, ec).lexically_normal().string();
  if (ec) candidate.inputFile = file;
  long schema = 0;
  if (!ParsePositiveLong(settings.at("schema_version")->value, &schema) || schema != 2)
    return InputError(*settings.at("schema_version"),
                      "only source schema_version=2 is supported");
  candidate.schemaVersion = 2;
  candidate.sourceCase = Unquote(settings.at("source_case")->value);
  candidate.frame = Unquote(settings.at("frame")->value);
  candidate.sourceLaw = Unquote(settings.at("source_law")->value);
  const std::string stripeMapping =
      Unquote(settings.at("stripe_surface_mapping")->value);
  if (stripeMapping != "local_vertical_to_reference_sphere")
    return InputError(*settings.at("stripe_surface_mapping"),
        "stripe_surface_mapping must be local_vertical_to_reference_sphere");
  candidate.stripeSurfaceMapping =
      StripeSurfaceMapping::LocalVerticalToReferenceSphere;
  if (!ParsePositive(settings.at("stripe_reference_radius_m")->value,
                     &candidate.stripeReferenceRadiusM))
    return InputError(*settings.at("stripe_reference_radius_m"),
                      "stripe_reference_radius_m must be finite and positive");
  const std::string mode = Unquote(settings.at("mode")->value);
  if (mode == "shape") candidate.mode = SourceMode::Shape;
  else if (mode == "absolute") candidate.mode = SourceMode::Absolute;
  else return InputError(*settings.at("mode"), "mode must be shape or absolute");
  if (candidate.sourceCase != "tenishev2014_hybrid_8")
    return InputError(*settings.at("source_case"),
        "this implementation supports only tenishev2014_hybrid_8");
  if (candidate.frame != "IAU_ENCELADUS")
    return InputError(*settings.at("frame"), "frame must be IAU_ENCELADUS");
  if (candidate.sourceLaw != "tenishev2010_point")
    return InputError(*settings.at("source_law"),
        "tenishev2014_hybrid_8 requires source_law=tenishev2010_point");

  std::set<std::string> ids;
  for (const AMPS::Input::Block& child : block.blocks) {
    if (child.kind != "source")
      return BlockError(child, "only source child blocks are allowed");
    Source source;
    status = ParseSourceBlock(child, &source);
    if (!status.success) return status;
    if (!ids.insert(source.id).second)
      return BlockError(child, "duplicate source id '" + source.id + "'");
    candidate.sources.push_back(std::move(source));
  }
  const std::vector<std::string>& required = HybridSourceIds();
  if (candidate.sources.size() != required.size())
    return BlockError(block,
        "tenishev2014_hybrid_8 requires exactly 8 jets and 7 stripes");
  for (const std::string& id : required)
    if (ids.count(id) == 0)
      return BlockError(block, "missing canonical hybrid source '" + id + "'");
  for (const Source& source : candidate.sources) {
    const auto found = std::find(required.begin(), required.end(), source.id);
    const bool requiresPointJet =
        std::distance(required.begin(), found) < 8;
    const SourceGeometry requiredGeometry = requiresPointJet
        ? SourceGeometry::PointJet : SourceGeometry::PolylineStripe;
    if (source.geometry != requiredGeometry)
      return BlockError(block, "canonical hybrid source '" + source.id +
          "' must use geometry=" + Name(requiredGeometry));
  }

  std::set<std::string> species;
  for (const Source& source : candidate.sources)
    for (const auto& entry : source.strength) species.insert(entry.first);
  for (const std::string& spec : species) {
    double sum = 0.0;
    for (const Source& source : candidate.sources) {
      const auto found = source.strength.find(spec);
      if (found == source.strength.end())
        return BlockError(block, "every source must define strength." + spec);
      sum += found->second;
    }
    if (!(sum > 0.0) || !std::isfinite(sum))
      return BlockError(block, "strength sum for " + spec + " is not positive");
    if (candidate.mode == SourceMode::Shape && std::abs(sum - 1.0) > 1.0e-12)
      return BlockError(block, "shape strengths for " + spec +
          " must sum to 1 within absolute tolerance 1e-12");
  }
  *output = std::move(candidate);
  return Status::Ok();
}

Status ParseConfiguration(const std::string& file, Configuration* output) {
  if (output == nullptr) return Status::Error("configuration output is null");
  AMPS::Input::Document document;
  const AMPS::Input::Status parsed = AMPS::Input::ParseFile(file, &document);
  if (!parsed.ok()) return Status::Error(parsed.message);
  const AMPS::Input::Block* block = nullptr;
  for (const AMPS::Input::Block& child : document.root.blocks) {
    if (child.kind == "application" && child.name == "enceladus") {
      if (block != nullptr) return BlockError(child, "duplicate Enceladus application block");
      block = &child;
    }
  }
  if (block == nullptr)
    return Status::Error(file + ": required '#begin application: enceladus' block not found");
  if (!block->blocks.empty())
    return BlockError(*block, "application: enceladus is a flat block");
  const auto settings = AssignmentMap(*block);
  const std::vector<std::string> required = {
      "schema_version", "source_distribution_file", "simulated_species",
      "injected_species", "chemistry_species", "selected_reactions",
      "reaction_products",
      "sample_length_steps", "total_iterations", "random_seed",
      "particle_splitting_merging", "solar_activity",
      "plasma_background_mode"};
  Status status = RequireKeys(*block, settings, required);
  if (!status.success) return status;

  Configuration candidate;
  std::error_code ec;
  const fs::path inputPath = fs::absolute(file, ec).lexically_normal();
  candidate.inputFile = ec ? file : inputPath.string();
  long schema = 0;
  if (!ParsePositiveLong(settings.at("schema_version")->value, &schema) || schema != 2)
    return InputError(*settings.at("schema_version"),
                      "only application schema_version=2 is supported");
  candidate.schemaVersion = 2;
  status = ResolveInputPath(inputPath, settings.at("source_distribution_file")->value,
                            &candidate.sourceDistributionFile);
  if (!status.success) return InputError(*settings.at("source_distribution_file"), status.message);
  candidate.simulatedSpecies = ParseList(settings.at("simulated_species")->value);
  candidate.injectedSpecies = ParseList(settings.at("injected_species")->value);
  candidate.chemistrySpecies =
      ParseOptionalList(settings.at("chemistry_species")->value);
  candidate.selectedReactions =
      ParseOptionalList(settings.at("selected_reactions")->value);
  const std::string reactionProducts =
      Unquote(settings.at("reaction_products")->value);
  if (reactionProducts != "discard")
    return InputError(*settings.at("reaction_products"),
        "reaction_products must be discard in application schema 2");
  candidate.reactionProducts = ReactionProducts::Discard;
  if (!ParsePositiveLong(settings.at("sample_length_steps")->value,
                         &candidate.sampleLengthSteps))
    return InputError(*settings.at("sample_length_steps"),
                      "sample_length_steps must be a positive integer");
  if (!ParsePositiveLong(settings.at("total_iterations")->value,
                         &candidate.totalIterations))
    return InputError(*settings.at("total_iterations"),
                      "total_iterations must be a positive integer");
  if (!ParseNonnegativeU64(settings.at("random_seed")->value,
                           &candidate.randomSeed))
    return InputError(*settings.at("random_seed"),
                      "random_seed must be a nonnegative integer");
  const std::string activity = Unquote(settings.at("solar_activity")->value);
  if (activity == "quiet") candidate.solarActivity = SolarActivity::Quiet;
  else if (activity == "active") candidate.solarActivity = SolarActivity::Active;
  else return InputError(*settings.at("solar_activity"),
                         "solar_activity must be quiet or active");
  const std::string control =
      Unquote(settings.at("particle_splitting_merging")->value);
  if (control == "off") candidate.populationControl = PopulationControl::Off;
  else if (control == "velocity_shift")
    candidate.populationControl = PopulationControl::VelocityShift;
  else if (control == "scatter")
    candidate.populationControl = PopulationControl::Scatter;
  else return InputError(*settings.at("particle_splitting_merging"),
      "particle_splitting_merging must be off, velocity_shift, or scatter");

  const std::string plasmaMode =
      Unquote(settings.at("plasma_background_mode")->value);
  if (plasmaMode == "none")
    candidate.plasmaBackgroundMode = PlasmaBackgroundMode::None;
  else if (plasmaMode == "cache_only")
    candidate.plasmaBackgroundMode = PlasmaBackgroundMode::CacheOnly;
  else if (plasmaMode == "cache_or_ices")
    candidate.plasmaBackgroundMode = PlasmaBackgroundMode::CacheOrIces;
  else
    return InputError(*settings.at("plasma_background_mode"),
        "plasma_background_mode must be none, cache_only, or cache_or_ices");
  const auto plasmaCache = settings.find("plasma_cache_directory");
  if (plasmaCache != settings.end()) {
    status = ResolveInputPath(
        inputPath, plasmaCache->second->value, &candidate.plasmaCacheDirectory);
    if (!status.success) return InputError(*plasmaCache->second, status.message);
  }

  for (const auto& entry : settings) {
    const std::string& key = entry.first;
    const AMPS::Input::Assignment& item = *entry.second;
    if (std::find(required.begin(), required.end(), key) != required.end()) continue;
    if (key == "plasma_cache_directory") continue;
    double scalar = 0.0;
    std::uint64_t count = 0;
    if (key.rfind("total_source_rate.", 0) == 0 &&
        key.size() > 20 && key.substr(key.size() - 4) == "_s-1") {
      const std::string species = key.substr(18, key.size() - 22);
      if (!ParsePositive(item.value, &scalar))
        return InputError(item, key + " must be finite and positive");
      candidate.totalSourceRatePerS[species] = scalar;
    } else if (key.rfind("characteristic_speed.", 0) == 0 &&
               key.size() > 27 && key.substr(key.size() - 6) == "_m_s-1") {
      const std::string species = key.substr(21, key.size() - 27);
      if (!ParsePositive(item.value, &scalar))
        return InputError(item, key + " must be finite and positive");
      candidate.characteristicSpeedMPerS[species] = scalar;
    } else if (key.rfind("courant.", 0) == 0) {
      const std::string species = key.substr(8);
      if (!ParsePositive(item.value, &scalar) || scalar > 1.0)
        return InputError(item, key + " must be in (0,1]");
      candidate.courantNumber[species] = scalar;
    } else if (key.rfind("primary_macroparticles_per_step.", 0) == 0) {
      const std::string species = key.substr(32);
      if (!ParsePositiveU64(item.value, &count))
        return InputError(item, key + " must be a positive integer");
      candidate.primaryMacroparticlesPerStep[species] = count;
    } else if (key.rfind("plasma_ion_fraction.", 0) == 0) {
      const std::string species = key.substr(20);
      if (!ParseNonnegative(item.value, &scalar) || scalar > 1.0)
        return InputError(item, key + " must be in [0,1]");
      candidate.plasmaIonFraction[species] = scalar;
    } else if (key == "chemistry_weight_reference_number_density_m-3") {
      if (!ParsePositive(item.value,
          &candidate.chemistryWeightReferenceNumberDensityM3))
        return InputError(item, key + " must be finite and positive");
    } else if (key == "chemistry_weight_reference_relative_speed_m_s-1") {
      if (!ParsePositive(item.value,
          &candidate.chemistryWeightReferenceRelativeSpeedMPerS))
        return InputError(item, key + " must be finite and positive");
    } else if (key == "split_position_shift_cell_diagonal_fraction") {
      if (!ParseNonnegative(
              item.value,
              &candidate.splitPositionShiftCellDiagonalFraction))
        return InputError(item, key + " must be finite and nonnegative");
    } else if (key == "split_velocity_shift_thermal_speed_fraction") {
      if (!ParseNonnegative(
              item.value,
              &candidate.splitVelocityShiftThermalSpeedFraction))
        return InputError(item, key + " must be finite and nonnegative");
    } else if (key == "min_particles_per_cell") {
      if (!ParsePositiveInt(item.value, &candidate.minParticlesPerCell))
        return InputError(item, key + " must be a positive integer");
    } else if (key == "max_particles_per_cell") {
      if (!ParsePositiveInt(item.value, &candidate.maxParticlesPerCell))
        return InputError(item, key + " must be a positive integer");
    } else if (key == "nonuniform_split_position_shift") {
      if (!ParseBoolean(item.value,
          &candidate.nonuniformSplitPositionShift))
        return InputError(item, key + " must be true or false");
    } else if (key == "flyby_density_output") {
      if (!ParseBoolean(item.value, &candidate.flybyDensityOutput))
        return InputError(item, key + " must be true or false");
    } else if (key == "flyby_id") {
      candidate.flybyId = Unquote(item.value);
      if (!IsPortableName(candidate.flybyId))
        return InputError(item,
            "flyby_id may contain only letters, digits, '_', '-', and '.'");
    } else if (key == "flyby_trajectory_file") {
      status = ResolveInputPath(
          inputPath, item.value, &candidate.flybyTrajectoryFile);
      if (!status.success) return InputError(item, status.message);
    } else if (key == "flyby_mesh_refinement") {
      if (!ParseBoolean(item.value, &candidate.flybyMeshRefinement))
        return InputError(item, key + " must be true or false");
    } else if (key == "flyby_refinement_target_cell_size_m") {
      if (!ParsePositive(item.value,
                         &candidate.flybyRefinementTargetCellSizeM))
        return InputError(item, key + " must be finite and positive");
    } else if (key == "flyby_refinement_radius_m") {
      if (!ParsePositive(item.value, &candidate.flybyRefinementRadiusM))
        return InputError(item, key + " must be finite and positive");
    } else if (key == "flyby_refinement_time_min_s") {
      if (!ParseFinite(item.value, &candidate.flybyRefinementTimeMinS))
        return InputError(item, key + " must be finite");
    } else if (key == "flyby_refinement_time_max_s") {
      if (!ParseFinite(item.value, &candidate.flybyRefinementTimeMaxS))
        return InputError(item, key + " must be finite");
    } else if (key == "ices_installation_directory") {
      status = ResolveInputPath(
          inputPath, item.value, &candidate.icesInstallationDirectory);
      if (!status.success) return InputError(item, status.message);
    } else if (key == "ices_case") {
      candidate.icesCase = Unquote(item.value);
      if (candidate.icesCase.empty())
        return InputError(item, "ices_case must not be empty");
    } else {
      return InputError(item, "unknown Enceladus application key '" + key + "'");
    }
  }
  const std::vector<std::string> flybyRefinementControls = {
      "flyby_refinement_target_cell_size_m",
      "flyby_refinement_radius_m",
      "flyby_refinement_time_min_s",
      "flyby_refinement_time_max_s"};
  if (candidate.flybyMeshRefinement) {
    for (const std::string& key : flybyRefinementControls)
      if (settings.find(key) == settings.end())
        return BlockError(*block,
            "flyby_mesh_refinement=true requires key '" + key + "'");
  } else {
    for (const std::string& key : flybyRefinementControls) {
      const auto found = settings.find(key);
      if (found != settings.end())
        return InputError(*found->second,
            "flyby refinement controls require flyby_mesh_refinement=true");
    }
  }
  *output = std::move(candidate);
  return Status::Ok();
}

const std::vector<Reaction>& Burger2007Reactions() {
  // Burger et al. (2007), JGR 112 A06219, Tables 1 and 2. Cross sections are
  // converted from A^2 to m^2 (1 A^2=1e-20 m^2); photo rates are already the
  // quiet/active endpoints at 9.5 AU. Excess energies are the corresponding
  // Huebner-1992/updated AMPS channel values and are used only for product
  // kinematics, never to alter Burger's frequencies.
  static const std::vector<Reaction> table = {
    {"burger2007.photo.h2o_to_oh_h", ReactionKind::Photo, "H2O", "",
      {{"OH",1},{"H",1}}, 1.1e-7,2.0e-7,3.41*kElectronVoltJ,4.04*kElectronVoltJ,0,0},
    {"burger2007.photo.h2o_to_h2_o", ReactionKind::Photo, "H2O", "",
      {{"H2",1},{"O",1}}, 6.6e-9,1.6e-8,3.84*kElectronVoltJ,3.94*kElectronVoltJ,0,0},
    {"burger2007.photo.h2o_to_h_h_o", ReactionKind::Photo, "H2O", "",
      {{"H",2},{"O",1}}, 8.3e-9,2.1e-8,0.7*kElectronVoltJ,0.7*kElectronVoltJ,0,0},
    {"burger2007.photo.h2o_to_h2op_e", ReactionKind::Photo, "H2O", "",
      {{"H2O_PLUS",1},{"ELECTRON",1}}, 3.6e-9,9.1e-9,12.4*kElectronVoltJ,15.2*kElectronVoltJ,0,0},
    {"burger2007.photo.h2o_to_h_ohp_e", ReactionKind::Photo, "H2O", "",
      {{"H",1},{"OH_PLUS",1},{"ELECTRON",1}}, 6.1e-10,1.6e-9,18.6*kElectronVoltJ,23.2*kElectronVoltJ,0,0},
    {"burger2007.photo.h2o_to_h2_op_e", ReactionKind::Photo, "H2O", "",
      {{"H2",1},{"O_PLUS",1},{"ELECTRON",1}}, 6.4e-11,2.4e-10,36.5*kElectronVoltJ,39.8*kElectronVoltJ,0,0},
    {"burger2007.photo.h2o_to_oh_hp_e", ReactionKind::Photo, "H2O", "",
      {{"OH",1},{"H_PLUS",1},{"ELECTRON",1}}, 1.4e-10,4.5e-10,25.0*kElectronVoltJ,30.5*kElectronVoltJ,0,0},
    {"burger2007.cx.hplus_h2o_to_h_h2op", ReactionKind::ChargeExchange,
      "H2O", "H_PLUS", {{"H",1},{"H2O_PLUS",1}}, 0,0,0,0,49e-20,42e-20},
    {"burger2007.cx.oplus_h2o_to_o_h2op", ReactionKind::ChargeExchange,
      "H2O", "O_PLUS", {{"O",1},{"H2O_PLUS",1}}, 0,0,0,0,41e-20,36e-20},
    {"burger2007.cx.h2op_h2o_resonant", ReactionKind::ChargeExchange,
      "H2O", "H2O_PLUS", {{"H2O",1},{"H2O_PLUS",1}},
      0,0,0,0,10e-20,8.1e-20},
    {"burger2007.cx.ohplus_h2o_to_oh_h2op", ReactionKind::ChargeExchange,
      "H2O", "OH_PLUS", {{"OH",1},{"H2O_PLUS",1}}, 0,0,0,0,10e-20,8.1e-20}
  };
  return table;
}

const Reaction* FindReaction(const std::string& id) {
  for (const Reaction& reaction : Burger2007Reactions())
    if (reaction.id == id) return &reaction;
  return nullptr;
}

Status EvaluateReactionRate(const Reaction& reaction,
                            SolarActivity activity,
                            double totalPlasmaNumberDensityM3,
                            double plasmaReactantFraction,
                            double relativeSpeedMPerS,
                            double* ratePerS) {
  if (ratePerS == nullptr) return Status::Error("reaction-rate output is null");
  if (reaction.kind == ReactionKind::Photo) {
    *ratePerS = activity == SolarActivity::Quiet
        ? reaction.quietRatePerS : reaction.activeRatePerS;
    return (*ratePerS > 0.0 && std::isfinite(*ratePerS))
        ? Status::Ok() : Status::Error("photo-reaction table contains an invalid rate");
  }
  if (!std::isfinite(totalPlasmaNumberDensityM3) ||
      totalPlasmaNumberDensityM3 < 0.0 ||
      !std::isfinite(plasmaReactantFraction) || plasmaReactantFraction < 0.0 ||
      plasmaReactantFraction > 1.0)
    return Status::Error("charge exchange requires nonnegative coupled density and ion fraction in [0,1]");
  if (totalPlasmaNumberDensityM3 == 0.0 || plasmaReactantFraction == 0.0) {
    *ratePerS = 0.0;
    return Status::Ok();
  }
  if (!std::isfinite(relativeSpeedMPerS) || relativeSpeedMPerS < 15000.0 ||
      relativeSpeedMPerS > 26000.0)
    return Status::Error("Burger 2007 charge-exchange cross sections are qualified only for relative speeds in [15000,26000] m/s");
  const double alpha = (relativeSpeedMPerS - 15000.0) / 11000.0;
  const double crossSection = reaction.crossSectionAt15KmS_M2 + alpha *
      (reaction.crossSectionAt26KmS_M2 - reaction.crossSectionAt15KmS_M2);
  *ratePerS = totalPlasmaNumberDensityM3 * plasmaReactantFraction *
              crossSection * relativeSpeedMPerS;
  return (*ratePerS >= 0.0 && std::isfinite(*ratePerS))
      ? Status::Ok() : Status::Error("computed charge-exchange rate is invalid");
}

Status ValidateConfiguration(const Configuration& configuration,
                             const SourceDistribution& distribution) {
  Status status = RequireUniqueNonempty(configuration.simulatedSpecies,
                                        "simulated_species");
  if (!status.success) return status;
  if (configuration.flybyDensityOutput) {
    if (!IsPortableName(configuration.flybyId) ||
        configuration.flybyTrajectoryFile.empty())
      return Status::Error(
          "flyby_density_output=true requires flyby_id and "
          "flyby_trajectory_file");
  } else if (!configuration.flybyId.empty() ||
             !configuration.flybyTrajectoryFile.empty()) {
    return Status::Error(
        "flyby_density_output=false forbids flyby_id and "
        "flyby_trajectory_file");
  }
  if (configuration.flybyMeshRefinement) {
    if (!configuration.flybyDensityOutput)
      return Status::Error(
          "flyby_mesh_refinement=true requires flyby_density_output=true");
    if (!(configuration.flybyRefinementTargetCellSizeM > 0.0) ||
        !(configuration.flybyRefinementRadiusM > 0.0) ||
        !(configuration.flybyRefinementTimeMaxS >
          configuration.flybyRefinementTimeMinS))
      return Status::Error(
          "flyby_mesh_refinement=true requires positive target cell size and "
          "radius plus time_max_s greater than time_min_s");
  } else if (configuration.flybyRefinementTargetCellSizeM != 0.0 ||
             configuration.flybyRefinementRadiusM != 0.0 ||
             configuration.flybyRefinementTimeMinS != 0.0 ||
             configuration.flybyRefinementTimeMaxS != 0.0) {
    return Status::Error(
        "flyby_mesh_refinement=false forbids flyby refinement controls");
  }
  status = RequireUniqueNonempty(configuration.injectedSpecies,
                                 "injected_species");
  if (!status.success) return status;
  std::set<std::string> chemistryUnique;
  for (const std::string& species : configuration.chemistrySpecies) {
    if (species.empty() || !chemistryUnique.insert(species).second)
      return Status::Error("chemistry_species contains an empty or duplicate value");
  }
  if (!configuration.chemistrySpecies.empty() &&
      configuration.selectedReactions.empty())
    return Status::Error(
        "chemistry_species must be none when selected_reactions is none");
  for (const std::string& species : configuration.injectedSpecies) {
    if (species != "H2O" && species != "CO2")
      return Status::Error("only H2O and CO2 can be primary injected species");
    if (!Contains(configuration.simulatedSpecies, species))
      return Status::Error("injected species " + species + " is not simulated");
    if (configuration.primaryMacroparticlesPerStep.count(species) == 0)
      return Status::Error("missing primary_macroparticles_per_step." + species);
  }
  for (const std::string& species : configuration.simulatedSpecies)
    if (!Contains(configuration.injectedSpecies, species))
      return Status::Error(
          "application schema 2 parent-loss-only mode permits only primary "
          "neutral H2O/CO2 species in simulated_species; " + species +
          " is not injected");
  for (const auto& item : configuration.primaryMacroparticlesPerStep)
    if (!Contains(configuration.injectedSpecies, item.first))
      return Status::Error("primary_macroparticles_per_step." + item.first +
                           " is set for a species that is not injected");
  for (const auto& item : configuration.totalSourceRatePerS)
    if (!Contains(configuration.injectedSpecies, item.first))
      return Status::Error("total_source_rate." + item.first +
                           "_s-1 is set for a species that is not injected");
  for (const std::string& species : configuration.chemistrySpecies)
    if (!Contains(configuration.simulatedSpecies, species))
      return Status::Error("chemistry species " + species + " is not simulated");
  if (Contains(configuration.simulatedSpecies, "H3O_PLUS") ||
      Contains(configuration.chemistrySpecies, "H3O_PLUS"))
    return Status::Error("H3O_PLUS is outside the authorized reaction network");
  if (configuration.plasmaBackgroundMode == PlasmaBackgroundMode::None) {
    if (!configuration.plasmaCacheDirectory.empty() ||
        !configuration.icesInstallationDirectory.empty() ||
        !configuration.icesCase.empty())
      return Status::Error(
          "plasma_background_mode=none forbids plasma_cache_directory and "
          "all ICES settings");
  } else {
    if (configuration.plasmaCacheDirectory.empty())
      return Status::Error(
          "cache_only and cache_or_ices require plasma_cache_directory");
    if (configuration.plasmaBackgroundMode ==
            PlasmaBackgroundMode::CacheOrIces &&
        (configuration.icesInstallationDirectory.empty() ||
         configuration.icesCase.empty()))
      return Status::Error(
          "cache_or_ices requires ices_installation_directory and ices_case");
    if (configuration.plasmaBackgroundMode == PlasmaBackgroundMode::CacheOnly &&
        (!configuration.icesInstallationDirectory.empty() ||
         !configuration.icesCase.empty()))
      return Status::Error(
          "cache_only forbids ices_installation_directory and ices_case");
  }
  if (configuration.plasmaBackgroundMode == PlasmaBackgroundMode::CacheOrIces &&
      !IsSafeIcesInstallationPath(configuration.icesInstallationDirectory))
    return Status::Error(
        "ices_installation_directory contains characters unsafe for the "
        "legacy ICES command adapter; use only letters, digits, '/', '_', "
        "'-', and '.'");
  if (configuration.plasmaBackgroundMode == PlasmaBackgroundMode::CacheOrIces &&
      !IsSafeIcesCase(configuration.icesCase))
    return Status::Error(
        "ices_case must be one path component containing only letters, "
        "digits, '_', '-', and '.' (but not '.' or '..')");

  std::set<std::string> selected;
  std::set<std::string> selectedPlasmaReactants;
  std::vector<const Reaction*> selectedReactionTable;
  bool needsChargeReference = false;
  for (const std::string& id : configuration.selectedReactions) {
    if (!selected.insert(id).second)
      return Status::Error("selected_reactions contains duplicate id " + id);
    const Reaction* reaction = FindReaction(id);
    if (reaction == nullptr)
      return Status::Error("unknown or unsupported reaction id " + id);
    if (reaction->kind != ReactionKind::Photo)
      return Status::Error(
          "application schema 2 parent-loss-only mode permits only "
          "photolytic reactions; charge exchange is not enabled");
    selectedReactionTable.push_back(reaction);
    if (!Contains(configuration.chemistrySpecies, reaction->parent))
      return Status::Error("reaction parent " + reaction->parent +
                           " is absent from chemistry_species");
    if (reaction->kind == ReactionKind::ChargeExchange) {
      if (configuration.plasmaIonFraction.count(reaction->plasmaReactant) == 0)
        return Status::Error("missing plasma_ion_fraction." +
                             reaction->plasmaReactant);
      selectedPlasmaReactants.insert(reaction->plasmaReactant);
      if (configuration.plasmaIonFraction.at(reaction->plasmaReactant) > 0.0)
        needsChargeReference = true;
    }
  }
  // Products are intentionally discarded in schema 2, so no selected channel
  // may depend on a parent that would have been supplied by another reaction.
  // Every reacting parent must be a primary neutral population.
  for (const Reaction* reaction : selectedReactionTable)
    if (!Contains(configuration.injectedSpecies, reaction->parent))
      return Status::Error("reaction parent " + reaction->parent +
          " is unreachable from the selected primary injected species");

  double fractionSum = 0.0;
  for (const auto& item : configuration.plasmaIonFraction) {
    if (selectedPlasmaReactants.count(item.first) == 0)
      return Status::Error("plasma_ion_fraction." + item.first +
          " is set for a reactant absent from selected_reactions");
    if (!std::isfinite(item.second) || item.second < 0.0 || item.second > 1.0)
      return Status::Error("plasma_ion_fraction." + item.first +
          " must be finite and in [0,1]");
    fractionSum += item.second;
  }
  if (fractionSum > 1.0 + 1.0e-12)
    return Status::Error("plasma ion fractions sum to more than one");
  if (needsChargeReference &&
      (!(configuration.chemistryWeightReferenceNumberDensityM3 > 0.0) ||
       !(configuration.chemistryWeightReferenceRelativeSpeedMPerS > 0.0)))
    return Status::Error("charge-exchange weight planning requires both chemistry_weight_reference_* values");
  if (needsChargeReference &&
      (configuration.chemistryWeightReferenceRelativeSpeedMPerS < 15000.0 ||
       configuration.chemistryWeightReferenceRelativeSpeedMPerS > 26000.0))
    return Status::Error(
        "chemistry_weight_reference_relative_speed_m_s-1 must be in the "
        "Burger 2007 bracket [15000,26000]");
  if (!needsChargeReference &&
      (configuration.chemistryWeightReferenceNumberDensityM3 != 0.0 ||
       configuration.chemistryWeightReferenceRelativeSpeedMPerS != 0.0))
    return Status::Error(
        "chemistry_weight_reference_* values are forbidden without a selected "
        "charge-exchange channel having nonzero plasma ion fraction");

  for (const std::string& species : configuration.simulatedSpecies) {
    if (configuration.characteristicSpeedMPerS.count(species) == 0)
      return Status::Error("missing characteristic_speed." + species + "_m_s-1");
    if (configuration.courantNumber.count(species) == 0)
      return Status::Error("missing courant." + species);
  }
  for (const auto& item : configuration.characteristicSpeedMPerS)
    if (!Contains(configuration.simulatedSpecies, item.first))
      return Status::Error("characteristic_speed." + item.first +
                           "_m_s-1 is set for a species that is not simulated");
  for (const auto& item : configuration.courantNumber)
    if (!Contains(configuration.simulatedSpecies, item.first))
      return Status::Error("courant." + item.first +
                           " is set for a species that is not simulated");
  if (configuration.populationControl != PopulationControl::Off) {
    if (configuration.minParticlesPerCell <= 0 ||
        configuration.maxParticlesPerCell <= configuration.minParticlesPerCell)
      return Status::Error("enabled splitting/merging requires 0 < min_particles_per_cell < max_particles_per_cell");
    if (configuration.populationControl == PopulationControl::VelocityShift &&
        (!(configuration.splitPositionShiftCellDiagonalFraction > 0.0) ||
         !(configuration.splitVelocityShiftThermalSpeedFraction > 0.0)))
      return Status::Error(
          "velocity_shift requires positive cell-diagonal and thermal-speed "
          "shift fractions");
    if (configuration.populationControl == PopulationControl::Scatter &&
        (configuration.splitPositionShiftCellDiagonalFraction != 0.0 ||
         configuration.splitVelocityShiftThermalSpeedFraction != 0.0 ||
         configuration.nonuniformSplitPositionShift))
      return Status::Error(
          "scatter forbids velocity-shift fractions and "
          "nonuniform_split_position_shift");
  } else if (configuration.minParticlesPerCell != 0 ||
             configuration.maxParticlesPerCell != 0 ||
             configuration.splitPositionShiftCellDiagonalFraction != 0.0 ||
             configuration.splitVelocityShiftThermalSpeedFraction != 0.0 ||
             configuration.nonuniformSplitPositionShift) {
    return Status::Error(
        "particle_splitting_merging=off forbids all split/merge parameter keys");
  }

  std::set<std::string> distributionSpecies;
  for (const Source& source : distribution.sources)
    for (const auto& strength : source.strength)
      distributionSpecies.insert(strength.first);
  for (const std::string& species : configuration.injectedSpecies)
    if (distributionSpecies.count(species) == 0)
      return Status::Error("source file has no strengths for injected species " + species);
  if (distribution.mode == SourceMode::Shape) {
    for (const std::string& species : configuration.injectedSpecies)
      if (configuration.totalSourceRatePerS.count(species) == 0)
        return Status::Error("shape mode requires total_source_rate." + species + "_s-1");
  } else {
    if (!configuration.totalSourceRatePerS.empty())
      return Status::Error("absolute source mode forbids all total_source_rate.* keys");
  }
  return Status::Ok();
}

Status BuildNumericalPlan(const Configuration& configuration,
                          const SourceDistribution& distribution,
                          double minimumAllocatedCellSizeM,
                          NumericalPlan* output) {
  if (output == nullptr) return Status::Error("numerical-plan output is null");
  Status status = ValidateConfiguration(configuration, distribution);
  if (!status.success) return status;
  if (!std::isfinite(minimumAllocatedCellSizeM) ||
      minimumAllocatedCellSizeM <= 0.0)
    return Status::Error("minimum allocated cell size must be finite and positive");
  NumericalPlan plan;
  plan.globalTimeStepS = std::numeric_limits<double>::infinity();
  for (const std::string& species : configuration.simulatedSpecies) {
    const double dt = configuration.courantNumber.at(species) *
        minimumAllocatedCellSizeM /
        configuration.characteristicSpeedMPerS.at(species);
    if (dt < plan.globalTimeStepS) {
      plan.globalTimeStepS = dt;
      plan.limitingSpecies = species;
    }
  }
  if (!(plan.globalTimeStepS > 0.0) || !std::isfinite(plan.globalTimeStepS))
    return Status::Error("global Courant time step is invalid");

  for (const std::string& species : configuration.injectedSpecies) {
    double rate = 0.0;
    if (distribution.mode == SourceMode::Shape) {
      rate = configuration.totalSourceRatePerS.at(species);
    } else {
      for (const Source& source : distribution.sources)
        rate += source.strength.at(species);
    }
    const double count = static_cast<double>(
        configuration.primaryMacroparticlesPerStep.at(species));
    plan.physicalSourceRatePerS[species] = rate;
    plan.expectedMacroparticlesPerStep[species] = count;
    plan.statisticalWeight[species] = rate * plan.globalTimeStepS / count;
  }

  // Schema 2 records every selected photo frequency for the runtime survival
  // calculation but intentionally plans no daughter population or weight.
  for (const std::string& id : configuration.selectedReactions) {
    const Reaction* reaction = FindReaction(id);
    double rate = 0.0;
    status = EvaluateReactionRate(*reaction, configuration.solarActivity,
                                  0.0, 0.0, 0.0, &rate);
    if (!status.success) return Status::Error(id + ": " + status.message);
    plan.reactionRatePerS[id] = rate;
  }
  *output = std::move(plan);
  return Status::Ok();
}

Vec3 PointJetNormal(const Source& source) {
  const double lat = source.latitudeDeg * kPi / 180.0;
  const double lon = -source.westLongitudeDeg * kPi / 180.0;
  return Vec3{std::cos(lat) * std::cos(lon),
              std::cos(lat) * std::sin(lon), std::sin(lat)};
}

Vec3 PointJetPositionM(const Source& source, double radiusM) {
  const Vec3 normal = PointJetNormal(source);
  return Vec3{radiusM * normal.x, radiusM * normal.y, radiusM * normal.z};
}

Vec3 PointJetAxis(const Source& source) {
  const double lat = source.latitudeDeg * kPi / 180.0;
  const double lon = -source.westLongitudeDeg * kPi / 180.0;
  const double tilt = source.tiltDeg * kPi / 180.0;
  const double azimuth = source.azimuthDeg * kPi / 180.0;
  const Vec3 normal = PointJetNormal(source);
  const Vec3 north{-std::sin(lat) * std::cos(lon),
                   -std::sin(lat) * std::sin(lon), std::cos(lat)};
  const Vec3 east{-std::sin(lon), std::cos(lon), 0.0};
  return Normalize(Vec3{
      std::cos(tilt) * normal.x + std::sin(tilt) *
          (std::cos(azimuth) * north.x + std::sin(azimuth) * east.x),
      std::cos(tilt) * normal.y + std::sin(tilt) *
          (std::cos(azimuth) * north.y + std::sin(azimuth) * east.y),
      std::cos(tilt) * normal.z + std::sin(tilt) *
          (std::cos(azimuth) * north.z + std::sin(azimuth) * east.z)});
}

Status StripeNodeWeights(const Source& source, std::vector<double>* weights) {
  if (weights == nullptr) return Status::Error("stripe-weight output is null");
  if (source.geometry != SourceGeometry::PolylineStripe ||
      source.verticesM.size() < 2)
    return Status::Error("stripe needs at least two vertices");
  weights->assign(source.verticesM.size(), 0.0);
  double total = 0.0;
  for (std::size_t i = 1; i < source.verticesM.size(); ++i) {
    const Vec3 delta{source.verticesM[i].x - source.verticesM[i-1].x,
                     source.verticesM[i].y - source.verticesM[i-1].y,
                     source.verticesM[i].z - source.verticesM[i-1].z};
    const double length = Norm(delta);
    if (!(length > 0.0) || !std::isfinite(length))
      return Status::Error("stripe contains a zero-length or invalid segment");
    (*weights)[i-1] += 0.5 * length;
    (*weights)[i] += 0.5 * length;
    total += length;
  }
  for (double& weight : *weights) weight /= total;
  return Status::Ok();
}

namespace {

// Return J_3(a)=integral_0^infinity t^3 exp[-(t-a)^2] dt.  The recurrence
// follows integration by parts.  Long-double intermediates retain accuracy
// down to a=-2, the documented recurrence boundary also used by the
// standalone model.  The current tenishev2014_hybrid_8 H2O/CO2 parameters
// have a>=-1.71 throughout their outward hemispheres.
long double ShiftedMoment3Total(long double a) {
  const long double sqrtPi = std::sqrt(std::acos(-1.0L));
  const long double j0 = 0.5L * sqrtPi * std::erfc(-a);
  const long double e0 = std::exp(-a*a);
  const long double j1 = a*j0 + 0.5L*e0;
  const long double j2 = a*j1 + 0.5L*j0;
  return a*j2 + j1;
}

long double ShiftedMoment3Partial(long double x, long double a) {
  const long double sqrtPi = std::sqrt(std::acos(-1.0L));
  const long double ex = std::exp(-(x-a)*(x-a));
  const long double e0 = std::exp(-a*a);
  const long double j0 = 0.5L*sqrtPi*(std::erf(x-a)+std::erf(a));
  const long double j1 = a*j0 - 0.5L*ex + 0.5L*e0;
  const long double j2 = a*j1 - 0.5L*x*ex + 0.5L*j0;
  return a*j2 - 0.5L*x*x*ex + j1;
}

}  // namespace

Status Tenishev2010FluxMoment3(double directionCosine, double speedRatio,
                              double* moment) {
  if (moment == nullptr) return Status::Error("Tenishev moment output is null");
  if (!std::isfinite(directionCosine) || directionCosine < -1.0 ||
      directionCosine > 1.0 || !std::isfinite(speedRatio) ||
      speedRatio < 0.0)
    return Status::Error("Tenishev moment requires c in [-1,1] and finite s>=0");
  const long double c = directionCosine;
  const long double s = speedRatio;
  const long double a = s*c;
  if (a < -2.0L)
    return Status::Error(
        "Tenishev recurrence is qualified only for speed_ratio*direction_cosine>=-2");
  const long double value =
      std::exp(-s*s*(1.0L-c*c))*ShiftedMoment3Total(a);
  *moment = static_cast<double>(value);
  if (!std::isfinite(*moment) || !(*moment > 0.0))
    return Status::Error("Tenishev flux moment is non-finite or non-positive");
  return Status::Ok();
}

Status Tenishev2010SpeedQuantile(double directionCosine, double speedRatio,
                                double probability,
                                double* speedOverThermal) {
  if (speedOverThermal == nullptr)
    return Status::Error("Tenishev speed-quantile output is null");
  if (!std::isfinite(directionCosine) || directionCosine < -1.0 ||
      directionCosine > 1.0 || !std::isfinite(speedRatio) ||
      speedRatio < 0.0 || !std::isfinite(probability) ||
      !(probability > 0.0 && probability < 1.0))
    return Status::Error(
        "Tenishev quantile requires c in [-1,1], finite s>=0, and p in (0,1)");
  const long double a = static_cast<long double>(speedRatio) *
                        static_cast<long double>(directionCosine);
  if (a < -2.0L)
    return Status::Error(
        "Tenishev recurrence is qualified only for speed_ratio*direction_cosine>=-2");
  const long double total = ShiftedMoment3Total(a);
  if (!(total > 0.0L) || !std::isfinite(total))
    return Status::Error("Tenishev conditional speed normalization is invalid");

  // Eight thermal speeds beyond a positive drift makes the omitted Gaussian
  // tail negligible. Expand defensively if a probability very near one still
  // lies beyond that bracket, then use monotone bisection on the analytic CDF.
  long double lower = 0.0L;
  long double upper = std::max(8.0L, a+8.0L);
  const long double target = static_cast<long double>(probability)*total;
  for (int expansion = 0;
       expansion < 8 && ShiftedMoment3Partial(upper, a) < target;
       ++expansion) upper *= 2.0L;
  if (ShiftedMoment3Partial(upper, a) < target)
    return Status::Error("Tenishev speed quantile could not bracket its CDF");
  for (int iteration = 0; iteration < 96; ++iteration) {
    const long double middle = 0.5L*(lower+upper);
    if (ShiftedMoment3Partial(middle, a) < target) lower = middle;
    else upper = middle;
  }
  *speedOverThermal = static_cast<double>(0.5L*(lower+upper));
  return std::isfinite(*speedOverThermal) && *speedOverThermal > 0.0
      ? Status::Ok()
      : Status::Error("Tenishev speed quantile is invalid");
}

}  // namespace ForwardModel
}  // namespace Enceladus
