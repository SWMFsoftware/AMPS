#include "application_input.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <utility>

namespace SEP3D {
namespace RuntimeModel {
namespace {

namespace fs = std::filesystem;

constexpr std::size_t kMaximumIncludeDepth = 64;

struct LogicalLine {
  std::string file;
  std::size_t line = 0;
  std::string text;
};

std::string Trim(const std::string& text) {
  const auto first = std::find_if_not(text.begin(), text.end(),
      [](unsigned char c) { return std::isspace(c) != 0; });
  if (first == text.end()) return std::string();
  const auto last = std::find_if_not(text.rbegin(), text.rend(),
      [](unsigned char c) { return std::isspace(c) != 0; }).base();
  return std::string(first, last);
}

std::string Lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

Core::Status Invalid(const LogicalLine& line, const std::string& reason) {
  std::ostringstream message;
  message << line.file << ':' << line.line
          << ": srcSEP3D input error: " << reason;
  if (!line.text.empty()) message << "\n  line: " << line.text;
  return Core::Status(Core::StatusCode::InvalidInput, message.str());
}

Core::Status InvalidFile(const std::string& file, const std::string& reason) {
  return Core::Status(Core::StatusCode::InvalidInput,
      file + ": srcSEP3D input error: " + reason);
}

// Parse the deliberately small include grammar.  Quoted and angle-bracket
// paths permit whitespace; an unquoted path consumes the complete trimmed
// remainder.  Rejecting trailing text after a quoted path catches misspelled
// comments (comments begin with !, not #) instead of guessing user intent.
Core::Status IncludeTarget(const LogicalLine& line, std::string* target) {
  const std::string trimmed = Trim(line.text);
  const std::string directive = "#include";
  if (trimmed.size() <= directive.size() ||
      std::isspace(static_cast<unsigned char>(trimmed[directive.size()])) == 0) {
    return Invalid(line, "#include requires a file name");
  }
  std::string rest = Trim(trimmed.substr(directive.size()));
  if (rest.empty()) return Invalid(line, "#include requires a file name");

  if (rest.front() == '"' || rest.front() == '<') {
    const char close = rest.front() == '"' ? '"' : '>';
    const std::size_t end = rest.find(close, 1);
    if (end == std::string::npos)
      return Invalid(line, "#include file name is missing its closing delimiter");
    if (!Trim(rest.substr(end + 1)).empty())
      return Invalid(line, "unrecognized text follows the #include file name");
    rest = rest.substr(1, end - 1);
  }
  if (rest.empty()) return Invalid(line, "#include file name is empty");
  *target = rest;
  return Core::Status::OK();
}

Core::Status ExpandFile(const fs::path& requested,
                        std::vector<std::string>* active,
                        std::vector<std::string>* files,
                        std::vector<LogicalLine>* lines,
                        const LogicalLine* includeSite) {
  if (active->size() >= kMaximumIncludeDepth) {
    return includeSite == nullptr
        ? InvalidFile(requested.string(), "maximum #include depth (64) exceeded")
        : Invalid(*includeSite, "maximum #include depth (64) exceeded");
  }

  std::error_code pathError;
  fs::path absolute = fs::absolute(requested, pathError);
  if (pathError) {
    const std::string reason = "cannot resolve input path '" +
        requested.string() + "': " + pathError.message();
    return includeSite == nullptr ? InvalidFile(requested.string(), reason)
                                  : Invalid(*includeSite, reason);
  }
  absolute = absolute.lexically_normal();
  fs::path canonical = fs::weakly_canonical(absolute, pathError);
  if (pathError) canonical = absolute;
  const std::string identity = canonical.string();
  if (std::find(active->begin(), active->end(), identity) != active->end()) {
    const std::string reason = "recursive #include cycle reaches '" + identity + "'";
    return includeSite == nullptr ? InvalidFile(identity, reason)
                                  : Invalid(*includeSite, reason);
  }

  std::ifstream stream(canonical);
  if (!stream) {
    const std::string reason = "cannot open input file '" + identity + "'";
    return includeSite == nullptr ? InvalidFile(identity, reason)
                                  : Invalid(*includeSite, reason);
  }

  active->push_back(identity);
  if (std::find(files->begin(), files->end(), identity) == files->end())
    files->push_back(identity);

  std::string accumulated;
  std::size_t logicalStart = 0;
  std::string physical;
  std::size_t physicalLine = 0;
  while (std::getline(stream, physical)) {
    ++physicalLine;
    if (!physical.empty() && physical.back() == '\r') physical.pop_back();

    // A comment is removed before inspecting the continuation marker.  Thus
    // ``value = 1 ! \\`` does not continue, while ``value = \\ ! note`` does.
    const std::size_t comment = physical.find('!');
    std::string fragment = comment == std::string::npos
        ? physical : physical.substr(0, comment);
    fragment = Trim(fragment);
    const bool continued = !fragment.empty() && fragment.back() == '\\';
    if (continued) fragment = Trim(fragment.substr(0, fragment.size() - 1));

    if (logicalStart == 0) logicalStart = physicalLine;
    if (!fragment.empty()) {
      if (!accumulated.empty()) accumulated.push_back(' ');
      accumulated += fragment;
    }
    if (continued) continue;

    LogicalLine line{identity, logicalStart, accumulated};
    accumulated.clear();
    logicalStart = 0;
    const std::string trimmed = Trim(line.text);
    if (trimmed.rfind("#include", 0) == 0) {
      std::string target;
      Core::Status status = IncludeTarget(line, &target);
      if (!status.ok()) { active->pop_back(); return status; }
      fs::path child(target);
      if (child.is_relative()) child = canonical.parent_path() / child;
      status = ExpandFile(child, active, files, lines, &line);
      if (!status.ok()) { active->pop_back(); return status; }
    } else if (!trimmed.empty()) {
      line.text = trimmed;
      lines->push_back(std::move(line));
    }
  }
  if (!stream.eof()) {
    active->pop_back();
    return InvalidFile(identity, "I/O failure while reading the input file");
  }
  if (logicalStart != 0) {
    LogicalLine dangling{identity, logicalStart, accumulated + " \\"};
    active->pop_back();
    return Invalid(dangling, "line continuation reaches end of file");
  }
  active->pop_back();
  return Core::Status::OK();
}

bool ParseUnsigned(const std::string& text, std::uint64_t* value) {
  if (text.empty() || text.front() == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      parsed > std::numeric_limits<std::uint64_t>::max()) return false;
  *value = static_cast<std::uint64_t>(parsed);
  return true;
}

bool ParseFiniteDouble(const std::string& text, double* value) {
  if (text.empty() || value == nullptr) return false;
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

}  // namespace

Core::Status ParseSep3dApplicationInput(
    const std::string& path, Sep3dApplicationInput* result) {
  if (result == nullptr)
    return Core::Status(Core::StatusCode::InvalidInput,
                        "srcSEP3D application-input output is null");
  if (path.empty())
    return Core::Status(Core::StatusCode::InvalidInput,
                        "srcSEP3D application-input path is empty");

  std::vector<std::string> active;
  std::vector<std::string> files;
  std::vector<LogicalLine> lines;
  Core::Status status = ExpandFile(path, &active, &files, &lines, nullptr);
  if (!status.ok()) return status;

  Sep3dApplicationInput candidate;
  candidate.rootFile = files.empty() ? path : files.front();
  candidate.expandedFiles = files;
  bool insideSection = false;
  bool insideSep3d = false;
  bool insideReducedShock = false;
  bool insideParticleInjection = false;
  bool sawSep3d = false;
  bool sawMaximumTimeSteps = false;
  bool sawValue = false;
  bool sawShockModel = false;
  bool sawBackgroundModel = false;
  bool sawSourceModel = false;
  bool sawMaximumSpeed = false;
  bool sawMargin = false;
  bool sawNormalizationRadius = false;
  bool sawReducedShock = false;
  bool sawParticleInjection = false;
  bool sawWeightingModel = false;
  bool sawPowerModel = false;
  bool sawMinimumEnergy = false;
  bool sawMaximumEnergy = false;
  bool sawFixedPowerIndex = false;
  bool sawMaximumInjectionEvents = false;
  std::vector<std::string> reducedKeys;
  LogicalLine sectionStart;
  LogicalLine subsectionStart;

  for (const LogicalLine& line : lines) {
    const std::string lower = Lower(line.text);
    const std::string begin = "#section begin:";
    if (lower.rfind(begin, 0) == 0) {
      if (insideSection)
        return Invalid(line, "nested #section begin is not allowed");
      const std::string name = Lower(Trim(line.text.substr(begin.size())));
      if (name.empty()) return Invalid(line, "#section begin is missing a section name");
      if (name.find_first_of(" \t") != std::string::npos)
        return Invalid(line, "unrecognized text follows the section name");
      insideSection = true;
      insideSep3d = name == "sep3d";
      sectionStart = line;
      if (insideSep3d) {
        if (sawSep3d)
          return Invalid(line, "duplicate sep3d section is not allowed");
        sawSep3d = true;
      }
      continue;
    }
    if (lower.rfind("#section begin", 0) == 0)
      return Invalid(line, "expected '#section begin: NAME'");
    if (lower.rfind("#section end", 0) == 0) {
      if (lower != "#section end")
        return Invalid(line, "unrecognized text follows '#section end'");
      if (!insideSection)
        return Invalid(line, "#section end has no matching #section begin");
      if (insideReducedShock || insideParticleInjection)
        return Invalid(line,
            "srcSEP3D subsection is missing '#subsection end'");
      insideSection = false;
      insideSep3d = false;
      continue;
    }
    if (lower.rfind("#section", 0) == 0)
      return Invalid(line, "unrecognized #section directive");
    if (!insideSep3d) continue;

    const std::string subsectionBegin = "#subsection begin:";
    if (lower.rfind(subsectionBegin, 0) == 0) {
      if (insideReducedShock || insideParticleInjection)
        return Invalid(line, "nested #subsection begin is not allowed");
      const std::string name = Lower(Trim(
          line.text.substr(subsectionBegin.size())));
      if (name != "reduced-shock-surface" &&
          name != "shock-particle-injection")
        return Invalid(line, "unsupported srcSEP3D subsection '" + name + "'");
      if (name == "reduced-shock-surface") {
        if (sawReducedShock)
          return Invalid(line,
              "duplicate reduced-shock-surface subsection is not allowed");
        insideReducedShock = true;
        sawReducedShock = true;
        candidate.reducedShockAssetDirectory =
            fs::path(line.file).parent_path().string();
      } else {
        if (sawParticleInjection)
          return Invalid(line,
              "duplicate shock-particle-injection subsection is not allowed");
        insideParticleInjection = true;
        sawParticleInjection = true;
      }
      subsectionStart = line;
      continue;
    }
    if (lower.rfind("#subsection begin", 0) == 0)
      return Invalid(line,
          "expected '#subsection begin: reduced-shock-surface' or "
          "'#subsection begin: shock-particle-injection'");
    if (lower.rfind("#subsection end", 0) == 0) {
      if (lower != "#subsection end")
        return Invalid(line,
            "unrecognized text follows '#subsection end'");
      if (!insideReducedShock && !insideParticleInjection)
        return Invalid(line,
            "#subsection end has no matching #subsection begin");
      if (insideReducedShock && candidate.reducedShockConfiguration.empty())
        return Invalid(line,
            "reduced-shock-surface subsection contains no model parameters");
      if (insideParticleInjection) {
        std::vector<std::string> missing;
        if (!sawWeightingModel) missing.push_back("statistical_weight_model");
        if (!sawPowerModel) missing.push_back("phase_space_power_model");
        if (!sawMinimumEnergy) missing.push_back("minimum_energy_j");
        if (!sawMaximumEnergy) missing.push_back("maximum_energy_j");
        if (!sawMaximumInjectionEvents)
          missing.push_back("maximum_events_per_species_per_step");
        if (candidate.momentumPowerLawModel == "constant" &&
            !sawFixedPowerIndex)
          missing.push_back("phase_space_power_index");
        if (!missing.empty()) {
          std::ostringstream reason;
          reason << "shock-particle-injection subsection is missing required ";
          for (std::size_t index = 0; index < missing.size(); ++index) {
            if (index != 0) reason << ", ";
            reason << missing[index];
          }
          return Invalid(line, reason.str());
        }
        if (candidate.momentumPowerLawModel == "compression-ratio" &&
            sawFixedPowerIndex)
          return Invalid(line, "phase_space_power_index is inactive and must "
              "be omitted for phase_space_power_model=compression-ratio");
      }
      insideReducedShock = false;
      insideParticleInjection = false;
      continue;
    }
    if (lower.rfind("#subsection", 0) == 0)
      return Invalid(line, "unrecognized #subsection directive");

    const std::size_t equal = line.text.find('=');
    if (equal == std::string::npos)
      return Invalid(line, "expected 'NAME = VALUE'");
    const std::string key = Lower(Trim(line.text.substr(0, equal)));
    const std::string value = Trim(line.text.substr(equal + 1));
    if (key.empty() || value.empty())
      return Invalid(line, "input assignment is missing its name or value");
    if (insideReducedShock) {
      if (std::find(reducedKeys.begin(), reducedKeys.end(), key) !=
          reducedKeys.end())
        return Invalid(line, "reduced shock parameter '" + key +
            "' is specified more than once");
      reducedKeys.push_back(key);
      candidate.reducedShockConfiguration += key + "=" + value + "\n";
      continue;
    }
    if (insideParticleInjection) {
      auto duplicate = [&](bool seen) {
        return seen ? Invalid(line, key + " is specified more than once")
                    : Core::Status::OK();
      };
      if (key == "statistical_weight_model") {
        Core::Status unique = duplicate(sawWeightingModel);
        if (!unique.ok()) return unique;
        candidate.particleWeightingModel = Lower(value);
        if (candidate.particleWeightingModel !=
                "constant-statistical-weight" &&
            candidate.particleWeightingModel !=
                "log-uniform-momentum-importance")
          return Invalid(line, "unsupported statistical_weight_model '" +
              candidate.particleWeightingModel + "'");
        sawWeightingModel = true;
      } else if (key == "phase_space_power_model") {
        Core::Status unique = duplicate(sawPowerModel);
        if (!unique.ok()) return unique;
        candidate.momentumPowerLawModel = Lower(value);
        if (candidate.momentumPowerLawModel != "constant" &&
            candidate.momentumPowerLawModel != "compression-ratio")
          return Invalid(line, "unsupported phase_space_power_model '" +
              candidate.momentumPowerLawModel + "'");
        sawPowerModel = true;
      } else if (key == "minimum_energy_j") {
        Core::Status unique = duplicate(sawMinimumEnergy);
        if (!unique.ok()) return unique;
        if (!ParseFiniteDouble(value, &candidate.minimumInjectionEnergyJ) ||
            candidate.minimumInjectionEnergyJ <= 0.0)
          return Invalid(line, "minimum_energy_j must be finite and positive");
        sawMinimumEnergy = true;
      } else if (key == "maximum_energy_j") {
        Core::Status unique = duplicate(sawMaximumEnergy);
        if (!unique.ok()) return unique;
        if (!ParseFiniteDouble(value, &candidate.maximumInjectionEnergyJ) ||
            candidate.maximumInjectionEnergyJ <= 0.0)
          return Invalid(line, "maximum_energy_j must be finite and positive");
        sawMaximumEnergy = true;
      } else if (key == "phase_space_power_index") {
        Core::Status unique = duplicate(sawFixedPowerIndex);
        if (!unique.ok()) return unique;
        if (!ParseFiniteDouble(value,
                &candidate.fixedPhaseSpacePowerIndex) ||
            candidate.fixedPhaseSpacePowerIndex <= 2.0)
          return Invalid(line, "phase_space_power_index must be finite and "
              "greater than two for f(p) proportional to p^(-q)");
        sawFixedPowerIndex = true;
      } else if (key == "maximum_events_per_species_per_step") {
        Core::Status unique = duplicate(sawMaximumInjectionEvents);
        if (!unique.ok()) return unique;
        if (!ParseUnsigned(value,
                &candidate.maximumInjectionEventsPerSpeciesPerStep) ||
            candidate.maximumInjectionEventsPerSpeciesPerStep == 0)
          return Invalid(line, "maximum_events_per_species_per_step must be "
              "a positive unsigned integer");
        sawMaximumInjectionEvents = true;
      } else {
        return Invalid(line, "unrecognized shock-particle-injection setting '" +
            key + "'");
      }
      continue;
    }

    auto duplicate = [&](bool seen) {
      return seen ? Invalid(line, key + " is specified more than once")
                  : Core::Status::OK();
    };
    if (key == "maximum_time_steps") {
      Core::Status unique = duplicate(sawMaximumTimeSteps);
      if (!unique.ok()) return unique;
      if (!ParseUnsigned(value, &candidate.maximumTimeSteps) ||
          candidate.maximumTimeSteps == 0)
        return Invalid(line,
            "maximum_time_steps must be a positive unsigned integer");
      sawMaximumTimeSteps = true;
    } else if (key == "particles_per_iteration") {
      Core::Status unique = duplicate(sawValue);
      if (!unique.ok()) return unique;
      if (!ParseUnsigned(value, &candidate.particlesPerIteration) ||
          candidate.particlesPerIteration == 0)
        return Invalid(line,
            "particles_per_iteration must be a positive unsigned integer");
      sawValue = true;
      candidate.valueFile = line.file;
      candidate.valueLine = line.line;
    } else if (key == "shock_model") {
      Core::Status unique = duplicate(sawShockModel);
      if (!unique.ok()) return unique;
      candidate.shockModel = Lower(value);
      if (candidate.shockModel != "reduced-shock-surface")
        return Invalid(line,
            "unsupported shock_model '" + candidate.shockModel + "'");
      sawShockModel = true;
    } else if (key == "background_plasma_model") {
      Core::Status unique = duplicate(sawBackgroundModel);
      if (!unique.ok()) return unique;
      candidate.backgroundPlasmaModel = Lower(value);
      if (candidate.backgroundPlasmaModel != "corona-swcme-ambient")
        return Invalid(line, "unsupported background_plasma_model '" +
            candidate.backgroundPlasmaModel + "'");
      sawBackgroundModel = true;
    } else if (key == "source_model") {
      Core::Status unique = duplicate(sawSourceModel);
      if (!unique.ok()) return unique;
      candidate.sourceModel = Lower(value);
      if (candidate.sourceModel != "accepted-shock-incident-flux")
        return Invalid(line,
            "unsupported source_model '" + candidate.sourceModel + "'");
      sawSourceModel = true;
    } else if (key == "maximum_particle_speed_m_s") {
      Core::Status unique = duplicate(sawMaximumSpeed);
      if (!unique.ok()) return unique;
      if (!ParseFiniteDouble(value, &candidate.maximumParticleSpeedMPerS) ||
          candidate.maximumParticleSpeedMPerS <= 0.0 ||
          candidate.maximumParticleSpeedMPerS > Core::Const::c)
        return Invalid(line,
            "maximum_particle_speed_m_s must be in (0,c]");
      sawMaximumSpeed = true;
    } else if (key == "time_step_margin_factor") {
      Core::Status unique = duplicate(sawMargin);
      if (!unique.ok()) return unique;
      if (!ParseFiniteDouble(value, &candidate.timeStepMarginFactor) ||
          candidate.timeStepMarginFactor <= 0.0 ||
          candidate.timeStepMarginFactor > 1.0)
        return Invalid(line,
            "time_step_margin_factor must be in (0,1]");
      sawMargin = true;
    } else if (key == "source_normalization_radius_m") {
      Core::Status unique = duplicate(sawNormalizationRadius);
      if (!unique.ok()) return unique;
      if (!ParseFiniteDouble(value,
              &candidate.sourceNormalizationRadiusM) ||
          candidate.sourceNormalizationRadiusM <= 0.0)
        return Invalid(line,
            "source_normalization_radius_m must be finite and positive");
      sawNormalizationRadius = true;
    } else {
      return Invalid(line, "unrecognized srcSEP3D setting '" + key + "'");
    }
  }

  if (insideReducedShock || insideParticleInjection)
    return Invalid(subsectionStart,
        "subsection reaches end of expanded input without '#subsection end'");
  if (insideSection)
    return Invalid(sectionStart, "section reaches end of expanded input without '#section end'");
  if (!sawSep3d)
    return InvalidFile(candidate.rootFile,
        "missing required '#section begin: sep3d' section");
  std::vector<std::string> missing;
  if (!sawMaximumTimeSteps) missing.push_back("maximum_time_steps");
  if (!sawValue) missing.push_back("particles_per_iteration");
  if (!sawShockModel) missing.push_back("shock_model");
  if (!sawBackgroundModel) missing.push_back("background_plasma_model");
  if (!sawSourceModel) missing.push_back("source_model");
  if (!sawMaximumSpeed) missing.push_back("maximum_particle_speed_m_s");
  if (!sawMargin) missing.push_back("time_step_margin_factor");
  if (!sawNormalizationRadius)
    missing.push_back("source_normalization_radius_m");
  if (!sawReducedShock) missing.push_back("reduced-shock-surface subsection");
  if (!sawParticleInjection)
    missing.push_back("shock-particle-injection subsection");
  if (!missing.empty()) {
    std::ostringstream reason;
    reason << "sep3d section is missing required ";
    for (std::size_t index=0;index<missing.size();++index) {
      if (index!=0) reason << ", ";
      reason << missing[index];
    }
    return InvalidFile(candidate.rootFile,reason.str());
  }
  if (candidate.maximumInjectionEnergyJ <=
      candidate.minimumInjectionEnergyJ)
    return InvalidFile(candidate.rootFile,
        "maximum_energy_j must be greater than minimum_energy_j");
  *result = std::move(candidate);
  return Core::Status::OK();
}

std::string Sep3dApplicationInputSummary(
    const Sep3dApplicationInput& input) {
  std::ostringstream summary;
  summary << "[srcSEP3D] input summary\n"
          << "  root_file=" << input.rootFile << '\n'
          << "  expanded_file_count=" << input.expandedFiles.size() << '\n';
  for (const std::string& file : input.expandedFiles)
    summary << "  expanded_file=" << file << '\n';
  summary << "  maximum_time_steps=" << input.maximumTimeSteps << '\n'
          << "  particles_per_iteration=" << input.particlesPerIteration
          << " (per compiled species)\n"
          << "  shock_model=" << input.shockModel << '\n'
          << "  background_plasma_model=" << input.backgroundPlasmaModel
          << '\n'
          << "  source_model=" << input.sourceModel << '\n'
          << "  maximum_particle_speed_m_s="
          << input.maximumParticleSpeedMPerS << '\n'
          << "  time_step_margin_factor=" << input.timeStepMarginFactor
          << '\n'
          << "  source_normalization_radius_m="
          << input.sourceNormalizationRadiusM << '\n'
          << "  statistical_weight_model="
          << input.particleWeightingModel << '\n'
          << "  phase_space_power_model="
          << input.momentumPowerLawModel << '\n'
          << "  minimum_injection_energy_j="
          << input.minimumInjectionEnergyJ << '\n'
          << "  maximum_injection_energy_j="
          << input.maximumInjectionEnergyJ << '\n'
          << "  phase_space_power_index="
          << (input.momentumPowerLawModel == "constant"
                  ? std::to_string(input.fixedPhaseSpacePowerIndex)
                  : std::string("derived from local compression")) << '\n'
          << "  maximum_events_per_species_per_step="
          << input.maximumInjectionEventsPerSpeciesPerStep << '\n'
          << "  reduced_shock_asset_directory="
          << input.reducedShockAssetDirectory << '\n'
          << "  value_source=" << input.valueFile << ':' << input.valueLine
          << '\n';
  return summary.str();
}

}  // namespace RuntimeModel
}  // namespace SEP3D
