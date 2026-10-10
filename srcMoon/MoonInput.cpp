#include "MoonInput.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <utility>

// Runtime parser for the application-owned portion of a shared AMPS input
// file.  Parsing happens before PIC/MPI initialization, performs no global
// mutation, and commits the completed Configuration only after every syntax,
// value, and path check succeeds.
namespace Moon {
namespace Runtime {
namespace {

namespace fs = std::filesystem;

// The executable has one immutable application configuration.  These objects
// are written only by InstallConfiguration(), before amps_init(), and are read
// thereafter by the production initialization path.
bool gInstalled = false;
Configuration gConfiguration;

struct Setting {
  // Preserve the unparsed value and its source line so conversion errors can
  // identify the exact input rather than reporting a generic startup failure.
  std::string value;
  std::size_t line = 0;
};

// Input syntax is ASCII-oriented, but the unsigned-char conversion avoids the
// undefined behavior of passing a negative signed char to std::isspace().
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

// Quotes are optional around scalar/path values.  Deliberately do not perform
// shell expansion, environment substitution, or escape processing: the run
// receipt must identify exactly the path supplied by the user.
std::string Unquote(const std::string& value) {
  if (value.size() >= 2 &&
      ((value.front() == '"' && value.back() == '"') ||
       (value.front() == '\'' && value.back() == '\''))) {
    return value.substr(1, value.size() - 2);
  }
  return value;
}

// Reject partial parses (for example "100 km"), overflow, NaN/infinity, zero,
// and negative lengths.  surface_mesh_resolution_m is an SI length.
bool ParsePositiveFinite(const std::string& text, double* result) {
  errno = 0;
  char* end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(value) || value <= 0.0) return false;
  *result = value;
  return true;
}

// Give all parser diagnostics a uniform file:line prefix.  Line zero denotes
// a whole-file error such as a missing required section or key.
std::string At(const fs::path& file, std::size_t line,
               const std::string& message) {
  std::ostringstream out;
  out << file.string();
  if (line != 0) out << ':' << line;
  out << ": srcMoon input error: " << message;
  return out.str();
}

// Resolve relative paths against the input file—not the launch directory—so
// rerunning the same receipt from a different output directory is invariant.
bool AbsolutePath(const fs::path& inputDirectory, const Setting& setting,
                  const std::string& key, fs::path* output,
                  std::string* error, const fs::path& inputFile) {
  fs::path value(Unquote(Trim(setting.value)));
  if (value.empty()) {
    *error = At(inputFile, setting.line, key + " cannot be empty");
    return false;
  }
  if (value.is_relative()) value = inputDirectory / value;
  std::error_code ec;
  value = fs::absolute(value, ec).lexically_normal();
  if (ec) {
    *error = At(inputFile, setting.line,
        "cannot resolve " + key + ": " + ec.message());
    return false;
  }
  *output = value;
  return true;
}

}  // namespace

const char* SurfaceGeometryName(SurfaceGeometry geometry) {
  return geometry == SurfaceGeometry::Sphere ? "sphere" : "lola";
}

bool ParseApplicationInput(const std::string& path, Configuration* result,
                           std::string* error) {
  if (result == nullptr || error == nullptr) return false;
  error->clear();

  std::error_code ec;
  fs::path inputFile = fs::absolute(fs::path(path), ec).lexically_normal();
  if (ec) {
    *error = path + ": srcMoon input error: cannot resolve input path: " +
        ec.message();
    return false;
  }

  std::ifstream stream(inputFile);
  if (!stream) {
    *error = At(inputFile, 0, "cannot open input file");
    return false;
  }

  // This is a deliberately small state machine.  Other named sections are
  // skipped without interpreting their contents, while malformed section
  // boundaries are rejected globally because they make ownership ambiguous.
  std::map<std::string, Setting> settings;
  bool insideSection = false;
  bool insideMoon = false;
  bool sawMoon = false;
  std::size_t lineNumber = 0;
  std::string line;
  while (std::getline(stream, line)) {
    ++lineNumber;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    // AMPS application inputs use '!' for comments.  '#' is reserved for the
    // section delimiters and is therefore not treated as an inline comment.
    const std::size_t comment = line.find('!');
    if (comment != std::string::npos) line.erase(comment);
    line = Trim(line);
    if (line.empty()) continue;

    const std::string lower = Lower(line);
    const std::string begin = "#section begin:";
    if (lower.rfind(begin, 0) == 0) {
      if (insideSection) {
        *error = At(inputFile, lineNumber,
                    "nested #section begin is not allowed");
        return false;
      }
      const std::string name = Lower(Trim(line.substr(begin.size())));
      if (name.empty() || name.find_first_of(" \t") != std::string::npos) {
        *error = At(inputFile, lineNumber,
                    "expected '#section begin: NAME'");
        return false;
      }
      insideSection = true;
      insideMoon = name == "moon";
      if (insideMoon) {
        if (sawMoon) {
          *error = At(inputFile, lineNumber,
                      "duplicate moon section is not allowed");
          return false;
        }
        sawMoon = true;
      }
      continue;
    }
    if (lower.rfind("#section begin", 0) == 0) {
      *error = At(inputFile, lineNumber,
                  "expected '#section begin: NAME'");
      return false;
    }
    if (lower.rfind("#section end", 0) == 0) {
      if (lower != "#section end") {
        *error = At(inputFile, lineNumber,
                    "unrecognized text follows '#section end'");
        return false;
      }
      if (!insideSection) {
        *error = At(inputFile, lineNumber,
                    "#section end has no matching #section begin");
        return false;
      }
      insideSection = false;
      insideMoon = false;
      continue;
    }
    if (lower.rfind("#section", 0) == 0) {
      *error = At(inputFile, lineNumber,
                  "unrecognized #section directive");
      return false;
    }
    if (!insideMoon) continue;

    // Only assignments inside the unique moon section belong to this parser.
    // Use the first '=' so paths or future opaque values may contain another.
    const std::size_t equal = line.find('=');
    if (equal == std::string::npos) {
      *error = At(inputFile, lineNumber, "expected 'NAME = VALUE'");
      return false;
    }
    const std::string key = Lower(Trim(line.substr(0, equal)));
    const std::string value = Trim(line.substr(equal + 1));
    if (key.empty() || value.empty()) {
      *error = At(inputFile, lineNumber,
                  "input assignment is missing its name or value");
      return false;
    }
    if (settings.find(key) != settings.end()) {
      *error = At(inputFile, lineNumber,
                  key + " is specified more than once");
      return false;
    }
    settings.emplace(key, Setting{value, lineNumber});
  }
  if (!stream.eof()) {
    *error = At(inputFile, lineNumber, "I/O failure while reading input");
    return false;
  }
  if (insideSection) {
    *error = At(inputFile, lineNumber,
                "input ends before the matching '#section end'");
    return false;
  }
  if (!sawMoon) {
    *error = At(inputFile, 0, "required moon section was not found");
    return false;
  }

  // Require a complete run receipt in both sphere and LOLA modes.  This makes
  // changing surface_geometry explicit and prevents a later LOLA run from
  // inheriting undeclared paths from compiled defaults or process state.
  const char* required[] = {
      "spice_path", "surface_geometry", "surface_mesh_resolution_m",
      "lola_product_id", "lola_image_file", "lola_label_file",
      "surface_cea_file", "surface_tecplot_file"};
  for (const char* key : required) {
    if (settings.find(key) == settings.end()) {
      *error = At(inputFile, 0,
                  std::string("moon section is missing required ") + key);
      return false;
    }
  }
  for (const auto& item : settings) {
    if (std::find(std::begin(required), std::end(required), item.first) ==
        std::end(required)) {
      *error = At(inputFile, item.second.line,
                  "unrecognized moon setting '" + item.first + "'");
      return false;
    }
  }

  // Build a local candidate and assign *result only at the end.  Callers can
  // safely retain their previous object if any validation below fails.
  Configuration candidate;
  candidate.inputFile = inputFile.string();
  const std::string geometry = Lower(Unquote(settings["surface_geometry"].value));
  if (geometry == "sphere") candidate.surfaceGeometry = SurfaceGeometry::Sphere;
  else if (geometry == "lola") candidate.surfaceGeometry = SurfaceGeometry::Lola;
  else {
    *error = At(inputFile, settings["surface_geometry"].line,
                "surface_geometry must be sphere or lola");
    return false;
  }

  if (!ParsePositiveFinite(settings["surface_mesh_resolution_m"].value,
                           &candidate.surfaceMeshResolutionM)) {
    *error = At(inputFile, settings["surface_mesh_resolution_m"].line,
                "surface_mesh_resolution_m must be finite and positive");
    return false;
  }
  candidate.lolaProductId = Unquote(settings["lola_product_id"].value);
  if (candidate.lolaProductId.empty()) {
    *error = At(inputFile, settings["lola_product_id"].line,
                "lola_product_id cannot be empty");
    return false;
  }

  const fs::path inputDirectory = inputFile.parent_path();
  fs::path spiceRoot, image, label, cea, tecplot;
  if (!AbsolutePath(inputDirectory, settings["spice_path"], "spice_path",
                    &spiceRoot, error, inputFile) ||
      !AbsolutePath(inputDirectory, settings["lola_image_file"],
                    "lola_image_file", &image, error, inputFile) ||
      !AbsolutePath(inputDirectory, settings["lola_label_file"],
                    "lola_label_file", &label, error, inputFile) ||
      !AbsolutePath(inputDirectory, settings["surface_cea_file"],
                    "surface_cea_file", &cea, error, inputFile) ||
      !AbsolutePath(inputDirectory, settings["surface_tecplot_file"],
                    "surface_tecplot_file", &tecplot, error, inputFile)) {
    return false;
  }

  // The requested root follows /home/vtenishe/SPICE: toolkit files live in
  // cspice/, while mission and generic kernels live in Kernels/.
  const fs::path toolkit = spiceRoot / "cspice";
  const fs::path kernels = spiceRoot / "Kernels";
  if (!fs::is_regular_file(toolkit / "include" / "SpiceUsr.h") ||
      !fs::is_regular_file(toolkit / "lib" / "cspice.a") ||
      !fs::is_directory(kernels)) {
    *error = At(inputFile, settings["spice_path"].line,
        "spice_path must contain cspice/include/SpiceUsr.h, "
        "cspice/lib/cspice.a, and Kernels/");
    return false;
  }
  // Validate inputs now; output parents are created later by the atomic mesh
  // writers.  std::ifstream still provides the final readable-file check when
  // the DEM is loaded, including permission and concurrent-removal failures.
  if (!fs::is_regular_file(image) || !fs::is_regular_file(label)) {
    *error = At(inputFile, 0,
        "lola_image_file and lola_label_file must name readable files");
    return false;
  }
  if (cea == tecplot) {
    *error = At(inputFile, 0,
                "surface_cea_file and surface_tecplot_file must differ");
    return false;
  }

  candidate.spiceRoot = spiceRoot.string();
  candidate.spiceToolkitDirectory = toolkit.string();
  candidate.spiceKernelDirectory = kernels.string();
  candidate.lolaImageFile = image.string();
  candidate.lolaLabelFile = label.string();
  candidate.surfaceCeaFile = cea.string();
  candidate.surfaceTecplotFile = tecplot.string();
  *result = std::move(candidate);
  return true;
}

bool InstallConfiguration(const Configuration& configuration,
                          std::string* error) {
  if (error == nullptr) return false;
  if (gInstalled) {
    *error = "srcMoon runtime configuration was already installed";
    return false;
  }
  gConfiguration = configuration;
  gInstalled = true;
  error->clear();
  return true;
}

bool HasConfiguration() { return gInstalled; }

const Configuration& GetConfiguration() {
  // The no-argument regression predates application input.  Preserve it by
  // returning an immutable empty configuration, while production callers use
  // HasConfiguration() to decide whether any fields are meaningful.
  if (!gInstalled) {
    static const Configuration legacy;
    return legacy;
  }
  return gConfiguration;
}

}  // namespace Runtime
}  // namespace Moon
