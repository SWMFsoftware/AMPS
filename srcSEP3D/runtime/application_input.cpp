#include "application_input.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
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
  bool sawSep3d = false;
  bool sawValue = false;
  LogicalLine sectionStart;

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
      if (insideSep3d && !sawValue)
        return Invalid(line, "sep3d section is missing required particles_per_iteration");
      insideSection = false;
      insideSep3d = false;
      continue;
    }
    if (lower.rfind("#section", 0) == 0)
      return Invalid(line, "unrecognized #section directive");
    if (!insideSep3d) continue;

    const std::size_t equal = line.text.find('=');
    if (equal == std::string::npos)
      return Invalid(line, "expected 'particles_per_iteration = INTEGER'");
    const std::string key = Lower(Trim(line.text.substr(0, equal)));
    const std::string value = Trim(line.text.substr(equal + 1));
    if (key != "particles_per_iteration")
      return Invalid(line, "unrecognized srcSEP3D setting '" + key + "'");
    if (sawValue)
      return Invalid(line, "particles_per_iteration is specified more than once");
    if (value.empty())
      return Invalid(line, "particles_per_iteration is missing its integer value");
    if (!ParseUnsigned(value, &candidate.particlesPerIteration))
      return Invalid(line, "particles_per_iteration must be an unsigned integer");
    sawValue = true;
    candidate.valueFile = line.file;
    candidate.valueLine = line.line;
  }

  if (insideSection)
    return Invalid(sectionStart, "section reaches end of expanded input without '#section end'");
  if (!sawSep3d)
    return InvalidFile(candidate.rootFile,
        "missing required '#section begin: sep3d' section");
  if (!sawValue)
    return InvalidFile(candidate.rootFile,
        "sep3d section is missing required particles_per_iteration");

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
  summary << "  particles_per_iteration=" << input.particlesPerIteration
          << " (per compiled species)\n"
          << "  value_source=" << input.valueFile << ':' << input.valueLine
          << '\n';
  return summary.str();
}

}  // namespace RuntimeModel
}  // namespace SEP3D
