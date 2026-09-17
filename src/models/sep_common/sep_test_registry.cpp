#include "sep_test_registry.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace SEP {
namespace Testing {

namespace {

// The standalone executable installs this record once, after command-line
// validation and before registry dispatch. Component callbacks receive a const
// view, preventing one case from redirecting a later case's input or evidence.
ExecutionContext executionContext;

}  // namespace

void SetExecutionContext(const ExecutionContext& context) {
  executionContext = context;
}

const ExecutionContext& GetExecutionContext() {
  return executionContext;
}
namespace {

std::string FoldCase(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  return value;
}

bool DescriptorIdLess(const Descriptor& left, const Descriptor& right) {
  return left.id < right.id;
}

bool DescriptorPointerIdLess(const Descriptor* left, const Descriptor* right) {
  return left->id < right->id;
}

std::string JsonEscape(const std::string& value) {
  std::ostringstream out;
  for (std::size_t i = 0; i < value.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(value[i]);
    switch (c) {
      case '\\': out << "\\\\"; break;
      case '"': out << "\\\""; break;
      case '\b': out << "\\b"; break;
      case '\f': out << "\\f"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (c < 0x20U) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<unsigned>(c) << std::dec << std::setfill(' ');
        }
        else out << static_cast<char>(c);
    }
  }
  return out.str();
}

std::string XmlEscape(const std::string& value) {
  std::string result;
  for (std::size_t i = 0; i < value.size(); ++i) {
    switch (value[i]) {
      case '&': result += "&amp;"; break;
      case '<': result += "&lt;"; break;
      case '>': result += "&gt;"; break;
      case '"': result += "&quot;"; break;
      case '\'': result += "&apos;"; break;
      default: result += value[i];
    }
  }
  return result;
}

bool OpenReport(const std::string& path, std::ofstream* output,
                std::string* error) {
  if (path.empty()) {
    if (error) *error = "structured report path is empty";
    return false;
  }
  output->open(path.c_str(), std::ios::out | std::ios::trunc);
  if (!output->good()) {
    if (error) *error = "cannot open structured report: " + path;
    return false;
  }
  return true;
}

void WriteJsonFloatingPoint(std::ostream& out, double value) {
  // JSON has no NaN or infinity tokens. Preserve non-finite failure evidence as
  // a quoted spelling instead of emitting invalid JSON or silently coercing it
  // to zero. Consumers should treat number-or-string as an explicit diagnostic
  // value, while ordinary successful metrics remain JSON numbers.
  if (std::isfinite(value)) out << std::setprecision(17) << value;
  else if (std::isnan(value)) out << "\"nan\"";
  else out << (value > 0.0 ? "\"+infinity\"" : "\"-infinity\"");
}

}  // namespace

int Summary::ExitCode() const {
  if (errors != 0) return 2;
  if (failed != 0) return 1;
  return 0;
}

void Summary::Add(const Result& result) {
  results.push_back(result);
  if (result.status == Status::Pass) ++passed;
  else if (result.status == Status::Fail) ++failed;
  else if (result.status == Status::Skip) ++skipped;
  else ++errors;
}

const char* StatusName(Status status) {
  switch (status) {
    case Status::Pass: return "PASS";
    case Status::Fail: return "FAIL";
    case Status::Skip: return "SKIP";
    case Status::Error: return "ERROR";
  }

  return "ERROR";
}

const char* InitializationName(InitializationLevel level) {
  switch (level) {
    case InitializationLevel::None: return "none";
    case InitializationLevel::FieldLineModel: return "field-line-model";
  }

  return "unknown";
}

const char* RuntimeClassName(RuntimeClass runtime) {
  return runtime == RuntimeClass::Routine ? "routine" : "extended";
}

RefinementOrderEstimate EstimateRefinementOrder(
    double coarseError, double fineError,
    double coarseResolution, double fineResolution) {
  RefinementOrderEstimate result;
  // Zero error cannot define a logarithmic order, even though it may indicate
  // an exact special case.  Refinement tests must choose a nontrivial solution
  // and positive norms so an accidentally exact/empty calculation is visible.
  if (!std::isfinite(coarseError) || !std::isfinite(fineError) ||
      !std::isfinite(coarseResolution) ||
      !std::isfinite(fineResolution) || coarseError <= 0.0 ||
      fineError <= 0.0 || coarseResolution <= fineResolution ||
      fineResolution <= 0.0) {
    result.message =
        "refinement order requires finite positive errors and h_coarse>h_fine>0";
    return result;
  }

  result.coarseToFineErrorRatio = coarseError / fineError;
  result.coarseToFineResolutionRatio =
      coarseResolution / fineResolution;
  result.observedOrder = std::log(result.coarseToFineErrorRatio) /
      std::log(result.coarseToFineResolutionRatio);
  if (!std::isfinite(result.observedOrder)) {
    result.message = "refinement order is not finite";
    return result;
  }
  result.valid = true;
  result.message = "observed refinement order calculated";
  return result;
}

Registry::Registry(std::vector<Descriptor> descriptors)
    : descriptors_(std::move(descriptors)) {
  // Validate once, before a selection can trigger expensive model setup.  A
  // duplicate ID is a programming error: accepting it would make --test ID
  // depend on construction/link order and therefore break reproducibility.
  std::set<std::string> foldedIds;
  for (const Descriptor& descriptor : descriptors_) {
    if (descriptor.id.empty() || descriptor.name.empty() ||
        descriptor.group.empty() || descriptor.description.empty() ||
        descriptor.supportedBuildModes.empty() ||
        descriptor.seedPolicy.empty() || descriptor.stateIsolation.empty() ||
        !descriptor.callback) {
      throw std::invalid_argument(
          "component-test descriptors require complete metadata and a callback");
    }

    const std::string foldedId = FoldCase(descriptor.id);
    if (!foldedIds.insert(foldedId).second) {
      throw std::invalid_argument("duplicate component-test ID: " + descriptor.id);
    }
  }

  std::sort(descriptors_.begin(), descriptors_.end(), DescriptorIdLess);
}

std::vector<const Descriptor*> Registry::Select(
    const std::vector<std::string>& ids,
    const std::vector<std::string>& groups,
    bool allRoutine) const {
  std::set<std::string> requestedIds;
  std::set<std::string> requestedGroups;
  for (const std::string& id : ids) requestedIds.insert(FoldCase(id));
  for (const std::string& group : groups) requestedGroups.insert(FoldCase(group));

  // Diagnose every unknown selector before initialization.  Silent empty group
  // matches are especially dangerous because they could otherwise turn a test
  // command into a successful no-op.
  for (const std::string& requestedId : requestedIds) {
    bool found = false;
    for (const Descriptor& descriptor : descriptors_) {
      if (FoldCase(descriptor.id) == requestedId) found = true;
    }
    if (!found) throw std::invalid_argument("unknown component-test ID: " + requestedId);
  }

  for (const std::string& requestedGroup : requestedGroups) {
    bool found = false;
    for (const Descriptor& descriptor : descriptors_) {
      if (FoldCase(descriptor.group) == requestedGroup) found = true;
    }
    if (!found) {
      throw std::invalid_argument("unknown component-test group: " + requestedGroup);
    }
  }

  std::vector<const Descriptor*> selected;
  for (const Descriptor& descriptor : descriptors_) {
    const bool selectedById =
        requestedIds.count(FoldCase(descriptor.id)) != 0;
    const bool selectedByGroup =
        requestedGroups.count(FoldCase(descriptor.group)) != 0;
    const bool selectedByAll =
        allRoutine && descriptor.runtime == RuntimeClass::Routine;

    if (selectedById || selectedByGroup || selectedByAll) {
      selected.push_back(&descriptor);
    }
  }

  std::sort(selected.begin(), selected.end(), DescriptorPointerIdLess);
  selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
  return selected;
}

Result Registry::RunOne(const Descriptor& descriptor) const {
  const std::chrono::steady_clock::time_point start =
      std::chrono::steady_clock::now();
  Result result;

  try {
    result = descriptor.callback();
  }
  catch (const std::exception& exception) {
    result.status = Status::Error;
    result.message = std::string("unexpected exception: ") + exception.what();
  }
  catch (...) {
    result.status = Status::Error;
    result.message = "unexpected non-standard exception";
  }

  result.id = descriptor.id;
  result.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  if (result.message.empty()) result.message = "test returned no diagnostic message";

  // A legacy callback may return PASS while retaining an old hidden failure
  // counter.  Treat any positive/non-finite assertion_failures metric as FAIL,
  // making it impossible for a converted diagnostic to hide an assertion from
  // the process exit code and structured reports.
  if (result.status == Status::Pass) {
    for (std::size_t i = 0; i < result.metrics.size(); ++i) {
      if (result.metrics[i].name == "assertion_failures" &&
          (!std::isfinite(result.metrics[i].value) ||
           result.metrics[i].value > 0.0)) {
        result.status = Status::Fail;
        result.message = "callback reported PASS with assertion failures: " +
                         result.message;
        break;
      }
    }
  }
  return result;
}

Summary Registry::Run(const std::vector<const Descriptor*>& selected,
                      std::ostream& out) const {
  Summary summary;
  for (const Descriptor* descriptor : selected) {
    if (!descriptor) {
      Result result;
      result.status = Status::Error;
      result.id = "<null>";
      result.message = "selection contains a null test descriptor";
      summary.Add(result);
      PrintResult(result, out);
      continue;
    }

    Result result = RunOne(*descriptor);
    summary.Add(result);
    PrintResult(result, out);
  }

  out << "Component-test summary: PASS=" << summary.passed
      << " FAIL=" << summary.failed << " SKIP=" << summary.skipped
      << " ERROR=" << summary.errors << '\n';
  return summary;
}

bool WriteJsonSummary(const Summary& summary, const std::string& path,
                      std::string* error) {
  std::ofstream out;
  if (!OpenReport(path, &out, error)) return false;
  out << "{\n  \"schema\": \"srcsep-component-tests-v1\",\n"
      << "  \"exit_code\": " << summary.ExitCode() << ",\n"
      << "  \"totals\": {\"passed\": " << summary.passed
      << ", \"failed\": " << summary.failed
      << ", \"skipped\": " << summary.skipped
      << ", \"errors\": " << summary.errors << "},\n"
      << "  \"results\": [\n";
  for (std::size_t i = 0; i < summary.results.size(); ++i) {
    const Result& result = summary.results[i];
    out << "    {\"id\": \"" << JsonEscape(result.id)
        << "\", \"status\": \"" << StatusName(result.status)
        << "\", \"message\": \"" << JsonEscape(result.message)
        << "\", \"elapsed_seconds\": ";
    WriteJsonFloatingPoint(out, result.elapsedSeconds);
    out << ", \"seed\": ";
    if (result.hasSeed) out << result.seed;
    else out << "null";
    out << ", \"configuration\": [";
    for (std::size_t j = 0; j < result.configuration.size(); ++j) {
      if (j != 0) out << ", ";
      out << '"' << JsonEscape(result.configuration[j]) << '"';
    }
    out << "], \"metrics\": [";
    for (std::size_t j = 0; j < result.metrics.size(); ++j) {
      if (j != 0) out << ", ";
      const Metric& metric = result.metrics[j];
      out << "{\"name\": \"" << JsonEscape(metric.name)
          << "\", \"value\": ";
      WriteJsonFloatingPoint(out, metric.value);
      out << ", \"tolerance\": ";
      WriteJsonFloatingPoint(out, metric.tolerance);
      out << ", \"comparison\": \"" << JsonEscape(metric.comparison)
          << "\", \"units\": \"" << JsonEscape(metric.units) << "\"}";
    }
    out << "], \"artifacts\": [";
    for (std::size_t j = 0; j < result.artifacts.size(); ++j) {
      if (j != 0) out << ", ";
      out << '"' << JsonEscape(result.artifacts[j]) << '"';
    }
    out << "]}" << (i + 1 == summary.results.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
  out.close();
  if (out.fail()) {
    if (error) *error = "failed while writing JSON report: " + path;
    return false;
  }
  return true;
}

bool WriteJUnitSummary(const Summary& summary, const std::string& path,
                       std::string* error) {
  std::ofstream out;
  if (!OpenReport(path, &out, error)) return false;
  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<testsuite name=\"srcSEP component tests\" tests=\""
      << summary.results.size() << "\" failures=\"" << summary.failed
      << "\" errors=\"" << summary.errors << "\" skipped=\""
      << summary.skipped << "\">\n";
  for (std::size_t i = 0; i < summary.results.size(); ++i) {
    const Result& result = summary.results[i];
    out << "  <testcase classname=\"srcSEP\" name=\""
        << XmlEscape(result.id) << "\" time=\"" << std::setprecision(17)
        << result.elapsedSeconds << "\">\n";
    if (result.status == Status::Fail)
      out << "    <failure message=\"" << XmlEscape(result.message)
          << "\"/>\n";
    else if (result.status == Status::Error)
      out << "    <error message=\"" << XmlEscape(result.message)
          << "\"/>\n";
    else if (result.status == Status::Skip)
      out << "    <skipped message=\"" << XmlEscape(result.message)
          << "\"/>\n";
    out << "    <system-out>";
    if (result.hasSeed) out << "seed=" << result.seed << "\n";
    for (std::size_t j = 0; j < result.configuration.size(); ++j)
      out << XmlEscape(result.configuration[j]) << "\n";
    for (std::size_t j = 0; j < result.metrics.size(); ++j)
      out << XmlEscape(result.metrics[j].name) << '='
          << result.metrics[j].value << "\n";
    for (std::size_t j = 0; j < result.artifacts.size(); ++j)
      out << "artifact=" << XmlEscape(result.artifacts[j]) << "\n";
    out << "</system-out>\n  </testcase>\n";
  }
  out << "</testsuite>\n";
  out.close();
  if (out.fail()) {
    if (error) *error = "failed while writing JUnit report: " + path;
    return false;
  }
  return true;
}

void Registry::PrintList(std::ostream& out) const {
  out << "Registered srcSEP standalone component tests\n";
  out << "ID | group | class | initialization | build modes | description\n";
  for (const Descriptor& descriptor : descriptors_) {
    out << descriptor.id << " | " << descriptor.group << " | "
        << RuntimeClassName(descriptor.runtime) << " | "
        << InitializationName(descriptor.initialization) << " | "
        << descriptor.supportedBuildModes << " | "
        << descriptor.description << '\n';
    out << "  seed: " << descriptor.seedPolicy << '\n'
        << "  state/isolation: " << descriptor.stateIsolation << '\n';
  }
}

void PrintResult(const Result& result, std::ostream& out) {
  out << '[' << result.id << "] " << StatusName(result.status)
      << " duration_s=" << std::fixed << std::setprecision(6)
      << result.elapsedSeconds << " message=\"" << result.message << "\"\n";

  if (result.hasSeed) out << "  seed=" << result.seed << '\n';
  for (const std::string& setting : result.configuration) {
    out << "  configuration " << setting << '\n';
  }
  for (const Metric& metric : result.metrics) {
    out << "  metric " << metric.name << '=' << std::setprecision(12)
        << metric.value << ' ' << metric.comparison << ' ' << metric.tolerance;
    if (!metric.units.empty()) out << ' ' << metric.units;
    out << '\n';
  }
  for (const std::string& artifact : result.artifacts) {
    out << "  artifact " << artifact << '\n';
  }
}

}  // namespace Testing
}  // namespace SEP
