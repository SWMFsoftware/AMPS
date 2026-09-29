#include "field_line_bundle_adapter.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace SEP1D { namespace Adapters {
namespace {

bool ParseFiniteDouble(const std::string& text, double* value) {
  if (value == nullptr || text.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ParsePointCount(const std::string& text, std::size_t* value) {
  if (value == nullptr || text.empty() || text.front() == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      parsed > static_cast<unsigned long long>(
          std::numeric_limits<std::size_t>::max())) return false;
  *value = static_cast<std::size_t>(parsed);
  return true;
}

SEP::Core::Result<std::vector<LineMeshRequest>> LineMeshFailure(
    const std::string& message) {
  return SEP::Core::Result<std::vector<LineMeshRequest>>::Failure(
      SEP::Core::StatusCode::InvalidConfiguration, message);
}

}  // namespace

SEP::Core::Result<std::vector<LineMeshRequest>> ParseLineMeshRequests(
    const std::map<std::string, std::string>& normalizedAssignments) {
  const std::string prefix = "line_mesh.";
  std::map<std::string, std::map<std::string, std::string>> records;

  for (const auto& assignment : normalizedAssignments) {
    if (assignment.first.compare(0, prefix.size(), prefix) != 0) continue;
    const std::string remainder = assignment.first.substr(prefix.size());
    const std::size_t separator = remainder.find('.');
    if (separator == std::string::npos || separator == 0 ||
        separator + 1 == remainder.size() ||
        remainder.find('.', separator + 1) != std::string::npos) {
      return LineMeshFailure("invalid repeated line_mesh key '" +
                             assignment.first + "'");
    }
    const std::string meshId = remainder.substr(0, separator);
    const std::string key = remainder.substr(separator + 1);
    const std::set<std::string> allowed = {
        "line_id", "start_s_m", "end_s_m", "point_count", "resampling"};
    if (allowed.count(key) == 0)
      return LineMeshFailure("unknown key '" + key + "' in [line_mesh." +
                             meshId + "]");
    records[meshId][key] = assignment.second;
  }

  if (records.empty())
    return LineMeshFailure("at least one [line_mesh.ID] section is required");

  std::vector<LineMeshRequest> requests;
  std::set<std::string> selectedLines;
  for (const auto& record : records) {
    const auto& values = record.second;
    const std::vector<std::string> required = {
        "line_id", "start_s_m", "end_s_m", "point_count", "resampling"};
    for (const std::string& key : required)
      if (values.count(key) == 0)
        return LineMeshFailure("incomplete [line_mesh." + record.first +
                               "]: missing " + key);

    LineMeshRequest request;
    request.stableMeshId = record.first;
    request.stableLineId = values.at("line_id");
    if (request.stableLineId.empty())
      return LineMeshFailure("line_mesh." + record.first +
                             ".line_id must not be empty");
    if (!ParseFiniteDouble(values.at("start_s_m"),
                           &request.startArcLengthM) ||
        !ParseFiniteDouble(values.at("end_s_m"), &request.endArcLengthM) ||
        !ParsePointCount(values.at("point_count"), &request.pointCount))
      return LineMeshFailure("line_mesh." + record.first +
                             " requires finite bare-SI bounds and an integer "
                             "point_count");
    if (!(request.startArcLengthM < request.endArcLengthM) ||
        request.pointCount < 2)
      return LineMeshFailure("line_mesh." + record.first +
                             " requires start_s_m < end_s_m and point_count >= 2");
    if (values.at("resampling") != "conservative-positive-one-sided")
      return LineMeshFailure("line_mesh." + record.first +
                             ".resampling must be conservative-positive-one-sided");
    if (!selectedLines.insert(request.stableLineId).second)
      return LineMeshFailure("physical line '" + request.stableLineId +
                             "' has more than one [line_mesh.ID] record");
    requests.push_back(std::move(request));
  }
  return SEP::Core::Result<std::vector<LineMeshRequest>>::Success(
      std::move(requests));
}

SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>> BuildLineMeshes(
    const SEP::FieldLine::FieldLineSet& bundle,
    const std::vector<LineMeshRequest>& requests) {
  const SEP::Core::Status valid = SEP::FieldLine::ValidateFieldLineSet(bundle);
  if (!valid.ok())
    return SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>>::Failure(
        valid.code, valid.message);
  if (requests.empty())
    return SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>>::Failure(
        SEP::Core::StatusCode::InvalidConfiguration,
        "at least one typed line mesh request is required");

  std::vector<SEP::FieldLine::LineRecord> meshes;
  meshes.reserve(requests.size());
  std::set<std::string> lineIds;
  for (const LineMeshRequest& request : requests) {
    if (!lineIds.insert(request.stableLineId).second)
      return SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>>::Failure(
          SEP::Core::StatusCode::InvalidConfiguration,
          "duplicate typed mesh request for line '" + request.stableLineId + "'");
    const SEP::FieldLine::LineRecord* selected = nullptr;
    for (const auto& line : bundle.lines)
      if (line.stableLineId == request.stableLineId) selected = &line;
    if (selected == nullptr)
      return SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>>::Failure(
          SEP::Core::StatusCode::OutOfDomain,
          "line mesh selects absent stable line '" + request.stableLineId + "'");
    auto mesh = SEP::FieldLine::ResampleLine(
        *selected, request.startArcLengthM, request.endArcLengthM,
        request.pointCount);
    if (!mesh.ok())
      return SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>>::Failure(
          mesh.status.code, "line_mesh." + request.stableMeshId + ": " +
                                mesh.status.message);
    meshes.push_back(std::move(mesh.value));
  }
  return SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>>::Success(
      std::move(meshes));
}

SEP::Core::Status FieldLineBundleAdapter::LoadAndPublish(
    const std::string& directory, const SnapshotPublisher& publisher) {
  if (!publisher)
    return SEP::Core::Status::Failure(
        SEP::Core::StatusCode::InvalidConfiguration,
        "field-line bundle adapter requires a snapshot publisher");
  auto candidate = SEP::FieldLine::ReadBundle(directory);
  if (!candidate.ok()) return candidate.status;
  // Publish before replacing the owning handle.  If application-specific
  // coefficient registration fails, the previously active bundle remains
  // valid and visible to particles.
  auto status = publisher(candidate.value);
  if (!status.ok()) return status;
  bundle_ = std::make_shared<const SEP::FieldLine::FieldLineSet>(
      std::move(candidate.value));
  identity_ = bundle_->bundleId;
  return SEP::Core::Status::Success();
}

SEP::Core::Result<SEP::FieldLine::NodeState>
FieldLineBundleAdapter::Evaluate(const std::string& stableLineId,
                                 double arcLengthM) const {
  if (!bundle_)
    return SEP::Core::Result<SEP::FieldLine::NodeState>::Failure(
        SEP::Core::StatusCode::InvalidState,
        "field-line bundle has not been published");
  for (const auto& line : bundle_->lines)
    if (line.stableLineId == stableLineId)
      return SEP::FieldLine::InterpolateBounded(line, arcLengthM);
  return SEP::Core::Result<SEP::FieldLine::NodeState>::Failure(
      SEP::Core::StatusCode::OutOfDomain,
      "stable field-line ID is not present in the bundle");
}

} }  // namespace SEP1D::Adapters
