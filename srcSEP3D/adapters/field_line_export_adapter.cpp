#include "field_line_export_adapter.h"

#include "sep_field_line_bundle_io.h"

#include <utility>

namespace SEP3D { namespace Adapters {

SEP::Core::Status ExportFieldLineBundle(
    const FieldLineExportRequest& request,
    const SEP::CoronalCME::FieldLineStateEvaluator& evaluator) {
  std::vector<SEP::FieldLine::LineRecord> lines;
  lines.reserve(request.lines.size());
  for (const auto& lineRequest : request.lines) {
    auto line = SEP::CoronalCME::TraceSunConnectedLine(lineRequest, evaluator);
    if (!line.ok()) return line.status;
    lines.push_back(line.value);
  }
  auto set = SEP::CoronalCME::BuildFieldLineSet(
      request.frame, request.epochS, request.historyBeginS,
      request.historyEndS, request.rotationProvenance, std::move(lines));
  if (!set.ok()) return set.status;
  return SEP::FieldLine::WriteBundleTransactional(
      set.value, request.outputDirectory);
}

} }  // namespace SEP3D::Adapters
