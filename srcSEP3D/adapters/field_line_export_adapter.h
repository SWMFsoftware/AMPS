#ifndef SRCSEP3D_FIELD_LINE_EXPORT_ADAPTER_H
#define SRCSEP3D_FIELD_LINE_EXPORT_ADAPTER_H

#include "sep_coronal_cme/field_line_reduction.h"
#include "sep_status.h"

#include <string>
#include <vector>

namespace SEP3D { namespace Adapters {

// This adapter contains no tracing or source physics.  The application fills
// requests from its parsed repeated [field_line.ID] sections and supplies an
// evaluator backed by owning immutable background snapshots.  All reduction,
// validation, stable ordering, checksums, and publication remain in the two
// shared model libraries.
struct FieldLineExportRequest {
  std::string frame;
  double epochS = 0.0;
  double historyBeginS = 0.0;
  double historyEndS = 0.0;
  std::string rotationProvenance;
  std::vector<SEP::CoronalCME::FieldLineTraceRequest> lines;
  std::string outputDirectory;
};

SEP::Core::Status ExportFieldLineBundle(
    const FieldLineExportRequest& request,
    const SEP::CoronalCME::FieldLineStateEvaluator& evaluator);

} }  // namespace SEP3D::Adapters

#endif  // SRCSEP3D_FIELD_LINE_EXPORT_ADAPTER_H
