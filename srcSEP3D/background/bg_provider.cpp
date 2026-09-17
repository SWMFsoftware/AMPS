#include "bg_provider.h"

namespace SEP3D {
namespace Background {

Core::Status BackgroundProvider::EvaluateBatchDetailed(
    const double* xM, const double* yM, const double* zM,
    std::size_t count, BackgroundSample* output,
    Core::Status* perSampleStatus) const {
  if ((count != 0) && (xM == nullptr || yM == nullptr || zM == nullptr ||
                       output == nullptr || perSampleStatus == nullptr)) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "batch background arrays must not be null");
  }

  Core::Status aggregate = Core::Status::OK();
  for (std::size_t i = 0; i < count; ++i) {
    const BackgroundSample candidate = Evaluate({xM[i], yM[i], zM[i]});
    Core::Status sampleStatus = candidate.status;
    if (sampleStatus.ok() && !candidate.valid) {
      sampleStatus = Core::Status(
          Core::StatusCode::BackgroundInvalid,
          "provider returned an invalid sample without an error status");
    }
    perSampleStatus[i] = sampleStatus;
    if (sampleStatus.ok() && candidate.valid) {
      output[i] = candidate;
    } else if (aggregate.ok()) {
      aggregate = sampleStatus;
    }
  }
  return aggregate;
}

}  // namespace Background
}  // namespace SEP3D
