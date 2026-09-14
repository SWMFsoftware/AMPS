#ifndef SRCSEP_VALIDATION_CASES_CONTROLLED_TRANSPORT_MODELS_H
#define SRCSEP_VALIDATION_CASES_CONTROLLED_TRANSPORT_MODELS_H

#include <string>
#include <vector>

namespace SEP {
namespace Validation {

// Execute one CV02-CV05 numerical model inside the linked srcSEP/AMPS process.
// The Python orchestration layer supplies reviewed SI-valued arguments and an
// isolated output path. This function calls only production transport kernels;
// analytical and finite-volume references deliberately live outside the linked
// application so the model cannot manufacture its own expected solution.
bool RunControlledTransportModel(const std::string& caseId,
                                 const std::vector<std::string>& arguments,
                                 const std::string& outputPath,
                                 std::string* error);

}  // namespace Validation
}  // namespace SEP

#endif  // SRCSEP_VALIDATION_CASES_CONTROLLED_TRANSPORT_MODELS_H
