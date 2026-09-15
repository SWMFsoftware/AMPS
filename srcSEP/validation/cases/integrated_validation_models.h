#ifndef SRCSEP_VALIDATION_CASES_INTEGRATED_VALIDATION_MODELS_H
#define SRCSEP_VALIDATION_CASES_INTEGRATED_VALIDATION_MODELS_H

#include <string>
#include <vector>

namespace SEP {
namespace Validation {

// Execute one IV01-IV06 numerical stage inside the linked srcSEP/AMPS process.
// The arguments are strict SI name/value tokens generated from the reviewed
// case JSON.  This function writes raw model evidence only; independent
// references, acceptance decisions, and figures remain outside the executable.
bool RunIntegratedValidationModel(const std::string& caseId,
                                  const std::vector<std::string>& arguments,
                                  const std::string& outputPath,
                                  std::string* error);

}  // namespace Validation
}  // namespace SEP

#endif
