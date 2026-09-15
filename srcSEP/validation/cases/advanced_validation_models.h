#ifndef SRCSEP_VALIDATION_CASES_ADVANCED_VALIDATION_MODELS_H
#define SRCSEP_VALIDATION_CASES_ADVANCED_VALIDATION_MODELS_H

#include <string>
#include <vector>

namespace SEP {
namespace Validation {

// Execute one CV06-CV12 controlled numerical stage inside the linked AMPS
// process. Arguments are reviewed SI-valued name/value tokens produced from
// the case JSON. References and scoring remain external so production output
// cannot generate its own expected answer.
bool RunAdvancedValidationModel(const std::string& caseId,
                                const std::vector<std::string>& arguments,
                                const std::string& outputPath,
                                std::string* error);

}  // namespace Validation
}  // namespace SEP

#endif  // SRCSEP_VALIDATION_CASES_ADVANCED_VALIDATION_MODELS_H
