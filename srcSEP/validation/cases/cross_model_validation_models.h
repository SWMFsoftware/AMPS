#ifndef SRCSEP_VALIDATION_CROSS_MODEL_VALIDATION_MODELS_H
#define SRCSEP_VALIDATION_CROSS_MODEL_VALIDATION_MODELS_H

#include <string>
#include <vector>

namespace SEP { namespace Validation {

// Execute the native half of XM01-XM03 inside the selected srcSEP/AMPS
// application. XM01 advances production focused-transport particles. XM02
// advances a publication-informed controlled first-passage ensemble with the
// same production core. XM03 advances an event-informed Parker-transport
// ensemble on an Earth-connected spiral and compares that linked result with
// the Earth observations digitized from Liu et al. Figure 12.  No case reads a
// reference table as model output.
bool RunCrossModelValidationModel(const std::string& caseId,
    const std::vector<std::string>& arguments, const std::string& outputPath,
    std::string* error);

}}

#endif
