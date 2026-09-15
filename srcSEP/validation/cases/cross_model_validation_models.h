#ifndef SRCSEP_VALIDATION_CROSS_MODEL_VALIDATION_MODELS_H
#define SRCSEP_VALIDATION_CROSS_MODEL_VALIDATION_MODELS_H

#include <string>
#include <vector>

namespace SEP { namespace Validation {

// Execute the native half of XM01-XM03 inside the selected srcSEP/AMPS
// application. XM01 advances production focused-transport particles. XM02
// advances a publication-informed controlled first-passage ensemble with the
// same production core. XM03 validates and transactionally normalizes a
// production output table; the Python layer never substitutes a reference as
// a model result.
bool RunCrossModelValidationModel(const std::string& caseId,
    const std::vector<std::string>& arguments, const std::string& outputPath,
    std::string* error);

}}

#endif
