#ifndef SRCSEP_VALIDATION_CROSS_MODEL_VALIDATION_MODELS_H
#define SRCSEP_VALIDATION_CROSS_MODEL_VALIDATION_MODELS_H

#include <string>
#include <vector>

namespace SEP { namespace Validation {

// Execute the native half of XM01-XM03, OV01-OV05, and EV01-EV02 inside the
// selected srcSEP/AMPS application. XM01-XM03 retain their independent-model
// and event-reconstruction roles. OV/EV calls advance the same production
// Parker core and emit model-only profiles; publication measurements and
// sealed campaign targets remain exclusively in the Python evidence layer.
// No native case reads a reference table as model output.
bool RunCrossModelValidationModel(const std::string& caseId,
    const std::vector<std::string>& arguments, const std::string& outputPath,
    std::string* error);

}}

#endif
