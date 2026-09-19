#ifndef SRCSEP_UTIL_SEP_INITIALIZATION_VALIDATION_H
#define SRCSEP_UTIL_SEP_INITIALIZATION_VALIDATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

#include <vector>

namespace SEP {
namespace Testing {

// Native descriptors for the file-driven Parker/domain/mesh initialization.
// The callbacks are AMPS-independent and are linked both into the production
// --all-tests catalog and the focused source-only gate.
std::vector<Descriptor> InitializationDescriptors();

}  // namespace Testing
}  // namespace SEP

#endif  // SRCSEP_UTIL_SEP_INITIALIZATION_VALIDATION_H
