#ifndef SEP_UTIL_SEP_SCIENTIFIC_VALIDATION_H
#define SEP_UTIL_SEP_SCIENTIFIC_VALIDATION_H

#include "sep_test_registry.h"

#include <vector>

namespace SEP {
namespace Testing {

// Step 15 numerical verification cases are ordinary registry descriptors so
// the production CLI, the source-only runner, JSON, and JUnit all observe the
// same callbacks and acceptance thresholds.  Coupled SWCME validation is kept
// in a C++17 runner because SWCME intentionally uses a newer language level
// than the C++11 srcSEP transport core.
std::vector<Descriptor> ScientificValidationDescriptors();

}  // namespace Testing
}  // namespace SEP

#endif  // SEP_UTIL_SEP_SCIENTIFIC_VALIDATION_H
