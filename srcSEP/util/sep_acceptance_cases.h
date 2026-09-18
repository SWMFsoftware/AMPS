#ifndef SEP_UTIL_SEP_ACCEPTANCE_CASES_H
#define SEP_UTIL_SEP_ACCEPTANCE_CASES_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

#include <vector>

namespace SEP {
namespace Testing {

// Return the dependency-light Step 13 acceptance fixtures that can execute in
// both the linked srcSEP CLI and a source-only test build.  Keeping the callback
// implementations outside component_tests.cpp prevents the production registry
// from acquiring a second, subtly different copy of the cross-mover rules.
std::vector<Descriptor> AcceptanceCaseDescriptors();

}  // namespace Testing
}  // namespace SEP

#endif  // SEP_UTIL_SEP_ACCEPTANCE_CASES_H
