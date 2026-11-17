#ifndef SEP_UTIL_SEP_RNG_VALIDATION_H
#define SEP_UTIL_SEP_RNG_VALIDATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

#include <vector>

namespace SEP {
namespace Testing {

// One descriptor catalog for the native CLI and dependency-light lifecycle
// runner. These exercise the exact shared resume/stage functions used by PIC;
// real PIC buffer, MPI migration, and native restart remain additional gates.
std::vector<Descriptor> RngLifecycleDescriptors();

}  // namespace Testing
}  // namespace SEP

#endif
