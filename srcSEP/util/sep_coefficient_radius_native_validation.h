#ifndef SRCSEP_UTIL_SEP_COEFFICIENT_RADIUS_NATIVE_VALIDATION_H
#define SRCSEP_UTIL_SEP_COEFFICIENT_RADIUS_NATIVE_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
namespace SEP { namespace Testing {
// Extended because this read-only gate requires its own real immutable
// --input initialization, rather than overriding configuration in-process.
std::vector<Descriptor> CoefficientRadiusNativeDescriptors();
} }
#endif
