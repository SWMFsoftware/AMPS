#ifndef SRCSEP_UTIL_SEP_ANALYTIC_RADIAL_NATIVE_VALIDATION_H
#define SRCSEP_UTIL_SEP_ANALYTIC_RADIAL_NATIVE_VALIDATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

namespace SEP { namespace Testing {
// These callbacks inspect real startup/PIC storage and actual list adapters.
// They require the native field-line initialization, unlike the pure tensor
// references. No source sampler, moving grid or production closure is implied.
std::vector<Descriptor> AnalyticRadialNativeDescriptors();
} }

#endif
