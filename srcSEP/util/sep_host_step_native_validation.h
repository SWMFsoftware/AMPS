#ifndef SRCSEP_UTIL_SEP_HOST_STEP_NATIVE_VALIDATION_H
#define SRCSEP_UTIL_SEP_HOST_STEP_NATIVE_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
namespace SEP { namespace Testing {
// Explicit-protocol, one-case-per-process host/PIC timestep diagnostics. PIC's
// private statistics/first-call flags cannot be restored by a value fixture;
// the application enforces process exclusivity before mesh initialization.
std::vector<Descriptor> HostStepNativeDescriptors();
} }
#endif
