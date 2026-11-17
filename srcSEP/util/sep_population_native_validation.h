#ifndef SRCSEP_UTIL_SEP_POPULATION_NATIVE_VALIDATION_H
#define SRCSEP_UTIL_SEP_POPULATION_NATIVE_VALIDATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
#include <vector>

namespace SEP { namespace Testing {

// Real PIC record/controller checks. These require the native field-line model;
// they do not substitute synthetic partition tests for MPI transfer or restart.
std::vector<Descriptor> PopulationNativeDescriptors();

}} // namespace SEP::Testing
#endif
