#ifndef SRCSEP_UTIL_SEP_ANALYTIC_RADIAL_VALIDATION_H
#define SRCSEP_UTIL_SEP_ANALYTIC_RADIAL_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
namespace SEP { namespace Testing {
// Independent fixed-radial/tensor and geometry-only API contracts. Native
// startup/list storage is tested separately by ARBN01/02, never inferred here.
std::vector<Descriptor> AnalyticRadialDescriptors();
} }
#endif
