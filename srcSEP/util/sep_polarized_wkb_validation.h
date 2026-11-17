#ifndef SRCSEP_UTIL_SEP_POLARIZED_WKB_VALIDATION_H
#define SRCSEP_UTIL_SEP_POLARIZED_WKB_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
// One callback catalog is compiled by both native registry and neutral runner.
// Every case declares controlled ideal SI inputs and independently integrated
// references; a PASS never changes the physical production capability gate.
namespace SEP { namespace Testing { std::vector<Descriptor> PolarizedWkbDescriptors(); } }
#endif
