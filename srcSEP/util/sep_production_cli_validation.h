#ifndef SRCSEP_PRODUCTION_CLI_VALIDATION_H
#define SRCSEP_PRODUCTION_CLI_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
namespace SEP { namespace Testing {
// Native linked software contracts. The PCLI11 application-state matrix is
// supplemented by separate fresh-process startup evidence, not a physics run.
std::vector<Descriptor> ProductionCliDescriptors();
} }
#endif
