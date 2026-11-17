#ifndef SRCSEP_UTIL_SEP_MAXWELL_JUTTNER_FLUX_VALIDATION_H
#define SRCSEP_UTIL_SEP_MAXWELL_JUTTNER_FLUX_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
namespace SEP { namespace Testing {
// The same callbacks execute in dependency-light sanitizers and the native
// registry. These gates qualify only the full-particle material-plane kernel.
std::vector<Descriptor> MaxwellJuttnerFluxDescriptors();
} }
#endif
