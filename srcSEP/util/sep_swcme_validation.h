#ifndef SEP_UTIL_SEP_SWCME_VALIDATION_H
#define SEP_UTIL_SEP_SWCME_VALIDATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

#include <vector>

namespace SEP {
namespace Testing {

// Return the native-registry tests introduced with improvements D01--D03.
//
// These descriptors are defined outside component_tests.cpp so the exact same
// callbacks can also be compiled by the dependency-light focused gates.  This
// prevents the Make-only checks and `amps --test <ID>` from evolving into two
// different definitions of the same acceptance criterion.
//
// D01 and D02 exercise the production srcSEP adapter and the canonical SWCME
// model/configuration implementation.  D03PRE is intentionally a *preflight*:
// it proves that the linked binary contains the three production mover
// contracts and the refresh diagnostics consumed by the D03 campaign.  The
// multi-process/restart campaign itself must remain outside the executable;
// recursively rebuilding and launching the running program from a component
// callback would be unsafe and would not be valid MPI evidence.
std::vector<Descriptor> SwcmeImprovementDescriptors();

}  // namespace Testing
}  // namespace SEP

#endif  // SEP_UTIL_SEP_SWCME_VALIDATION_H
