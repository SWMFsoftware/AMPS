#ifndef SEP_UTIL_SEP_MOVER_VALIDATION_H
#define SEP_UTIL_SEP_MOVER_VALIDATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

#include <vector>

namespace SEP {
namespace Testing {

// ControlledMoverDescriptors is the single catalog for the dependency-light
// Parker, coefficient-driven focused-transport, and event-driven MFP tests.
// Both the linked CLI and the sanitizer runners consume these descriptors, so
// a stable ID cannot silently execute different physics in the two contexts.
std::vector<Descriptor> ControlledMoverDescriptors();

}  // namespace Testing
}  // namespace SEP

#endif  // SEP_UTIL_SEP_MOVER_VALIDATION_H
