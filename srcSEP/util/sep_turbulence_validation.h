#ifndef SEP_UTIL_SEP_TURBULENCE_VALIDATION_H
#define SEP_UTIL_SEP_TURBULENCE_VALIDATION_H

#include "sep_test_registry.h"

#include <vector>

namespace SEP {
namespace Testing {

// These descriptors exercise the authoritative turbulence core with small,
// fully owned states.  They are shared by the linked component-test CLI and
// the source-only sanitizer runner; neither path substitutes a second physics
// implementation for the production Advance() function.
std::vector<Descriptor> ControlledTurbulenceDescriptors();

}  // namespace Testing
}  // namespace SEP

#endif  // SEP_UTIL_SEP_TURBULENCE_VALIDATION_H
