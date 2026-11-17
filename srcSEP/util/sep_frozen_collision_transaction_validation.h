#ifndef SRCSEP_UTIL_SEP_FROZEN_COLLISION_TRANSACTION_VALIDATION_H
#define SRCSEP_UTIL_SEP_FROZEN_COLLISION_TRANSACTION_VALIDATION_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)
namespace SEP { namespace Testing {
// All fixtures are neutral, deterministic and independently selectable. No
// external domain limits, native PIC state or accepted-event RNG is invented.
std::vector<Descriptor> FrozenCollisionTransactionDescriptors();
} }
#endif
