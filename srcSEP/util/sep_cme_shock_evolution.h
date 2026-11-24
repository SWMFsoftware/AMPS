#ifndef SEP_UTIL_CME_SHOCK_EVOLUTION_H
#define SEP_UTIL_CME_SHOCK_EVOLUTION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

namespace SEP::Testing {
// Run the prepared, fixed-cadence shock regression on the designated writer.
// The registry owns MPI and reporting. This callback owns only test state.
Result RunCmeShockEvolution();
}
#endif
