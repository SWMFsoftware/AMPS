#ifndef SRCSEP_CLI_CAMPAIGN_VALIDATION_H
#define SRCSEP_CLI_CAMPAIGN_VALIDATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_test_registry.h)

namespace SEP { namespace Testing {
// Exercises the native campaign harness only; this is not a physics gate.
std::vector<Descriptor> CliCampaignDescriptors();
} }

#endif
