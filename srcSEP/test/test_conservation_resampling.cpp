#include "../util/sep_conservation_resampling_validation.h"
#include "seed_campaign_checks.h"
#include <iostream>
int main() {
  const SEP::Testing::Registry registry(SEP::Testing::ConservationResamplingDescriptors());
  const auto summary=registry.Run(registry.Select({}, {"conservative-resampling"},false),std::cout);
  std::string error;
  if(!SEP::Testing::WriteJsonSummary(summary,"conservative-resampling.json",&error)||
     !SEP::Testing::WriteJUnitSummary(summary,"conservative-resampling.xml",&error)) {
    std::cerr<<error<<'\n';return 2;
  }
  const bool seeded=SEP::Testing::CheckAppliedCampaignSeeds(registry,
      registry.Select({"CRS04","CRS09"},{},false),std::cout);
  return seeded?summary.ExitCode():1;
}
