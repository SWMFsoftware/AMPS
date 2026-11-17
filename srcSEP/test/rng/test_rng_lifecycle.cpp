#include "../../util/sep_rng_validation.h"
#include "../seed_campaign_checks.h"

#include <iostream>
#include <string>

int main() {
  const SEP::Testing::Registry registry(SEP::Testing::RngLifecycleDescriptors());
  const auto selected = registry.Select({}, {"rng-lifecycle"}, false);
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  std::string error;
  if (!SEP::Testing::WriteJsonSummary(summary, "rng-lifecycle-results.json", &error) ||
      !SEP::Testing::WriteJUnitSummary(summary, "rng-lifecycle-results.xml", &error)) {
    std::cerr << "report error: " << error << '\n';
    return 2;
  }
  const bool seeded=SEP::Testing::CheckAppliedCampaignSeeds(registry,
      registry.Select({"RNG01","RNG04","RNG10"},{},false),std::cout);
  return seeded?summary.ExitCode():1;
}
