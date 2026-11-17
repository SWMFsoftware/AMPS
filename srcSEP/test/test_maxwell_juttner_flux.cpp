#include "sep_maxwell_juttner_flux_validation.h"
#include "seed_campaign_checks.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::MaxwellJuttnerFluxDescriptors());
  const auto summary=registry.Run(registry.Select({}, {}, true),std::cout);std::string error;
  if(!SEP::Testing::WriteJsonSummary(summary,"maxwell-juttner-results.json",&error)||
      !SEP::Testing::WriteJUnitSummary(summary,"maxwell-juttner-results.xml",&error)){
    std::cerr<<error<<'\n';return 2;
  }
  // Verify actual observable sensitivity and order independence as well as the
  // per-callback distribution gates; a metadata-only seed override would fail.
  const bool seedChecks=SEP::Testing::CheckAppliedCampaignSeeds(registry,
      registry.Select({"MJ03","MJ04","MJ06"}, {}, false),std::cout);
  return summary.ExitCode()!=0?summary.ExitCode():(seedChecks?0:1);
}
