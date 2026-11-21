#include "sep_maxwell_juttner_volume_validation.h"
#include "seed_campaign_checks.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::MaxwellJuttnerVolumeDescriptors());
  const auto report=registry.Run(registry.Select({}, {}, true),std::cout);std::string error;
  if(!SEP::Testing::WriteJsonSummary(report,"maxwell-juttner-volume-results.json",&error)||
     !SEP::Testing::WriteJUnitSummary(report,"maxwell-juttner-volume-results.xml",&error)) {
    std::cerr<<error<<'\n';return 2;
  }
  const bool seeds=SEP::Testing::CheckAppliedCampaignSeeds(registry,
      registry.Select({"MJV01","MJV02"},{},false),std::cout);
  return report.ExitCode()!=0?report.ExitCode():(seeds?0:1);
}
