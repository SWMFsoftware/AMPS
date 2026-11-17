#include "sep_flux_tube_state_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::FluxTubeStateDescriptors());
  const auto summary=registry.Run(registry.Select({}, {}, true),std::cout);std::string error;
  if(!SEP::Testing::WriteJsonSummary(summary,"flux-state-results.json",&error)||
      !SEP::Testing::WriteJUnitSummary(summary,"flux-state-results.xml",&error)){
    std::cerr<<error<<'\n';return 2;
  }
  return summary.ExitCode();
}
