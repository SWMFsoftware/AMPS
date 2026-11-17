#include "sep_standard_batsrus_wave_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::StandardBatsrusWaveDescriptors());
  const auto summary=registry.Run(registry.Select({}, {}, true),std::cout);std::string error;
  if(!SEP::Testing::WriteJsonSummary(summary,"standard-batsrus-results.json",&error)||
      !SEP::Testing::WriteJUnitSummary(summary,"standard-batsrus-results.xml",&error)){
    std::cerr<<error<<'\n';return 2;
  }
  return summary.ExitCode();
}
