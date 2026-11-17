#include "sep_analytic_shock_epoch_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::AnalyticShockEpochDescriptors());
  const auto summary=registry.Run(registry.Select({}, {}, true),std::cout);
  std::string error;
  if(!SEP::Testing::WriteJsonSummary(summary,"analytic-shock-epoch-results.json",&error)||
      !SEP::Testing::WriteJUnitSummary(summary,"analytic-shock-epoch-results.xml",&error)) {
    std::cerr<<error<<'\n';return 2;
  }
  return summary.ExitCode();
}
