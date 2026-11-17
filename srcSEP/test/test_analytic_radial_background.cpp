#include "sep_analytic_radial_validation.h"
#include <iostream>

int main() {
  SEP::Testing::Registry registry(SEP::Testing::AnalyticRadialDescriptors());
  const auto summary=registry.Run(registry.Select({}, {}, true),std::cout);
  std::string error;
  if(!SEP::Testing::WriteJsonSummary(summary,"analytic-radial-results.json",&error)||
      !SEP::Testing::WriteJUnitSummary(summary,"analytic-radial-results.xml",&error)) {
    std::cerr<<error<<'\n';return 2;
  }
  return summary.ExitCode();
}
