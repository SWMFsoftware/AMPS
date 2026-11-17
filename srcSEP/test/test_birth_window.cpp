#include "sep_birth_window_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::BirthWindowDescriptors());
  const auto report=registry.Run(registry.Select({}, {}, true),std::cout);
  std::string error;
  if(!SEP::Testing::WriteJsonSummary(report,"birth-window-results.json",&error)||
      !SEP::Testing::WriteJUnitSummary(report,"birth-window-results.xml",&error)) {
    std::cerr<<error<<'\n';return 2;
  }
  return report.ExitCode();
}
