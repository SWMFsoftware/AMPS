#include "sep_swcme_validation.h"

#include <iostream>
#include <string>
#include <vector>

int main() {
  // D02's source-only gate is a thin launcher for the registry callback.  This
  // preserves its fast compiler/configuration feedback while ensuring that the
  // linked AMPS executable and the standalone target use identical assertions,
  // units, failure messages, and structured-result semantics.
  SEP::Testing::Registry registry(
      SEP::Testing::SwcmeImprovementDescriptors());
  const std::vector<const SEP::Testing::Descriptor*> selected =
      registry.Select(std::vector<std::string>(1, "D02"),
                      std::vector<std::string>(), false);
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  return summary.ExitCode();
}
