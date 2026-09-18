#include "sep_swcme_validation.h"

#include <iostream>
#include <string>
#include <vector>

int main() {
  // Run the exact extended registry callback that `test/run_tests.py --all`
  // discovers as D03PRE.  This source-only executable verifies the preflight
  // contract; it intentionally cannot substitute for the external campaign,
  // whose MPI decompositions and restart launches require a configured AMPS
  // executable and site-owned inputs.
  SEP::Testing::Registry registry(
      SEP::Testing::SwcmeImprovementDescriptors());
  const std::vector<const SEP::Testing::Descriptor*> selected =
      registry.Select(std::vector<std::string>(1, "D03PRE"),
                      std::vector<std::string>(), false);
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  return summary.ExitCode();
}
