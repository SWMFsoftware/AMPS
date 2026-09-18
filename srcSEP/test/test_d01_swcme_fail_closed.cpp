#include "sep_swcme_validation.h"

#include <iostream>
#include <string>
#include <vector>

int main() {
  // The dependency-light executable intentionally consumes the same descriptor
  // factory linked by ComponentTestRegistry().  It therefore tests the exact
  // D01 callback, metadata validation, result normalization, and exit-code
  // contract used by `amps --test D01`; no second Make-only implementation can
  // silently diverge from the native `test/run_tests.py --all` path.
  SEP::Testing::Registry registry(
      SEP::Testing::SwcmeImprovementDescriptors());
  const std::vector<const SEP::Testing::Descriptor*> selected =
      registry.Select(std::vector<std::string>(1, "D01"),
                      std::vector<std::string>(), false);
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  return summary.ExitCode();
}
