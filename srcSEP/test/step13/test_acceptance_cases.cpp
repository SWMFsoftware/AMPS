#include "../../util/sep_acceptance_cases.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main() {
  // Exercise the same descriptors linked into ComponentTestRegistry.  The
  // source-only runner is deliberately a thin registry client so CI without an
  // AMPS checkout cannot drift from the production --all-tests callbacks.
  const SEP::Testing::Registry registry(
      SEP::Testing::AcceptanceCaseDescriptors());
  const std::vector<const SEP::Testing::Descriptor*> selected =
      registry.Select({}, {}, true);
  std::ostringstream output;
  const SEP::Testing::Summary summary = registry.Run(selected, output);
  std::cout << output.str();

  std::string error;
  if (!SEP::Testing::WriteJsonSummary(summary, "acceptance-results.json", &error)) {
    std::cerr << "REPORT ERROR: " << error << '\n';
    return 2;
  }
  if (!SEP::Testing::WriteJUnitSummary(summary, "acceptance-results.xml", &error)) {
    std::cerr << "REPORT ERROR: " << error << '\n';
    return 2;
  }
  if (selected.size() != 4U) {
    std::cerr << "REGISTRY ERROR: expected four Step 13 acceptance fixtures\n";
    return 2;
  }
  return summary.ExitCode();
}
