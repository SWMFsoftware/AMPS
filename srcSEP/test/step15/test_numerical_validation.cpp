#include "../../util/sep_scientific_validation.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  const SEP::Testing::Registry registry(
      SEP::Testing::ScientificValidationDescriptors());
  std::vector<std::string> ids;
  if (argc == 2) {
    ids.push_back(argv[1]);
  }
  else if (argc != 1) {
    std::cerr << "Usage: test_numerical_validation [VAL01|VAL02|VAL03]\n";
    return 2;
  }

  // An omitted ID means the explicit Step 15 source-only campaign, including
  // extended cases.  We do not use allRoutine because VAL02/VAL03 are
  // intentionally excluded from production --all-tests due to their ensemble
  // size, but the focused validation target must execute all three.
  std::vector<const SEP::Testing::Descriptor*> selected;
  if (ids.empty()) {
    selected = registry.Select({"VAL01", "VAL02", "VAL03"}, {}, false);
  }
  else {
    try {
      selected = registry.Select(ids, {}, false);
    }
    catch (const std::exception& exception) {
      std::cerr << "Selection error: " << exception.what() << '\n';
      return 2;
    }
  }
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  std::string error;
  if (!SEP::Testing::WriteJsonSummary(
          summary, "step15-numerical-results.json", &error) ||
      !SEP::Testing::WriteJUnitSummary(
          summary, "step15-numerical-results.xml", &error)) {
    std::cerr << "Report error: " << error << '\n';
    return 2;
  }
  return summary.ExitCode();
}
