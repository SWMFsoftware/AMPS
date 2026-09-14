#include "../../util/sep_mover_validation.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: test_registered_movers parker|fte-dmumu|fte-mfp\n";
    return 2;
  }
  const std::string group = argv[1];
  if (group != "parker" && group != "fte-dmumu" && group != "fte-mfp") {
    std::cerr << "unknown controlled mover group: " << group << '\n';
    return 2;
  }

  const SEP::Testing::Registry registry(
      SEP::Testing::ControlledMoverDescriptors());
  std::vector<const SEP::Testing::Descriptor*> selected;
  try {
    // Explicit group selection includes extended statistical cases.  This is
    // intentionally different from --all-tests, whose bounded policy excludes
    // extended ensembles from routine production regression runs.
    selected = registry.Select({}, {group}, false);
  }
  catch (const std::exception& exception) {
    std::cerr << "selection error: " << exception.what() << '\n';
    return 2;
  }

  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  std::string error;
  const std::string stem = "controlled-" + group + "-results";
  if (!SEP::Testing::WriteJsonSummary(summary, stem + ".json", &error) ||
      !SEP::Testing::WriteJUnitSummary(summary, stem + ".xml", &error)) {
    std::cerr << "report error: " << error << '\n';
    return 2;
  }
  return summary.ExitCode();
}
