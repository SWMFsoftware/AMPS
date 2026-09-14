#include "../../util/sep_test_registry.h"
#include "../../util/sep_turbulence_validation.h"

#include <iostream>
#include <string>
#include <vector>

int main() {
  const SEP::Testing::Registry registry(
      SEP::Testing::ControlledTurbulenceDescriptors());
  // The focused runner intentionally selects the explicit group so it executes
  // any future extended turbulence cases as well as the bounded routine set.
  const std::vector<const SEP::Testing::Descriptor*> selected =
      registry.Select({}, {"turbulence"}, false);
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  std::string error;
  if (!SEP::Testing::WriteJsonSummary(
          summary, "controlled-turbulence-results.json", &error) ||
      !SEP::Testing::WriteJUnitSummary(
          summary, "controlled-turbulence-results.xml", &error)) {
    std::cerr << "report error: " << error << '\n';
    return 2;
  }
  return summary.ExitCode();
}
