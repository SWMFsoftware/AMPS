#include "../../util/sep_wave_normalization_validation.h"

#include <iostream>
#include <string>

int main() {
  const SEP::Testing::Registry registry(SEP::Testing::WaveNormalizationDescriptors());
  const auto selected = registry.Select({}, {"wave-normalization"}, false);
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  std::string error;
  if (!SEP::Testing::WriteJsonSummary(summary, "wave-normalization-results.json", &error) ||
      !SEP::Testing::WriteJUnitSummary(summary, "wave-normalization-results.xml", &error)) {
    std::cerr << "report error: " << error << '\n';
    return 2;
  }
  return summary.ExitCode();
}
