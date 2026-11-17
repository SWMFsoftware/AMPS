#include "sep_real_wave_normalization_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::RealWaveNormalizationDescriptors());
  return registry.Run(registry.Select({"RWN01","RWN02"}, {}, false),std::cout).ExitCode();
}
