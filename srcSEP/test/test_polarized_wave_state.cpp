#include "sep_polarized_wave_state_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::PolarizedWaveStateDescriptors());
  return registry.Run(registry.Select({}, {}, true),std::cout).ExitCode();
}
