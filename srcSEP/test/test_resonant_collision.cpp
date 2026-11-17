#include "sep_resonant_collision_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::ResonantCollisionDescriptors());
  return registry.Run(registry.Select({}, {}, true),std::cout).ExitCode();
}
