#include "sep_background_derivative_validation.h"
#include <iostream>
int main() {
  SEP::Testing::Registry registry(SEP::Testing::BackgroundDerivativeDescriptors());
  return registry.Run(registry.Select({}, {}, true),std::cout).ExitCode();
}
