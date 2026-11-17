#include "sep_polarized_wkb_validation.h"
#include <iostream>
// Execute the shared native callbacks directly: there is no second test model
// or private success flag, and an authoritative FAIL/ERROR reaches process exit.
int main() {
  SEP::Testing::Registry registry(SEP::Testing::PolarizedWkbDescriptors());
  return registry.Run(registry.Select({}, {}, true),std::cout).ExitCode();
}
