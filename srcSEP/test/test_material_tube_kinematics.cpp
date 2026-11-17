#include "sep_material_tube_kinematics_validation.h"
#include <iostream>
// Native registry and neutral runner execute exactly the same deterministic
// callbacks; authoritative statuses propagate directly to the process exit.
int main() {
  SEP::Testing::Registry registry(SEP::Testing::MaterialTubeKinematicsDescriptors());
  return registry.Run(registry.Select({}, {}, true),std::cout).ExitCode();
}
