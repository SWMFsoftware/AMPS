#include "sep_frozen_collision_transaction_validation.h"
#include <iostream>
// Exercise the same individually selectable descriptors as native AMPS. This
// executable owns no MPI/PIC initialization and propagates every registry error.
int main() {
  SEP::Testing::Registry registry(SEP::Testing::FrozenCollisionTransactionDescriptors());
  return registry.Run(registry.Select({}, {}, true),std::cout).ExitCode();
}
