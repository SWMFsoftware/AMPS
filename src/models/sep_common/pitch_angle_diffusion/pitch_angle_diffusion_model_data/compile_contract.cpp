#include "pitch_angle_contract.hpp"
#include <type_traits>
static_assert(std::has_virtual_destructor<amps_scattering::ScalarCoefficient>::value);
static_assert(std::has_virtual_destructor<amps_scattering::TensorCoefficient>::value);
static_assert(!std::is_base_of<amps_scattering::ScalarCoefficient,
                              amps_scattering::TensorCoefficient>::value);
int main() { return 0; }
