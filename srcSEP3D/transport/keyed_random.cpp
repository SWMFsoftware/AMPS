#include "keyed_random.h"

#include "../core/sep3d_types.h"

#include <cmath>

namespace SEP3D {
namespace Transport {
namespace {

std::uint64_t Mix(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

void Fold(std::uint64_t value, std::uint64_t* hash) {
  *hash = Mix(*hash ^ Mix(value));
}

}  // namespace

std::uint64_t KeyedRandomStream::Hash(const RandomKey& key,
                                      std::uint64_t drawCounter) {
  std::uint64_t hash = UINT64_C(0x6a09e667f3bcc909);
  Fold(key.campaignSeed, &hash);
  Fold(key.particleId, &hash);
  Fold(key.step, &hash);
  Fold(key.substep, &hash);
  Fold(static_cast<std::uint64_t>(key.purpose), &hash);
  Fold(drawCounter, &hash);
  return hash;
}

double KeyedRandomStream::UniformOpen01() {
  const std::uint64_t bits = Hash(state_.key, state_.drawCounter++);
  // Retain the high 53 bits and center the integer in its 2^-53-wide bin.
  // The largest result is below one and the smallest is strictly positive.
  const std::uint64_t mantissa = bits >> 11;
  return (static_cast<double>(mantissa) + 0.5) *
         (1.0 / 9007199254740992.0);
}

double KeyedRandomStream::Normal01() {
  // No cached second deviate is used: one normal always consumes exactly two
  // counters, which keeps restart and purpose-isolation bookkeeping obvious.
  const double radius = std::sqrt(-2.0 * std::log(UniformOpen01()));
  const double angle = 2.0 * Core::Const::kPi * UniformOpen01();
  return radius * std::cos(angle);
}

}  // namespace Transport
}  // namespace SEP3D
