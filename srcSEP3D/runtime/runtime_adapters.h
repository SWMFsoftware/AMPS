// AMPS-independent front ends for the common Runtime lifecycle.  Phase B will
// populate complete physical fields; R2 publishes validated snapshot metadata
// through exactly the same Runtime methods for analytic and SWMF authorities.

#ifndef SEP3D_RUNTIME_RUNTIME_ADAPTERS_H
#define SEP3D_RUNTIME_RUNTIME_ADAPTERS_H

#include "runtime.h"

namespace SEP3D { namespace Background { class BackgroundSnapshot; } }

namespace SEP3D {
namespace RuntimeModel {

class StandaloneAdapter {
 public:
  Core::Status Initialize(Runtime* runtime) const;
  Core::Status PublishFrozenParker(Runtime* runtime, double epochS,
                                   double validUntilS,
                                   std::uint64_t generation) const;
  Core::Status PublishSnapshot(
      Runtime* runtime,
      const Background::BackgroundSnapshot& snapshot) const;
};

class SwmfAdapter {
 public:
  Core::Status Initialize(Runtime* runtime) const;
  Core::Status PublishImported(Runtime* runtime, double epochS,
                               double validUntilS,
                               std::uint64_t generation,
                               bool complete) const;
  Core::Status PublishSnapshot(
      Runtime* runtime,
      const Background::BackgroundSnapshot& snapshot) const;
};

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_RUNTIME_ADAPTERS_H
