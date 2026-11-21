#ifndef SRCSEP_MAXWELL_JUTTNER_VOLUME_H
#define SRCSEP_MAXWELL_JUTTNER_VOLUME_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <array>

namespace SEP { namespace Sources { namespace MaxwellJuttnerVolume {
struct Input {
  // SWMF supplies the fluid thermodynamic state; the caller must explicitly
  // establish its density/temperature basis. No temperature partition or
  // composition is invented by this momentum kernel. All quantities are SI.
  double properTemperatureK=0.0,restMassKg=0.0;
  std::array<double,3> bulkVelocityMPerS{{0,0,0}};
};
struct Sample {
  // Guiding-center state is stored in the local plasma frame. Keep gamma'-1
  // separately so cold kinetic energy is not erased by adding one. The phase
  // measure is that of an instantaneous simulation-frame spatial volume.
  std::array<double,3> restMomentumKgMPerS{{0,0,0}};
  double restKineticGamma=0.0,bulkGamma=1.0,simulationGamma=1.0;
};
Transport::Status BulkGamma(const std::array<double,3>& velocity,double* gamma);
Transport::Status Draw(const Input& input,Transport::KeyedRandomStream* random,
                       Sample* sample);
} } }
#endif
