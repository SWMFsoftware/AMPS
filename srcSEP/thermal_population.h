#ifndef SRCSEP_THERMAL_POPULATION_H
#define SRCSEP_THERMAL_POPULATION_H
#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
namespace SEP { namespace ParticleSource {
// Fresh imported proton initial condition only. This is a volume distribution,
// not a boundary flux or a repair of existing/loaded particle extension bytes.
Transport::Status PopulateImportedThermalInitialCondition();
} }
#endif
