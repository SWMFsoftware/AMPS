#ifndef SRCSEP_PHYSICS_CHECKPOINT_H
#define SRCSEP_PHYSICS_CHECKPOINT_H

#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <cstdint>
#include <string>

namespace SEP { namespace PhysicsCheckpoint {

// Collective, corrected-model physics checkpoint. First supported boundary:
// fixed native lines, standalone analytic provider, explicit typed coefficient
// models, MPI-only execution, species-global weights and single/species-global
// timesteps, tracker off,
// identical writer/reader rank count and segment ownership. All particle/wave
// workers must have joined. Unsupported capabilities fail before publication.
// Audited native physics fields are captured; diagnostic/output continuation
// is false. Schema four adds explicit fieldwise vertex XYZ, validates particle
// Cartesian views against those decoded endpoints, and requires strictly
// positive drawn MFP residual hazard. Earlier corrected schemas are rejected.
// Schema-three same/fresh-process fixtures passed; this hardened schema still
// needs their linked rerun. Full driver resume, material geometry, live four-
// polarization authority and output/sample continuation remain separate gates.
// The run purpose is persisted and must match on load. Retained-operator
// diagnostic checkpoints cannot bypass the production physical-closure gate.
Transport::Status SaveCollective(const std::string& path,
                                  std::uint64_t completedSteps);
Transport::Status LoadCollective(const std::string& path,
                                  std::uint64_t* completedSteps);

// Deterministic local physics image, with allocation/list handles normalized.
// Useful for exact continuation fixtures; this is not a persistence substitute.
Transport::Status LocalPhysicsImage(std::uint64_t completedSteps,
                                    std::string* image);

} }
#endif
