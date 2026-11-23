#ifndef SRCSEP_PHYSICS_CHECKPOINT_H
#define SRCSEP_PHYSICS_CHECKPOINT_H

#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <cstdint>
#include <functional>
#include <string>

namespace SEP { namespace PhysicsCheckpoint {

// Collective, corrected-model physics checkpoint. Retained fixed boundary:
// fixed native lines, standalone analytic provider, explicit typed coefficient
// models, MPI-only execution, species-global weights and single/species-global
// timesteps, tracker off,
// identical writer/reader rank count and segment ownership. All particle/wave
// workers must have joined. Unsupported capabilities fail before publication.
// Schema five additionally encodes the one-line prescribed reduced-core
// material owner: topology, SI positions/epoch, generation, stable labels,
// local collecting and replicated completed samples, per-segment exposure,
// output counters/buffer phase, source ledgers and complete particle history.
// Material restore allocates a detached native line/table and reconstructs
// core connectivity/cadence before collective publication. A changed MPI
// count is admitted only for the keyed reduced-source consumer: stable labels
// determine new owners, additive measures are summed once and carried on rank
// zero, while completed MPI-reduced numerators/exposure remain single copies.
// This does not repartition arbitrary legacy rank-local RNG source models.
// Fixed-line output continuation and live four-polarization wave authority
// remain unsupported. Earlier corrected schemas are explicitly rejected.
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
// Diagnostic view of fields that MUST be identical on all material replicas.
// This omits additive local measures and particles and cannot be loaded as a
// checkpoint. Native tests retain it to locate an MPI agreement failure without
// relaxing agreement or exposing implementation allocation addresses.
Transport::Status LocalReplicatedPhysicsImage(std::uint64_t completedSteps,
                                              std::string* image);

namespace Diagnostics {
// A rejection-only observation point AFTER detached geometry, ownership,
// source/flux plans and replacement particles have been prepared, but BEFORE
// collective readiness and publication. Counts are values, never mutable
// handles to staged/live physics. This proves rollback of late rank-local
// resource failures, which an early CRC rejection cannot exercise.
struct PreparedReplacement {
  bool material=false;
  std::uint64_t lines=0,particles=0;
};
using PreparationProbe=std::function<Transport::Status(const PreparedReplacement&)>;
// Deliberately separate from the production API: only ControlledDiagnostic
// processes may install a probe. The ordinary LoadCollective has no hook,
// environment switch or configuration option that can bypass validation.
Transport::Status LoadWithPreparationProbe(const std::string& path,
    std::uint64_t* completedSteps,const PreparationProbe& probe);
}

} }
#endif
