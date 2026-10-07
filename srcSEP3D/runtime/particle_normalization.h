// ============================================================================
// Mesh/shock-derived global particle numerics for srcSEP3D
//
// This dependency-light layer contains the physical equations used at the
// AMPS boundary but no PIC or MPI symbols. Tests can therefore distinguish a
// defect in the rate/weight algebra from a defect in native mesh reduction.
// ============================================================================

#ifndef SEP3D_RUNTIME_PARTICLE_NORMALIZATION_H
#define SEP3D_RUNTIME_PARTICLE_NORMALIZATION_H

#include "run_configuration.h"

#include "diagnostics.h"

#include <cstdint>
#include <vector>

namespace SEP3D { namespace RuntimeModel {

// One global step limits a model particle moving at v_max to the configured
// fraction of the smallest allocated characteristic cell length:
//   dt = margin * h_min / v_max.
// All values are SI. No species-dependent step is introduced by this model.
Core::Status CalculateMeshGlobalTimeStep(
    double minimumCellSizeM,double maximumParticleSpeedMPerS,
    double marginFactor,double* timeStepS);

// AMPS schedules observations on integer global ticks. Preserve the exact CFL
// step above and move each requested physical cadence to the first tick at or
// after it. Thus output is never made more frequent than requested, and the
// existing exact integer-cadence validation remains fully active.
Core::Status AlignObserverCadencesToGlobalStep(
    double timeStepS,std::vector<ObserverOptions>* observers);

// Match each generated AMPS species to an explicitly available upstream
// plasma population (electron, proton, or alpha), then apply
//   W_s = Ndot_s * dt / N_model,s.
// Unknown ions and zero-abundance populations fail closed: silently assigning
// the proton rate would change composition and represented charge/mass.
Core::Status CalculateSpeciesParticleNormalizations(
    const std::vector<CompiledSpeciesRecord>& species,
    const SEP::CoronaSwcme::ShockFront::IncidentParticleFlux& flux,
    double timeStepS,std::uint64_t particlesPerIteration,
    std::vector<SpeciesParticleNormalization>* normalizations);

} } // namespace SEP3D::RuntimeModel

#endif
