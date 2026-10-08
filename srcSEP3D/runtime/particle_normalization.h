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
#include "provider.h"

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

// One entry retains the physical triangle identity and its cumulative gross
// incident population rate.  Only SolvedFastShock records contribute:
//   Ndot_f,s = n_1,s (V_sh,n-U_1.n) A_f.
// A_f is the provider's exact curved quadrature area, while particle positions
// are sampled on the corresponding planar triangle carried by the epoch.
struct SurfaceFaceParticleRate {
  std::size_t triangleIndex = 0;
  std::uint64_t stableId = 0;
  double physicalRatePerS = 0.0;
  double cumulativeRatePerS = 0.0;
  double compressionRatio = 0.0;
};

struct SurfaceParticleRateDistribution {
  std::uint64_t generation = 0;
  double epochS = 0.0;
  double acceptedAreaM2 = 0.0;
  double physicalRatePerS = 0.0;
  std::vector<SurfaceFaceParticleRate> faces;
};

// A dependency-light candidate event. All MPI ranks build the same ordered
// list from semantic random keys. Native code subsequently allocates an event
// only on the rank owning its sampled AMR point.
struct SurfaceInjectionEvent {
  std::uint64_t stableId = 0;
  std::uint64_t triangleStableId = 0;
  std::size_t triangleIndex = 0;
  double eventTimeS = 0.0;
  double remainingStepFraction = 0.0;
  double momentumKgMPerS = 0.0;
  double compressionRatio = 0.0;
  Core::Vec3 positionM;
};

struct SurfaceInjectionBatch {
  SurfaceParticleRateDistribution distribution;
  std::vector<SurfaceInjectionEvent> events;
};

Core::Status BuildSurfaceParticleRateDistribution(
    const SEP::CoronaSwcme::ShockFront::Provider& provider,
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const CompiledSpeciesRecord& species,
    SurfaceParticleRateDistribution* distribution);

// Constant-W production sampler. Waiting times are exponential with
// lambda=Ndot_total/W_s; conditional face probability is Ndot_f/Ndot_total.
// Thus each face has the correct independent Poisson intensity without
// changing a particle's species base weight. The finite maximum is a fatal
// guard, not a truncation or renormalization.
Core::Status GenerateConstantWeightSurfaceInjectionBatch(
    const SEP::CoronaSwcme::ShockFront::Provider& provider,
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const CompiledSpeciesRecord& species,
    const SourceOptions& source,double macroparticleWeight,double intervalS,
    std::uint64_t campaignSeed,std::uint64_t step,
    SurfaceInjectionBatch* batch);

// Reserved second representation. It is intentionally callable so the input
// selection reaches a typed boundary, but it must not return particles until
// its individual-weight correction and conservation tests are implemented.
Core::Status GenerateLogUniformMomentumImportanceBatch(
    const SEP::CoronaSwcme::ShockFront::Provider& provider,
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const CompiledSpeciesRecord& species,
    const SourceOptions& source,double macroparticleWeight,double intervalS,
    std::uint64_t campaignSeed,std::uint64_t step,
    SurfaceInjectionBatch* batch);

} } // namespace SEP3D::RuntimeModel

#endif
