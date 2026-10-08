// ============================================================================
// Phase-A AMPS particle-buffer boundary.
//
// This is the only srcSEP3D header (besides the application umbrella) that
// translates the Cartesian AMPS record and linked-cell ABI into the portable
// Adapters::ParticleRecord consumed by the Phase-P cores. Numerical equations
// remain below this boundary and therefore remain unit-testable without AMPS.
// ============================================================================

#ifndef SEP3D_AMPS_PARTICLE_ADAPTER_H
#define SEP3D_AMPS_PARTICLE_ADAPTER_H

#include "pic.h"

#include "../adapters/particle_ledger.h"

// InjectionPlan is intentionally incomplete here.  This header is reachable
// from AMPS-facing declarations, while its concrete definition eventually
// includes canonical SWCME/sep_common headers.  Only the implementation that
// materializes source particles needs that dependency and include path.
namespace SEP3D { namespace Adapters { struct InjectionPlan; } }

namespace SEP3D {
namespace AMPS {
namespace Movers {

// A resolver supplies the complete frozen local environment. It is invoked
// before every accepted transport substep and must sample only the Runtime's
// pinned background/turbulence generations. Keeping this callback explicit
// prevents the particle adapter from reaching into mutable SWMF arrays while
// still allowing a particle to cross cells during one requested AMPS step.
using LocalRecordResolver = Core::Status (*)(
    const Core::Vec3& positionM, int species, double momentumKgMPerS,
    double mu, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    Adapters::LocalTransportRecord* local);

// Population resampling needs only the magnetic direction at a newly
// constructed phase point.  Keeping that lookup separate from the full
// transport resolver avoids evaluating an unrelated D_mumu or mean-free-path
// model (which may legitimately be singular at the provisional pitch angle).
using MagneticDirectionResolver = Core::Status (*)(
    const Core::Vec3& positionM,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node, Core::Vec3* bHat);

struct Context {
  LocalRecordResolver resolveLocal = nullptr;
  MagneticDirectionResolver resolveMagneticDirection = nullptr;
  Adapters::ParticleLedger* ledger = nullptr;
  Adapters::ExpandingShock shock;
  // Hard failure guard for a malformed local limiter.  This is deliberately
  // part of the frozen context instead of a file-scope magic number so it can
  // be fingerprinted and restored with the run configuration.
  std::uint64_t maximumSubsteps = 100000;
};

// Install is legal before particle motion begins. The pointed-to ledger is
// host-owned and must outlive the AMPS run; the other context data are copied.
Core::Status InstallContext(const Context& context);
bool ContextInstalled();
Core::Status UpdateShock(const Adapters::ExpandingShock& shock);

// Request exactly one packed persistent record before AMPS freezes its
// particle-buffer layout. AMPS checkpoint/restart then carries stochastic
// identifiers and pitch state with every particle automatically.
Core::Status RequestParticleStorage();
long int ParticleStateOffset();

// Initialize a newly allocated AMPS particle from the common SWCME source
// adapter. Existing uninitialized particles are rejected by MoveParticle;
// synthesizing a stable ID from an allocation slot would break reproducibility.
Core::Status InitializeParticle(long int ptr,
                                const Adapters::ParticleRecord& particle,
                                double remainingFirstStepFraction = 1.0);
Core::Status ReadParticle(long int ptr, Adapters::ParticleRecord* particle);

struct InjectionOutcome {
  Core::Status status;
  std::uint64_t allocated = 0;
  std::uint64_t rejected = 0;
};

// Materialize a validated R05 plan through AMPS's canonical particle-buffer
// initializer, then append srcSEP3D's persistent stochastic state.  The plan
// is complete before this call, so no source decision depends on allocation
// order or MPI traversal.
InjectionOutcome InjectParticles(const Adapters::InjectionPlan& plan);

// Convert a requested Cartesian unit direction into the (mu,gyrophase)
// coordinates consumed by the maintained guiding-centre mover.  This is the
// inverse of the adapter's deterministic basis construction and is used by
// the reduced-front source to impose anti-sunward launch without bypassing
// persistent particle state.
Core::Status GyrotropicCoordinatesForDirection(
    const Core::Vec3& direction,const Core::Vec3& bHat,
    double* mu,double* gyrophaseRad);

struct PopulationControlRequest {
  std::uint64_t step = 0;
  std::uint64_t campaignSeed = 0;
  unsigned minimumParticlesPerCellPerSpecies = 0;
  unsigned targetParticlesPerCellPerSpecies = 0;
  unsigned maximumParticlesPerCellPerSpecies = 0;
};

struct PopulationControlReport {
  Core::Status status;
  std::uint64_t occupiedCellSpecies = 0;
  std::uint64_t splitOperations = 0;
  std::uint64_t mergeOperations = 0;
  std::uint64_t particlesBefore = 0;
  std::uint64_t particlesAfter = 0;
  double maximumRelativeWeightResidual = 0.0;
  double maximumRelativeMomentumResidual = 0.0;
  double maximumRelativeEnergyResidual = 0.0;
};

// Apply population control to owner-local AMPS cell lists at a joined
// timestep boundary.  Splitting is an exact equal-weight clone; merging is a
// relativistic 3-to-2 reconstruction.  Both paths rewrite srcSEP3D persistent
// identity/gyrotropic/event state so no cloned stable ID or stale momentum can
// enter the next mover call.
PopulationControlReport ApplyPopulationControl(
    const PopulationControlRequest& request);

// The single validating production dispatcher selected by the AMPS mover
// macro. Runtime configuration chooses one of exactly three registered Phase-P
// cores; there are no parallel legacy mover entry points.
int MoveParticle(long int ptr, double dtTotal,
                 cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

}  // namespace Movers
}  // namespace AMPS
}  // namespace SEP3D

#endif  // SEP3D_AMPS_PARTICLE_ADAPTER_H
