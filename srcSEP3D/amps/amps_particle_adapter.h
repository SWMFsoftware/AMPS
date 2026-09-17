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

namespace SEP3D {
namespace AMPS {
namespace Movers {

// A resolver supplies the complete frozen local environment. It is invoked
// once at the start of an AMPS mover call and must sample only the Runtime's
// pinned background/turbulence generations. Keeping this callback explicit
// prevents the particle adapter from reaching into mutable SWMF arrays.
using LocalRecordResolver = Core::Status (*)(
    const Core::Vec3& positionM, int species, double momentumKgMPerS,
    double mu, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    Adapters::LocalTransportRecord* local);

struct Context {
  LocalRecordResolver resolveLocal = nullptr;
  Adapters::ParticleLedger* ledger = nullptr;
  Adapters::ExpandingSphericalShock shock;
};

// Install is legal before particle motion begins. The pointed-to ledger is
// host-owned and must outlive the AMPS run; the other context data are copied.
Core::Status InstallContext(const Context& context);

// Request exactly one packed persistent record before AMPS freezes its
// particle-buffer layout. AMPS checkpoint/restart then carries stochastic
// identifiers and pitch state with every particle automatically.
Core::Status RequestParticleStorage();
long int ParticleStateOffset();

// Initialize a newly allocated AMPS particle from the common SWCME source
// adapter. Existing uninitialized particles are rejected by MoveParticle;
// synthesizing a stable ID from an allocation slot would break reproducibility.
Core::Status InitializeParticle(long int ptr,
                                const Adapters::ParticleRecord& particle);

// The single validating production dispatcher selected by the AMPS mover
// macro. Runtime configuration chooses one of exactly two registered Phase-P
// cores; there are no parallel legacy mover entry points.
int MoveParticle(long int ptr, double dtTotal,
                 cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

}  // namespace Movers
}  // namespace AMPS
}  // namespace SEP3D

#endif  // SEP3D_AMPS_PARTICLE_ADAPTER_H
