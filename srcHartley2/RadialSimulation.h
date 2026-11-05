#ifndef _GENERIC_COMET_RADIAL_SIMULATION_H_
#define _GENERIC_COMET_RADIAL_SIMULATION_H_

#include <string>

template <class T> class cTreeNodeAMR;
namespace PIC { namespace Mesh { class cDataBlockAMR; } }

namespace Comet {
namespace RadialSimulation {

// These are numerical controls, not comet properties. CharacteristicSpeed is
// the maximum resolved transport speed used by the cell-crossing timestep
// criterion; it is indexed by the compiled AMPS species number.
extern bool Enabled;
extern double CourantNumber;
extern double CharacteristicSpeed[];
extern double InnerBoundaryTemperature;
// MaxIterations is a hard safety ceiling. SamplingLength defines the number
// of particle steps in each statistically independent convergence window.
extern long int MaxIterations;
extern long int SamplingLength;
// LoadBalanceInterval redistributes complete AMPS blocks independently of the
// longer science-sampling interval.  This prevents a radially advancing
// population front from exhausting one MPI rank while the global population
// is still well within the distributed particle-buffer capacity.
extern long int LoadBalanceInterval;
// When enabled, the population at each completed sampling break is measured
// independently for every species and owner block.  The block-local weight is
// changed by an integral power of two until the expected mean population per
// cell lies in this inclusive range.  Empty species/block populations remain
// empty because adaptation only resamples particles that already exist.  Once
// populated blocks have been adapted, their W/dt carrier constraints are
// propagated by power-of-two weight increases through face-connected empty
// neighbors.  This prevents the first particle entering a previously empty
// block from being multiplied by an obsolete, anomalously small local weight;
// the propagation itself never creates a particle.
extern bool ParticleWeightAdaptationEnabled;
// A positive value requests regular adaptation at this many completed model
// iterations, independently of SamplingLength.  The source default -1 is the
// legacy contract: adapt only at complete sampling breaks.  Both paths run
// after PIC::TimeStep has completed transport, chemistry, collisions, and all
// continuation exchanges, so no in-flight particle is ever split or removed.
extern long int ParticleWeightAdaptationInterval;
extern long int MinimumParticlesPerCell;
extern long int MaximumParticlesPerCell;
// Stage-1 input resolution produces one immutable policy row per compiled
// species.  The legacy scalar range remains the global fallback and keeps old
// input behavior unchanged.  Locked species ignore regular population and
// photochemical weight adaptation; emergency common-factor coarsening remains
// a safety override and is reported separately when it occurs.
extern long int SpeciesMinimumParticlesPerCell[];
extern long int SpeciesMaximumParticlesPerCell[];
extern int ParticleWeightLocked[];
extern int PopulationPolicySource[];
extern int PopulationClass[];

namespace PolicySource {
enum Type {Global=0,DirectSource=1,ChemicalOnly=2,ExplicitSpecies=3};
}
namespace PolicyClass {
enum Type {Other=0,DirectSource=1,ChemicalOnly=2};
}
// Emergency roulette runs immediately before the regular particle exchange
// when either (a) any MPI rank reaches this fraction of its local particle-
// buffer capacity or (b) the preceding sampling break found an owner block
// whose all-species mean population per cell exceeded
// 10*nTotalSpecies*MaximumParticlesPerCell.  The input value controls trigger
// (a) and is a fraction in (0,1), e.g. 0.75 for 75%.
extern double EmergencyParticleBufferUtilization;
// Transactional allocation control stops ordinary work at the high-water
// value above, coarsens collectively until the maximum rank is at or below
// AllocationLowWatermark, and keeps AllocationReserveFraction of every rank's
// particle slots unavailable to application-level chemistry transactions.
extern double AllocationLowWatermark;
extern double AllocationReserveFraction;
// Maximum actual destination multiplicity at any directed block face.  A
// value of one removes mover-side clone allocation entirely by equalizing the
// W/dt carrier rate over each face-connected mesh component.
extern long int MaximumFaceTransferMultiplicity;
// A chemistry event may create several fixed-weight model daughters when the
// parent and product carrier rates W/dt differ.  At a sampling break, every
// reaction transitively reachable from any represented parent is inspected.
// Product weights—whether that product is present or absent—are raised by
// powers of two until the sum of worst-case residual-resampling allocations
// for one selected reaction event does not exceed this hard limit.  The final
// check uses face-propagated carrier rates, so a parent arriving from another
// block cannot bypass the bound.  Empty cells remain empty.
extern long int MaximumPhotochemicalProductsPerReaction;
// ConvergenceTolerance applies to the volume-weighted relative L1 change of
// each species' complete radial number-density profile. The run stops only
// after every species passes for RequiredConvergedWindows successive windows.
extern double ConvergenceTolerance;
extern int RequiredConvergedWindows;

bool Validate(std::string *error=NULL);

// Complete the AMPS initialization that follows radial mesh construction,
// install the spherical source and boundary callbacks, and advance the run.
void Initialize();
void Run();

// Run an integrated split/roulette regression on initialized radial blocks.
// This is a focused validation mode: it replaces the live particle lists and
// must therefore be invoked instead of, never during, a science run.
bool ValidateParticleWeightAdaptation();

// Collective callback installed in PIC::Parallel.  It is public only because
// the core stores its address; application code must not call it directly.
void EmergencyParticleWeightAdaptation();
// Called collectively by transactional allocation stages.  It records a
// one-shot request and invokes the same common-factor coarsening used by the
// pre-MPI-exchange hook, but continues to AllocationLowWatermark so the
// deferred stage can resume with deterministic headroom.
void RequestAllocationBurstEmergency();

// The following functions match AMPS callback signatures.  They are public so
// the core can store function pointers to them; case code should normally call
// only Validate, Initialize, and Run.
double LocalTimeStep(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
double LocalBlockInjectionRate(int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
bool InjectionBlockIndicator(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
long int InjectBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
int ProcessDomainBoundary(long int ptr,double *position,double *velocity,
    int intersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

}
}

#endif
