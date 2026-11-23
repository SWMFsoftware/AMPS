#ifndef SRCSEP_REDUCED_SHOCK_PARTICLES_H
#define SRCSEP_REDUCED_SHOCK_PARTICLES_H
#include "transport_common.h"
#include "util/sep_reduced_source.h"
#include <array>
namespace SEP { namespace ReducedParticles { namespace Native {
// Pointer-free history follows native migration and exact splitting. Schema
// zero explicitly identifies a legacy/non-reduced source; it is never filled
// in later by guessing a birth energy from an evolved particle's momentum.
struct SourceRecord {std::uint64_t schema=0;double birthKineticEnergyJ=0;};
// Additive rank-local physical accounting. Six (number,birth J,event J)
// triplets are proposed, excluded, born, shock return, domain exit, geometry
// retirement; final entries are signed background/wave work [J]. This owner
// state is checkpointed only after workers and geometry accounting have joined.
// It is not reconstructed from survivors: absorbed particles are already gone.
using LedgerState=std::array<double,20>;
Transport::Status CaptureLedger(LedgerState*);
Transport::Status ValidateLedger(const LedgerState&);
void RestoreLedger(const LedgerState&) noexcept;
bool Enabled();
Transport::Status Inject(const Transport::Birth::DriverWindow&);
long int InjectForNativeDriver();
// Joined geometry retirement is a discrete owner operation. Inventory capture
// and comparison audit its actual removed particles, including birth energy,
// separately from shock returns and particle-mover domain exits.
Transport::Status BeforeGeometry();
Transport::Status AfterGeometry(double epochS);
void TransportWork(const Transport::PICAdapter::ParticleContext&,double beforeMomentum,
    double afterMomentum,double waveWorkPerPhysicalParticleJ=0);
enum class Terminal {ShockReturn,DomainExit,GeometryRetirement};
void RecordTerminal(const Transport::PICAdapter::ParticleContext&,Terminal,double eventMomentum);
} } }
#endif
