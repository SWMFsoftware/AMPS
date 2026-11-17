#ifndef SRCSEP_UTIL_SEP_RESONANT_COLLISION_H
#define SRCSEP_UTIL_SEP_RESONANT_COLLISION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP { namespace Transport { namespace ResonantCollision {

enum class Branch { AlongB=0, AgainstB=1 };
// Fourier helicity convention: omega-k*v_parallel+eta*Omega_signed=0,
// eta=+1 for Right and -1 for Left. Signed particle charge selects resonance.
enum class Helicity { Left=0, Right=1 };
enum class DomainAuthority { Unspecified, ExplicitDeclaredLimits, ValidationStudy };

struct DomainContract {
  DomainAuthority authority=DomainAuthority::Unspecified;
  // No numerical physical-domain limits are assumed. A qualified production
  // contract supplies both limits and the provenance of their qualification.
  double maximumDeltaBOverB=0.0;
  double maximumKThermalIonGyroradius=0.0;
  double thermalIonGyroradiusM=0.0;
  double maximumAlfvenSpeedOverC=0.0;
  double ionInertialLengthM=0.0;
  double maximumKIonInertialLength=0.0;
  double ionCyclotronFrequencyRadPerS=0.0;
  double maximumWaveFrequencyOverIonCyclotron=0.0;
  std::string plasmaCompositionIdentity;
  std::string thermalStateIdentity;
  std::string provenance;
};

struct Spectrum {
  // Explicit positive bin edges [1/m], N+1 entries for N bins. Arrays contain
  // MAGNETIC cell-integrated joules, indexed 2*branch+helicity. Unpolarized
  // two-array storage cannot supply this capability without initial helicity.
  std::vector<double> kEdgesPerM;
  std::array<std::vector<double>,4> magneticEnergyJ;
  bool polarizationResolved=false;
  double volumeM3=0.0;
  double meanMagneticFieldT=0.0;
  double alfvenSpeedMPerS=0.0; // Positive magnitude, branch supplies sign.
  double vacuumPermeabilityHPerM=0.0;
  double epochS=0.0;
  std::uint64_t generation=0;
  std::string identity;
  DomainContract domain;
};

Status ValidateSpectrum(const Spectrum& spectrum);
Status ValidateDomain(const Spectrum& spectrum);
// Current implementation capability, independent of user-provided numbers.
// Aggregate-only legacy production paths cannot claim the approved polarized
// resonant closure merely because a neutral validation kernel is available.
Status RequireImplementedProductionClosure(bool particleWaveExchangeRequested,
                                           bool localWaveEvolutionRequested);

struct Particle {
  double massKg=0.0;
  double chargeC=0.0;
  // Both momentum components are in the plasma frame, relative to local B.
  double parallelMomentumKgMPerS=0.0;
  double perpendicularMomentumKgMPerS=0.0;
};

struct Orbit {
  Branch branch=Branch::AlongB;
  double massKg=0.0;
  double chargeC=0.0;
  double waveSpeedMPerS=0.0;
  double waveGamma=1.0;
  double energyInWaveFrameJ=0.0;
  double momentumInWaveFrameKgMPerS=0.0;
  double invariantJ=0.0; // C=E_plasma-u*p_parallel=E_wave/Gamma.
  double initialWavePitch=0.0;
  bool stationaryInWaveFrame=false;
};

struct OrbitResult {Status status;Orbit orbit;};
OrbitResult PrepareOrbit(const Spectrum& spectrum,Branch branch,const Particle& particle);
struct ParticleResult {Status status;Particle particle;double totalEnergyJ=0.0;};
ParticleResult ToPlasmaFrame(const Orbit& orbit,double wavePitch);

struct Resonance {
  Status status;
  bool active=false;
  std::string zeroReason;
  Helicity helicity=Helicity::Left;
  std::size_t bin=0;
  double kPerM=0.0;
  double magneticPowerT2M=0.0; // Integral P_h dk = real variance of this helicity.
  // Plasma (p,mu_B) FP parameter. D_zetazeta in plasma time is
  // (1-zeta^2)*nu/(2*Gamma_wave^2), not (1-zeta^2)*nu/2.
  double frequencyPerS=0.0;
};
Resonance EvaluateResonance(const Spectrum& spectrum,const Orbit& orbit,double wavePitch);

// On one branch the diffusion tensor in (mu_B,p) is rank one and tangent to
// C=E-u*p*mu_B. In (C,zeta=mu_wave), its particle measure is J dC dzeta,
// J=Gamma*E_plasma*p_wave/c^2, and J*D_zetazeta=
// (pi/2)*q^2*P_h(k_res)*(1-zeta^2)/abs(zeta), where integral P_h dk
// is the actual real magnetic variance of one circular helicity. The 1/2
// comes from random-phase force autocorrelation, not magnetic/kinetic energy.
// This form is regular for v<v_A;
// (C,mu_plasma) instead folds there and is not a global coordinate chart.
struct Generator {
  Spectrum spectrum;
  Orbit orbit;
  std::vector<double> edges;
  std::vector<double> centers;
  std::vector<double> measure;
  std::vector<double> conductance; // N-1 shared interior-face conductances.
  std::vector<double> leftRatePerS,rightRatePerS;
  bool qualifiedForProduction=false;
  bool withinDeclaredLimits=false;
};
struct GeneratorResult {Status status;Generator generator;};

// Add every mapped physical k edge and both finite-band gap boundaries to a
// uniform pitch mesh. Conductances integrate the exact resistance between
// cell centers. A positive-width region of zero resonant power disconnects
// the mesh; no pitch-angle bridge or numerical coefficient floor is inserted.
GeneratorResult BuildGenerator(const Spectrum& spectrum,const Orbit& orbit,
                              std::size_t uniformPitchCells);

// Forward probability equation on cell masses. Its shared-face flux is
// G*(prob_i/measure_i-prob_j/measure_j), so mass is conserved and equilibrium
// is proportional to the exact phase-space measure. Rates are nonnegative.
Status ProbabilityDerivative(const Generator& generator,
    const std::vector<double>& probability,std::vector<double>* derivative);

struct Jump {
  Status status;
  bool active=false;
  double waitingTimeS=0.0;
  std::size_t targetCell=0;
  double targetWavePitch=0.0;
};
Jump SampleJump(const Generator& generator,std::size_t currentCell,RandomStream& random);
ScalarResult TotalRate(const Generator& generator,std::size_t currentCell);
// Draw marks only after an event's integrated hazard has expired. Drivers
// carry residual optical depth across geometry, snapshot and outer-step cuts.
Jump SampleJumpAtEvent(const Generator& generator,std::size_t currentCell,
                       RandomStream& random);
std::size_t CellAt(const Generator& generator,double wavePitch);

struct BinExchange {
  Branch branch=Branch::AlongB;
  Helicity helicity=Helicity::Left;
  std::size_t bin=0;
  double particleEnergyChangeJ=0.0;
  double magneticEnergyChangeJ=0.0;
};
struct Proposal {
  Status status;
  Particle before,after;
  double beforeWavePitch=0.0,afterWavePitch=0.0;
  double particleEnergyChangeJ=0.0;
  double parallelMomentumChangeKgMPerS=0.0;
  std::vector<BinExchange> exchange;
  std::string spectrumIdentity;
  std::uint64_t generation=0;
  double epochS=0.0;
};
// An accepted neighboring-cell jump is only a private proposal. Macro weight
// scales both sides; total wave energy is 2*magnetic energy. The coupled
// transaction must approve budgets before publishing either particle or waves.
Proposal ProposeJump(const Generator& generator,double beforeWavePitch,
                    double afterWavePitch,double statisticalWeight);
// Revalidate a mutable wire/scratch proposal against its immutable spectrum,
// selected branch and macro weight. This performs no random draw, phase edit
// or budget publication. The deterministic elastic/resonance partition is
// shared with ProposeJump. It validates physical bookkeeping, not the event's
// stochastic selection or membership in a particular pitch-grid generator;
// those identities belong to the enclosing event/particle transaction.
Status ValidateProposal(const Spectrum& spectrum,Branch branch,
                        const Proposal& proposal,double statisticalWeight);

} } }
#endif
