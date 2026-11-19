#ifndef SRCMOON_SPECIES_OXYGEN_CASCADE_H
#define SRCMOON_SPECIES_OXYGEN_CASCADE_H
#include "../../src/models/chemical_database/Runtime.h"
#include <array>
#include <optional>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace OxygenCascade {
using Real=long double;
using Vector3=std::array<Real,3>;
using Status=AMPS::ChemicalDatabase::Status;
using Catalogue=AMPS::ChemicalDatabase::Catalogue;
enum class EnergyBasis { Observed, Theoretical };
struct Level {
  unsigned index=0;
  std::string configuration,multiplicity,L,J,observedNative,theoreticalNative,rowSha256;
  Real observedCmInverse=0,theoreticalCmInverse=0;
};
struct Link {
  unsigned lower=0,upper=0,nativeLine=0;
  Real rateS=0,signedVacuumWavelengthAngstrom=0,branchProbability=0;
  std::string rateNative,wavelengthNative,rowSha256,uncertaintyStatus;
};
struct DecayLaw {
  bool stableWithinSelectedNetwork=false;
  Real totalRateS=0;
  std::optional<Real> lifetimeS;
  std::vector<std::size_t> links;
  std::vector<Real> branchProbabilities;
};
struct PhotonConstants {
  // Explicit SI constants and source identities. No implicit units conversion
  // or missing-constant fallback is selected by this service.
  Real planckJS=0,speedOfLightMS=0;
  std::string planckSource,lightSource;
};
struct StateGap {
  unsigned lower=0,upper=0;
  EnergyBasis basis=EnergyBasis::Observed;
  Real lowerEnergyJ=0,upperEnergyJ=0,releasedStateEnergyJ=0;
  Real signedNativeVacuumWavelengthAngstrom=0;
  Real planckJS=0,speedOfLightMS=0;
  std::string lowerNativeCmInverse,upperNativeCmInverse;
  std::string levelTableSha256,linkTableSha256,planckSource,lightSource;
  // releasedStateEnergyJ is NOT automatically the emitted photon energy:
  // finite emitter recoil requires independently supplied physical masses.
};
struct UniformPair {Real waiting=0,branch=0;};
struct DecayEvent {Real offsetS=0;std::size_t link=0;};
struct History {
  unsigned initialLevel=0,finalLevel=0;
  Real durationS=0;
  std::size_t waitingDrawsUsed=0,branchDrawsUsed=0;
  std::array<unsigned,9> photonCounts{};
  std::vector<DecayEvent> events;
};
class Network {
 public:
  // Only the reviewed complete selected five-level/nine-link source is
  // accepted. The complete native originals are also checked. Transactional;
  // source access does not admit an isotope/state/quenching physical profile.
  Status LoadFrozenSelectedChiantI(const Catalogue& catalogue);
  bool loaded()const{return loaded_;}
  const std::vector<Level>& levels()const{return levels_;}
  const std::vector<Link>& links()const{return links_;}
  const AMPS::ChemicalDatabase::BundleIdentity& identity()const{return identity_;}
  Status GetDecayLaw(unsigned level,DecayLaw* output)const;
  Status SelectDecay(unsigned level,Real unitUniform,std::size_t* link)const;
  Status WaitingTime(unsigned level,Real openUniform,Real* seconds)const;
  Status GetStateGap(std::size_t link,EnergyBasis basis,const PhotonConstants&,StateGap*)const;
  // Mathematical radiative-only evolution with supplied uniforms, not a
  // claim that collisional/upward/quenching processes vanish on the Moon.
  // Native lower endpoints strictly decrease: at most four decays occur.
  Status AdvanceRadiativeOnly(unsigned initialLevel,Real durationS,
      const std::vector<UniformPair>& draws,History* output)const;
 private:
  bool loaded_=false;
  std::vector<Level> levels_;
  std::vector<Link> links_;
  std::array<DecayLaw,5> laws_{};
  AMPS::ChemicalDatabase::BundleIdentity identity_;
};
struct Bound {Real absolute=-1,relative=-1;};
struct RecoilInput {
  StateGap stateGap;
  PhotonConstants constants;
  Real physicalParentMassKg=0,physicalDaughterMassKg=0;
  std::string parentMassSource,daughterMassSource,isotopeIdentity,frameIdentity;
  Vector3 photonDirection{};
  Bound massEnergyConsistencyJ,energyJ,momentumKgMS;
  Real unitDirectionTolerance=-1;
};
struct RecoilLedger {
  StateGap stateGap;
  PhotonConstants constants;
  Real physicalParentMassKg=0,physicalDaughterMassKg=0;
  std::string frameConvention="parent-rest-inertial";
  Real photonEnergyJ=0,photonMomentumMagnitudeKgMS=0,atomicKineticEnergyJ=0;
  Real stateEnergyReleaseJ=0,massEnergyReleaseJ=0,energyResidualJ=0,massEnergyResidualJ=0;
  Vector3 photonMomentumKgMS{},atomicMomentumKgMS{},atomicVelocityMS{};
  std::string isotopeIdentity,frameIdentity,parentMassSource,daughterMassSource;
};
// Exact two-body relativistic rest-emitter relation, with stable recoil KE.
// Caller-supplied physical masses must match the declared state gap within
// explicit arithmetic bounds. All outputs remain unchanged on every refusal.
// This mathematical evaluator does not approve a physical application profile.
Status SolveRestFramePhotonRecoil(const RecoilInput&,RecoilLedger* output);
}}}
#endif
