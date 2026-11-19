#ifndef SRCMOON_SPECIES_REDUCED_KINEMATICS_H
#define SRCMOON_SPECIES_REDUCED_KINEMATICS_H

#include "ReactionExecution.h" // identity structs only; no executable dependency
#include <array>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace ReducedKinematics {
using Real=long double;
using Vector3=std::array<Real,3>;
constexpr Real SpeedOfLightMS=299792458.0L; // exact SI definition, not a cutoff
struct MassValue {Real kg=0;std::string sourceIdentity;};
struct IsotopeMass {std::string identity;MassValue neutralAtom;};
struct EnergyValue {Real joules=0;std::string sourceIdentity;};
struct StateDefinition {
  std::string identity,state;
  std::vector<ReactionExecution::NucleusCount> nuclei;
  std::int64_t elementaryCharge=0;
  EnergyValue formation,excitation;
};

// Only source-explicit resolution constructs a usable state. No mutable
// mass/energy/state setters exist, and no native species-name inference occurs.
class ResolvedState {
 public:
  bool valid()const{return valid_;}
  Real referenceMassKg()const{return mass_;}
  const StateDefinition& definition()const{return definition_;}
  const std::vector<std::string>& sourceIdentities()const{return sources_;}
  const std::string& referenceInventoryIdentity()const{return inventory_;}
 private:
  friend bool ResolveState(const StateDefinition&,const std::vector<IsotopeMass>&,
      const MassValue&,ResolvedState*,std::string*);
  bool valid_=false;Real mass_=0;StateDefinition definition_;
  std::vector<std::string> sources_;
  std::string inventory_;
};
bool ResolveState(const StateDefinition& definition,
    const std::vector<IsotopeMass>& isotopes,const MassValue& electron,
    ResolvedState* output,std::string* error);
struct EnergyOccurrence {
  bool incoming=false;std::size_t occurrence=0;
  std::string identity,state;
  EnergyValue formation,excitation;
};
struct FormationBalance {
  Real incomingStateEnergyJ=0,outgoingStateEnergyJ=0,qJ=0;
  std::vector<std::string> sourceIdentities;
  std::vector<EnergyOccurrence> terms;
};
// Include each occurrence explicitly; photons/baths are separate terms.
bool FormationQ(const std::vector<ResolvedState>& incoming,
    const std::vector<ResolvedState>& outgoing,FormationBalance* output,std::string* error);

struct ErrorLimits {
  Real maximumRelativeInertialMassError=-1;
  Real maximumRelativeAccelerationError=-1;
  Real maximumRelativeRelativisticKineticError=-1;
};
struct ModelErrorReport {
  Real relativeInertialMassError=0,relativeAccelerationError=0;
  Real relativeRelativisticKineticError=0,speedMS=0;
  std::vector<std::string> sourceIdentities;
};
bool EvaluateModelErrors(const ResolvedState& reference,const MassValue& physicalRestMass,
    Real speedMS,const ErrorLimits& limits,ModelErrorReport* output,std::string* error);

// Numerical bounds are caller inputs, not physical uncertainty budgets. All
// fields have invalid defaults so omission cannot silently choose a tolerance.
struct Bound {Real absolute=-1,relative=-1;};
struct NumericControls {
  // The unit-vector bound must lie in[0,1). Accepted encodings are normalized;
  // this bound never changes
  // the magnitude of photon momentum or prescribed relative kinetic energy.
  Real maximumClassicalSpeedMS=-1,unitVectorTolerance=-1;
  Bound referenceMassKg,momentumKgMS,energyJ,thresholdJ;
};
struct Context {std::string frameIdentity,sourceIdentity;};
struct Product {
  // occurrence is the input-array position. A later stoichiometric adapter
  // must retain its own explicit term/occurrence mapping (including equal H).
  ResolvedState state;std::size_t occurrence=0;
  Vector3 momentumKgMS{{0,0,0}},velocityMS{{0,0,0}};
};
// Shared production isotropic sampler. Both photon kernels and event adapters
// use exactly these two supplied open uniforms; no internal RNG or retry draw
// exists. The stable sine retains transverse components for tiny positive u.
bool IsotropicDirectionFromOpenUniforms(Real polarOpenUniform,
    Real azimuthOpenUniform,Vector3* output,std::string* error);
struct FinalState {
  std::string modelIdentity,frameIdentity;
  std::vector<Product> products;
  Vector3 totalMomentumKgMS{{0,0,0}},momentumResidualKgMS{{0,0,0}};
  Real availableCOMKineticJ=0,actualCOMKineticJ=0,totalKineticJ=0;
  Real energyResidualJ=0,maximumProductSpeedMS=0;
  bool photonIncluded=false;
  Real photonEnergyJ=0,stateQJ=0;
  // Retain the source energy components separately: a small supplied excess
  // must not vanish when displayed beside a much larger threshold.
  Real channelThresholdJ=0,sourceMeanExcessJ=0,stateAvailableKineticJ=0;
  Vector3 photonMomentumKgMS{{0,0,0}};
  std::vector<std::string> sourceIdentities;
  FormationBalance formationBalance;
};
struct TwoBodyInput {
  std::array<ResolvedState,2> products;
  Vector3 totalMomentumKgMS{{0,0,0}},relativeDirection{{0,0,0}};
  Real availableCOMKineticJ=-1;
  Context context;NumericControls controls;
};
// General SAME production solver used by PhotonTwoBody and future explicit
// CEX wrappers. It chooses no incoming law, Q, scattering angle or reservoir.
bool GeneralTwoBodyFinalState(const TwoBodyInput& input,FinalState* output,std::string* error);

struct PhotonTwoBodyInput {
  ResolvedState parent;
  std::array<ResolvedState,2> products;
  Vector3 parentMomentumKgMS{{0,0,0}},photonDirection{{0,0,0}};
  EnergyValue channelThreshold,meanExcess;
  Real polarOpenUniform=0,azimuthOpenUniform=0;
  Context context;NumericControls controls;
};
// PARENT-REST inertial frame only: nonzero parent momentum is refused.
// Complete products, reference material mass/charge, source-state threshold,
// photon momentum and all energy terms are checked before output publication.
bool PhotonTwoBody(const PhotonTwoBodyInput& input,FinalState* output,std::string* error);

struct ThreeBodyInput {
  std::array<ResolvedState,3> products;
  Vector3 totalMomentumKgMS{{0,0,0}};
  Real availableCOMKineticJ=-1;
  std::array<Vector3,3> standardNormals{};
  Context context;NumericControls controls;
};
// Uniform classical microcanonical measure in mass-weighted coordinates.
// Exactly three neutral occurrences; finite supplied draws, no internal RNG.
bool ThreeNeutralMicrocanonical(const ThreeBodyInput& input,FinalState* output,std::string* error);
}}}
#endif
