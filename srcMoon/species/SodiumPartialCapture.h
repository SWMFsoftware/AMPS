#ifndef SRCMOON_SPECIES_SODIUM_PARTIAL_CAPTURE_H
#define SRCMOON_SPECIES_SODIUM_PARTIAL_CAPTURE_H

#include "ReducedKinematics.h"
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace SodiumPartialCapture {
namespace K=ReducedKinematics;
using Real=long double;
struct Constants {
  Real atomicMassUnitKg=0,electronVoltJ=0;
  std::string sourceIdentity;
};
struct Node {Real energyEvPerU=0,crossSectionM2=0;};
class Rate {
 public:
  bool valid()const{return valid_;}
  Real rateS()const{return rateS_;}
  Real relativeSpeedMS()const{return relativeSpeedMS_;}
  Real crossSectionM2()const{return crossSectionM2_;}
  Real energyEvPerU()const{return energyEvPerU_;}
  Real protonDensityM3()const{return protonDensityM3_;}
  const std::string& sourceIdentity()const{return sourceIdentity_;}
  const Constants& constants()const{return constants_;}
 private:
  friend class Registry;
  // Only a loaded source registry can create a valid receipt. Immutable
  // values bind density, constants and native node together; an evidence
  // string alone cannot authorize edited cross sections, speeds or rates.
  bool valid_=false;
  Real rateS_=0,relativeSpeedMS_=0,crossSectionM2_=0,energyEvPerU_=0,protonDensityM3_=0;
  std::string sourceIdentity_;
  Constants constants_;
};

// Measured GROUND-PRODUCT PARTIAL capture channel, not the Morgan total:
// H+ + Na -> H(1s) + Na+(2p6). The native label leaves the initial Na state
// unspecified; Na(3s) below is the explicitly approved REFERENCE cohort,
// not measured initial-state metadata. That crosswalk remains unqualified.
// The exact acquired CSV must match its
// frozen numerical hash. Only the 12 native 4--25keV/u nodes are admitted;
// interpolation/extrapolation and thermal-distribution tail handling are not
// supplied by this deliberately bounded first P07 reference experiment.
class Registry {
 public:
  bool LoadNativeNodes(const std::string& path,std::string* error);
  bool RateAtNativeNode(Real energyEvPerU,Real protonDensityM3,
      const Constants&,Rate* output,std::string* error)const;
  const std::vector<Node>& nodes()const{return nodes_;}
  bool loaded()const{return !nodes_.empty();}
  bool physicalTargetStateQualified()const{return false;}
  const std::string& sourceIdentity()const{return source_;}
 private:
  std::vector<Node> nodes_;
  std::string source_;
};

struct ForwardInput {
  K::ResolvedState incomingProton,incomingSodium,outgoingHydrogen,outgoingSodiumIon;
  K::Vector3 relativeIncidentDirection{{0,0,0}};
  Rate sourceRate;Constants constants;
  K::Context context;K::NumericControls controls;
};

// SOURCE TARGET-REST frame: incoming Na is at rest. Relative speed is the
// rate receipt's source energy/speed mapping, not a separately guessed wind
// speed. Use source-explicit state energies to derive physical Q. All heavy
// products are returned by the same general production two-body solver as
// photon recoil. No projectile reservoir or PIC particle is mutated here.
// Forward COM scattering is the approved EXPERIMENTAL closure, not an angular
// measurement. The caller must account for the prescribed incoming proton.
bool ForwardProducts(const Registry&,const ForwardInput&,K::FinalState* output,std::string* error);
} } }
#endif
