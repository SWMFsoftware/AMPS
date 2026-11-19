#ifndef SRCMOON_SPECIES_ION_COLLECTOR_H
#define SRCMOON_SPECIES_ION_COLLECTOR_H

#include "ReactionExecution.h"
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace IonCollector {
namespace R=ReactionExecution;
using Vector3=std::array<double,3>;
using ExtendedVector3=std::array<long double,3>;

// Explicit model-state energies [J per incident physical particle]. Formation
// and excitation use the SAME named energy zero. NaN defaults deliberately
// refuse an unavailable value; a declared zero must actually be supplied.
// This boundary invents no source state and does not certify the evidence.
struct StateDefinition {
  R::SpeciesIdentity material;
  double formationEnergyJ=std::numeric_limits<double>::quiet_NaN();
  double excitationEnergyJ=std::numeric_limits<double>::quiet_NaN();
  std::string energyReferenceIdentity,sourceEvidenceIdentity;
};

// Actual selected impact geometry, not a ray-preview or a source-time pose.
// LSO velocity coordinates and surface velocity must use the same convention.
// In frozen-environment runs geometryEpochEtS differs from advancing transport
// epoch; retain both rather than labelling the frozen pose as evolving.
struct ImpactFrame {
  Vector3 outwardNormalLso{{0,0,0}},surfaceVelocityLsoMS{{0,0,0}};
  std::array<double,36> lsoToBody36{},bodyToLso36{};
  bool stateTransformsSupplied=false;
  std::string geometryIdentity;
  double geometryEpochEtS=std::numeric_limits<double>::quiet_NaN();
};
struct Incident {
  StateDefinition state;
  double representedNumber=0; // actual W_node,s * individual correction, no W/dt
  double elementaryChargeC=0; // explicit production constant [C]
  Vector3 positionLsoM{{0,0,0}},velocityLsoMS{{0,0,0}};
  double impactEpochEtS=std::numeric_limits<double>::quiet_NaN();
  ImpactFrame frame;
};
struct RepresentedNucleus {
  std::string identity;
  std::int64_t multiplicity=0; // exact integer identity per incident particle
  long double representedNumber=0;
};
struct Record {
  Incident incident; // retains exact species/state/frame and physical epoch
  std::vector<RepresentedNucleus> nuclei;
  long double representedChargeElementary=0,chargeC=0,modelMassKg=0;
  ExtendedVector3 momentumFacetKgMS{{0,0,0}},momentumLsoKgMS{{0,0,0}};
  long double kineticEnergyFacetJ=0,kineticEnergyLsoJ=0;
  long double formationEnergyJ=0,excitationEnergyJ=0,totalEnergyFacetJ=0;
};
struct Ledger {
  // Append-only process-local records. No MPI collective occurs in an impact
  // callback. The owner aggregates/writes these artificial-boundary transfers.
  // Energy is transferred to an ideal collector ledger; it is NOT surface heat
  // and establishes no implantation/neutralization/residence/release outcome.
  std::vector<Record> records;
};

// Pure deterministic services, no PIC/SPICE/RNG/I/O dependencies. Failure
// preserves output/ledger. Collect reserves all storage BEFORE publication;
// failure never consumes or deletes an incident model packet.
bool ValidateStates(const std::vector<StateDefinition>&,std::string* error);
bool Prepare(const Incident&,Record* output,std::string* error);
bool Collect(const Incident&,Ledger*,std::string* error);
} } }
#endif
