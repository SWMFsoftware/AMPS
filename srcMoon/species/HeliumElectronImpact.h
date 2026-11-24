#ifndef SRCMOON_HELIUM_ELECTRON_IMPACT_H
#define SRCMOON_HELIUM_ELECTRON_IMPACT_H

#include "ElectronImpact.h"
#include "PhysicalStateProfile.h"

#if __has_include("../../models/chemical_database/ReservoirProducts.h")
#include "../../models/chemical_database/ReservoirProducts.h"
#else
#include "../../src/models/chemical_database/ReservoirProducts.h"
#endif

namespace Moon { namespace Species { namespace HeliumElectronImpact {
namespace R=ReactionExecution;
namespace B=::AMPS::ChemicalDatabase::ReservoirProducts;
using ::AMPS::ChemicalDatabase::Code;
using ::AMPS::ChemicalDatabase::Status;

// This pilot is an explicitly prescribed Maxwellian for stationary 4He at
// 10 eV. All vectors below use the same instantaneous inertial axes.
struct Bath {
  double thermalEnergyEv=0,densityM3=-1;
  B::Vector3 bulkVelocityMS{},targetVelocityMS{};
  std::string frameIdentity;
};
struct Evaluation {
  double coefficientM3S=0,rateS=0,logSlope=0,fitThresholdEv=0;
  B::Real meanIncomingEnergyJ=0,physicalThresholdJ=0;
};
struct Admission {
  B::Evidence evidence;
  B::NumericalBounds numerical;
  R::KinematicPolicy planner;
  std::string crosswalkPath,modelErrorPath;
  B::Real thresholdCrosswalkMaximumDifferenceJ=-1;
  // Reviewed bound for the actual conditional energy/slope crosswalk in
  // this bath. The original rate-fit error by itself supplies no such bound.
  B::Real maximumAbsoluteLogSlopeError=-1;
};
struct Composition {
  R::KinematicPolicy policy;
  R::EventPlan plan;
  B::Result reservoir;
};

// Explicit, immutable normal-controller composition service. Prepare alone
// grants no activation. Compose requires a reviewed source/observable
// admission and supplies the policy consumed by the usual particle planner.
class Provider {
 public:
  static Status Prepare(const PhysicalStateProfile::Registry& profile,
      const std::string& nativeD05Path,const std::string& nativeD05Sha256,
      int neutralCompiledSpecies,int ionCompiledSpecies,Provider* output);
  Status Evaluate(const Bath&,Evaluation* output) const;
  Status Compose(const Bath&,const Admission&,const R::ParentState&,
      const R::TraceInterval&,double acceptedEventEpochS,double ionReferenceWeight,
      const R::OpticalDepthClock& ionClock,Composition* output) const;
  const std::vector<R::SpeciesIdentity>& catalog() const {return catalog_;}
  const R::Reaction& reaction() const {return reaction_;}
 private:
  ElectronImpact::MaxwellianVoronovTable rates_;
  ElectronImpact::Coefficient source_;
  std::string d05Identity_,profileIdentity_;
  B::Real electronVoltJ_=0;
  std::array<B::State,3> states_;
  std::vector<R::SpeciesIdentity> catalog_;
  R::Reaction reaction_;
  B::MaxwellianInput Incoming(const Bath&,const Evaluation&) const;
};
}}}
#endif
