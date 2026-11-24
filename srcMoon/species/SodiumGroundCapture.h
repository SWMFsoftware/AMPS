// One source-selected reference experiment for the collision of a monoenergetic
// proton beam with 23Na. All vectors in one call use the same instantaneous
// Cartesian axes. There is no interpolation in energy or implied branching to
// excited states. Production supplies the parent, interval, weights and RNG.
#ifndef SRCMOON_SPECIES_SODIUM_GROUND_CAPTURE_H
#define SRCMOON_SPECIES_SODIUM_GROUND_CAPTURE_H
#include "ReducedKinematics.h"
#if __has_include("../../models/chemical_database/Runtime.h")
#include "../../models/chemical_database/Runtime.h"
#else
#include "../../src/models/chemical_database/Runtime.h"
#endif

namespace Moon { namespace Species { namespace SodiumGroundCapture {
namespace R=ReactionExecution;
namespace K=ReducedKinematics;
namespace DB=::AMPS::ChemicalDatabase;
using Real=K::Real;
extern const char ModelId[];

struct Definition {
  // In this order: Na(3s), H+(bare), Na+(2p6), H(1s); all physical
  // formation energies refer to the same separated neutral-atom zero.
  std::array<K::ResolvedState,4> states;
  int sodiumSpecies=-1,sodiumIonSpecies=-1,hydrogenSpecies=-1;
  Real atomicMassUnitKg=0,electronVoltJ=0;
  std::string constantsIdentity;
  K::NumericControls numerical;
  R::ConservationBounds ledgerBounds;
  std::string frame;
};
struct Beam {
  Real protonDensityM3=-1;
  K::Vector3 protonVelocityMS{};
};
struct Evaluation { Real energyEvPerU=0,relativeSpeedMS=0,crossSectionM2=0,rateS=0; };
struct Composition {
  Evaluation evaluation;
  K::FinalState products; // [Na+, H] with momentum and energy in frame
  R::KinematicPolicy policy;
  R::EventPlan plan;
};
class Provider {
 public:
  // The point CSV comes from srcMoon/tools/export_ground_capture.py. Both the
  // installed database selection and the derived lexical projection are pinned.
  static DB::Status Prepare(const DB::Catalogue&,const std::string& pointCsv,
      const Definition&,Provider* output);
  DB::Status Evaluate(const Beam&,const R::Vector3& sodiumVelocityMS,Evaluation* output) const;
  // term 0 = Na+, term 1 = H. Each inherits the full physical weight of the
  // parent via its own reference weight and receives a separately fresh clock.
  DB::Status Compose(const Beam&,const R::ParentState&,const R::TraceInterval&,
      Real polarOpenUniform,Real azimuthOpenUniform,const std::array<double,2>& referenceWeights,
      const std::array<R::OpticalDepthClock,2>& freshClocks,Composition* output) const;
  const std::vector<R::SpeciesIdentity>& catalog() const {return catalog_;}
  const R::Reaction& reaction() const {return reaction_;}
  const std::string& frame() const {return definition_.frame;}
 private:
  Definition definition_;
  std::vector<R::SpeciesIdentity> catalog_;
  std::array<double,4> internal_{};
  R::Reaction reaction_;
  Real energyEvPerU_=0,crossSectionM2_=0;
  K::FormationBalance formation_;
};
}}}
#endif
