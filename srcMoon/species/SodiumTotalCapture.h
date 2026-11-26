// Total one-electron capture of a monoenergetic proton beam on ground-state
// 23Na, approved 2026-11-25 for the first I24 transport case:
//
//   H+ + Na(3s) -> Na+(2p6) + H(n),  n = 1, 2, 3
//
// * Rate nu = n_p sigma_tot(E) g. sigma_tot is the installed Morgan et al.
//   (1985) Table 2 Na law (ChargeExchange::ProtonNaMorgan1985LogLog). E is the
//   projectile energy per atomic mass unit in the target rest frame,
//   E = (1/2) u g^2 [eV/u]. Only 1000 <= E <= 5000 eV/u is admitted.
// * Final shell n is selected with one supplied open uniform from the
//   Shingal et al. (1986) Table 4 A34 shell partials, interpolated log-log
//   between adjacent printed rows and normalized (no extrapolation).
// * Products follow zero-angle exchange in the pair centre-of-mass frame: H
//   continues along the incident relative velocity, Na+ opposite. The
//   available COM energy is (1/2) mu g^2 + IE(H) - IE(Na) - E_exc(n).
// * The prompt H(n) decays to H(1s) before transport. E_exc(n) is booked at
//   the event as radiated energy leaving into the untracked environment;
//   the photon momentum is omitted and its H recoil is reported as a bound.
// * The omitted transverse impulse is bounded by the approved maximum COM
//   angle of 0.4 degrees (Courbin et al. 1990 extent at 0.5 keV/u).
//
// E_exc(n) is the statistical-weight centroid of the NIST ASD H I nlj levels.
// Both tables come from srcMoon/test/integration-tests/I24/prepare_total_capture.py;
// their SHA-256s and byte counts are compiled below, so other bytes are refused.
// Units are SI unless a name says otherwise. All vectors in one call use the
// same instantaneous Cartesian axes of the declared event frame.
// Contract: srcMoon/test/integration-tests/I24/TOTAL_CAPTURE_CONTRACT.md
#ifndef SRCMOON_SPECIES_SODIUM_TOTAL_CAPTURE_H
#define SRCMOON_SPECIES_SODIUM_TOTAL_CAPTURE_H
#include "ReducedKinematics.h"
#include "ChargeExchange.h"
#if __has_include("../../models/chemical_database/Runtime.h")
#include "../../models/chemical_database/Runtime.h"
#else
#include "../../src/models/chemical_database/Runtime.h"
#endif
#include <array>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace SodiumTotalCapture {
namespace R=ReactionExecution;
namespace K=ReducedKinematics;
namespace DB=::AMPS::ChemicalDatabase;
using Real=K::Real;
extern const char ModelId[];
constexpr int Shells=3;
// Approved envelope [eV/u] and maximum pair-COM scattering angle [deg].
constexpr Real MinimumEnergyEvPerU=1000.0L,MaximumEnergyEvPerU=5000.0L;
constexpr Real MaximumComAngleDeg=0.4L;

struct Definition {
  // In this order: Na(3s), H+(bare), Na+(2p6), H(1s). All formation energies
  // refer to the separated neutral ground atoms (as in SodiumGroundCapture).
  std::array<K::ResolvedState,4> states;
  int sodiumSpecies=-1,sodiumIonSpecies=-1,hydrogenSpecies=-1;
  Real atomicMassUnitKg=0,electronVoltJ=0;
  std::string constantsIdentity;
  // Directory holding the two total-capture v1 publications. The shared
  // registry checks their SHA-256s before exposing the Morgan law.
  std::string morganPublicationDirectory;
  // Outputs of prepare_total_capture.py (digests compiled in the provider).
  std::string shingalTablePath,hydrogenLevelsPath;
  K::NumericControls numerical;
  R::ConservationBounds ledgerBounds;
  std::string frame;
};
struct Beam {
  Real protonDensityM3=-1;
  K::Vector3 protonVelocityMS{};
};
struct Evaluation {
  Real energyEvPerU=0,relativeSpeedMS=0,totalCrossSectionM2=0,rateS=0;
  std::array<Real,Shells> shellFraction{}; // n = 1, 2, 3; sums to one
};
struct Composition {
  Evaluation evaluation;
  int shell=0;                     // selected principal quantum number of prompt H(n)
  Real shellQJ=0;                  // IE(H) - IE(Na) - E_exc(n) [J per event]
  Real radiatedEnergyJ=0;          // E_exc(n) [J], booked to the untracked environment
  Real finalComMomentumKgMS=0;     // p_f, magnitude of each product's COM momentum
  Real omittedImpulseBoundKgMS=0;  // 2 p_f sin(theta_max/2)
  Real sodiumIonRecoilBoundJ=0;    // impulse^2 / (2 m_Na+)
  Real photonRecoilSpeedBoundMS=0; // E_exc(n) / (m_H c), omitted
  K::FinalState products;          // [Na+, H(1s)]
  R::KinematicPolicy policy;
  R::EventPlan plan;
};
class Provider {
 public:
  static DB::Status Prepare(const Definition&,Provider* output);
  // Rate, Morgan sigma and shell fractions for the beam and a target velocity.
  DB::Status Evaluate(const Beam&,const R::Vector3& sodiumVelocityMS,Evaluation* output) const;
  // term 0 = Na+, term 1 = H. Each product inherits the parent's physical
  // weight through its own reference weight and receives a fresh clock.
  // shellOpenUniform must lie strictly inside (0,1).
  DB::Status Compose(const Beam&,const R::ParentState&,const R::TraceInterval&,
      Real shellOpenUniform,const std::array<double,2>& referenceWeights,
      const std::array<R::OpticalDepthClock,2>& freshClocks,Composition* output) const;
  const std::vector<R::SpeciesIdentity>& catalog() const {return catalog_;}
  const R::Reaction& reaction() const {return reaction_;}
  const std::string& frame() const {return definition_.frame;}
  // Shell excitation E_exc(n) [J] for n = 1, 2, 3 (n = 1 is zero).
  Real shellExcitationJ(int n) const {return (n>=1&&n<=Shells)?excitationJ_[n-1]:-1;}
 private:
  struct Row {Real energyEvPerU=0;std::array<Real,Shells> sigmaM2{};};
  Definition definition_;
  std::vector<R::SpeciesIdentity> catalog_;
  std::array<double,4> internal_{};
  R::Reaction reaction_;
  K::FormationBalance formation_;
  ChargeExchange::Registry morgan_;
  std::vector<Row> shingal_;
  std::array<Real,Shells> excitationJ_{};
};
}}}
#endif
