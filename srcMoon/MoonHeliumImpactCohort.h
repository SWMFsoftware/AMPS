// Explicit normal-startup owner for the I23 prescribed-electron-bath coupling:
//
//   e(bath) + 4He(ground) -> 4He+(1s ground) + 2 e(bath)
//
// The electron bath is a declared, spatially uniform and time-constant,
// isotropic zero-drift Maxwellian at rest in the frozen LSO frame (user
// approval 2026-11-25, srcMoon/test/integration-tests/I23/HELIUM_ADMISSION_20261125.md).
// It is NOT a lunar plasma model: no wake/shadow structure, no solar-wind
// drift, no measured D06 forcing. The hazard is n_e * k_Voronov(T_e) [s^-1],
// valid only for a stationary, force-free 4He target in that frame; every
// accepted trajectory segment is checked against that domain and a moving or
// accelerated parent stops the run rather than silently using the wrong rate.
//
// Products come exclusively from the production HeliumElectronImpact::Compose
// (inherited-parent-velocity reduced reservoir). The heavy ion is tracked; the
// incoming and both outgoing electrons are signed bath terms. Compose refuses
// without the reviewed admission (threshold/state crosswalk, conditional-mean
// slope error, recoil/observable bound); this owner never supplies defaults.
//
// Units: SI throughout (m, s, m/s, kg, J, C); bath thermal energy k_B*T_e in
// eV as declared by D05. Time: particle physical epochs are SPICE ET seconds.
// Registration is explicit and programmatic; no .in key installs it.
#ifndef SRCMOON_HELIUM_IMPACT_COHORT_H
#define SRCMOON_HELIUM_IMPACT_COHORT_H
#include "MoonReactionController.h"
#include "MoonSemanticRandom.h"
#include "species/HeliumElectronImpact.h"
#include <memory>

namespace Moon { namespace Chemistry { namespace HeliumImpactCohort {
namespace H=Species::HeliumElectronImpact;
namespace R=Species::ReactionExecution;
extern const char OwnerIdentity[];

struct Request {
  H::Provider provider;       // prepared from the authenticated profile and frozen D05 table
  H::Bath bath;               // frameIdentity must be the LSO event frame
  H::Admission admission;     // reviewed evidence; planner.identity names the kinematic provider
  SemanticRandom::Configuration random;
  double elementaryChargeC=0; // exact SI value, checked against the PIC charge unit
  std::string elementaryChargeIdentity;
  std::size_t maximumRootRefinements=0; // explicit numerical effort for the exposure clock
};
struct Receipt {
  SemanticRandom::State key;
  R::ParentState parent;
  R::TraceInterval interval;
  double epochEt=0;
  double rateS=0;             // actual-event n_e*k(T_e) [s^-1]
  H::Composition composition; // policy, plan and signed reservoir moments
};
class Prepared {
 public:
  // Register before normal amps_init/particle storage. The context survives
  // the caller's scope and is installed after the generated masses, weights
  // and timesteps exist. A second owner, a non-LSO bath/event frame or an
  // incomplete admission is refused before any particle can be allocated.
  static bool RegisterForNormalStartup(const Request&,std::shared_ptr<Prepared>* owner,std::string* error);
  const ReactionLedger& ledger()const{return ledger_;}
  // Composition occurs only at an accepted event. Publication, weights and
  // actual PIC readback must be taken from the paired COMMITTED ledger entry.
  const std::vector<Receipt>& receipts()const{return receipts_;}
  const Request& request()const{return request_;}
 private:
  Request request_;
  ReactionLedger ledger_;
  std::vector<Receipt> receipts_;
  static bool Install(void*,std::string*);
  static bool Resolve(long int,Controller::Node*,void*,std::size_t*,double*,std::string*);
  static bool History(const Controller::ExposureRequest&,void*,Controller::ExposureHistory*,std::string*);
  static bool Kinematics(const Controller::EventContext&,R::KinematicPolicy*,std::string*);
};
}}}
#endif
