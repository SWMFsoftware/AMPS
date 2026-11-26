// Normal-startup owner of the approved I24 total-capture case (2026-11-25):
// H+ + Na(3s) -> Na+ + H(n), n = 1, 2, 3, for stationary force-free 23Na in a
// prescribed monoenergetic proton beam (species/SodiumTotalCapture.h).
// Register before amps_init/particle storage; the controller then executes
// the channel through the normal reaction path in the frozen LSO frame.
#ifndef SRCMOON_SODIUM_TOTAL_CAPTURE_COHORT_H
#define SRCMOON_SODIUM_TOTAL_CAPTURE_COHORT_H
#include "MoonReactionController.h"
#include "MoonSemanticRandom.h"
#include "species/SodiumTotalCapture.h"
#include <memory>

namespace Moon { namespace Chemistry { namespace SodiumTotalCaptureCohort {
namespace S=Species::SodiumTotalCapture;
namespace R=Species::ReactionExecution;
struct Request {
  S::Provider provider;
  S::Beam beam;                     // density [m^-3] and velocity [m/s] in LSO
  SemanticRandom::Configuration random;
  double elementaryChargeC=0;
  std::string elementaryChargeIdentity;
  std::size_t maximumRootRefinements=0; // explicitly supplied numerical effort
};
struct Receipt {
  SemanticRandom::State key;
  R::ParentState parent;
  R::TraceInterval interval;
  double epochEt=0;
  double shellUniform=0;            // the keyed Kinematics draw that selected H(n)
  S::Composition composition;
};
class Prepared {
 public:
  // Refuses a second owner, a non-LSO provider, an inadmissible beam (energy
  // outside 1-5 keV/u, negative density), zero solver effort or a wrong e.
  static bool RegisterForNormalStartup(const Request&,std::shared_ptr<Prepared>* owner,std::string* error);
  const ReactionLedger& ledger()const{return ledger_;}
  // Composition occurs only at an accepted event. Publication, weights and
  // PIC readback must be taken from the paired COMMITTED ledger entry.
  const std::vector<Receipt>& receipts()const{return receipts_;}
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
