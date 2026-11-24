#ifndef SRCMOON_SODIUM_CAPTURE_COHORT_H
#define SRCMOON_SODIUM_CAPTURE_COHORT_H
#include "MoonReactionController.h"
#include "MoonSemanticRandom.h"
#include "species/SodiumGroundCapture.h"
#include <memory>

namespace Moon { namespace Chemistry { namespace SodiumCaptureCohort {
namespace S=Species::SodiumGroundCapture;
namespace R=Species::ReactionExecution;
struct Request {
  S::Provider provider;
  S::Beam beam;
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
  S::Composition composition;
};
class Prepared {
 public:
  // Register before normal amps_init/particle storage. The immutable context
  // survives caller scope and is installed after generated masses and weights
  // exist. Another owner or a different state-to-species binding is refused.
  static bool RegisterForNormalStartup(const Request&,std::shared_ptr<Prepared>* owner,std::string* error);
  const ReactionLedger& ledger()const{return ledger_;}
  // Composition occurs only at an accepted event. Publication/weights/actual
  // PIC setter readback must be taken from the paired COMMITTED ledger entry.
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
