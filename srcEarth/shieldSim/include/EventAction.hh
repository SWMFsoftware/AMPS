#ifndef SHIELDSIM_EVENT_ACTION_HH
#define SHIELDSIM_EVENT_ACTION_HH

/* ============================================================================
 * EventAction.hh
 *
 * Defines event boundaries for statistical estimation. RunAction resets its
 * per-event ledgers at BeginOfEventAction and commits one independent history
 * sample at EndOfEventAction. This prevents step-to-step correlations inside a
 * primary history from being mistaken for independent Monte Carlo samples.
 * ========================================================================== */

#include <G4UserEventAction.hh>

class G4Event;
class RunAction;

class EventAction : public G4UserEventAction {
public:
  explicit EventAction(RunAction* runAction);
  ~EventAction() override;
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction  (const G4Event*) override;
private:
  RunAction* fRunAction=nullptr;
};

#endif // SHIELDSIM_EVENT_ACTION_HH
