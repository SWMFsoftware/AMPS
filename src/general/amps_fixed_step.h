#ifndef AMPS_FIXED_STEP_H
#define AMPS_FIXED_STEP_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace AMPS { namespace Coupling {

// A component's completed physical time remains in its coupler ledger. This
// value binds one upcoming interval [s]; it is not another elapsed-time clock.
// PIC's counter-OFF convention presents the end epoch during particle work,
// while the coupler may commit that endpoint only after successful execution.
class FixedStepContract {
 public:
  // Numerical epochs and duration are SI seconds from the component's declared
  // simulation origin. Keep duration independently: subtracting two rounded
  // late-time epochs need not recover the configured physical duration. The
  // sequence identifies accepted windows, never host calls or particle ages.
  struct Window {
    double startEpochS=0.0,durationS=0.0,endEpochS=0.0;
    std::uint64_t sequence=0;
  };
  // These observations never change clock state. An active interval spans
  // planning through verified particle completion, and cannot be overwritten.
  // NominalDurationS is meaningful only after successful Configure().
  bool Configured() const {return configured_;}
  bool Active() const {return phase_!=Phase::Idle;}
  bool Planned() const {return phase_==Phase::Planned;}
  double NominalDurationS() const {return nominalDurationS_;}
  // Receive-time metadata can distinguish bootstrap (no accepted interval)
  // from a background epoch that may already have been consumed. This value
  // never counts skipped host calls and is not an elapsed-time accumulator.
  std::uint64_t Sequence() const {return sequence_;}
  // Capture dt once after model initialization. Subsequent registrations may
  // confirm the same value, but cannot change sampling exposure or the native
  // mover duration. Validate before writing so invalid input leaves the
  // existing contract intact; no tolerance-based replacement is permitted.
  void Configure(double durationS) {
    if(!std::isfinite(durationS)||durationS<=0.0)
      throw std::logic_error("coupled AMPS timestep must be finite and positive [s]");
    if(configured_) {CheckDuration(durationS);return;}
    nominalDurationS_=durationS;configured_=true;
  }
  // All callers compare the actual production dt, not a cached candidate or
  // a host-limit remainder. Exact equality is intentional: both values come
  // from the initialized timestep storage, not independent measurements.
  // Throws are translated to the normal AMPS fatal path at the host boundary.
  void CheckDuration(double actualDurationS) const {
    if(!configured_||actualDurationS!=nominalDurationS_)
      throw std::logic_error("coupled AMPS timestep cannot change after initialization; sampling requires the configured constant dt");
  }
  // Publish a const view only while a real host interval is bound. Returning
  // the last stale interval after commit could give a newly born particle the
  // previous step's age, so inactive access is an explicit contract error.
  const Window& Current() const {
    if(!Active())throw std::logic_error("externally owned AMPS clock has no bound full-step window");
    return window_;
  }
  // completedEpochS belongs to the coupler's successful-execution ledger.
  // Compute one endpoint with the unchanged nominal dt and keep that exact
  // stored result for PIC and commit. An endpoint that overflows or rounds
  // back to departure cannot describe a numerical step and must reject.
  // Host comparisons are exact; nextafter/epsilon adjustments would invent
  // execution time. A false result changes no interval, sequence or clock.
  bool Plan(double completedEpochS,double hostLimitS,double actualDurationS) {
    CheckDuration(actualDurationS);
    if(Active())throw std::logic_error("cannot replace an active coupled AMPS step");
    const double end=completedEpochS+nominalDurationS_;
    if(!std::isfinite(completedEpochS)||!std::isfinite(hostLimitS)||
       !std::isfinite(end)||end<=completedEpochS)
      throw std::logic_error("coupled AMPS step has invalid or unrepresentable physical endpoints [s]");
    // Never shorten dt, nudge an endpoint or count host acknowledgements as
    // particle evolution. A skipped call leaves both the ledger and PIC alone.
    if(end>hostLimitS)return false;
    if(sequence_==std::numeric_limits<std::uint64_t>::max())
      throw std::logic_error("coupled AMPS step sequence overflow");
    window_={completedEpochS,nominalDurationS_,end,++sequence_};
    phase_=Phase::Planned;return true;
  }
  // Called after the departure background/source window has been frozen and
  // before PIC executes. This reproduces the audited counter-OFF presentation
  // of the end epoch without enabling PIC's own increment. Particle departure
  // and newborn remaining duration continue to come from the frozen window,
  // not subtraction from this mutable counter. This is a phase publication,
  // not yet a claim that physical evolution has completed.
  void EnterParticlePhase(double& picEpochS,double actualDurationS) {
    CheckDuration(actualDurationS);
    if(phase_!=Phase::Planned||picEpochS!=window_.startEpochS)
      throw std::logic_error("coupled AMPS particle phase requires its original physical departure");
    picEpochS=window_.endEpochS;phase_=Phase::Executing;
  }
  // Only the application, after PIC returns, readers/workers join, and native
  // source-window closure succeeds, calls this boundary. Verify the end and
  // unchanged dt without repairing either one. If work failed after mutation,
  // retain the active interval for diagnosis; there is no rollback/retry claim.
  void CompleteParticlePhase(double picEpochS,double actualDurationS) {
    CheckDuration(actualDurationS);
    if(phase_!=Phase::Executing||picEpochS!=window_.endEpochS)
      throw std::logic_error("coupled AMPS particle phase did not reach its declared full-step endpoint");
    phase_=Phase::Complete;
  }
  // The wrapper invokes this only after the entire application step returns,
  // including its post-particle work. Return the already verified endpoint
  // for assignment to the one completed-time ledger; never increment it a
  // second time or derive it from a scheduler acknowledgement.
  double Commit(double picEpochS,double actualDurationS) {
    CheckDuration(actualDurationS);
    if(phase_!=Phase::Complete||picEpochS!=window_.endEpochS)
      throw std::logic_error("cannot commit an uncompleted coupled AMPS physical step");
    phase_=Phase::Idle;return window_.endEpochS;
  }
 private:
  // This state machine prevents premature completion, duplicate PIC entry and
  // replacement of a live interval. None of these fields is a second elapsed
  // physical clock, and no field is written by particle worker threads.
  enum class Phase {Idle,Planned,Executing,Complete};
  bool configured_=false;
  double nominalDurationS_=0.0;
  std::uint64_t sequence_=0;
  Window window_;
  Phase phase_=Phase::Idle;
};
} }
#endif
