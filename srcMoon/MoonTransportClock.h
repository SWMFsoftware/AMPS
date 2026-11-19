#ifndef SRCMOON_MOON_TRANSPORT_CLOCK_H
#define SRCMOON_MOON_TRANSPORT_CLOCK_H

#include <cstdint>
#include <string>

namespace Moon { namespace TransportTime {

// A transport interval is relative SI seconds from an EXPLICIT physical time
// origin. Absolute ET is a projection for existing PIC/SPICE-facing records,
// never the repeatedly incremented authoritative accumulator. No SPICE call,
// PIC pointer, RNG, source coefficient or observation belongs in this class.
struct Step {
  std::uint64_t sequence=0;
  double originEtS=0,durationS=0,beginEtS=0,endEtS=0;
  long double relativeBeginS=0,relativeEndS=0;
};
struct BirthTime {
  double offsetS=0,remainingS=0,epochEtS=0;
  long double relativeS=0;
};

class Clock {
 public:
  // All bool functions leave both state and caller output unchanged on error.
  // Initialize is one-shot: restarting an application requires a separately
  // designed full clock/particle/source checkpoint, not a hidden reinitialize.
  bool Initialize(double explicitOriginEtS,std::string* error);
  bool Begin(double physicalDurationS,Step* output,std::string* error);
  bool Commit(const Step& token,std::string* error);
  // Abort cancels a PRE-WORK begin only. This pure class owns no particles;
  // callers must never use it as rollback after publishing physical work.
  bool AbortUnaccepted(const Step& token,std::string* error);
  bool BirthFromOffset(double offsetS,BirthTime* output,std::string* error) const;
  bool BirthFromRemaining(double remainingS,BirthTime* output,std::string* error) const;
  bool Project(long double relativeS,double* epochEtS,std::string* error) const;
  bool initialized() const {return initialized_;}
  bool active() const {return active_;}
  double originEtS() const {return origin_;}
  long double elapsedS() const {return elapsed_-compensation_;}
  std::uint64_t completedSteps() const {return completed_;}
  const Step* step() const {return active_?&step_:nullptr;}
 private:
  bool Matches(const Step& token) const;
  bool initialized_=false,active_=false;
  double origin_=0;
  long double elapsed_=0,compensation_=0,nextElapsed_=0,nextCompensation_=0;
  std::uint64_t completed_=0;
  Step step_;
};

// Two rounded absolute epoch additions can differ from a single origin-relative
// projection by at most two adjacent binary64 values for the supported native
// interval construction. This is an arithmetic envelope, NOT a timestep,
// chemical optical-depth tolerance or observational acceptance threshold.
// Runtime completion refuses larger discrepancies rather than expanding it.
bool EpochProjectionAgrees(double actualEtS,double projectedEtS);

} // namespace TransportTime

// Normal Moon adapter. Default independent histories are exact pass-through.
// Initialize receives the UTC->ET value from ordinary SPICE initialization;
// Begin/Commit run ONLY at the normal outer-step boundary, where collectives
// are safe. Source/mover callbacks read the immutable active step, no MPI.
bool CommonTransportTimeSelected();
bool InitializeCommonTransportTime(double explicitOriginEtS,std::string* error);
bool BeginCommonTransportStep(double frozenExternalEpochEtS,std::string* error);
bool CommitCommonTransportStep(std::string* error);
const TransportTime::Clock* GetCommonTransportClock();

} // namespace Moon
#endif
