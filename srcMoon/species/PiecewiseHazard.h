#ifndef SRCMOON_SPECIES_PIECEWISE_HAZARD_H
#define SRCMOON_SPECIES_PIECEWISE_HAZARD_H

#include <cstddef>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace PiecewiseHazard {

// A caller DECLARES this affine law over one explicit physical-time segment.
// These SI endpoint rates are not an observational interpolation policy.
struct Channel {
  std::size_t reaction=0;
  long double beginRateS=0,endRateS=0; // s^-1, nonnegative and finite
};
struct Segment {
  long double beginS=0,endS=0; // nonnegative relative physical seconds
  std::vector<Channel> channels;
};
struct InstantaneousRate { std::size_t reaction=0; long double rateS=0; };
struct ChannelExposure { std::size_t reaction=0; long double integratedHazard=0; };
struct Exposure {
  long double integratedHazard=0; // dimensionless analytic integral
  std::vector<ChannelExposure> channels;
};
// Exact native clock value is high+low. Both components are binary64; low
// may be signed and separated from high by the entire binary64 exponent
// range. Never collapse this pair before a threshold comparison.
struct CompensatedDepth {double high=0,low=0;};
struct NumericControls;
struct EventLocation;
enum class Status { Complete,Event,IncompleteSupport,InvalidInput,NumericalError };
class ClockUpdate {
 public:
  bool valid()const{return valid_;}
  bool event()const{return event_;}
  const CompensatedDepth& origin()const{return origin_;}
  const CompensatedDepth& residual()const{return residual_;}
  // Rigorous absolute bound for rounding the EXACT rational no-event
  // residual back to native two-double storage. Zero means exact. This is
  // representation diagnostics, never a physical acceptance tolerance.
  long double residualRoundingErrorBound()const{return roundingBound_;}
 private:
  friend Status LocateEvent(const std::vector<Segment>&,long double,long double,
      const CompensatedDepth&,const NumericControls&,EventLocation*,std::string*);
  bool valid_=false,event_=false;
  CompensatedDepth origin_,residual_;
  long double roundingBound_=0;
};
struct EventLocation {
  // Complete/no-event: NaN time, no rates, full-window integratedHazard.
  // Event: located time and instantaneous rates, with integratedHazard equal
  // to the requested optical crossing depth. Representational upward time
  // rounding does not authorize an extra physical optical-depth debit.
  long double eventTimeS=0;
  std::size_t segmentIndex=0;
  long double integratedHazard=0,totalInstantaneousRateS=0;
  std::vector<InstantaneousRate> instantaneousRates;
  // Set ONLY by the exact compensated overload. For that overload an event's
  // integratedHazard is a display approximation; this immutable receipt,
  // rather than a rounded scalar, owns the exact event/no-event decision.
  ClockUpdate clockUpdate;
};
struct NumericControls {
  // REQUIRED explicit work bound. There is deliberately no physical time or
  // optical-depth acceptance epsilon and no implicit maximum-step policy.
  std::size_t maximumRootRefinements;
};

// Gaps are preserved and permitted; overlap, missing/changed channel entries,
// nonfinite endpoints or negative rates are invalid. No table sorting occurs.
bool ValidateHistory(const std::vector<Segment>& history,std::string* error);

// All result objects remain unchanged on incomplete/invalid/numerical returns.
// Complete from LocateEvent means the ENTIRE window was covered but no event
// occurred. Event uses LEFT-limit rates at an exact first-hitting segment end;
// total rate zero remains explicit for downstream refusal, never a fitted CDF.
// Both functions are pure: no I/O, RNG, optical-clock or particle mutation.
Status Integrate(const std::vector<Segment>& history,long double beginS,
    long double endS,Exposure* output,std::string* error);
Status LocateEvent(const std::vector<Segment>& history,long double beginS,
    long double endS,long double opticalDepth,const NumericControls& controls,
    EventLocation* output,std::string* error);
// Same support/left-limit/first-representable-time contract, now using ALL
// exact binary components. No fixed-precision approximation decides crossing.
Status LocateEvent(const std::vector<Segment>& history,long double beginS,
    long double endS,const CompensatedDepth& opticalDepth,const NumericControls& controls,
    EventLocation* output,std::string* error);

}}}
#endif
