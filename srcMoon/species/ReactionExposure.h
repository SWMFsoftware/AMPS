#ifndef SRCMOON_SPECIES_REACTION_EXPOSURE_H
#define SRCMOON_SPECIES_REACTION_EXPOSURE_H

#include "PiecewiseHazard.h"
#include "ReactionExecution.h"

namespace Moon { namespace Species { namespace ReactionExposure {

// A pure bridge between an explicitly declared time-dependent hazard history
// and the SAME compensated optical clock/CDF used by the production constant
// controller. Times are relative physical seconds, rates s^-1 and depths
// dimensionless. This does not construct rates from a particle trajectory.
using Status=PiecewiseHazard::Status;
struct Preview {
  PiecewiseHazard::EventLocation location;
  bool event=false;
};
struct Decision {
  Preview preview;
  ReactionExecution::ExposureClockAdvance advance;
};

// Preview does not consume a random number, clock or particle. Complete means
// a fully supported interval with no crossing; Event contains its first
// representable crossing and the left-limit rates at an exact segment end.
// Unsupported future portions are refused even if an earlier event exists:
// support preflight must be independent of a particular random draw.
Status Locate(const std::vector<PiecewiseHazard::Segment>& history,
    long double beginS,long double endS,
    const ReactionExecution::OpticalDepthClock& origin,
    const PiecewiseHazard::NumericControls& controls,
    Preview* output,std::string* error);

// Reconstruct the decision from the immutable interval-start clock/history.
// No-event decisions ignore channelUniform entirely. Crossing decisions use
// lambda_i(t_event)/sum(lambda_i(t_event)), NOT integrated channel fractions.
// A zero left-limit total, an unrepresentable positive binary64 channel rate
// or incomplete support is an explicit refusal, never time shifting,
// renormalization, extrapolation or an optical-depth acceptance tolerance.
// Every exact binary origin component participates in arbitrary-size rational
// comparison/location. The no-event residual is then split into the native
// two-double clock, with its explicit residualDepthRoundingErrorBound receipt;
// positive storage underflow is refused. No low component is projected away.
// Caller must commit this result only after accepting the same trajectory.
// Every error leaves output unchanged. Physical provider/controller wiring
// remains a separate qualification task; this service enables no process.
Status Evaluate(const std::vector<PiecewiseHazard::Segment>& history,
    long double beginS,long double endS,
    const ReactionExecution::OpticalDepthClock& origin,
    const PiecewiseHazard::NumericControls& controls,double channelUniform,
    Decision* output,std::string* error);

}}}
#endif
