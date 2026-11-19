// ANALYTICAL CONTROLS ONLY: fictitious channel IDs/rates are an independent
// mathematical oracle for the compiled production exposure/clock bridge.
// They are never lunar reaction data, observations or daughter kinematics.
#include "srcMoon/species/ReactionExposure.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace P=Moon::Species::PiecewiseHazard;
namespace R=Moon::Species::ReactionExecution;
namespace E=Moon::Species::ReactionExposure;
int main() {
  int checks=0,failed=0;
  auto check=[&](bool value,const char* name) {
    ++checks;if(!value)++failed;
    std::cout<<(value?"PASS ":"FAIL ")<<name<<'\n';
  };
  std::string error;
  R::OpticalDepthClock clock;
  check(R::OpticalDepthClock::FromOpenUniform(std::exp(-1.),&clock,&error),"origin_clock");
  const long double depth=clock.remainingExtended();
  const P::NumericControls controls{256}; // explicit numerical work bound
  // Channel 7 has lambda=1; channel 9 has lambda=2t. Thus H=t+t^2,
  // t=(-1+sqrt(1+4H))/2. Instantaneous and integrated fractions differ.
  const std::vector<P::Segment> history{{0,2,{{7,1,1},{9,0,4}}}};
  const long double exact=2*depth/(1+std::sqrt(1+4*depth));
  E::Preview preview;
  check(E::Locate(history,0,2,clock,controls,&preview,&error)==P::Status::Event,"pure_first_crossing");
  check(std::abs(preview.location.eventTimeS-exact)<32*std::numeric_limits<long double>::epsilon(),"quadratic_analytic_event_time");
  check(!clock.consumed() && clock.fresh() && clock.remainingExtended()==depth,"preview_preserves_origin");
  E::Decision event;
  // u=.5 falls in channel 9's instantaneous CDF; integrated fractions
  // would incorrectly choose channel 7. This is the primary regression.
  check(E::Evaluate(history,0,2,clock,controls,.5,&event,&error)==P::Status::Event && event.advance.reaction==9,"instantaneous_not_integrated_channel_CDF");
  const long double observedTime=event.preview.location.eventTimeS;
  const long double observedProbability=event.preview.location.instantaneousRates.front().rateS/event.preview.location.totalInstantaneousRateS;
  check(event.advance.clock.consumed() && !clock.consumed(),"accepted_copy_consumed_only");
  E::Decision retry;
  check(E::Evaluate(history,0,2,clock,controls,.5,&retry,&error)==P::Status::Event && retry.preview.location.eventTimeS==event.preview.location.eventTimeS && retry.advance.reaction==event.advance.reaction,"immutable_retry_same_event");
  E::Decision noevent;
  check(E::Evaluate(history,0,.25L,clock,controls,std::numeric_limits<double>::quiet_NaN(),&noevent,&error)==P::Status::Complete && !noevent.advance.event,"no_event_ignores_channel_draw");
  check(std::abs(noevent.advance.clock.remainingExtended()-(depth-.3125L))<4*std::numeric_limits<double>::epsilon(),"analytic_integrated_debit");
  E::Decision continuation;
  check(E::Evaluate(history,.25L,2,noevent.advance.clock,controls,.5,&continuation,&error)==P::Status::Event && std::abs(continuation.preview.location.eventTimeS-exact)<16*std::numeric_limits<double>::epsilon(),"partition_preserves_absolute_crossing");
  const std::vector<P::Segment> gap{{0,.5L,{{7,1,1}}},{.6L,2,{{7,1,1}}}};
  E::Decision sentinel;sentinel.advance.reaction=12345;
  check(E::Evaluate(gap,0,2,clock,controls,.5,&sentinel,&error)==P::Status::IncompleteSupport && sentinel.advance.reaction==12345,"gap_refusal_transactional");
  const std::vector<P::Segment> futureGap{{0,2,{{7,10,10}}},{3,4,{{7,10,10}}}};
  check(E::Evaluate(futureGap,0,4,clock,controls,.5,&sentinel,&error)==P::Status::IncompleteSupport,"support_preflight_independent_of_early_draw");
  check(E::Evaluate(history,0,2,clock,controls,1.,&sentinel,&error)==P::Status::NumericalError && sentinel.advance.reaction==12345,"invalid_event_draw_transactional");
  // Set total rate H on [0,1], so crossing is EXACTLY the discontinuity.
  const std::vector<P::Segment> jump{{0,1,{{7,depth,depth},{9,0,0}}},{1,2,{{7,0,0},{9,depth,depth}}}};
  check(E::Evaluate(jump,0,2,clock,controls,.5,&event,&error)==P::Status::Event && event.preview.location.eventTimeS==1 && event.advance.reaction==7,"endpoint_uses_declared_left_limit");
  const std::vector<P::Segment> zeroAtEnd{{0,1,{{7,2*depth,0}}}};
  check(E::Locate(zeroAtEnd,0,1,clock,controls,&preview,&error)==P::Status::Event && preview.location.totalInstantaneousRateS==0,"zero_rate_endpoint_preview_explicit");
  check(E::Evaluate(zeroAtEnd,0,1,clock,controls,.5,&sentinel,&error)==P::Status::NumericalError && sentinel.advance.reaction==12345,"zero_rate_endpoint_refuses_CDF_no_time_shift");
  const std::vector<P::Segment> shadow{{0,2,{{7,0,0}}}};
  check(E::Evaluate(shadow,0,2,clock,controls,std::numeric_limits<double>::quiet_NaN(),&noevent,&error)==P::Status::Complete && noevent.advance.clock.remainingExtended()==depth,"zero_exposure_keeps_pending_depth");
  const std::vector<P::Segment> duplicate{{0,2,{{7,1,1},{7,1,1}}}};
  check(E::Evaluate(duplicate,0,2,clock,controls,.5,&sentinel,&error)==P::Status::InvalidInput,"duplicate_channel_refused");
  check(E::Evaluate(history,0,2,event.advance.clock,controls,.5,&sentinel,&error)==P::Status::InvalidInput,"consumed_origin_refused");
  R::ExposureClockAdvance tiny;R::OpticalDepthClock pending=clock;
  const long double debit=std::ldexp(1.L,-60);
  bool good=true;
  for(int i=0;i<1024;++i) {
    good=good && R::AdvanceIntegratedHazardClock(debit,pending,{},std::numeric_limits<double>::quiet_NaN(),&tiny,&error);
    pending=tiny.clock;
  }
  check(good && pending.remainingExtended()<depth && std::abs(pending.remainingExtended()-(depth-1024*debit))<2*std::numeric_limits<long double>::epsilon(),"sub_ULP_exposure_compensation");
  std::size_t selected=12345;
  check(!R::SelectHazardChannel({{7,0}},.5,&selected,&error) && selected==12345,"shared_selector_zero_total_transactional");
  check(R::SelectHazardChannel({{7,1},{9,3}},.5,&selected,&error) && selected==9,"shared_selector_constant_control");
  // A tiny nonzero compensation is lost by collapsing the exact double pair
  // into one long double. The exact pending depth exceeds predecessor(1),
  // but a rounded scalar would locate an event there. The full-component
  // bridge must retain the positive tail; only the legacy scalar API refuses.
  R::OpticalDepthClock one,rounded;R::ExposureClockAdvance debited;
  R::OpticalDepthClock::FromOpenUniform(std::exp(-1.),&one,&error);
  const long double predecessor=std::nextafter(1.L,0.L);
  const double lostDebit=static_cast<double>(.75L*(1.L-predecessor));
  check(R::ConsumeExposure(one,lostDebit,&rounded,&error),"precision_counterexample_clock_constructed");
  const std::vector<P::Segment> precise{{0,predecessor,{{7,1,1}}}};
  check(E::Evaluate(precise,0,predecessor,rounded,controls,.5,&sentinel,&error)==P::Status::Complete && !sentinel.advance.event && sentinel.advance.clock.valid() && !sentinel.advance.clock.consumed(),"compensated_pending_depth_no_early_event");
  check(!R::AdvanceIntegratedHazardClock(predecessor,rounded,{{7,1}},.5,&debited,&error) && !debited.event,"public_integrated_clock_refuses_lossy_pending_depth");
  std::cout<<std::setprecision(21)<<"METRIC actual_event_time="<<observedTime<<" actual_instantaneous_channel_7="<<observedProbability<<" reference_event_time="<<exact<<'\n';
  std::cout<<"SUMMARY checks="<<checks<<" failed="<<failed<<'\n';
  return failed?1:0;
}
