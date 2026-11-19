// Mathematical production-path controls, not physical reaction/observation
// data. Python Fraction independently checks every serialized binary depth.
#include "srcMoon/species/ReactionExposure.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
namespace R=Moon::Species::ReactionExecution;
namespace P=Moon::Species::PiecewiseHazard;
namespace E=Moon::Species::ReactionExposure;
int Checks=0,Failures=0;
void Check(const std::string& name,bool pass){++Checks;if(!pass)++Failures;std::cout<<(pass?"PASS ":"FAIL ")<<name<<'\n';}
void Metric(const char* name,const R::OpticalDepthClock& before,const R::OpticalDepthClock& after,long double errorBound) {
  double a=0,b=0,c=0,d=0;std::string error;
  if(!before.PendingDepthComponents(&a,&b,&error)||!after.PendingDepthComponents(&c,&d,&error))return;
  std::cout<<std::hexfloat<<"CLOCK "<<name<<' '<<a<<' '<<b<<' '<<c<<' '<<d<<' '<<errorBound<<std::defaultfloat<<'\n';
}
int main(){
  std::string error;const P::NumericControls controls{256};
  R::OpticalDepthClock origin;
  Check("origin",R::OpticalDepthClock::FromOpenUniform(std::exp(-1.),&origin,&error));
  const long double base=origin.remainingExtended();
  const double rate=7e-6;
  const std::vector<P::Segment> ordinary{{0,1,{{7,rate,rate}}}};
  R::OpticalDepthClock pending=origin,afterThousand;bool good=true;long double totalBound=0;
  for(int i=0;i<10000;++i){E::Decision decision;
    if(E::Evaluate(ordinary,0,1,pending,controls,NAN,&decision,&error)!=P::Status::Complete||decision.advance.event){good=false;break;}
    totalBound+=decision.advance.residualDepthRoundingErrorBound;pending=decision.advance.clock;
    if(i==999)afterThousand=pending;
  }
  Check("ten_thousand_realistic_no_event_outer_intervals",good);
  Check("dyadic_ordinary_debits_exact_in_native_pair",good&&totalBound==0);
  Metric("ordinary10000",origin,pending,totalBound);
  R::OpticalDepthClock partitioned=origin;const std::vector<P::Segment> shortStep{{0,.125L,{{7,rate,rate}}}};
  for(int i=0;i<8000&&good;++i){E::Decision decision;
    if(E::Evaluate(shortStep,0,.125L,partitioned,controls,NAN,&decision,&error)!=P::Status::Complete){good=false;break;}
    partitioned=decision.advance.clock;
  }
  Check("equivalent_dyadic_partition_same_full_clock",good&&partitioned.SamePendingDepth(afterThousand));
  const long double predecessor=std::nextafter(base,0.L),gap=base-predecessor;
  R::OpticalDepthClock precise;const double debit=static_cast<double>(.75L*gap);
  Check("near_endpoint_pair_constructed",R::ConsumeExposure(origin,debit,&precise,&error));
  long double unchanged=-7;
  Check("legacy_exact_scalar_query_still_refuses_loss",!precise.ExactPendingDepth(&unchanged,&error)&&unchanged==-7);
  const std::vector<P::Segment> near{{0,predecessor,{{7,1,1}}}};
  E::Decision noevent;
  Check("full_components_no_early_event_at_predecessor",E::Evaluate(near,0,predecessor,precise,controls,.5,&noevent,&error)==P::Status::Complete&&!noevent.advance.event);
  Metric("near_endpoint",precise,noevent.advance.clock,noevent.advance.residualDepthRoundingErrorBound);
  const std::vector<P::Segment> tail{{0,gap,{{7,1,1}}}};E::Decision event;
  Check("retained_positive_tail_has_its_real_event",E::Evaluate(tail,0,gap,noevent.advance.clock,controls,.5,&event,&error)==P::Status::Event&&event.advance.clock.consumed());
  Check("tail_event_exact_time",event.preview.location.eventTimeS==gap/4);
  Check("full_pair_first_representable_event_not_projected_early",E::Evaluate({{0,base,{{7,1,1}}}},0,base,precise,controls,.5,&event,&error)==P::Status::Event&&event.preview.location.eventTimeS==base);
  R::ExposureClockAdvance direct;
  Check("legacy_scalar_integrated_api_retains_refusal",!R::AdvanceIntegratedHazardClock(predecessor,precise,{{7,1}},.5,&direct,&error));
  R::OpticalDepthClock fullSpan;
  Check("full_double_exponent_span_pair_constructed",R::ConsumeExposure(origin,std::numeric_limits<double>::denorm_min(),&fullSpan,&error));
  Check("cache_detects_low_component_hidden_in_scalar",!fullSpan.SamePendingDepth(origin)&&fullSpan.remainingExtended()==origin.remainingExtended());
  Check("cache_identical_full_pair",fullSpan.SamePendingDepth(fullSpan));
  Check("full_exponent_span_no_event",E::Evaluate(near,0,predecessor,fullSpan,controls,.5,&noevent,&error)==P::Status::Complete);
  Metric("full_exponent_span",fullSpan,noevent.advance.clock,noevent.advance.residualDepthRoundingErrorBound);
  Check("full_exponent_span_tail_event",E::Evaluate(tail,0,gap,noevent.advance.clock,controls,.5,&event,&error)==P::Status::Event&&event.preview.location.eventTimeS==gap);
  // Integral over0..1 of a declared ramp0..2 on0..3 is EXACTLY1/3.
  const std::vector<P::Segment> thirds{{0,3,{{7,0,2}}}};
  Check("non_dyadic_integral_supported",E::Evaluate(thirds,0,1,origin,controls,NAN,&noevent,&error)==P::Status::Complete);
  Check("non_dyadic_storage_bound_explicit",noevent.advance.residualDepthRoundingErrorBound>0);
  Metric("third",origin,noevent.advance.clock,noevent.advance.residualDepthRoundingErrorBound);
  double high=0,low=0;origin.PendingDepthComponents(&high,&low,&error);
  P::EventLocation untouched;untouched.eventTimeS=-99;
  Check("invalid_component_pair_transactional",P::LocateEvent(ordinary,0,1,P::CompensatedDepth{1,-2},controls,&untouched,&error)==P::Status::InvalidInput&&untouched.eventTimeS==-99);
  Check("nonfinite_component_pair_transactional",P::LocateEvent(ordinary,0,1,P::CompensatedDepth{1,NAN},controls,&untouched,&error)==P::Status::InvalidInput&&untouched.eventTimeS==-99);
  Check("positive_native_storage_underflow_refused",P::LocateEvent({{0,.5L,{{7,std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::denorm_min()}}}},0,.5L,
      P::CompensatedDepth{std::numeric_limits<double>::denorm_min(),0},controls,&untouched,&error)==P::Status::NumericalError&&untouched.eventTimeS==-99);
  R::ExposureClockAdvance sentinel;sentinel.reaction=12345;
  Check("unproven_default_receipt_refused",!R::ApplyLocatedHazardClock(origin,P::ClockUpdate{}, {},NAN,&sentinel,&error)&&sentinel.reaction==12345);
  E::Preview preview;
  Check("real_receipt_available",E::Locate(ordinary,0,1,origin,controls,&preview,&error)==P::Status::Complete);
  Check("receipt_matches_all_origin_bits",!R::ApplyLocatedHazardClock(fullSpan,preview.location.clockUpdate,{},NAN,&sentinel,&error)&&sentinel.reaction==12345);
  const std::vector<P::Segment> channels{{0,2,{{7,1,1},{9,0,4}}}};
  Check("instantaneous_CDF_retained",E::Evaluate(channels,0,2,origin,controls,.5,&event,&error)==P::Status::Event&&event.advance.reaction==9);
  const auto saved=origin;
  Check("query_and_commit_do_not_mutate_origin",origin.SamePendingDepth(saved)&&origin.fresh());
  double h=-7,l=-9;
  Check("consumed_components_query_transactional",!event.advance.clock.PendingDepthComponents(&h,&l,&error)&&h==-7&&l==-9);
  std::cout<<std::hexfloat<<"INPUT rate "<<rate<<" predecessor "<<predecessor<<" gap "<<gap<<std::defaultfloat<<'\n';
  std::cout<<"SUMMARY checks="<<Checks<<" failed="<<Failures<<'\n';return Failures?1:0;
}
