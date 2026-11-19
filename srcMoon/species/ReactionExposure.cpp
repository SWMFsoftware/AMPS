#include "ReactionExposure.h"

#include <cmath>
#include <limits>
#include <utility>

namespace Moon { namespace Species { namespace ReactionExposure {
namespace {
Status Invalid(const char* message,std::string* error) {
  if(error)*error=message;
  return Status::InvalidInput;
}
}

Status Locate(const std::vector<PiecewiseHazard::Segment>& history,
    long double beginS,long double endS,
    const ReactionExecution::OpticalDepthClock& origin,
    const PiecewiseHazard::NumericControls& controls,
    Preview* output,std::string* error) {
  if(!output || !origin.valid() || origin.consumed())
    return Invalid("exposure preview requires output and a pending origin optical clock",error);
  PiecewiseHazard::CompensatedDepth depth;
  if(!origin.PendingDepthComponents(&depth.high,&depth.low,error))return Status::NumericalError;
  Preview value;
  const auto status=PiecewiseHazard::LocateEvent(history,beginS,endS,
      depth,controls,&value.location,error);
  if(status!=Status::Complete && status!=Status::Event)return status;
  value.event=status==Status::Event;
  *output=std::move(value);
  return status;
}

Status Evaluate(const std::vector<PiecewiseHazard::Segment>& history,
    long double beginS,long double endS,
    const ReactionExecution::OpticalDepthClock& origin,
    const PiecewiseHazard::NumericControls& controls,double uniform,
    Decision* output,std::string* error) {
  if(!output)return Invalid("exposure decision requires output storage",error);
  Decision value;
  const auto status=Locate(history,beginS,endS,origin,controls,&value.preview,error);
  if(status!=Status::Complete && status!=Status::Event)return status;
  std::vector<ReactionExecution::Hazard> rates;
  if(value.preview.event) {
    for(const auto& rate:value.preview.location.instantaneousRates) {
      const double represented=static_cast<double>(rate.rateS);
      if(!std::isfinite(represented) || (rate.rateS>0 && represented==0)) {
        if(error)*error="instantaneous channel rate is outside binary64 clock storage [s^-1]";
        return Status::NumericalError;
      }
      rates.push_back({rate.reaction,represented});
    }
  }
  // Only the immutable exact-location receipt owns crossing and the native
  // residual pair. A displayed scalar area/depth may have lost low bits and
  // must never be used to reclassify or debit this accepted prefix.
  if(!ReactionExecution::ApplyLocatedHazardClock(
      origin,value.preview.location.clockUpdate,rates,uniform,
      &value.advance,error))return Status::NumericalError;
  if(value.advance.event!=value.preview.event) {
    if(error)*error="exposure location and optical-clock crossing disagree";
    return Status::NumericalError;
  }
  *output=std::move(value);
  return status;
}

}}}
