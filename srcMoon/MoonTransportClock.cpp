#include "MoonTransportClock.h"

#include <cmath>
#include <limits>

namespace Moon { namespace TransportTime {
namespace {
bool Fail(std::string* error,const char* message) {
  if(error)*error=message;
  return false;
}
}
bool Clock::Initialize(double origin,std::string* error) {
  if(!error)return false;
  if(initialized_ || !std::isfinite(origin))
    return Fail(error,"common transport clock needs one explicit finite origin ET; reinitialization/restart is unsupported");
  origin_=origin;initialized_=true;error->clear();return true;
}
bool Clock::Project(long double relative,double* output,std::string* error) const {
  if(!error)return false;
  if(!initialized_ || !output || !std::isfinite(relative) || relative<0)
    return Fail(error,"transport projection needs an initialized explicit origin and nonnegative finite relative seconds");
  const long double exact=static_cast<long double>(origin_)+relative;
  if(!std::isfinite(exact) || std::abs(exact)>std::numeric_limits<double>::max())
    return Fail(error,"transport epoch projection exceeds binary64 ET storage");
  const double projected=static_cast<double>(exact);
  if(!std::isfinite(projected))return Fail(error,"nonfinite transport epoch projection");
  *output=projected;error->clear();return true;
}
bool Clock::Begin(double duration,Step* output,std::string* error) {
  if(!error)return false;
  if(!initialized_ || active_ || !output || !std::isfinite(duration) || !(duration>0) ||
      completed_==std::numeric_limits<std::uint64_t>::max())
    return Fail(error,"common transport Begin requires initialized idle clock, positive finite seconds and an unused sequence");
  // Compensated elapsed-relative addition. In contrast, absoluteET += dt can
  // accumulate a large epoch-scale rounding error even for exactly equal dt.
  const long double increment=static_cast<long double>(duration)-compensation_;
  const long double next=elapsed_+increment;
  const long double nextCompensation=(next-elapsed_)-increment;
  Step candidate;candidate.sequence=completed_+1;candidate.originEtS=origin_;
  candidate.durationS=duration;candidate.relativeBeginS=elapsedS();
  candidate.relativeEndS=next-nextCompensation;
  if(!std::isfinite(candidate.relativeEndS) || !(candidate.relativeEndS>candidate.relativeBeginS) ||
      !Project(candidate.relativeBeginS,&candidate.beginEtS,error) ||
      !Project(candidate.relativeEndS,&candidate.endEtS,error))
    return Fail(error,"positive common step is unrepresentable in relative seconds / absolute ET; review physical time resolution");
  if(!(candidate.endEtS>candidate.beginEtS))
    return Fail(error,"positive common outer step cannot advance binary64 ET at this origin; increase step only after numerical review or qualify a finer epoch representation");
  nextElapsed_=next;nextCompensation_=nextCompensation;step_=candidate;
  active_=true;*output=candidate;error->clear();return true;
}
bool Clock::Matches(const Step& token) const {
  return active_ && token.sequence==step_.sequence && token.originEtS==step_.originEtS &&
    token.durationS==step_.durationS && token.beginEtS==step_.beginEtS && token.endEtS==step_.endEtS &&
    token.relativeBeginS==step_.relativeBeginS && token.relativeEndS==step_.relativeEndS;
}
bool Clock::Commit(const Step& token,std::string* error) {
  if(!error)return false;
  if(!Matches(token))return Fail(error,"common transport Commit has no matching immutable active-step token; duplicate/stale/forged commits are refused");
  elapsed_=nextElapsed_;compensation_=nextCompensation_;completed_=step_.sequence;
  active_=false;error->clear();return true;
}
bool Clock::AbortUnaccepted(const Step& token,std::string* error) {
  if(!error)return false;
  if(!Matches(token))return Fail(error,"common transport abort has no matching PRE-WORK step token");
  active_=false;error->clear();return true;
}
bool Clock::BirthFromOffset(double offset,BirthTime* output,std::string* error) const {
  if(!error)return false;
  if(!active_ || !output || !std::isfinite(offset) || offset<0 || offset>step_.durationS)
    return Fail(error,"source departure offset is outside the active common physical interval [0,dt]");
  BirthTime candidate;candidate.offsetS=offset;candidate.remainingS=step_.durationS-offset;
  candidate.relativeS=offset==step_.durationS?step_.relativeEndS:
    offset==0?step_.relativeBeginS:step_.relativeBeginS+static_cast<long double>(offset);
  if(candidate.relativeS<step_.relativeBeginS || candidate.relativeS>step_.relativeEndS ||
      (offset>0 && candidate.relativeS==step_.relativeBeginS) ||
      (offset<step_.durationS && candidate.relativeS==step_.relativeEndS) ||
      !Project(candidate.relativeS,&candidate.epochEtS,error))
    return Fail(error,"positive source offset cannot be represented in the declared relative time coordinate");
  *output=candidate;error->clear();return true;
}
bool Clock::BirthFromRemaining(double remaining,BirthTime* output,std::string* error) const {
  if(!error)return false;
  if(!active_ || !output || !std::isfinite(remaining) || remaining<0 || remaining>step_.durationS)
    return Fail(error,"source remaining seconds are outside the active common interval [0,dt]");
  BirthTime candidate;candidate.remainingS=remaining;candidate.offsetS=step_.durationS-remaining;
  // The supplied unused physical interval is authoritative. Subtract it from
  // the relative END directly, avoiding two rounded offset/epoch additions.
  candidate.relativeS=remaining==step_.durationS?step_.relativeBeginS:
    remaining==0?step_.relativeEndS:step_.relativeEndS-static_cast<long double>(remaining);
  if(candidate.relativeS<step_.relativeBeginS || candidate.relativeS>step_.relativeEndS ||
      (remaining>0 && candidate.relativeS==step_.relativeEndS) ||
      (remaining<step_.durationS && candidate.relativeS==step_.relativeBeginS) ||
      !Project(candidate.relativeS,&candidate.epochEtS,error))
    return Fail(error,"source remaining interval has an unrepresentable relative birth coordinate");
  *output=candidate;error->clear();return true;
}
bool EpochProjectionAgrees(double actual,double projected) {
  if(!std::isfinite(actual) || !std::isfinite(projected))return false;
  double lower=projected,upper=projected;
  for(int addition=0;addition<2;++addition) {
    lower=std::nextafter(lower,-std::numeric_limits<double>::infinity());
    upper=std::nextafter(upper,std::numeric_limits<double>::infinity());
  }
  return actual>=lower && actual<=upper;
}
} } // namespace Moon::TransportTime
