#include "PiecewiseHazard.h"
#include "ExactArithmetic.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace Moon { namespace Species { namespace PiecewiseHazard {
namespace {
using Real=long double;
static_assert(std::numeric_limits<Real>::radix==2,
    "exact affine threshold comparisons require a binary long-double representation");
Status Refuse(Status status,const std::string& message,std::string* error) {
  if(error)*error=message;
  return status;
}
bool Finite(Real x){return std::isfinite(x);}

// Neumaier accumulation retains short exposures beside a much larger one in
// returned approximate integrals and analytical seeds. Threshold decisions
// below use exact rational values, never this rounded display accumulator.
struct Sum {
  Real high=0,low=0;
  bool Add(Real x) {
    const Real next=high+x;
    if(!Finite(next))return false;
    low+=std::fabs(high)>=std::fabs(x)?(high-next)+x:(x-next)+high;
    high=next;return Finite(low);
  }
  Real Value()const{return high+low;}
};

// Exact threshold decisions require more than a rounded cumulative area. In
// particular, integral(2-2t,0,1-epsilon)=1-epsilon^2 can round to1 well before
// t=1. Every finite binary input is rational, so comparisons can instead be
// exact without an acceptance epsilon or a third-party arithmetic dependency.
// This private unsigned base-2^32 integer implements ONLY the positive
// operations needed by affine integrals. It has no application/state role.
struct Integer {
  std::vector<std::uint32_t> words;
  explicit Integer(std::uint64_t value=0) {
    if(value){words.push_back(static_cast<std::uint32_t>(value));if(value>>32)words.push_back(value>>32);}
  }
  bool Zero()const{return words.empty();}
  void Trim(){while(!words.empty()&&words.back()==0)words.pop_back();}
  int Compare(const Integer& b)const {
    if(words.size()!=b.words.size())return words.size()<b.words.size()?-1:1;
    for(std::size_t i=words.size();i>0;--i)if(words[i-1]!=b.words[i-1])return words[i-1]<b.words[i-1]?-1:1;
    return 0;
  }
  Integer ShiftLeft(unsigned bits)const {
    if(Zero())return Integer();
    Integer out;out.words.assign(bits/32,0);std::uint64_t carry=0;
    for(auto word:words) {const std::uint64_t part=(static_cast<std::uint64_t>(word)<<(bits%32))|carry;
      out.words.push_back(static_cast<std::uint32_t>(part));carry=part>>32;}
    if(carry)out.words.push_back(static_cast<std::uint32_t>(carry));
    return out;
  }
  unsigned TrailingZeros()const {
    if(Zero())return 0;
    unsigned bits=0;for(auto word:words) {if(word==0){bits+=32;continue;}
      while((word&1)==0){++bits;word>>=1;}break;}return bits;
  }
  void ShiftRight(unsigned bits) {
    const std::size_t whole=bits/32;
    if(whole>=words.size()){words.clear();return;}
    words.erase(words.begin(),words.begin()+whole);const unsigned part=bits%32;
    if(part)for(std::size_t i=0;i<words.size();++i)
      words[i]=static_cast<std::uint32_t>((static_cast<std::uint64_t>(words[i])>>part)|
          (i+1<words.size()?static_cast<std::uint64_t>(words[i+1])<<(32-part):0));
    Trim();
  }
};
Integer Add(const Integer& a,const Integer& b) {
  Integer out;std::uint64_t carry=0;const auto size=std::max(a.words.size(),b.words.size());
  for(std::size_t i=0;i<size;++i) {const std::uint64_t part=carry+
      (i<a.words.size()?static_cast<std::uint64_t>(a.words[i]):0)+(i<b.words.size()?static_cast<std::uint64_t>(b.words[i]):0);
    out.words.push_back(static_cast<std::uint32_t>(part));carry=part>>32;}
  if(carry)out.words.push_back(static_cast<std::uint32_t>(carry));
  return out;
}
Integer Subtract(const Integer& a,const Integer& b) {
  if(a.Compare(b)<0)throw std::logic_error("negative exact positive-arithmetic difference");
  Integer out;std::uint64_t borrow=0;
  for(std::size_t i=0;i<a.words.size();++i) {
    const std::uint64_t sub=borrow+(i<b.words.size()?static_cast<std::uint64_t>(b.words[i]):0);
    const std::uint64_t word=a.words[i];out.words.push_back(static_cast<std::uint32_t>(word-sub));borrow=word<sub;
  }
  out.Trim();return out;
}
Integer Multiply(const Integer& a,const Integer& b) {
  if(a.Zero()||b.Zero())return Integer();
  Integer out;out.words.assign(a.words.size()+b.words.size(),0);
  for(std::size_t i=0;i<a.words.size();++i) {
    std::uint64_t carry=0;
    for(std::size_t j=0;j<b.words.size();++j) {
      // The largest product+stored-word+carry is2^64-1, so this arithmetic
      // is exact in uint64_t; no compiler-specific wider integer is used.
      const std::uint64_t part=static_cast<std::uint64_t>(a.words[i])*b.words[j]+out.words[i+j]+carry;
      out.words[i+j]=static_cast<std::uint32_t>(part);carry=part>>32;
    }
    out.words[i+b.words.size()]=static_cast<std::uint32_t>(carry);
  }
  out.Trim();return out;
}
struct Rational {
  Integer numerator,denominator{1};
  explicit Rational(Real value=0) {
    if(value==0)return;
    int exponent=0;Real fraction=std::frexp(value,&exponent);
    for(int bit=0;bit<std::numeric_limits<Real>::digits;++bit) {
      numerator=numerator.ShiftLeft(1);fraction*=2;
      if(fraction>=1){numerator=Add(numerator,Integer(1));fraction-=1;}
    }
    const int shift=exponent-std::numeric_limits<Real>::digits;
    if(shift>=0)numerator=numerator.ShiftLeft(static_cast<unsigned>(shift));
    else denominator=denominator.ShiftLeft(static_cast<unsigned>(-shift));
    Normalize();
  }
  void Normalize() {
    if(numerator.Zero()){denominator=Integer(1);return;}
    const unsigned common=std::min(numerator.TrailingZeros(),denominator.TrailingZeros());
    numerator.ShiftRight(common);denominator.ShiftRight(common);
  }
};
Rational Plus(const Rational& a,const Rational& b) {
  Rational out;out.numerator=Add(Multiply(a.numerator,b.denominator),Multiply(b.numerator,a.denominator));
  out.denominator=Multiply(a.denominator,b.denominator);out.Normalize();return out;
}
Rational Minus(const Rational& a,const Rational& b) {
  Rational out;out.numerator=Subtract(Multiply(a.numerator,b.denominator),Multiply(b.numerator,a.denominator));
  out.denominator=Multiply(a.denominator,b.denominator);out.Normalize();return out;
}
Rational Times(const Rational& a,const Rational& b) {
  Rational out;out.numerator=Multiply(a.numerator,b.numerator);out.denominator=Multiply(a.denominator,b.denominator);
  out.Normalize();return out;
}
Rational Divide(const Rational& a,const Rational& b) {
  if(b.numerator.Zero())throw std::logic_error("zero exact affine duration");
  Rational out;out.numerator=Multiply(a.numerator,b.denominator);out.denominator=Multiply(a.denominator,b.numerator);
  out.Normalize();return out;
}
int Compare(const Rational& a,const Rational& b) {
  return Multiply(a.numerator,b.denominator).Compare(Multiply(b.numerator,a.denominator));
}
bool ExactComponents(const CompensatedDepth& parts,Rational* output) {
  if(!std::isfinite(parts.high)||!std::isfinite(parts.low)||!(parts.high>0))return false;
  const Rational high(parts.high),low(std::fabs(parts.low));
  if(parts.low<0) {
    if(Compare(high,low)<=0)return false;
    *output=Minus(high,low);
  }else *output=Plus(high,low);
  return true;
}

double Binary64(std::uint64_t bits) {
  double value=0;
  if constexpr(sizeof(double)==sizeof(bits))std::memcpy(&value,&bits,sizeof(value));
  return value;
}
// Round an arbitrary-size positive rational to nearest binary64, ties-even.
// Comparisons use exact integers, not a fixed-width floating approximation.
// This bounded64-step search also handles the entire subnormal exponent span.
// Values above max_double are explicitly clamped ONLY for the HIGH splitting
// component; the exact leftover is then handled separately, never discarded.
double RoundBinary64(const Rational& value) {
  constexpr std::uint64_t maximum=UINT64_C(0x7fefffffffffffff);
  std::uint64_t lo=0,hi=maximum;
  while(lo<hi) {
    const auto middle=lo+(hi-lo+1)/2;
    if(Compare(Rational(Binary64(middle)),value)<=0)lo=middle;else hi=middle-1;
  }
  const double lower=Binary64(lo);const Rational exactLower(lower);
  if(lo==maximum||Compare(exactLower,value)==0)return lower;
  const double upper=Binary64(lo+1);
  const int side=Compare(Minus(value,exactLower),Minus(Rational(upper),value));
  return side<0||(side==0&&(lo&1)==0)?lower:upper;
}
bool StoreResidual(const Rational& exact,CompensatedDepth* output,Real* errorBound) {
  // Unsupported non-IEEE hosts retain the legacy scalar API; this new native
  // two-binary64 representation fails explicitly rather than misreading bits.
  if(sizeof(double)!=sizeof(std::uint64_t)||!std::numeric_limits<double>::is_iec559||
      std::numeric_limits<double>::digits!=53||Binary64(UINT64_C(0x3ff0000000000000))!=1)
    return false;
  const double high=RoundBinary64(exact);const Rational exactHigh(high);
  const bool negative=Compare(exact,exactHigh)<0;
  const Rational leftover=negative?Minus(exactHigh,exact):Minus(exact,exactHigh);
  if(Compare(leftover,Rational(std::numeric_limits<double>::max()))>0)return false;
  const double magnitude=RoundBinary64(leftover);
  CompensatedDepth parts{high,negative?-magnitude:magnitude};Rational encoded;
  if(!ExactComponents(parts,&encoded))return false; // positive underflow is never an event
  Real bound=0;
  if(Compare(encoded,exact)!=0) {
    // The second nearest-even rounding is the ONLY residual storage error.
    // A whole adjacent spacing safely bounds its half-ulp error, including
    // subnormal/zero/binade endpoints. This is an exact power-of-two bound,
    // not an added exposure tolerance. Exact residuals retain bound0.
    const double up=std::nextafter(magnitude,std::numeric_limits<double>::infinity());
    const double spacing=std::isfinite(up)?up-magnitude:
        magnitude-std::nextafter(magnitude,0.0);
    bound=spacing;
    if(!(bound>0)||!Finite(bound))return false;
  }
  *output=parts;*errorBound=bound;return true;
}
Rational ExactArea(const Segment& segment,Real begin,Real end) {
  const Rational a(begin),b(end),s(segment.beginS),e(segment.endS);
  const auto left=Plus(Minus(e,a),Minus(e,b)),right=Plus(Minus(a,s),Minus(b,s));
  Rational endpoints;
  for(const auto& channel:segment.channels)
    endpoints=Plus(endpoints,Plus(Times(Rational(channel.beginRateS),left),Times(Rational(channel.endRateS),right)));
  return Divide(Times(Minus(b,a),endpoints),Times(Rational(2),Minus(e,s)));
}

// Products of positive factors are evaluated in exponent space. For example,
// r=1e4000 and dt=1e-4000 have a perfectly ordinary exposure; neither an
// overflowed intermediate nor a vanished normalized time fraction may replace
// that exposure. Every term below is an EXACT analytical affine contribution,
// subject only to declared long-double floating arithmetic, not quadrature.
struct Scaled {Real fraction=0;int exponent=0;};
Scaled Product(Real a,Real b,Real c,Real denominator,bool half) {
  if(a==0||b==0||c==0)return {};
  int ea=0,eb=0,ec=0,ed=0,adjust=0;
  const Real fa=std::frexp(a,&ea),fb=std::frexp(b,&eb),fc=std::frexp(c,&ec);
  const Real fd=std::frexp(denominator,&ed);
  const Real f=std::frexp((fa*fb*fc)/fd,&adjust);
  return {f,ea+eb+ec-ed+adjust-(half?1:0)};
}
enum class Range {Finite,Underflow,Overflow};
Range PositiveSum(const std::vector<Scaled>& terms,Real* value) {
  int maximum=std::numeric_limits<int>::min();
  for(const auto& t:terms)if(t.fraction>0)maximum=std::max(maximum,t.exponent);
  if(maximum==std::numeric_limits<int>::min()){*value=0;return Range::Finite;}
  Sum normalized;
  for(const auto& t:terms)if(t.fraction>0)
    if(!normalized.Add(std::scalbn(t.fraction,t.exponent-maximum)))return Range::Overflow;
  *value=std::scalbn(normalized.Value(),maximum);
  if(!Finite(*value))return Range::Overflow;
  return *value==0?Range::Underflow:Range::Finite;
}
void AreaTerms(const Segment& s,const Channel& c,Real a,Real b,
               std::vector<Scaled>* terms) {
  const Real duration=s.endS-s.beginS,length=b-a;
  // Integral = length/2 * [rate(a)+rate(b)]. Keeping the four endpoint
  // contributions separate avoids cancellation at a zero-rate endpoint and
  // underflow in (a-segmentBegin)/duration before multiplication by a rate.
  terms->push_back(Product(c.beginRateS,length,s.endS-a,duration,true));
  terms->push_back(Product(c.beginRateS,length,s.endS-b,duration,true));
  terms->push_back(Product(c.endRateS,length,a-s.beginS,duration,true));
  terms->push_back(Product(c.endRateS,length,b-s.beginS,duration,true));
}
Range Area(const Segment& s,Real a,Real b,Real* value) {
  std::vector<Scaled> terms;terms.reserve(4*s.channels.size());
  for(const auto& c:s.channels)AreaTerms(s,c,a,b,&terms);
  return PositiveSum(terms,value);
}
Range Rate(const Segment& s,const Channel& c,Real time,Real* value) {
  // Equal endpoints describe the same instantaneous physical rate throughout
  // the segment. Evaluating two weighted terms can round their sum away
  // from that value at an interior event. The reaction owner independently
  // re-evaluates the beam at that event and checks its rate identity.
  if (c.beginRateS==c.endRateS) {
    *value=c.beginRateS;
    return Range::Finite;
  }
  return PositiveSum({Product(c.beginRateS,1,s.endS-time,s.endS-s.beginS,false),
      Product(c.endRateS,1,time-s.beginS,s.endS-s.beginS,false)},value);
}
struct Window {std::size_t segment;Real begin,end;};
Status CheckHistory(const std::vector<Segment>& history,std::string* error);
Status Windows(const std::vector<Segment>& history,Real begin,Real end,
               std::vector<Window>* windows,std::string* error) {
  if(!Finite(begin)||!Finite(end)||begin<0||end<begin)
    return Refuse(Status::InvalidInput,"invalid nonnegative relative-time query window [s]",error);
  if(history.empty())return Refuse(Status::IncompleteSupport,"hazard history has no declared support",error);
  const auto validity=CheckHistory(history,error);
  if(validity!=Status::Complete)return validity;
  if(begin==end) {
    for(const auto& s:history)if(s.beginS<=begin&&begin<=s.endS)return Status::Complete;
    return Refuse(Status::IncompleteSupport,"zero-duration query lies outside declared support",error);
  }
  Real cursor=begin;
  for(std::size_t index=0;index<history.size();++index) {
    const auto& s=history[index];
    if(s.endS<=cursor)continue;
    if(s.beginS>cursor)break;
    const Real stop=std::min(s.endS,end);
    windows->push_back({index,cursor,stop});cursor=stop;
    if(cursor==end)return Status::Complete;
  }
  return Refuse(Status::IncompleteSupport,"query includes a gap or extends beyond explicit hazard support; no rate was inferred",error);
}
Status Arithmetic(const char* message,std::string* error) {
  return Refuse(Status::NumericalError,message,error);
}

// The analytical quadratic supplies a fast seed; strict analytic-exposure
// comparisons, not a root-tolerance window, finish the location. Extreme
// exponents can make this normalized seed unavailable even when the actual
// time is representable; bracket certification remains valid in that case.
Real Seed(const Segment& segment,Real begin,Real end,Real depth,Real area) {
  Sum first,last;
  for(const auto& c:segment.channels) {
    Real a=0,b=0;
    if(Rate(segment,c,begin,&a)!=Range::Finite||Rate(segment,c,end,&b)!=Range::Finite)
      return begin+(end-begin)/2;
    if(!first.Add(a)||!last.Add(b))return begin+(end-begin)/2;
  }
  const Real scale=std::max(first.Value(),last.Value());
  if(first.Value()>0&&first.Value()==last.Value()) {
    const Real constant=begin+depth/first.Value();
    if(constant>begin&&constant<end)return constant;
  }
  if(!(scale>0)||!Finite(area)||!(area>0))return begin+(end-begin)/2;
  const Real a=first.Value()/scale,b=last.Value()/scale;
  const Real target=(depth/area)*(a/2+b/2);
  const Real discriminant=std::fma(2*(b-a),target,a*a);
  if(!(target>0)||discriminant<0||!Finite(discriminant))return begin+(end-begin)/2;
  const Real denominator=a+std::sqrt(discriminant);
  const Real fraction=denominator>0?2*target/denominator:0;
  const Real seed=begin+fraction*(end-begin);
  return seed>begin&&seed<end?seed:begin+(end-begin)/2;
}

Status CheckHistory(const std::vector<Segment>& history,std::string* error) {
  if(history.empty())return Refuse(Status::InvalidInput,"history must declare at least one segment",error);
  for(std::size_t index=0;index<history.size();++index) {
    const auto& s=history[index];
    if(!Finite(s.beginS)||!Finite(s.endS)||s.beginS<0||!(s.endS>s.beginS)||
        !Finite(s.endS-s.beginS)||s.channels.empty()||
        (index>0&&s.beginS<history[index-1].endS)) {
      return Refuse(Status::InvalidInput,"segments need finite ordered nonoverlapping positive durations and explicit channels",error);
    }
    std::set<std::size_t> seen;
    for(std::size_t k=0;k<s.channels.size();++k) {
      const auto& c=s.channels[k];
      if(!seen.insert(c.reaction).second||!Finite(c.beginRateS)||!Finite(c.endRateS)||c.beginRateS<0||c.endRateS<0||
          s.channels.size()!=history.front().channels.size()||c.reaction!=history.front().channels[k].reaction) {
        return Refuse(Status::InvalidInput,"channel IDs/order must remain explicit and unique; endpoint rates must be finite nonnegative [s^-1]",error);
      }
    }
  }
  if(error)error->clear();
  return Status::Complete;
}
}

bool ValidateHistory(const std::vector<Segment>& history,std::string* error) {
  try {return CheckHistory(history,error)==Status::Complete;}
  catch(const std::exception& exception) {
    Refuse(Status::NumericalError,std::string("history validation staging failed: ")+exception.what(),error);
    return false;
  }
}

Status Integrate(const std::vector<Segment>& history,Real begin,Real end,
                 Exposure* output,std::string* error) {
  if(!output)return Refuse(Status::InvalidInput,"exposure output storage is required",error);
  try {
    std::vector<Window> windows;const auto status=Windows(history,begin,end,&windows,error);
    if(status!=Status::Complete)return status;
    Exposure value;std::vector<Sum> channelSums(history.front().channels.size());Sum total;
    for(const auto& w:windows) {
      const auto& s=history[w.segment];
      for(std::size_t k=0;k<s.channels.size();++k) {
        std::vector<Scaled> terms;terms.reserve(4);AreaTerms(s,s.channels[k],w.begin,w.end,&terms);Real area=0;
        if(PositiveSum(terms,&area)!=Range::Finite||!channelSums[k].Add(area)||!total.Add(area))
          return Arithmetic("requested affine exposure is outside long-double arithmetic range",error);
      }
    }
    value.integratedHazard=total.Value();
    if(!Finite(value.integratedHazard))return Arithmetic("combined exposure is unrepresentable",error);
    for(std::size_t k=0;k<channelSums.size();++k)
      value.channels.push_back({history.front().channels[k].reaction,channelSums[k].Value()});
    *output=std::move(value);if(error)error->clear();return Status::Complete;
  }
  catch(const std::exception& exception) {
    return Refuse(Status::NumericalError,std::string("exposure staging failed before publication: ")+exception.what(),error);
  }
}

static Status LocateExact(const std::vector<Segment>& history,Real begin,Real end,
    Real depth,const Rational& exactDepth,const NumericControls& controls,
    EventLocation* output,Rational* finalResidual,std::string* error) {
  if(!output||!Finite(depth)||!(depth>0)||controls.maximumRootRefinements==0)
    return Refuse(Status::InvalidInput,"event location needs output, positive finite optical depth and an explicit positive refinement limit",error);
  try {
    std::vector<Window> windows;const auto status=Windows(history,begin,end,&windows,error);
    if(status!=Status::Complete)return status;
    Sum remaining;remaining.high=depth;Sum accepted;
    Rational exactAccepted;
    EventLocation value;value.eventTimeS=std::numeric_limits<Real>::quiet_NaN();
    for(const auto& w:windows) {
      const auto& segment=history[w.segment];Real area=0;const auto range=Area(segment,w.begin,w.end,&area);
      const auto exactArea=ExactArea(segment,w.begin,w.end);
      const int comparison=Compare(Plus(exactAccepted,exactArea),exactDepth);
      if(comparison<0) {
        if(range==Range::Overflow)
          return Arithmetic("sub-threshold segment exposure overflowed its numerical representation",error);
        // A rounded prefix can equal the input depth although its exact
        // rational exposure is still smaller. That approximate remainder is
        // a SEED only; it must never override the exact crossing decision.
        // An individually underflowed segment remains in exactAccepted and
        // can therefore influence the first representable event timestamp.
        if(!remaining.Add(-area)||!accepted.Add(area))
          return Arithmetic("accumulated optical exposure is unrepresentable",error);
        exactAccepted=Plus(exactAccepted,exactArea);
        continue;
      }
      Real lower=w.begin,upper=w.end;
      // Exact endpoint ownership is explicit. In particular, a falling rate
      // ending at zero must not be declared early merely because its tiny
      // positive tail rounds away in a near-total cumulative integral.
      if(comparison!=0) {
        Real seed=Seed(segment,lower,upper,remaining.Value(),area);
        bool located=false;
        for(std::size_t attempt=0;attempt<controls.maximumRootRefinements;++attempt) {
          if(std::nextafter(lower,upper)==upper){located=true;break;}
          Real trial=attempt==0?seed:lower+(upper-lower)/2;
          if(!(trial>lower&&trial<upper))trial=std::nextafter(lower,upper);
          const int hit=Compare(Plus(exactAccepted,ExactArea(segment,w.begin,trial)),exactDepth);
          if(hit==0){upper=trial;located=true;break;}
          if(hit>=0)upper=trial;else lower=trial;
        }
        if(!located)return Arithmetic("affine event certification exhausted the explicitly supplied refinement limit",error);
      }
      value.eventTimeS=upper;value.segmentIndex=w.segment;value.integratedHazard=depth;
      Sum total;
      for(const auto& c:segment.channels) {
        Real rate=0;
        if(Rate(segment,c,upper,&rate)!=Range::Finite||!total.Add(rate))
          return Arithmetic("instantaneous event rate is outside long-double arithmetic range",error);
        value.instantaneousRates.push_back({c.reaction,rate});
      }
      value.totalInstantaneousRateS=total.Value();
      if(finalResidual)*finalResidual=Rational();
      *output=std::move(value);if(error)error->clear();return Status::Event;
    }
    value.integratedHazard=accepted.Value();
    if(value.integratedHazard==0&&Compare(exactAccepted,Rational())>0)
      return Arithmetic("positive full-window exposure underflows long-double output storage",error);
    if(finalResidual)*finalResidual=Minus(exactDepth,exactAccepted);
    *output=std::move(value);if(error)error->clear();return Status::Complete;
  }
  catch(const std::exception& exception) {
    return Refuse(Status::NumericalError,std::string("event staging failed before publication: ")+exception.what(),error);
  }
}

Status LocateEvent(const std::vector<Segment>& history,Real begin,Real end,
    Real depth,const NumericControls& controls,EventLocation* output,std::string* error) {
  if(!output||!Finite(depth)||!(depth>0)||controls.maximumRootRefinements==0)
    return Refuse(Status::InvalidInput,"event location needs output, positive finite optical depth and an explicit positive refinement limit",error);
  try{return LocateExact(history,begin,end,depth,Rational(depth),controls,output,nullptr,error);}
  catch(const std::exception& exception){return Refuse(Status::NumericalError,std::string("scalar depth staging failed: ")+exception.what(),error);}
}

Status LocateEvent(const std::vector<Segment>& history,Real begin,Real end,
    const CompensatedDepth& parts,const NumericControls& controls,
    EventLocation* output,std::string* error) {
  if(!output||controls.maximumRootRefinements==0)
    return Refuse(Status::InvalidInput,"compensated event location needs output and an explicit positive refinement limit",error);
  try {
    Rational depth;
    if(!ExactComponents(parts,&depth))
      return Refuse(Status::InvalidInput,"pending binary components must define an exact positive finite native optical depth",error);
    // This scalar is a seed/display ONLY. All threshold/root comparisons in
    // LocateExact use the full arbitrary-size rational above.
    const Real display=static_cast<Real>(parts.high)+parts.low;
    EventLocation value;Rational remaining;
    const auto status=LocateExact(history,begin,end,display,depth,controls,&value,&remaining,error);
    if(status!=Status::Complete&&status!=Status::Event)return status;
    CompensatedDepth residual;Real bound=0;
    if(status==Status::Complete&&!StoreResidual(remaining,&residual,&bound))
      return Arithmetic("positive exact residual cannot be represented by the native two-double optical clock; no zero/event substitution was made",error);
    value.clockUpdate.valid_=true;value.clockUpdate.event_=status==Status::Event;
    value.clockUpdate.origin_=parts;value.clockUpdate.residual_=residual;
    value.clockUpdate.roundingBound_=bound;
    *output=std::move(value);if(error)error->clear();return status;
  }catch(const std::exception& exception){return Refuse(Status::NumericalError,std::string("compensated depth staging failed: ")+exception.what(),error);}
}
}}}

namespace Moon { namespace Species { namespace ExactArithmetic {
bool SignedSumMatchesDecimalScale(const std::vector<long double>& terms,
    long double sourceValue,std::uint64_t numerator,unsigned decimalPlaces,
    long double absoluteRepresentationEnvelope,bool* matches,std::string* error) {
  // Reuse this translation unit's private exact-rational implementation. In
  // particular, a large common energy zero must not enlarge a small gap's
  // permitted mismatch through cancellation in floating addition.
  using namespace PiecewiseHazard;
  if(!matches||terms.empty()||!std::isfinite(sourceValue)||sourceValue<0||
      !std::isfinite(absoluteRepresentationEnvelope)||absoluteRepresentationEnvelope<0||
      decimalPlaces>10000) {
    if(error)*error="exact scaled-sum comparison requires finite inputs/output and decimal exponent at most10000";
    return false;
  }
  try {
    Rational positive,negative;
    for(long double term:terms) {
      if(!std::isfinite(term)){if(error)*error="nonfinite exact signed-sum term";return false;}
      if(term<0)negative=Plus(negative,Rational(-term));
      else positive=Plus(positive,Rational(term));
    }
    Rational scale;scale.numerator=Integer(numerator);
    for(unsigned i=0;i<decimalPlaces;++i)scale.denominator=Multiply(scale.denominator,Integer(10));
    scale.Normalize();
    const auto target=Times(Rational(sourceValue),scale);
    const auto right=Plus(negative,target);
    const auto difference=Compare(positive,right)<0?Minus(right,positive):Minus(positive,right);
    const bool candidate=Compare(difference,Rational(absoluteRepresentationEnvelope))<=0;
    *matches=candidate;if(error)error->clear();return true;
  }catch(const std::exception& exception) {
    if(error)*error=std::string("exact scaled-sum staging failed: ")+exception.what();
    return false;
  }
}
}}}
