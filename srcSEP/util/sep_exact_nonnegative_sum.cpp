#include "sep_exact_nonnegative_sum.h"
#include <cmath>
#include <cstring>
#include <limits>

namespace SEP { namespace Transport { namespace {
static_assert(sizeof(double)==sizeof(std::uint64_t)&&std::numeric_limits<double>::is_iec559&&
    std::numeric_limits<double>::radix==2&&std::numeric_limits<double>::digits==53&&
    std::numeric_limits<double>::max_exponent==1024&&std::numeric_limits<double>::min_exponent==-1021,
    "exact account requires IEEE754 binary64");
Status Error(const char* text){return Status::Error(StatusCode::UnresolvedCoefficient,text);}
int Compare(const ExactNonnegativeSum& a,const ExactNonnegativeSum& b) {
  for(std::size_t i=a.words.size();i-->0;)if(a.words[i]!=b.words[i])return a.words[i]<b.words[i]?-1:1;
  return 0;
}
bool AddWord(ExactNonnegativeSum* sum,std::size_t word,std::uint64_t value) {
  // Unsigned overflow is defined modulo2^64; the comparison detects its
  // single carry. All arithmetic stays in integer space, including increments
  // far smaller than the diagnostic floating-point total's ULP.
  while(value) {
    if(word>=sum->words.size())return false;
    const auto old=sum->words[word];sum->words[word]+=value;
    value=sum->words[word]<old?1:0;++word;
  }
  return true;
}
ExactNonnegativeSum Encode(double value) {
  // A normal binary64 x has integer significand m and value
  // m*2^(exponentBits-1075), hence shift exponentBits-1 in our2^-1074 units.
  // Subnormals already use those units. Callers validate finite nonnegative
  // inputs before entering here; no floating scaling can underflow a bit.
  ExactNonnegativeSum result;std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));
  const auto exponent=static_cast<unsigned>((bits>>52)&2047);
  const std::uint64_t fraction=bits&UINT64_C(0xfffffffffffff);
  const std::uint64_t significand=exponent?fraction|UINT64_C(0x10000000000000):fraction;
  const std::size_t shift=exponent?exponent-1:0,word=shift/64,offset=shift%64;
  AddWord(&result,word,significand<<offset);
  if(offset)AddWord(&result,word+1,significand>>(64-offset));
  return result;
}
ExactNonnegativeSum Difference(const ExactNonnegativeSum& larger,const ExactNonnegativeSum& smaller) {
  // Precondition larger>=smaller follows from Compare at the call site. The
  // exact difference is used only for display-error diagnostics, never to
  // subtract an escaped physical increment from the authoritative account.
  ExactNonnegativeSum result;std::uint64_t borrow=0;
  for(std::size_t i=0;i<result.words.size();++i) {
    const auto first=larger.words[i]-smaller.words[i];
    const auto firstBorrow=larger.words[i]<smaller.words[i];
    result.words[i]=first-borrow;const bool next=firstBorrow||first<borrow;borrow=next?1:0;
  }
  return result;
}
bool Bit(const ExactNonnegativeSum& sum,int bit) {
  return bit>=0&&(sum.words[static_cast<std::size_t>(bit)/64]>>(bit%64)&1)!=0;
}
double Nearest(const ExactNonnegativeSum& sum) {
  // Construct the IEEE bits explicitly, so the current floating rounding mode
  // cannot affect this conversion. The first omitted bit plus all remaining
  // sticky bits implements round-to-nearest, ties-to-even. Validate has already
  // excluded values above DBL_MAX, so rounding cannot produce infinity.
  int top=-1;
  for(std::size_t word=sum.words.size();word-->0;)if(sum.words[word]) {
    int bit=63;while((sum.words[word]>>bit&1)==0)--bit;top=static_cast<int>(word*64)+bit;break;
  }
  if(top<0)return 0;
  std::uint64_t bits=0;
  if(top<52)bits=sum.words[0]; // Exact subnormal integer significand.
  else {
    std::uint64_t significand=0;
    for(int bit=top;bit>=top-52;--bit)significand=(significand<<1)|(Bit(sum,bit)?1:0);
    const int shift=top-52;bool sticky=false;
    for(int bit=0;bit<shift-1;++bit)sticky=sticky||Bit(sum,bit);
    if(shift>0&&Bit(sum,shift-1)&&(sticky||(significand&1)))++significand;
    if(significand==UINT64_C(0x20000000000000)){significand>>=1;++top;}
    const auto exponent=static_cast<std::uint64_t>(top-51);
    bits=(exponent<<52)|(significand&UINT64_C(0xfffffffffffff));
  }
  double result=0;std::memcpy(&result,&bits,sizeof(result));return result;
}
} // namespace

Status ValidateExactNonnegativeSum(const ExactNonnegativeSum& sum) {
  const auto maximum=Encode(std::numeric_limits<double>::max());
  return Compare(sum,maximum)<=0?Status::Ok():Error("exact external account exceeds its declared finite binary64 total bound");
}
Status AddExactNonnegative(const ExactNonnegativeSum& before,double increment,ExactNonnegativeSum* after) {
  auto status=ValidateExactNonnegativeSum(before);if(!status.ok())return status;
  if(!after||!std::isfinite(increment)||increment<0)return Error("exact external account needs a finite nonnegative increment/destination");
  // Work on a value copy: null/nonfinite/negative/overflow rejection leaves
  // both input and caller destination unchanged, including aliasing after=&before.
  auto candidate=before;const auto encoded=Encode(increment);
  for(std::size_t i=0;i<candidate.words.size();++i)
    if(!AddWord(&candidate,i,encoded.words[i]))return Error("exact external account integer carry overflow");
  status=ValidateExactNonnegativeSum(candidate);if(!status.ok())return status;
  *after=candidate;return Status::Ok();
}
Status ApproximateExactNonnegative(const ExactNonnegativeSum& sum,ExactSumApproximation* output) {
  auto status=ValidateExactNonnegativeSum(sum);if(!status.ok())return status;
  if(!output)return Error("exact account diagnostic destination is absent");
  ExactSumApproximation result;result.nearest=Nearest(sum);const auto rounded=Encode(result.nearest);
  const int sign=Compare(sum,rounded);
  const auto residual=sign<0?Difference(rounded,sum):Difference(sum,rounded);
  const double magnitude=Nearest(residual);result.signedRemainder=sign<0?-magnitude:magnitude;
  // A second nearest rounding may itself lose low bits (for example1+2^-54
  // +2^-108). Keep those bits in the integer authority and bound the display
  // error by one residual ULP; the two displayed terms are NOT called exact.
  if(Compare(residual,Encode(magnitude))!=0)
    result.absoluteRoundingBound=std::nextafter(magnitude,std::numeric_limits<double>::infinity())-magnitude;
  if(!std::isfinite(result.nearest)||!std::isfinite(result.signedRemainder)||!std::isfinite(result.absoluteRoundingBound))
    return Error("exact account diagnostics exceed finite representation");
  *output=result;return Status::Ok();
}
} }
