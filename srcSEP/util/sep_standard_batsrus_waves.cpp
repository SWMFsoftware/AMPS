#include "sep_standard_batsrus_waves.h"
#include <cmath>
#include <cstring>
#include <limits>

namespace SEP { namespace Waves { namespace StandardBatsrus { namespace {
namespace T=Transport;
namespace CP=Transport::CoefficientPhysics;
void HashWord(std::uint64_t* hash,std::uint64_t value) {
  // Canonical little-endian scalar hashing makes the source identity
  // independent of host byte order. This is a view digest, not a wire format.
  for(unsigned i=0;i<8;++i){*hash^=static_cast<unsigned char>(value>>(8*i));*hash*=UINT64_C(1099511628211);}
}
void HashDouble(std::uint64_t* hash,double value) {
  static_assert(sizeof(double)==sizeof(std::uint64_t),"standard BATSRUS view requires binary64-sized doubles");
  std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));HashWord(hash,bits);
}
T::Status ConvertEndpoint(const RawTotalAlfvenEndpoint& input,std::array<double,2>* output) {
  std::array<double,2> magnetic;
  for(unsigned branch=0;branch<2;++branch) {
    const double raw=input.waveJPerM3[branch];
    if(!std::isfinite(raw)||raw<0.0)return T::Status::Error(T::StatusCode::InvalidCoefficient,
        "standard BATSRUS raw total Alfven energy endpoints must be finite and nonnegative [J/m^3]");
    magnetic[branch]=0.5*raw;
    if(raw>0.0&&magnetic[branch]==0.0)return T::Status::Error(T::StatusCode::StepUnderflow,
        "positive standard BATSRUS total energy has no representable half-energy magnetic view");
  }
  const auto variance=CP::MagneticVarianceFromEnergyDensity(magnetic[0],magnetic[1]);
  if(!variance.status.ok())return variance.status;
  *output=magnetic;return T::Status::Ok();
}
} // namespace
const char* ContractIdentity() {
  return "standard-batsrus-total-alfven-energy-v1;Wave1=along-B;Wave2=opposite-B;magnetic-fraction=1/2;L/R=unspecified";
}
const char* Representation() {
  return "swmf-standard-batsrus-total-to-magnetic-half;Wave1=E+along-B;Wave2=E-opposite-B;L/R=unspecified";
}
T::Status InterpolateMagneticView(const RawTotalAlfvenEndpoint& begin,
    const RawTotalAlfvenEndpoint& end,double fraction,std::uint64_t generation,
    MagneticBranchView* output) {
  if(!output||!std::isfinite(fraction)||fraction<0.0||fraction>1.0||generation==0)
    return T::Status::Error(T::StatusCode::InvalidArgument,"standard BATSRUS view needs output, closed-unit interpolation fraction and nonzero import generation");
  std::array<double,2> first,last;auto status=ConvertEndpoint(begin,&first);if(!status.ok())return status;
  status=ConvertEndpoint(end,&last);if(!status.ok())return status;
  MagneticBranchView candidate;double interpolated[2];
  for(unsigned branch=0;branch<2;++branch) {
    // Wider convex arithmetic avoids overflowing raw endpoint sums and exposes
    // a positive value that would otherwise silently vanish on narrowing.
    const long double value=(1.0L-fraction)*first[branch]+static_cast<long double>(fraction)*last[branch];
    if(!std::isfinite(value)||value>std::numeric_limits<double>::max())
      return T::Status::Error(T::StatusCode::InvalidCoefficient,"standard BATSRUS magnetic interpolation overflowed");
    interpolated[branch]=static_cast<double>(value);
    if(value>0.0L&&interpolated[branch]==0.0)
      return T::Status::Error(T::StatusCode::StepUnderflow,"positive standard BATSRUS interpolated magnetic branch is unrepresentable");
  }
  const auto variance=CP::MagneticVarianceFromEnergyDensity(interpolated[0],interpolated[1]);
  if(!variance.status.ok())return variance.status;
  candidate.plusMagneticJPerM3=interpolated[0];candidate.minusMagneticJPerM3=interpolated[1];
  candidate.plusVarianceT2=variance.plusT2;candidate.minusVarianceT2=variance.minusT2;candidate.totalVarianceT2=variance.totalT2;
  candidate.generation=generation;candidate.branchDirectionsKnown=true;candidate.circularPolarizationsKnown=false;
  std::uint64_t hash=UINT64_C(1469598103934665603);
  for(const char* letter=ContractIdentity();*letter;++letter){hash^=static_cast<unsigned char>(*letter);hash*=UINT64_C(1099511628211);}
  HashWord(&hash,generation);HashDouble(&hash,fraction);
  for(double value:begin.waveJPerM3)HashDouble(&hash,value);
  for(double value:end.waveJPerM3)HashDouble(&hash,value);
  candidate.checksum=hash?hash:UINT64_C(1);
  *output=candidate;return T::Status::Ok();
}
} } }
