#include "sep_injection_spectrum.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP {
namespace Injection {
namespace {

Transport::Status Error(const std::string& text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,text);
}

std::uint64_t Mix(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

std::uint64_t Append(std::uint64_t state, std::uint64_t tag,
                     std::uint64_t value) {
  // Hash-combine the tagged value into the previous state rather than XORing
  // independent hashes.  XOR is commutative and made distinct semantic tuples
  // collide whenever two fields were swapped.
  return Mix(state ^ Mix(tag) ^ (Mix(value)+UINT64_C(0x517cc1b727220a95)+
      (state<<6)+(state>>2)));
}

Transport::Status ValidateSpectrum(const Spectrum& s) {
  if (!std::isfinite(s.minimum) || !std::isfinite(s.maximum) ||
      !std::isfinite(s.powerIndex) || s.minimum<=0.0 || s.maximum<=s.minimum)
    return Error("injection spectrum requires 0 < minimum < maximum and finite index");
  return Transport::Status::Ok();
}

// Stable integral of x^q between positive bounds.  The q=-1 logarithmic limit
// is evaluated explicitly and expm1 avoids cancellation near that limit.
double PowerIntegral(double lo, double hi, double q) {
  const double exponent=q+1.0;
  if (std::fabs(exponent)<1.0e-12) return std::log(hi/lo);
  return std::pow(lo,exponent)*std::expm1(exponent*std::log(hi/lo))/exponent;
}

double EffectiveCoordinateIndex(const Spectrum& s) {
  // A density per logarithmic coordinate acquires one factor of x relative to
  // the corresponding density per linear coordinate.
  return (s.measure==Measure::LogMomentum ||
          s.measure==Measure::LogKineticEnergy)
      ? s.powerIndex-1.0 : s.powerIndex;
}

}  // namespace

Transport::Status Validate(const Configuration& c) {
  const Transport::Status spectrum=ValidateSpectrum(c.spectrum);
  if (!spectrum.ok()) return spectrum;
  if (c.macroparticlesPerEvent==0 || !std::isfinite(c.injectionEfficiency) ||
      c.injectionEfficiency<0.0 || c.injectionEfficiency>1.0)
    return Error("injection count must be positive and efficiency must be in [0,1]");
  return Transport::Status::Ok();
}

Transport::ScalarResult ProbabilityDensity(const Spectrum& s, double x) {
  Transport::ScalarResult result;
  result.status=ValidateSpectrum(s);
  if (!result.status.ok() || !std::isfinite(x) || x<s.minimum || x>s.maximum) {
    if (result.status.ok())
      result.status=Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                             "spectrum coordinate is outside bounds");
    return result;
  }
  const double index=EffectiveCoordinateIndex(s);
  const double normalization=PowerIntegral(s.minimum,s.maximum,-index);
  result.value=std::pow(x,-index)/normalization;
  result.status=std::isfinite(result.value) && result.value>=0.0
      ? Transport::Status::Ok() : Error("spectrum density is not finite");
  return result;
}

Transport::ScalarResult InverseCdf(const Spectrum& s, double u) {
  Transport::ScalarResult result;
  result.status=ValidateSpectrum(s);
  if (!result.status.ok() || !std::isfinite(u) || !(u>0.0 && u<1.0)) {
    if (result.status.ok()) result.status=Error("inverse CDF requires u in (0,1)");
    return result;
  }
  const double exponent=1.0-EffectiveCoordinateIndex(s);
  if (std::fabs(exponent)<1.0e-12)
    result.value=s.minimum*std::exp(u*std::log(s.maximum/s.minimum));
  else {
    const double lo=std::pow(s.minimum,exponent);
    const double hi=std::pow(s.maximum,exponent);
    result.value=std::pow(lo+u*(hi-lo),1.0/exponent);
  }
  result.status=std::isfinite(result.value) ? Transport::Status::Ok()
                                            : Error("inverse CDF overflowed");
  return result;
}

Transport::ScalarResult ImportanceWeightForLogUniformProposal(
    const Spectrum& target, double x) {
  Transport::ScalarResult result=ProbabilityDensity(target,x);
  if (!result.status.ok()) return result;
  // q_log(x)=1/[x log(max/min)].  The ratio below is dimensionless and positive;
  // it is valid for any of the explicitly named target measures.
  result.value*=x*std::log(target.maximum/target.minimum);
  return result;
}

std::uint64_t HashRandomKey(const RandomKey& key) {
  std::uint64_t state=UINT64_C(0x5352435345504b59); // "SRCSEPKY"
  state=Append(state,1,key.campaign);
  state=Append(state,2,key.event);
  state=Append(state,3,key.fieldLine);
  state=Append(state,4,key.species);
  state=Append(state,5,key.macroparticle);
  state=Append(state,6,static_cast<std::uint64_t>(key.purpose));
  return state;
}

Transport::KeyedRandomStream MakeRandomStream(const RandomKey& key) {
  const std::uint64_t hash=HashRandomKey(key);
  // The semantic tuple is already reduced with tagged, ordered combination.
  // A fixed operator tag separates source streams from mover streams.
  return Transport::KeyedRandomStream(hash,key.macroparticle,
                                      UINT64_C(0x534f55524345),key.event);
}

std::string Fingerprint(const Configuration& c) {
  std::ostringstream canonical;
  canonical << std::setprecision(17) << c.campaignSeed << '|'
      << c.macroparticlesPerEvent << '|' << c.injectionEfficiency << '|'
      << static_cast<int>(c.spectrum.measure) << '|' << c.spectrum.minimum << '|'
      << c.spectrum.maximum << '|' << c.spectrum.powerIndex << '|'
      << static_cast<int>(c.angular);
  std::ostringstream out;
  std::uint64_t hash=UINT64_C(0x535045435452554d);
  const std::string bytes=canonical.str();
  for (std::size_t i=0;i<bytes.size();++i)
    hash=Append(hash,UINT64_C(7),static_cast<unsigned char>(bytes[i]));
  out << std::hex << std::setw(16) << std::setfill('0')
      << Mix(hash);
  return out.str();
}

}  // namespace Injection
}  // namespace SEP
