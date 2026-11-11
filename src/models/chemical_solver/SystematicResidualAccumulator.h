#ifndef _CHEMICAL_SOLVER_SYSTEMATIC_RESIDUAL_ACCUMULATOR_H_
#define _CHEMICAL_SOLVER_SYSTEMATIC_RESIDUAL_ACCUMULATOR_H_

#include <cmath>
#include <cstdint>
#include <limits>

namespace ChemicalSolver {
namespace SystematicResidualAccumulator {

// A residual update is planned before a caller reserves any particle slots.
// The caller commits NextPhase only after the complete physical transaction
// has passed its capacity preflight.  Keeping this small value type free of
// AMPS mesh dependencies makes the no-mutation-on-deferral contract directly
// testable and reusable by other chemistry applications.
struct cPlan {
  long int Multiplicity;
  double NextPhase;

  cPlan(long int multiplicity=0,double nextPhase=0.0) :
      Multiplicity(multiplicity),NextPhase(nextPhase) {}
};

inline bool IsValidPhase(double phase) {
  return std::isfinite(phase) && phase>=0.0 && phase<1.0;
}

inline cPlan Plan(double phase,double mean) {
  if (!IsValidPhase(phase) || !(mean>=0.0) || !std::isfinite(mean))
    return cPlan(-1,std::numeric_limits<double>::quiet_NaN());

  const double accumulated=phase+mean;
  const double integral=std::floor(accumulated);
  const double integerLimit=static_cast<double>(
      std::numeric_limits<long int>::max());
  // On LP64, double(LONG_MAX) rounds upward to 2^63.  Reject that boundary
  // before converting to long; otherwise the cast itself is out of range.
  const bool integralOutOfRange=
      std::numeric_limits<long int>::digits>
          std::numeric_limits<double>::digits ?
      integral>=integerLimit : integral>integerLimit;
  if (!std::isfinite(accumulated) || integralOutOfRange)
    return cPlan(-1,std::numeric_limits<double>::quiet_NaN());
  double next=accumulated-integral;

  // Roundoff can put a value that is mathematically zero at the upper edge of
  // the half-open interval. Normalize it without changing the integer count.
  if (next>=1.0) next=0.0;
  if (next<0.0 && next>
      -64.0*std::numeric_limits<double>::epsilon()) next=0.0;
  if (!IsValidPhase(next))
    return cPlan(-1,std::numeric_limits<double>::quiet_NaN());

  return cPlan(static_cast<long int>(integral),next);
}

// SplitMix64 is used as a counter-based phase generator, not as a physical-
// event RNG. A phase keyed by a science seed and stable model identifiers is
// independent of MPI owner and first-use order. Enabling systematic residual
// sampling therefore need not consume or reorder unrelated random draws.
inline std::uint64_t Mix(std::uint64_t value) {
  value+=UINT64_C(0x9e3779b97f4a7c15);
  value=(value^(value>>30))*UINT64_C(0xbf58476d1ce4e5b9);
  value=(value^(value>>27))*UINT64_C(0x94d049bb133111eb);
  return value^(value>>31);
}

inline std::uint64_t Combine(std::uint64_t state,std::uint64_t value) {
  return Mix(state^Mix(value));
}

inline double UnitPhase(std::uint64_t key) {
  // Use the 53 high-quality mantissa bits. The result is exactly in [0,1)
  // and is reproducible across IEEE-754 hosts.
  return static_cast<double>(Mix(key)>>11)*
      (1.0/static_cast<double>(UINT64_C(1)<<53));
}

} // namespace SystematicResidualAccumulator
} // namespace ChemicalSolver

#endif
