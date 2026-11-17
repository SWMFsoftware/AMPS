#ifndef SRCSEP_UTIL_SEP_EXACT_NONNEGATIVE_SUM_H
#define SRCSEP_UTIL_SEP_EXACT_NONNEGATIVE_SUM_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <array>
#include <cstdint>

namespace SEP { namespace Transport {
// Exact sum of finite nonnegative IEEE754 binary64 inputs. Words are an
// unsigned integer in units2^-1074 of the caller's declared SI quantity,
// least-significant word first. 34 words contain all2098 meaningful bits plus
// carry headroom. The supported exact total is <=largest finite binary64;
// exceeding that representation contract rejects before replacing the input.
// No data-dependent cutoff, floating compensation loss or physical tolerance
// enters addition. This is an account/diagnostic accumulator, not wave dynamics.
struct ExactNonnegativeSum {std::array<std::uint64_t,34> words{{0}};};
struct ExactSumApproximation {
  double nearest=0,signedRemainder=0;
  // Bound on |exact-(nearest+signedRemainder)|. Zero means both diagnostic
  // terms are exact. Otherwise one residual ULP bounds its nearest rounding;
  // this conservative bound remains representable even at subnormal scales.
  double absoluteRoundingBound=0;
};
Status ValidateExactNonnegativeSum(const ExactNonnegativeSum& sum);
Status AddExactNonnegative(const ExactNonnegativeSum& before,double increment,ExactNonnegativeSum* after);
Status ApproximateExactNonnegative(const ExactNonnegativeSum& sum,ExactSumApproximation* approximation);
} }
#endif
