#ifndef SRCMOON_SPECIES_EXACT_ARITHMETIC_H
#define SRCMOON_SPECIES_EXACT_ARITHMETIC_H

#include <cstdint>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace ExactArithmetic {

// Compare an EXACT signed sum of finite native binary inputs with
// sourceValue*numerator/10^decimalPlaces. No intermediate sum/product is
// rounded. The caller owns the dimensional meaning and a separately justified
// nonnegative absolute representation envelope. This numerical service does
// not establish a physical uncertainty or authorize an acceptance tolerance.
// A false return is an arithmetic/input failure; a successful false `matches`
// is an evaluated disagreement. Outputs are unchanged on failure. Allocation
// failures are explicit, and no third-party arbitrary-integer library is used.
// Implementation shares the exact integer engine in PiecewiseHazard.cpp.
bool SignedSumMatchesDecimalScale(const std::vector<long double>& terms,
    long double sourceValue,std::uint64_t numerator,unsigned decimalPlaces,
    long double absoluteRepresentationEnvelope,bool* matches,std::string* error);

}}}
#endif
