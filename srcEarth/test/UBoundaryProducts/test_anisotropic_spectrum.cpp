#include "../../gridless/AnisotropicSpectrum.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

// AnisotropicSpectrum.cpp only needs this small parser utility.  Defining it in
// the focused unit-test executable lets us compile and link the production source
// without pulling in the complete AMPS/SPICE parameter-parser dependency graph.
namespace EarthUtil {
std::string ToUpper(std::string value) {
  std::transform(value.begin(),value.end(),value.begin(),
      [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return value;
}
}

namespace {

int failures=0;

bool Near(double actual,double expected,double relative=1.0e-13,
          double absolute=1.0e-14) {
  return std::fabs(actual-expected)<=absolute+
      relative*std::max(std::fabs(actual),std::fabs(expected));
}

void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

} // namespace

int main() {
  const double dayside[3]={1.0,0.0,0.0};
  const double nightside[3]={-1.0,0.0,0.0};
  EarthUtil::AnisotropyParam par;

  // sin^2(alpha) has an exact full-sphere mean of 2/3.  This comparison
  // exercises the production EvalAnisotropyFactor ordering: first calculate
  // the raw PAD value and only then divide by its analytic mean.
  par.padModel="SINALPHA_N";
  par.padExponent=2.0;
  par.padNormalization="RAW";
  Check(Near(EvalAnisotropyFactor(par,0.6,dayside),0.64),
        "raw sin^2 PAD matches 1-cos^2(alpha)");
  par.padNormalization="UNIT_MEAN";
  Check(Near(EvalAnisotropyFactor(par,0.6,dayside),0.96),
        "unit-mean sin^2 PAD divides raw value by 2/3");

  // cos^2(alpha) has mean 1/3.  This covers the other analytic normalizer
  // branch and catches accidental normalization before f_pad is assigned.
  par.padModel="COSALPHA_N";
  Check(Near(EvalAnisotropyFactor(par,0.5,dayside),0.75),
        "unit-mean cos^2 PAD divides raw value by 1/3");

  // PAD and spatial normalizations are deliberately independent and multiply.
  // For day:night factors 3:1 the hemispheric mean is 2, hence normalized
  // weights 1.5 and 0.5.
  par.padModel="ISOTROPIC";
  par.spatialModel="DAYSIDE_NIGHTSIDE";
  par.daysideFactor=3.0;
  par.nightsideFactor=1.0;
  par.spatialNormalization="UNIT_MEAN";
  Check(Near(EvalAnisotropyFactor(par,0.0,dayside),1.5) &&
        Near(EvalAnisotropyFactor(par,0.0,nightside),0.5),
        "normalized day/night weights retain unit hemispheric mean");

  // Preserve the documented defensive behavior for an undefined pitch angle.
  par.padModel="COSALPHA_N";
  par.padNormalization="RAW";
  par.spatialModel="UNIFORM";
  par.spatialNormalization="RAW";
  Check(Near(EvalAnisotropyFactor(
                 par,std::numeric_limits<double>::quiet_NaN(),dayside),0.0),
        "NaN pitch-angle cosine uses the perpendicular fallback");

  bool rejected=false;
  try {
    par.padModel="NOT_A_PAD";
    (void)EvalAnisotropyFactor(par,0.0,dayside);
  }
  catch (const std::runtime_error&) {
    rejected=true;
  }
  Check(rejected,"unknown PAD model is rejected");

  if (failures!=0) {
    std::cerr << failures << " AnisotropicSpectrum test(s) failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: production AnisotropicSpectrum raw/normalized PAD and spatial weights\n";
  return EXIT_SUCCESS;
}
