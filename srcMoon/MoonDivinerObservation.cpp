#include "MoonDivinerObservation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace Moon {
namespace DivinerObservation {
namespace {

constexpr double kPlanckJs=6.62607015e-34;
constexpr double kLightSpeedMs=299792458.0;
constexpr double kBoltzmannJK=1.380649e-23;
constexpr double kStefanBoltzmannWm2K4=5.670374419e-8;
constexpr double kPi=3.141592653589793238462643383279502884;

struct ChannelDefinition {
  int channel;
  // Paige et al. spectral partition used to reconstruct bolometric exitance.
  double partitionMinimumMicrometre;
  double partitionMaximumMicrometre;
  // Measured half-maximum detector passband used only for diagnostics.
  double passbandMinimumMicrometre;
  double passbandMaximumMicrometre;
  // A negative value denotes the always-required longest-wave channel C9.
  double minimumDetectableTemperatureK;
};

// Paige et al. (2010) supplementary methods, Table in lines 54--63.  The C3
// partition begins at zero even though its detector responds near 8 um: the
// measured C3 brightness temperature represents all flux shortward of the
// C3/C4 boundary in the bolometric reconstruction.  C9 terminates at 1000 um,
// exactly as published; flux beyond 1000 um is therefore not manufactured.
constexpr std::array<ChannelDefinition,InfraredChannelCount> kChannels{{
    {3,0.0,8.075,7.55,8.05,190.0},
    {4,8.075,8.40,8.10,8.40,190.0},
    {5,8.40,13.0,8.38,8.68,180.0},
    {6,13.0,25.0,13.0,23.0,95.0},
    {7,25.0,50.0,25.0,41.0,60.0},
    {8,50.0,100.0,50.0,100.0,40.0},
    {9,100.0,1000.0,100.0,400.0,-1.0}}};

const ChannelDefinition* Definition(int channel) {
  if (channel<FirstInfraredChannel || channel>LastInfraredChannel) return nullptr;
  return &kChannels[static_cast<std::size_t>(channel-FirstInfraredChannel)];
}

// Planck spectral radiance per metre of wavelength.  expm1 retains precision
// in the long-wave/small-exponent limit.  Very large exponents have negligible
// radiance and are returned as zero before exp() could overflow.
double PlanckRadiancePerM(double wavelengthM,double temperatureK) {
  if (!(wavelengthM>0.0) || !(temperatureK>0.0)) return 0.0;
  const double exponent=kPlanckJs*kLightSpeedMs/
      (wavelengthM*kBoltzmannJK*temperatureK);
  if (exponent>700.0) return 0.0;
  return 2.0*kPlanckJs*kLightSpeedMs*kLightSpeedMs/
      (std::pow(wavelengthM,5)*std::expm1(exponent));
}

// Deterministic composite Simpson integration.  The wavelength bands are
// smooth, positive Planck functions.  4096 even panels resolve the narrowest
// 0.3-um channel and the 900-um C9 partition well below the U10 1e-9 relative
// verification tolerance without an adaptive, platform-dependent branch tree.
double IntegratePlanck(double minimumMicrometre,double maximumMicrometre,
                       double temperatureK) {
  if (!(maximumMicrometre>minimumMicrometre) || !(temperatureK>0.0) ||
      !std::isfinite(temperatureK)) return std::numeric_limits<double>::quiet_NaN();
  constexpr int panels=4096;
  // The zero lower limit is physically integrable but cannot be evaluated
  // directly.  At 1e-6 um the exponent makes the omitted contribution exactly
  // negligible at all lunar temperatures while avoiding division by zero.
  const double lowerM=std::max(minimumMicrometre,1.0e-6)*1.0e-6;
  const double upperM=maximumMicrometre*1.0e-6;
  const double step=(upperM-lowerM)/panels;
  double sum=PlanckRadiancePerM(lowerM,temperatureK)+
      PlanckRadiancePerM(upperM,temperatureK);
  for (int panel=1;panel<panels;++panel) {
    const double value=PlanckRadiancePerM(lowerM+panel*step,temperatureK);
    sum+=(panel%2==0 ? 2.0 : 4.0)*value;
  }
  return sum*step/3.0;
}

bool ValidTemperature(double value) {
  return std::isfinite(value) && value>0.0;
}

}  // namespace

double BolometricPartitionFraction(int channel,double temperatureK) {
  const ChannelDefinition* definition=Definition(channel);
  if (definition==nullptr || !ValidTemperature(temperatureK))
    return std::numeric_limits<double>::quiet_NaN();
  const double bandRadiance=IntegratePlanck(
      definition->partitionMinimumMicrometre,
      definition->partitionMaximumMicrometre,temperatureK);
  // Integral of radiance over wavelength is sigma*T^4/pi.  The pi factors in
  // numerator and denominator cancel, leaving this dimensionless fraction.
  return bandRadiance/(kStefanBoltzmannWm2K4*
      std::pow(temperatureK,4)/kPi);
}

bool BolometricBrightnessTemperatureK(
    const std::array<double,InfraredChannelCount>& input,
    double* bolometricTemperatureK,
    std::array<double,InfraredChannelCount>* temperaturesUsedK,
    std::string* error) {
  if (bolometricTemperatureK==nullptr || temperaturesUsedK==nullptr ||
      error==nullptr) return false;
  error->clear();
  if (!ValidTemperature(input.back())) {
    *error="Diviner C9 brightness temperature is missing or non-positive";
    return false;
  }
  *temperaturesUsedK=input;
  // Apply the publication's sensitivity rule from long to short wavelength.
  // A channel below Tmin is replaced by the *next longer* channel that is
  // itself above its threshold.  Walking backward retains that already chosen
  // longer-wave value, including a chain of several insensitive channels.
  for (int index=InfraredChannelCount-2;index>=0;--index) {
    const double threshold=kChannels[static_cast<std::size_t>(index)].
        minimumDetectableTemperatureK;
    if (!ValidTemperature(input[static_cast<std::size_t>(index)]) ||
        input[static_cast<std::size_t>(index)]<threshold) {
      (*temperaturesUsedK)[static_cast<std::size_t>(index)]=
          (*temperaturesUsedK)[static_cast<std::size_t>(index+1)];
    }
  }
  double fourthPower=0.0;
  for (std::size_t index=0;index<kChannels.size();++index) {
    const double temperature=(*temperaturesUsedK)[index];
    const double fraction=BolometricPartitionFraction(
        kChannels[index].channel,temperature);
    if (!std::isfinite(fraction) || fraction<0.0) {
      *error="failed to integrate a Diviner bolometric wavelength partition";
      return false;
    }
    fourthPower+=std::pow(temperature,4)*fraction;
  }
  if (!(fourthPower>0.0) || !std::isfinite(fourthPower)) {
    *error="Diviner bolometric reconstruction produced invalid exitance";
    return false;
  }
  *bolometricTemperatureK=std::pow(fourthPower,0.25);
  return true;
}

double GreySurfaceExitanceWm2(double temperatureK,double emissivity) {
  if (!ValidTemperature(temperatureK) || !std::isfinite(emissivity) ||
      emissivity<0.0 || emissivity>1.0)
    return std::numeric_limits<double>::quiet_NaN();
  return emissivity*kStefanBoltzmannWm2K4*std::pow(temperatureK,4);
}

double GreyBolometricBrightnessTemperatureK(double temperatureK,
                                             double emissivity) {
  const double exitance=GreySurfaceExitanceWm2(temperatureK,emissivity);
  if (!std::isfinite(exitance)) return std::numeric_limits<double>::quiet_NaN();
  return std::pow(exitance/kStefanBoltzmannWm2K4,0.25);
}

double RectangularBandRadianceWm2Sr(int channel,double temperatureK,
                                    double emissivity) {
  const ChannelDefinition* definition=Definition(channel);
  if (definition==nullptr || !ValidTemperature(temperatureK) ||
      !std::isfinite(emissivity) || emissivity<0.0 || emissivity>1.0)
    return std::numeric_limits<double>::quiet_NaN();
  return emissivity*IntegratePlanck(definition->passbandMinimumMicrometre,
      definition->passbandMaximumMicrometre,temperatureK);
}

bool RectangularBandBrightnessTemperatureK(int channel,double radiance,
    double* temperatureK,std::string* error) {
  if (temperatureK==nullptr || error==nullptr) return false;
  error->clear();
  if (Definition(channel)==nullptr) {
    *error="Diviner channel must be in [3,9]";
    return false;
  }
  if (!(radiance>0.0) || !std::isfinite(radiance)) {
    *error="Diviner band radiance must be finite and positive";
    return false;
  }
  // The integrated Planck radiance is strictly increasing in temperature.
  // Bracket every physically meaningful lunar value, expanding the upper end
  // if needed, then bisect to a relative temperature interval below 1e-12.
  double lower=1.0,upper=500.0;
  while (RectangularBandRadianceWm2Sr(channel,upper,1.0)<radiance &&
         upper<4000.0) upper*=2.0;
  if (RectangularBandRadianceWm2Sr(channel,upper,1.0)<radiance) {
    *error="Diviner band radiance lies outside the supported temperature bracket";
    return false;
  }
  for (int iteration=0;iteration<80;++iteration) {
    const double middle=0.5*(lower+upper);
    if (RectangularBandRadianceWm2Sr(channel,middle,1.0)<radiance)
      lower=middle;
    else upper=middle;
  }
  *temperatureK=0.5*(lower+upper);
  return true;
}

}  // namespace DivinerObservation
}  // namespace Moon
