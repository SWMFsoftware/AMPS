#ifndef SRCMOON_MOON_DIVINER_OBSERVATION_H
#define SRCMOON_MOON_DIVINER_OBSERVATION_H

#include <array>
#include <string>

namespace Moon {
namespace DivinerObservation {

// Diviner infrared channels 3--9 are stored in ascending channel order in
// every array accepted by this interface.  Keeping the count in the public
// contract prevents a caller from silently supplying only the four broadband
// thermal channels to the Paige et al. bolometric algorithm.
constexpr int FirstInfraredChannel = 3;
constexpr int LastInfraredChannel = 9;
constexpr int InfraredChannelCount = 7;

// Published Paige et al. (2010) bolometric-brightness-temperature operator.
// The input is one brightness temperature per Diviner channel [K], ordered
// C3,...,C9.  C3--C8 values below their published sensitivity thresholds are
// replaced by the next longer-wavelength channel that is above its own
// threshold; C9 must be finite and positive.  The returned temperature is the
// blackbody temperature with the same spectrally integrated exitance [K].
// Invalid/missing input returns false and an explanatory message; no channel
// is interpolated or filled with a model value.
bool BolometricBrightnessTemperatureK(
    const std::array<double,InfraredChannelCount>& brightnessTemperatureK,
    double* bolometricBrightnessTemperatureK,
    std::array<double,InfraredChannelCount>* temperaturesUsedK,
    std::string* error);

// Fraction of a blackbody's total hemispheric radiative exitance contained in
// one of the wavelength partitions used by the Paige et al. T_bol operator.
// Wavelength limits are the published *partition* limits, not the narrower
// detector passbands.  The result is dimensionless and the partitions sum to
// the 0--1000 micrometre fraction represented by Diviner channels 3--9.
double BolometricPartitionFraction(int channel,double temperatureK);

// Production grey-surface observation operator.  The dynamic thermal model
// predicts kinetic surface temperature T and uses a bolometric infrared
// emissivity epsilon in its boundary condition.  Comparing T directly with a
// blackbody-equivalent Diviner T_bol would be inconsistent; emitted exitance
// is epsilon*sigma*T^4 and its equivalent blackbody temperature is
// epsilon^(1/4)*T.  Units are W m^-2 and K, respectively.
double GreySurfaceExitanceWm2(double kineticTemperatureK,
                              double bolometricEmissivity);
double GreyBolometricBrightnessTemperatureK(double kineticTemperatureK,
                                             double bolometricEmissivity);

// Diagnostic forward operator for C6--C9.  The calibration report explicitly
// permits rectangular half-maximum passbands for uncertainty calculations.
// We use that documented approximation only for channel-residual diagnostics,
// never as the primary I17 bolometric score.  Radiance is band-integrated
// spectral radiance [W m^-2 sr^-1]; the inverse returns brightness temperature
// [K].  Grey emissivity is explicit because the production thermal state does
// not currently carry a wavelength-dependent emissivity spectrum.
double RectangularBandRadianceWm2Sr(int channel,double kineticTemperatureK,
                                    double greyEmissivity);
bool RectangularBandBrightnessTemperatureK(int channel,
    double radianceWm2Sr,double* brightnessTemperatureK,std::string* error);

}  // namespace DivinerObservation
}  // namespace Moon

#endif  // SRCMOON_MOON_DIVINER_OBSERVATION_H
