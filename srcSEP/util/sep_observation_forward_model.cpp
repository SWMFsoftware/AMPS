#include "sep_observation_forward_model.h"

#include <algorithm>
#include <cmath>

namespace SEP {
namespace Observation {
namespace {

Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                  message);
}

}  // namespace

Transport::Status ValidateChannel(const InstrumentChannel& c,
                                  std::size_t bins) {
  if (c.id.empty() || !(c.minimumEnergyJ >= 0.0) ||
      !(c.maximumEnergyJ > c.minimumEnergyJ) ||
      !(c.geometricFactorM2Sr > 0.0) || !(c.cadenceS > 0.0) ||
      !(c.nonparalyzableDeadTimeS >= 0.0) ||
      !(c.saturationCounts > 0.0) || !std::isfinite(c.maximumEnergyJ) ||
      !std::isfinite(c.geometricFactorM2Sr) || !std::isfinite(c.cadenceS) ||
      c.response.size() != bins)
    return Error("instrument channel metadata or response length is invalid");
  for (std::size_t i = 0; i < c.response.size(); ++i)
    if (!(c.response[i] >= 0.0 && c.response[i] <= 1.0) ||
        !std::isfinite(c.response[i]))
      return Error("instrument response must be finite and in [0,1]");
  return Transport::Status::Ok();
}

ChannelPrediction ForwardModel(const PopulationSpectrum& spectrum,
                               const InstrumentChannel& channel,
                               const BackgroundEstimate& background) {
  ChannelPrediction result;
  result.channelId = channel.id;
  const std::size_t bins = spectrum.differentialIntensity.size();
  if (spectrum.species.empty() || spectrum.energyEdgesJ.size() != bins + 1 ||
      spectrum.variance.size() != bins || !(background.expectedCounts >= 0.0) ||
      !(background.variance >= 0.0) || !std::isfinite(background.expectedCounts) ||
      !std::isfinite(background.variance)) {
    result.status = Error("population spectrum or background is malformed");
    return result;
  }
  result.status = ValidateChannel(channel, bins);
  if (!result.status.ok()) return result;

  long double counts = 0.0L;
  long double modelVariance = 0.0L;
  for (std::size_t i = 0; i < bins; ++i) {
    const double lower = spectrum.energyEdgesJ[i];
    const double upper = spectrum.energyEdgesJ[i + 1];
    if (!(upper > lower) || !std::isfinite(lower) || !std::isfinite(upper) ||
        !(spectrum.differentialIntensity[i] >= 0.0) ||
        !(spectrum.variance[i] >= 0.0) ||
        !std::isfinite(spectrum.differentialIntensity[i]) ||
        !std::isfinite(spectrum.variance[i])) {
      result.status = Error("spectrum bins must be ordered, finite, and non-negative");
      return result;
    }
    const double overlap = std::max(0.0,
        std::min(upper, channel.maximumEnergyJ) -
        std::max(lower, channel.minimumEnergyJ));
    const double exposure = overlap * channel.geometricFactorM2Sr *
                            channel.cadenceS * channel.response[i];
    counts += spectrum.differentialIntensity[i] * exposure;
    modelVariance += spectrum.variance[i] * exposure * exposure;
  }
  result.incidentCounts = static_cast<double>(counts);
  // A nonparalyzable counter records N/(1+N*tau/T).  The model is applied
  // before saturation and is explicit even when tau=0, avoiding an implicit
  // instrument-specific correction in later comparison scripts.
  result.recordedCounts = result.incidentCounts /
      (1.0 + result.incidentCounts * channel.nonparalyzableDeadTimeS /
             channel.cadenceS);
  if (result.recordedCounts > channel.saturationCounts) {
    result.recordedCounts = channel.saturationCounts;
    result.saturated = true;
  }
  result.backgroundSubtractedCounts =
      result.recordedCounts - background.expectedCounts;
  // Poisson counting variance and independently estimated background/model
  // variances add.  Saturated channels remain flagged because this symmetric
  // uncertainty is not a license to treat the clipped value as an unbiased datum.
  result.standardUncertainty = std::sqrt(std::max(0.0,
      result.recordedCounts + background.variance +
      static_cast<double>(modelVariance)));
  result.status = Transport::Status::Ok();
  return result;
}

}  // namespace Observation
}  // namespace SEP
