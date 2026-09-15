#ifndef SEP_UTIL_SEP_OBSERVATION_FORWARD_MODEL_H
#define SEP_UTIL_SEP_OBSERVATION_FORWARD_MODEL_H

#include "sep_transport_common.h"

#include <string>
#include <vector>

namespace SEP {
namespace Observation {

struct InstrumentChannel {
  std::string id;
  double minimumEnergyJ = 0.0;
  double maximumEnergyJ = 0.0;
  double geometricFactorM2Sr = 0.0;
  double cadenceS = 0.0;
  double nonparalyzableDeadTimeS = 0.0;
  double saturationCounts = 0.0;
  // Response is sampled on the simulation energy-bin grid and is a
  // dimensionless probability in [0,1] for a particle in that bin.
  std::vector<double> response;
};

struct PopulationSpectrum {
  std::string species;
  std::vector<double> energyEdgesJ;
  // Omnidirectional differential intensity [m^-2 s^-1 sr^-1 J^-1].
  std::vector<double> differentialIntensity;
  std::vector<double> variance;
};

struct BackgroundEstimate {
  double expectedCounts = 0.0;
  double variance = 0.0;
};

struct ChannelPrediction {
  Transport::Status status;
  std::string channelId;
  double incidentCounts = 0.0;
  double recordedCounts = 0.0;
  double backgroundSubtractedCounts = 0.0;
  double standardUncertainty = 0.0;
  bool saturated = false;
};

Transport::Status ValidateChannel(const InstrumentChannel& channel,
                                  std::size_t spectrumBins);
ChannelPrediction ForwardModel(const PopulationSpectrum& spectrum,
                               const InstrumentChannel& channel,
                               const BackgroundEstimate& background);

}  // namespace Observation
}  // namespace SEP

#endif  // SEP_UTIL_SEP_OBSERVATION_FORWARD_MODEL_H
