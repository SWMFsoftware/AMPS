#ifndef SHIELDSIM_PUBLISHED_SOURCE_ADAPTERS_HH
#define SHIELDSIM_PUBLISHED_SOURCE_ADAPTERS_HH

/*
 * Lossless Stage-2 adapters for pinned published source products.
 *
 * These adapters stop at the canonical SpectrumDataset boundary.  In
 * particular, neither an ISS AMS-02 measurement nor a directional GOES SGPS
 * record is silently converted into an isotropic deep-space boundary source.
 * Transport/folding is a separate Stage-3 gate.
 */

#include "SpectrumDataset.hh"

#include <cstdint>
#include <string>
#include <vector>

namespace shield_sim {

class AMS02HapiAdapter {
public:
  // interpretationPath pins the primary-publication rigidity edges and the
  // exact relativistic transformation used by the CDAWeb HAPI product.  It is
  // mandatory because the HAPI JSON calls its coordinate array "centers"
  // even though its Bin*Low name and the primary tables identify lower edges.
  static SpectrumDataset Load(const std::string& infoPath,
                              const std::string& dataPath,
                              const std::string& interpretationPath);
};

struct GoesSgpsRecord {
  std::int64_t startUnixSeconds = 0;
  std::int64_t midpointUnixSeconds = 0;
  std::int64_t endUnixSeconds = 0;
  std::string startUtc;
  std::string midpointUtc;
  std::string endUtc;
  int l1bRecordsInAverage = 0;
  int sensorUnit = -1;
  std::string spacecraftAxis;
  std::string nominalLook;
  int yawFlipFlag = 0;
  std::string channel;
  double energyLowMeV = 0.0;
  double energyHighMeV = 0.0;
  bool hasFlux = false;
  double fluxPerCm2SrSMeV = 0.0;
  bool hasUncertainty = false;
  double uncertaintyPerCm2SrSMeV = 0.0;
  int validOneSecondSamples = 0;
  int deadTimeFlaggedSamples = 0;
  int outOfBandFlaggedSamples = 0;
  int dynamicErrorFlaggedSamples = 0;

  bool IsUsable() const;
};

struct GoesEventWindow {
  std::string id;
  std::string role;
  std::string startUtc;
  std::string endUtc;
  std::int64_t startUnixSeconds = 0;
  std::int64_t endUnixSeconds = 0;
  std::string definition;
  std::string sourceUrl;
};

struct GoesFluenceChannel {
  std::string windowId;
  std::string sensorId;
  std::string spacecraftAxis;
  std::string nominalLook;
  std::string channel;
  double energyLowMeV = 0.0;
  double energyHighMeV = 0.0;
  double partialFluencePerCm2SrMeV = 0.0;
  double expectedSeconds = 0.0;
  double validSeconds = 0.0;
  double rejectedSeconds = 0.0;
  double missingSeconds = 0.0;
  double coverageFraction = 0.0;
  bool complete = false;
  std::string uncertaintyStatus;
};

struct GoesFluenceResult {
  std::vector<GoesEventWindow> windows;
  std::vector<GoesFluenceChannel> channels;
  Json canonicalJson;
  std::string contentHash;
};

class GOESRSGPSAdapter {
public:
  static GOESRSGPSAdapter Load(const std::string& csvPath,
                               const std::string& eventDefinitionPath);

  const std::vector<GoesSgpsRecord>& Records() const { return records_; }
  const std::vector<GoesEventWindow>& Windows() const { return windows_; }
  const Json& CanonicalMetadata() const { return canonicalMetadata_; }
  const std::string& ContentHash() const { return contentHash_; }

  // Integrates each sensor/channel independently. Values are time integrals of
  // published five-minute averages. Total-flux uncertainties are not combined
  // because their temporal covariance is unpublished; the result states that
  // limitation instead of assuming independent samples.
  GoesFluenceResult Integrate() const;

  // Public pure integration kernel used by analytic tests and by Integrate().
  // Exposing it also allows future source readers with the identical temporal
  // contract to reuse the audited boundary/gap logic.
  static std::vector<GoesFluenceChannel> IntegrateRecords(
      const std::vector<GoesSgpsRecord>& records,
      const std::vector<GoesEventWindow>& windows);

private:
  std::vector<GoesSgpsRecord> records_;
  std::vector<GoesEventWindow> windows_;
  Json canonicalMetadata_;
  std::string contentHash_;
};

} // namespace shield_sim

#endif
