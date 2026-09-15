#ifndef SEP_UTIL_SEP_VALIDATION_TOOLS_H
#define SEP_UTIL_SEP_VALIDATION_TOOLS_H

#include "sep_transport_common.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace SEP {
namespace ValidationTools {

struct SeedPanel {
  std::string specificationVersion;
  std::string observableId;
  std::vector<std::uint64_t> seeds;
};

// Seed panels are domain separated by observable and panel version.  Their
// deterministic generation makes ensemble correctness reproducible while
// remaining a distinct gate from byte-identical fixed-seed regression.
SeedPanel MakeSeedPanel(std::uint64_t campaignSeed,
                        const std::string& specificationVersion,
                        const std::string& observableId,
                        std::size_t count);

struct RunningStatistics {
  std::uint64_t count = 0;
  double mean = 0.0;
  double secondCentralMoment = 0.0;
  void Add(double value);
  double SampleVariance() const;
  double StandardError() const;
};

struct StatisticalGate {
  Transport::Status status;
  double estimate = 0.0;
  double standardError = 0.0;
  double zScore = 0.0;
  double confidenceLevel = 0.95;
  std::uint64_t sampleCount = 0;
};

StatisticalGate CompareMean(const RunningStatistics& statistics,
                            double analyticalMean,
                            double maximumAbsoluteZ,
                            double confidenceLevel);

struct GeneratedCase {
  std::uint64_t seed = 0;
  std::uint64_t index = 0;
  double value = 0.0;
};

// The bounded generator always includes IEEE and model boundaries before
// pseudo-random finite values.  A failure therefore logs a stable seed/index
// pair that can be promoted directly into a named regression test.
std::vector<GeneratedCase> GenerateBoundaryCases(std::uint64_t seed,
                                                 std::size_t randomCount);

struct PropertyFailure {
  bool failed = false;
  GeneratedCase counterexample;
  std::string property;
  std::string message;
};

typedef std::function<Transport::Status(double)> Property;
PropertyFailure CheckProperty(const std::string& name,
                              const std::vector<GeneratedCase>& cases,
                              const Property& property);

enum class FaultPoint {
  Allocation,
  FileWrite,
  Checksum,
  Snapshot,
  Coupler,
  Communication
};

class FaultInjector {
 public:
  FaultInjector(FaultPoint point, std::uint64_t failOnHit);
  bool ShouldFail(FaultPoint point);
  std::uint64_t hits() const { return hits_; }
 private:
  FaultPoint point_;
  std::uint64_t failOnHit_;
  std::uint64_t hits_;
};

struct WorkCounters {
  std::uint64_t moverSubsteps = 0;
  std::uint64_t rejectedSteps = 0;
  std::uint64_t hazardEvaluations = 0;
  std::uint64_t coefficientEvaluations = 0;
  std::uint64_t peakQueuedContributions = 0;
  std::uint64_t reductionBytes = 0;
  std::uint64_t spectralBinUpdates = 0;
  std::uint64_t outputBytes = 0;
  std::uint64_t checkpointBytes = 0;
  std::uint64_t peakResidentBytes = 0;
};

struct PerformanceSample {
  WorkCounters work;
  double wallSeconds = 0.0;
  std::string compiler;
  std::string hardware;
  int threads = 1;
  int ranks = 1;
};

struct PerformanceBaseline {
  WorkCounters expectedWork;
  double medianWallSeconds = 0.0;
  double maximumWallRatio = 1.25;
  std::string environmentId;
};

// Normalized work is deterministic and always enforced.  Wall time is enforced
// only when the caller supplies the same nonempty environment identity; this
// prevents a laptop or CI host change from masquerading as an algorithmic
// regression while retaining an actionable hardware-specific gate.
Transport::Status ComparePerformance(const PerformanceSample& sample,
                                     const PerformanceBaseline& baseline,
                                     const std::string& environmentId);

}  // namespace ValidationTools
}  // namespace SEP

#endif  // SEP_UTIL_SEP_VALIDATION_TOOLS_H
