#include "sep_validation_tools.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace ValidationTools {
namespace {

std::uint64_t HashString(const std::string& value, std::uint64_t seed) {
  std::uint64_t hash = seed;
  for (std::size_t i = 0; i < value.size(); ++i) {
    hash ^= static_cast<unsigned char>(value[i]);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

std::uint64_t Mix(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

bool SameWork(const WorkCounters& a, const WorkCounters& b) {
  return a.moverSubsteps == b.moverSubsteps &&
      a.rejectedSteps == b.rejectedSteps &&
      a.hazardEvaluations == b.hazardEvaluations &&
      a.coefficientEvaluations == b.coefficientEvaluations &&
      a.peakQueuedContributions == b.peakQueuedContributions &&
      a.reductionBytes == b.reductionBytes &&
      a.spectralBinUpdates == b.spectralBinUpdates &&
      a.outputBytes == b.outputBytes &&
      a.checkpointBytes == b.checkpointBytes;
}

}  // namespace

SeedPanel MakeSeedPanel(std::uint64_t campaignSeed,
                        const std::string& version,
                        const std::string& observable,
                        std::size_t count) {
  SeedPanel panel;
  panel.specificationVersion = version;
  panel.observableId = observable;
  std::uint64_t state = HashString(version,
      HashString(observable, campaignSeed ^ UINT64_C(14695981039346656037)));
  for (std::size_t i = 0; i < count; ++i) {
    state = Mix(state + static_cast<std::uint64_t>(i));
    panel.seeds.push_back(state);
  }
  return panel;
}

void RunningStatistics::Add(double value) {
  ++count;
  const double delta = value - mean;
  mean += delta / static_cast<double>(count);
  secondCentralMoment += delta * (value - mean);
}

double RunningStatistics::SampleVariance() const {
  return count > 1 ? secondCentralMoment / static_cast<double>(count - 1) : 0.0;
}

double RunningStatistics::StandardError() const {
  return count > 1 ? std::sqrt(SampleVariance() / static_cast<double>(count)) : 0.0;
}

StatisticalGate CompareMean(const RunningStatistics& statistics,
                            double analyticalMean, double maximumAbsoluteZ,
                            double confidenceLevel) {
  StatisticalGate result;
  result.estimate = statistics.mean;
  result.standardError = statistics.StandardError();
  result.confidenceLevel = confidenceLevel;
  result.sampleCount = statistics.count;
  if (statistics.count < 2 || !std::isfinite(analyticalMean) ||
      !(maximumAbsoluteZ > 0.0) || !(confidenceLevel > 0.0) ||
      !(confidenceLevel < 1.0) || !(result.standardError > 0.0)) {
    result.status = Transport::Status::Error(
        Transport::StatusCode::InvalidArgument,
        "statistical mean gate requires variance, sample count, and declared confidence");
    return result;
  }
  result.zScore = (statistics.mean - analyticalMean) / result.standardError;
  result.status = std::fabs(result.zScore) <= maximumAbsoluteZ
      ? Transport::Status::Ok()
      : Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                 "ensemble mean exceeds its declared z threshold");
  return result;
}

std::vector<GeneratedCase> GenerateBoundaryCases(std::uint64_t seed,
                                                 std::size_t randomCount) {
  std::vector<GeneratedCase> result;
  const double values[] = {0.0, -0.0, 1.0, -1.0,
      std::nextafter(1.0, 0.0), std::nextafter(-1.0, 0.0),
      std::numeric_limits<double>::min(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::quiet_NaN()};
  for (std::size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
    GeneratedCase generated;
    generated.seed = seed;
    generated.index = static_cast<std::uint64_t>(i);
    generated.value = values[i];
    result.push_back(generated);
  }
  std::uint64_t state = seed;
  for (std::size_t i = 0; i < randomCount; ++i) {
    state = Mix(state + static_cast<std::uint64_t>(i));
    const double unit = static_cast<double>(state >> 11) *
                        (1.0 / 9007199254740992.0);
    GeneratedCase generated;
    generated.seed = seed;
    generated.index = static_cast<std::uint64_t>(result.size());
    generated.value = 4.0 * unit - 2.0;
    result.push_back(generated);
  }
  return result;
}

PropertyFailure CheckProperty(const std::string& name,
                              const std::vector<GeneratedCase>& cases,
                              const Property& property) {
  PropertyFailure result;
  result.property = name;
  if (name.empty() || !property) {
    result.failed = true;
    result.message = "property name or callback is missing";
    return result;
  }
  for (std::size_t i = 0; i < cases.size(); ++i) {
    const Transport::Status status = property(cases[i].value);
    if (!status.ok()) {
      result.failed = true;
      result.counterexample = cases[i];
      result.message = status.message;
      return result;
    }
  }
  return result;
}

FaultInjector::FaultInjector(FaultPoint point, std::uint64_t failOnHit)
    : point_(point), failOnHit_(failOnHit), hits_(0) {}

bool FaultInjector::ShouldFail(FaultPoint point) {
  if (point != point_) return false;
  ++hits_;
  return failOnHit_ != 0 && hits_ == failOnHit_;
}

Transport::Status ComparePerformance(const PerformanceSample& sample,
                                     const PerformanceBaseline& baseline,
                                     const std::string& environmentId) {
  if (!SameWork(sample.work, baseline.expectedWork))
    return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                    "normalized algorithmic work changed from baseline");
  if (!(sample.wallSeconds >= 0.0) || !std::isfinite(sample.wallSeconds) ||
      !(baseline.medianWallSeconds >= 0.0) ||
      !std::isfinite(baseline.medianWallSeconds) ||
      !(baseline.maximumWallRatio >= 1.0))
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "performance sample or baseline is invalid");
  if (!environmentId.empty() && environmentId == baseline.environmentId &&
      baseline.medianWallSeconds > 0.0 &&
      sample.wallSeconds > baseline.maximumWallRatio * baseline.medianWallSeconds)
    return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                    "wall time exceeds the environment-specific budget");
  return Transport::Status::Ok();
}

}  // namespace ValidationTools
}  // namespace SEP
