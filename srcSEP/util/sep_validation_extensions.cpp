#include "sep_validation_extensions.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace SEP {
namespace ValidationExtensions {
namespace {

Transport::Status Error(Transport::StatusCode code,
                        const std::string& message) {
  return Transport::Status::Error(code, message);
}

bool Finite(double value) { return std::isfinite(value); }

double PoissonLogLikelihood(double observed, double expected) {
  if (expected == 0.0) return observed == 0.0
      ? 0.0 : -std::numeric_limits<double>::infinity();
  return observed * std::log(expected) - expected - std::lgamma(observed + 1.0);
}

// Acklam's rational approximation gives more than enough accuracy for sample-
// size planning.  It is kept private so scientific likelihood calculations do
// not accidentally treat this planning approximation as an exact normal CDF.
double InverseNormal(double p) {
  const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
      -2.759285104469687e+02, 1.383577518672690e+02,
      -3.066479806614716e+01, 2.506628277459239e+00};
  const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
      -1.556989798598866e+02, 6.680131188771972e+01,
      -1.328068155288572e+01};
  const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
      -2.400758277161838e+00, -2.549732539343734e+00,
       4.374664141464968e+00, 2.938163982698783e+00};
  const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
      2.445134137142996e+00, 3.754408661907416e+00};
  const double low = 0.02425;
  const double high = 1.0 - low;
  if (p < low) {
    const double q = std::sqrt(-2.0 * std::log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
  }
  if (p > high) return -InverseNormal(1.0 - p);
  const double q = p - 0.5;
  const double r = q * q;
  return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
         (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
}

double RelativeDifference(double a, double b) {
  return std::fabs(a - b) / std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

}  // namespace

LineageEstimate EstimateByLineage(const std::vector<LineageSample>& samples,
                                  std::size_t bins,
                                  EstimatorKind kind) {
  LineageEstimate result;
  result.kind = kind;
  if (bins == 0 || samples.empty()) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "lineage estimator requires samples and bins");
    return result;
  }
  typedef std::pair<std::uint64_t, std::uint64_t> RootWindow;
  std::map<RootWindow, std::vector<double> > clusters;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    const LineageSample& s = samples[i];
    if (s.rootLineage == 0 || s.samplingWindow == 0 || s.bin >= bins ||
        !Finite(s.statisticalWeight) || s.statisticalWeight < 0.0 ||
        !Finite(s.physicalValue) || !Finite(s.speedMPerS) ||
        s.speedMPerS < 0.0 || !Finite(s.jacobian) || s.jacobian < 0.0) {
      result.status = Error(Transport::StatusCode::InvalidArgument,
                            "lineage sample metadata or value is invalid");
      return result;
    }
    std::vector<double>& values = clusters[RootWindow(
        s.rootLineage, s.samplingWindow)];
    if (values.empty()) values.assign(bins, 0.0);
    double contribution = s.statisticalWeight * s.physicalValue * s.jacobian;
    if (kind == EstimatorKind::CrossingFlux) contribution *= s.speedMPerS;
    values[s.bin] += contribution;
  }
  result.independentRoots = clusters.size();
  result.estimate.assign(bins, 0.0);
  result.variance.assign(bins, 0.0);
  result.effectiveSampleSize.assign(bins, 0.0);
  result.covariance.assign(bins * bins, 0.0);
  for (std::map<RootWindow, std::vector<double> >::const_iterator i =
       clusters.begin(); i != clusters.end(); ++i)
    for (std::size_t b = 0; b < bins; ++b)
      result.estimate[b] += i->second[b];
  if (clusters.size() > 1) {
    const double n = static_cast<double>(clusters.size());
    for (std::map<RootWindow, std::vector<double> >::const_iterator i =
         clusters.begin(); i != clusters.end(); ++i) {
      for (std::size_t a = 0; a < bins; ++a) {
        const double da = i->second[a] - result.estimate[a] / n;
        for (std::size_t b = 0; b < bins; ++b) {
          const double db = i->second[b] - result.estimate[b] / n;
          result.covariance[a * bins + b] += da * db / (n - 1.0);
        }
      }
    }
  }
  for (std::size_t b = 0; b < bins; ++b) {
    result.variance[b] = result.covariance[b * bins + b];
    long double sum = 0.0L, sumSquare = 0.0L;
    for (std::map<RootWindow, std::vector<double> >::const_iterator i =
         clusters.begin(); i != clusters.end(); ++i) {
      sum += i->second[b];
      sumSquare += static_cast<long double>(i->second[b]) * i->second[b];
    }
    result.effectiveSampleSize[b] = sumSquare > 0.0L
        ? static_cast<double>(sum * sum / sumSquare) : 0.0;
  }
  result.status = Transport::Status::Ok();
  return result;
}

Transport::Status ValidateResponseMatrix(
    const InstrumentResponseMatrix& matrix) {
  const std::size_t expected = matrix.channels * matrix.species *
      matrix.directions * matrix.trueEnergyBins;
  if (matrix.instrumentId.empty() || matrix.calibrationVersion.empty() ||
      matrix.checksum.empty() || matrix.dataLevel.empty() ||
      !Finite(matrix.validFromS) || !Finite(matrix.validUntilS) ||
      matrix.validUntilS < matrix.validFromS || matrix.channels == 0 ||
      matrix.species == 0 || matrix.directions == 0 ||
      matrix.trueEnergyBins == 0 || matrix.probability.size() != expected)
    return Error(Transport::StatusCode::InvalidArgument,
                 "instrument response manifest or dimensions are invalid");
  for (std::size_t i = 0; i < matrix.probability.size(); ++i)
    if (!Finite(matrix.probability[i]) || matrix.probability[i] < 0.0 ||
        matrix.probability[i] > 1.0)
      return Error(Transport::StatusCode::InvalidCoefficient,
                   "instrument response probability is outside [0,1]");
  return Transport::Status::Ok();
}

ResponseFoldResult FoldInstrumentCounts(
    const InstrumentResponseMatrix& matrix,
    const std::vector<double>& trueCounts,
    const std::vector<CountObservation>& observations,
    double exposureS, double deadTimeS, DeadTimeModel deadTimeModel) {
  ResponseFoldResult result;
  result.status = ValidateResponseMatrix(matrix);
  const std::size_t inputSize = matrix.species * matrix.directions *
                                matrix.trueEnergyBins;
  if (!result.status.ok()) return result;
  if (trueCounts.size() != inputSize || observations.size() != matrix.channels ||
      !Finite(exposureS) || exposureS <= 0.0 || !Finite(deadTimeS) ||
      deadTimeS < 0.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "instrument fold inputs do not match response dimensions");
    return result;
  }
  result.incidentCounts.assign(matrix.channels, 0.0);
  result.recordedCounts.assign(matrix.channels, 0.0);
  result.poissonLogLikelihood.assign(matrix.channels, 0.0);
  result.censored.assign(matrix.channels, false);
  for (std::size_t c = 0; c < matrix.channels; ++c) {
    long double incident = 0.0L;
    for (std::size_t i = 0; i < inputSize; ++i) {
      if (!Finite(trueCounts[i]) || trueCounts[i] < 0.0) {
        result.status = Error(Transport::StatusCode::InvalidArgument,
                              "true instrument counts are invalid");
        return result;
      }
      incident += matrix.probability[c * inputSize + i] * trueCounts[i];
    }
    result.incidentCounts[c] = static_cast<double>(incident);
    const double rate = result.incidentCounts[c] / exposureS;
    if (deadTimeModel == DeadTimeModel::Nonparalyzable)
      result.recordedCounts[c] = exposureS * rate / (1.0 + rate * deadTimeS);
    else if (deadTimeModel == DeadTimeModel::Paralyzable)
      result.recordedCounts[c] = exposureS * rate * std::exp(-rate * deadTimeS);
    else result.recordedCounts[c] = result.incidentCounts[c];
    const CountObservation& observed = observations[c];
    if (!Finite(observed.counts) || observed.counts < 0.0 ||
        !Finite(observed.backgroundCounts) || observed.backgroundCounts < 0.0 ||
        (observed.saturated && (!Finite(observed.saturationThreshold) ||
                                observed.saturationThreshold <= 0.0))) {
      result.status = Error(Transport::StatusCode::InvalidArgument,
                            "count observation is invalid");
      return result;
    }
    const double expected = result.recordedCounts[c] + observed.backgroundCounts;
    if (observed.saturated) {
      // A saturated datum is right-censored.  The compact core records the
      // censoring and omits a false Gaussian/point-Poisson residual; mission
      // adapters may supply an exact survival likelihood for their electronics.
      result.censored[c] = true;
      result.poissonLogLikelihood[c] = 0.0;
    } else {
      result.poissonLogLikelihood[c] =
          PoissonLogLikelihood(observed.counts, expected);
    }
  }
  result.status = Transport::Status::Ok();
  return result;
}

Transport::Status ValidateNativeMatrixRow(
    const NativeMatrixTrace& trace, bool couplingExpected) {
  if (trace.observedLevel != Evidence::Level::NativeAmps &&
      trace.observedLevel != Evidence::Level::SwmfReplay &&
      trace.observedLevel != Evidence::Level::ObservationalValidation)
    return Error(Transport::StatusCode::UnsupportedConfiguration,
                 "source-only evidence cannot satisfy the native matrix gate");
  if (trace.configurationFingerprint.empty() ||
      trace.executableChecksum.empty() || trace.compiler.empty() ||
      trace.flags.empty() || trace.moverId.empty() ||
      trace.coefficientProviderId.empty() || trace.backgroundOwnerId.empty() ||
      trace.turbulenceOwnerId.empty() || !trace.moverEntered ||
      !trace.coefficientEntered || !trace.turbulenceEntered ||
      (couplingExpected && !trace.couplingEntered))
    return Error(Transport::StatusCode::InvalidParticleState,
                 "native matrix trace is incomplete");
  return Transport::Status::Ok();
}

Transport::Status CompareDecompositions(
    const DecompositionSignature& reference,
    const DecompositionSignature& candidate,
    double tolerance) {
  if (reference.ranks <= 0 || reference.threads <= 0 || candidate.ranks <= 0 ||
      candidate.threads <= 0 || !Finite(tolerance) || tolerance < 0.0 ||
      reference.configurationFingerprint.empty() ||
      reference.configurationFingerprint != candidate.configurationFingerprint)
    return Error(Transport::StatusCode::InvalidArgument,
                 "decomposition comparison metadata is invalid");
  if (reference.particleIdentityHash != candidate.particleIdentityHash ||
      reference.transactionOrderHash != candidate.transactionOrderHash ||
      reference.restartHash != candidate.restartHash ||
      reference.sourceEvents != candidate.sourceEvents ||
      reference.transactions != candidate.transactions)
    return Error(Transport::StatusCode::InvalidParticleState,
                 "discrete identity or transaction order changed with decomposition");
  if (reference.ledger.size() != candidate.ledger.size())
    return Error(Transport::StatusCode::InvalidArgument,
                 "decomposition ledgers have different schemas");
  for (std::size_t i = 0; i < reference.ledger.size(); ++i)
    if (!Finite(reference.ledger[i]) || !Finite(candidate.ledger[i]) ||
        RelativeDifference(reference.ledger[i], candidate.ledger[i]) > tolerance)
      return Error(Transport::StatusCode::OutOfDomain,
                   "floating ledger exceeds decomposition tolerance");
  return Transport::Status::Ok();
}

RefinementFit FitRefinementOrder(const ManufacturedCase& testCase) {
  RefinementFit result;
  const std::size_t n = testCase.resolution.size();
  if (testCase.id.empty() || testCase.activeTerms.empty() || n < 3 ||
      testCase.error.size() != n) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "manufactured refinement case is incomplete");
    return result;
  }
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, syy = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    if (!Finite(testCase.resolution[i]) || testCase.resolution[i] <= 0.0 ||
        !Finite(testCase.error[i]) || testCase.error[i] <= 0.0) {
      result.status = Error(Transport::StatusCode::InvalidArgument,
                            "refinement resolutions and errors must be positive");
      return result;
    }
    const double x = std::log(testCase.resolution[i]);
    const double y = std::log(testCase.error[i]);
    sx += x; sy += y; sxx += x*x; sxy += x*y; syy += y*y;
  }
  const double denominator = n*sxx - sx*sx;
  if (denominator == 0.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "refinement resolutions are degenerate");
    return result;
  }
  const double slope = (n*sxy - sx*sy) / denominator;
  // resolution is the physical mesh spacing h (not the number of cells), so
  // error=C*h^p has a positive log-log slope p.
  result.observedOrder = slope;
  result.intercept = (sy - slope*sx) / n;
  const double correlationDenominator =
      (n*sxx - sx*sx) * (n*syy - sy*sy);
  result.rSquared = correlationDenominator > 0.0
      ? (n*sxy - sx*sy) * (n*sxy - sx*sy) / correlationDenominator
      : 1.0;
  result.status = Transport::Status::Ok();
  return result;
}

PowerResult PlanNormalMeanPower(const PowerSpecification& specification) {
  PowerResult result;
  if (specification.observable.empty() ||
      !Finite(specification.standardizedEffect) ||
      specification.standardizedEffect <= 0.0 || !Finite(specification.alpha) ||
      specification.alpha <= 0.0 || specification.alpha >= 1.0 ||
      !Finite(specification.power) || specification.power <= 0.5 ||
      specification.power >= 1.0 || specification.comparisons == 0 ||
      specification.seedPanelVersion.empty()) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "statistical power specification is incomplete");
    return result;
  }
  result.adjustedAlpha = specification.alpha /
                         static_cast<double>(specification.comparisons);
  const double zAlpha = InverseNormal(1.0 - 0.5 * result.adjustedAlpha);
  const double zPower = InverseNormal(specification.power);
  result.requiredSamples = static_cast<std::size_t>(std::ceil(
      (zAlpha + zPower) * (zAlpha + zPower) /
      (specification.standardizedEffect * specification.standardizedEffect)));
  result.status = Transport::Status::Ok();
  return result;
}

Transport::ScalarResult KolmogorovSmirnovStatistic(
    const std::vector<double>& samples,
    const std::vector<double>& expected) {
  Transport::ScalarResult result;
  if (samples.empty() || samples.size() != expected.size()) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "KS statistic requires matched nonempty arrays");
    return result;
  }
  double statistic = 0.0;
  const double n = static_cast<double>(samples.size());
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (!Finite(samples[i]) || (i && samples[i] < samples[i - 1]) ||
        !Finite(expected[i]) || expected[i] < 0.0 || expected[i] > 1.0) {
      result.status = Error(Transport::StatusCode::InvalidArgument,
                            "KS samples or CDF values are invalid");
      return result;
    }
    const double lower = static_cast<double>(i) / n;
    const double upper = static_cast<double>(i + 1) / n;
    statistic = std::max(statistic, std::fabs(expected[i] - lower));
    statistic = std::max(statistic, std::fabs(upper - expected[i]));
  }
  result.value = statistic;
  result.status = Transport::Status::Ok();
  return result;
}

Transport::Status ValidateExternalEventManifest(
    const ExternalEventManifest& manifest,
    Evidence::Level claimedLevel) {
  if (manifest.eventId.empty() || manifest.backgroundChecksum.empty() ||
      manifest.geometryChecksum.empty() || manifest.shockChecksum.empty() ||
      manifest.preprocessing.empty() || manifest.exclusions.empty() ||
      manifest.configurationPrior.empty() || manifest.metrics.empty())
    return Error(Transport::StatusCode::InvalidArgument,
                 "external event manifest is incomplete");
  if (claimedLevel == Evidence::Level::ObservationalValidation &&
      (manifest.synthetic || manifest.role != EventRole::HeldOut ||
       manifest.calibrationChecksum.empty()))
    return Error(Transport::StatusCode::UnsupportedConfiguration,
                 "observational evidence requires real held-out calibrated data");
  if (claimedLevel == Evidence::Level::SwmfReplay && manifest.synthetic)
    return Error(Transport::StatusCode::UnsupportedConfiguration,
                 "synthetic background cannot satisfy the SWMF replay gate");
  return Transport::Status::Ok();
}

Transport::ScalarResult NormalizedMetricDistance(
    const EventMetrics& model, const EventMetrics& reference,
    const EventMetrics& sigma) {
  Transport::ScalarResult result;
  const double m[] = {model.onsetS, model.peak, model.fluence,
      model.anisotropy, model.spectralIndex, model.profileRmse};
  const double r[] = {reference.onsetS, reference.peak, reference.fluence,
      reference.anisotropy, reference.spectralIndex, reference.profileRmse};
  const double s[] = {sigma.onsetS, sigma.peak, sigma.fluence,
      sigma.anisotropy, sigma.spectralIndex, sigma.profileRmse};
  long double sum = 0.0L;
  for (std::size_t i = 0; i < 6; ++i) {
    if (!Finite(m[i]) || !Finite(r[i]) || !Finite(s[i]) || s[i] <= 0.0) {
      result.status = Error(Transport::StatusCode::InvalidArgument,
                            "event metric or uncertainty is invalid");
      return result;
    }
    const double z = (m[i] - r[i]) / s[i];
    sum += z * z;
  }
  result.value = std::sqrt(static_cast<double>(sum / 6.0L));
  result.status = Transport::Status::Ok();
  return result;
}

ScalingGate CompareScaling(const ScalingSample& reference,
                           const ScalingSample& candidate,
                           double minimumEfficiency,
                           double maximumWorkRatio,
                           double maximumMemoryRatio) {
  ScalingGate result;
  if (reference.workloadId.empty() ||
      reference.workloadId != candidate.workloadId ||
      reference.ranks <= 0 || reference.threadsPerRank <= 0 ||
      candidate.ranks <= 0 || candidate.threadsPerRank <= 0 ||
      reference.problemItems == 0 || candidate.problemItems == 0 ||
      !Finite(reference.wallSeconds) || reference.wallSeconds <= 0.0 ||
      !Finite(candidate.wallSeconds) || candidate.wallSeconds <= 0.0 ||
      !Finite(minimumEfficiency) || minimumEfficiency <= 0.0 ||
      minimumEfficiency > 1.0 || !Finite(maximumWorkRatio) ||
      maximumWorkRatio < 1.0 || !Finite(maximumMemoryRatio) ||
      maximumMemoryRatio < 1.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "scaling comparison metadata is invalid");
    return result;
  }
  const double refWorkers = reference.ranks * reference.threadsPerRank;
  const double candidateWorkers = candidate.ranks * candidate.threadsPerRank;
  const double sizeRatio = static_cast<double>(candidate.problemItems) /
                           reference.problemItems;
  const double idealWall = reference.wallSeconds * sizeRatio *
                           refWorkers / candidateWorkers;
  result.parallelEfficiency = idealWall / candidate.wallSeconds;
  result.operationsPerItem = static_cast<double>(candidate.operations) /
                             candidate.problemItems;
  result.bytesPerItem = static_cast<double>(candidate.communicationBytes) /
                        candidate.problemItems;
  result.memoryPerItem = static_cast<double>(candidate.peakResidentBytes) /
                         candidate.problemItems;
  result.queuePerItem = static_cast<double>(candidate.peakQueueRecords) /
                        candidate.problemItems;
  const double refWorkPerItem = static_cast<double>(reference.operations) /
                                reference.problemItems;
  const double refMemoryPerItem = static_cast<double>(reference.peakResidentBytes) /
                                  reference.problemItems;
  if (result.operationsPerItem > maximumWorkRatio * refWorkPerItem ||
      result.memoryPerItem > maximumMemoryRatio * refMemoryPerItem) {
    result.status = Error(Transport::StatusCode::OutOfDomain,
                          "algorithmic work or memory complexity regressed");
    return result;
  }
  if (!reference.environmentFingerprint.empty() &&
      reference.environmentFingerprint == candidate.environmentFingerprint &&
      result.parallelEfficiency < minimumEfficiency) {
    result.status = Error(Transport::StatusCode::OutOfDomain,
                          "parallel efficiency is below its compatible-environment gate");
    return result;
  }
  result.status = Transport::Status::Ok();
  return result;
}

}  // namespace ValidationExtensions
}  // namespace SEP
