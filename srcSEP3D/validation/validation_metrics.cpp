#include "validation_metrics.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <utility>

namespace SEP3D {
namespace Validation {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

bool ValidSeries(const std::vector<SeriesPoint>& series) {
  if (series.size() < 2) return false;
  for (std::size_t i = 0; i < series.size(); ++i) {
    if (!std::isfinite(series[i].coordinate) ||
        !std::isfinite(series[i].value) || series[i].value <= 0.0 ||
        (i != 0 && !(series[i].coordinate > series[i - 1].coordinate)))
      return false;
  }
  return true;
}

bool ValidPolicy(const ComparisonPolicy& policy) {
  return std::isfinite(policy.minimumCoverage) &&
      policy.minimumCoverage >= 0.0 && policy.minimumCoverage <= 1.0 &&
      std::isfinite(policy.onsetFractionOfPeak) &&
      policy.onsetFractionOfPeak > 0.0 &&
      policy.onsetFractionOfPeak <= 1.0 &&
      !std::isnan(policy.maximumLog10Rmse) &&
      policy.maximumLog10Rmse >= 0.0 &&
      !std::isnan(policy.maximumMedianAbsoluteLog10Error) &&
      policy.maximumMedianAbsoluteLog10Error >= 0.0 &&
      std::isfinite(policy.minimumLog10Correlation) &&
      policy.minimumLog10Correlation >= -1.0 &&
      policy.minimumLog10Correlation <= 1.0 &&
      !std::isnan(policy.maximumOnsetCoordinateError) &&
      policy.maximumOnsetCoordinateError >= 0.0 &&
      !std::isnan(policy.maximumPeakCoordinateError) &&
      policy.maximumPeakCoordinateError >= 0.0 &&
      !std::isnan(policy.maximumAbsoluteLog10PeakRatio) &&
      policy.maximumAbsoluteLog10PeakRatio >= 0.0 &&
      !std::isnan(policy.maximumAbsoluteLog10FluenceRatio) &&
      policy.maximumAbsoluteLog10FluenceRatio >= 0.0;
}

// Log-linear interpolation is appropriate for positive SEP intensities that
// commonly span many decades.  It also guarantees a positive interpolated
// value without a numerical floor that would contaminate log-space errors.
double InterpolateLogValue(const std::vector<SeriesPoint>& series, double x) {
  if (x < series.front().coordinate || x > series.back().coordinate)
    return std::numeric_limits<double>::quiet_NaN();
  auto high = std::lower_bound(
      series.begin(), series.end(), x,
      [](const SeriesPoint& point, double coordinate) {
        return point.coordinate < coordinate;
      });
  if (high == series.begin()) return high->value;
  if (high == series.end()) return series.back().value;
  if (high->coordinate == x) return high->value;
  const SeriesPoint& right = *high;
  const SeriesPoint& left = *(high - 1);
  const double fraction =
      (x - left.coordinate) / (right.coordinate - left.coordinate);
  return std::exp(std::log(left.value) +
                  fraction * (std::log(right.value) - std::log(left.value)));
}

double Median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  const std::size_t middle = values.size() / 2;
  if (values.size() % 2 != 0) return values[middle];
  return 0.5 * (values[middle - 1] + values[middle]);
}

double Correlation(const std::vector<double>& left,
                   const std::vector<double>& right) {
  const double leftMean =
      std::accumulate(left.begin(), left.end(), 0.0) / left.size();
  const double rightMean =
      std::accumulate(right.begin(), right.end(), 0.0) / right.size();
  double covariance = 0.0, leftVariance = 0.0, rightVariance = 0.0;
  for (std::size_t i = 0; i < left.size(); ++i) {
    const double dl = left[i] - leftMean;
    const double dr = right[i] - rightMean;
    covariance += dl * dr;
    leftVariance += dl * dl;
    rightVariance += dr * dr;
  }
  if (leftVariance == 0.0 || rightVariance == 0.0) {
    return left == right ? 1.0 : 0.0;
  }
  return covariance / std::sqrt(leftVariance * rightVariance);
}

std::size_t PeakIndex(const std::vector<double>& values) {
  return static_cast<std::size_t>(
      std::max_element(values.begin(), values.end()) - values.begin());
}

double OnsetCoordinate(const std::vector<double>& coordinates,
                       const std::vector<double>& values, double fraction) {
  const double threshold =
      *std::max_element(values.begin(), values.end()) * fraction;
  for (std::size_t i = 0; i < values.size(); ++i)
    if (values[i] >= threshold) return coordinates[i];
  return coordinates.back();  // unreachable for a positive finite series
}

double Trapezoid(const std::vector<double>& coordinates,
                 const std::vector<double>& values) {
  double integral = 0.0;
  for (std::size_t i = 1; i < values.size(); ++i) {
    integral += 0.5 * (values[i - 1] + values[i]) *
        (coordinates[i] - coordinates[i - 1]);
  }
  return integral;
}

bool AddOverflow(std::uint64_t left, std::uint64_t right) {
  return right > std::numeric_limits<std::uint64_t>::max() - left;
}

bool AddSizeOverflow(std::size_t left, std::size_t right) {
  return right > std::numeric_limits<std::size_t>::max() - left;
}

bool AddField(std::uint64_t value, std::uint64_t* destination) {
  if (AddOverflow(*destination, value)) return false;
  *destination += value;
  return true;
}

}  // namespace

ComparisonMetrics ComparePositiveSeries(
    const std::vector<SeriesPoint>& model,
    const std::vector<SeriesPoint>& reference,
    const ComparisonPolicy& policy) {
  ComparisonMetrics result;
  result.referencePoints = reference.size();
  if (!ValidSeries(model) || !ValidSeries(reference)) {
    result.status = Invalid(
        "scientific curves require at least two strictly ordered, positive, finite points");
    return result;
  }
  if (!ValidPolicy(policy)) {
    result.status = Invalid("scientific comparison policy is invalid");
    return result;
  }

  std::vector<double> coordinates, modelValues, referenceValues;
  for (const SeriesPoint& point : reference) {
    const double interpolated = InterpolateLogValue(model, point.coordinate);
    if (!std::isfinite(interpolated)) continue;
    coordinates.push_back(point.coordinate);
    modelValues.push_back(interpolated);
    referenceValues.push_back(point.value);
  }
  result.matchedPoints = modelValues.size();
  result.coverage = static_cast<double>(result.matchedPoints) /
      static_cast<double>(result.referencePoints);
  if (modelValues.size() < 2) {
    result.status = Invalid(
        "model/reference overlap contains fewer than two scientific points");
    return result;
  }

  if (policy.normalization == Normalization::OneGlobalAmplitude) {
    double meanLogRatio = 0.0;
    for (std::size_t i = 0; i < modelValues.size(); ++i)
      meanLogRatio += std::log(referenceValues[i] / modelValues[i]);
    result.amplitudeScale = std::exp(meanLogRatio / modelValues.size());
  } else if (policy.normalization == Normalization::UnitPeak) {
    result.amplitudeScale =
        *std::max_element(referenceValues.begin(), referenceValues.end()) /
        *std::max_element(modelValues.begin(), modelValues.end());
  }
  for (double& value : modelValues) value *= result.amplitudeScale;

  std::vector<double> modelLogs, referenceLogs, absoluteErrors;
  double squaredError = 0.0;
  for (std::size_t i = 0; i < modelValues.size(); ++i) {
    const double modelLog = std::log10(modelValues[i]);
    const double referenceLog = std::log10(referenceValues[i]);
    const double error = modelLog - referenceLog;
    modelLogs.push_back(modelLog);
    referenceLogs.push_back(referenceLog);
    absoluteErrors.push_back(std::fabs(error));
    squaredError += error * error;
  }
  result.log10Rmse = std::sqrt(squaredError / modelValues.size());
  result.medianAbsoluteLog10Error = Median(absoluteErrors);
  result.log10Correlation = Correlation(modelLogs, referenceLogs);

  const std::size_t modelPeak = PeakIndex(modelValues);
  const std::size_t referencePeak = PeakIndex(referenceValues);
  result.onsetCoordinateError = std::fabs(
      OnsetCoordinate(coordinates, modelValues, policy.onsetFractionOfPeak) -
      OnsetCoordinate(coordinates, referenceValues, policy.onsetFractionOfPeak));
  result.peakCoordinateError = std::fabs(
      coordinates[modelPeak] - coordinates[referencePeak]);
  result.log10PeakRatio = std::log10(
      modelValues[modelPeak] / referenceValues[referencePeak]);
  const double modelFluence = Trapezoid(coordinates, modelValues);
  const double referenceFluence = Trapezoid(coordinates, referenceValues);
  if (!(modelFluence > 0.0) || !(referenceFluence > 0.0)) {
    result.status = Invalid("covered scientific curve has zero fluence");
    return result;
  }
  result.log10FluenceRatio = std::log10(modelFluence / referenceFluence);

  auto requireMaximum = [&](double value, double maximum,
                            const char* name) {
    if (value > maximum) {
      std::ostringstream message;
      message << name << '=' << value << " exceeds " << maximum;
      result.violations.push_back(message.str());
    }
  };
  if (result.coverage < policy.minimumCoverage) {
    std::ostringstream message;
    message << "coverage=" << result.coverage << " is below "
            << policy.minimumCoverage;
    result.violations.push_back(message.str());
  }
  requireMaximum(result.log10Rmse, policy.maximumLog10Rmse, "log10_rmse");
  requireMaximum(result.medianAbsoluteLog10Error,
                 policy.maximumMedianAbsoluteLog10Error,
                 "median_absolute_log10_error");
  if (result.log10Correlation < policy.minimumLog10Correlation) {
    std::ostringstream message;
    message << "log10_correlation=" << result.log10Correlation
            << " is below " << policy.minimumLog10Correlation;
    result.violations.push_back(message.str());
  }
  requireMaximum(result.onsetCoordinateError,
                 policy.maximumOnsetCoordinateError,
                 "onset_coordinate_error");
  requireMaximum(result.peakCoordinateError,
                 policy.maximumPeakCoordinateError,
                 "peak_coordinate_error");
  requireMaximum(std::fabs(result.log10PeakRatio),
                 policy.maximumAbsoluteLog10PeakRatio,
                 "absolute_log10_peak_ratio");
  requireMaximum(std::fabs(result.log10FluenceRatio),
                 policy.maximumAbsoluteLog10FluenceRatio,
                 "absolute_log10_fluence_ratio");
  result.accepted = result.violations.empty();
  result.status = Core::Status::OK();
  return result;
}

IntegrationAudit AuditIntegration(
    const std::vector<RankPartition>& partitions,
    const IntegrationBudget& budget) {
  IntegrationAudit result;
  if (partitions.empty() ||
      !std::isfinite(budget.maximumParticleMaxToMean) ||
      budget.maximumParticleMaxToMean < 1.0 ||
      std::isnan(budget.maximumWallSeconds) || budget.maximumWallSeconds < 0.0) {
    result.status = Invalid("integration audit requires partitions and valid budgets");
    return result;
  }

  std::set<unsigned> ranks;
  std::map<Adapters::LedgerKey, Adapters::LedgerRow> ledgers;
  std::size_t maximumParticles = 0;
  for (const RankPartition& partition : partitions) {
    if (!ranks.insert(partition.rank).second ||
        !std::isfinite(partition.elapsedSeconds) ||
        partition.elapsedSeconds < 0.0) {
      result.status = Invalid("rank identity or resource measurement is invalid");
      return result;
    }
    maximumParticles = std::max(maximumParticles, partition.observations.size());
    result.maximumRankWallSeconds =
        std::max(result.maximumRankWallSeconds, partition.elapsedSeconds);
    if (AddSizeOverflow(result.totalPeakResidentBytes,
                        partition.peakResidentBytes)) {
      result.status = Invalid("rank memory sum overflowed");
      return result;
    }
    result.totalPeakResidentBytes += partition.peakResidentBytes;
    result.globalObservations.insert(result.globalObservations.end(),
                                     partition.observations.begin(),
                                     partition.observations.end());

    for (const Adapters::LedgerRow& row : partition.ledgerRows) {
      if (!row.closed || row.key.species < 0 ||
          row.advanced != row.activeEnd) {
        result.status = Invalid(
            "integration audit requires closed, internally consistent rank ledgers");
        return result;
      }
      Adapters::LedgerRow& total = ledgers[row.key];
      total.key = row.key;
      total.closed = true;
      if (!AddField(row.activeStart, &total.activeStart) ||
          !AddField(row.injected, &total.injected) ||
          !AddField(row.advanced, &total.advanced) ||
          !AddField(row.escaped, &total.escaped) ||
          !AddField(row.absorbed, &total.absorbed) ||
          !AddField(row.failed, &total.failed) ||
          !AddField(row.shockCrossings, &total.shockCrossings) ||
          !AddField(row.activeEnd, &total.activeEnd)) {
        result.status = Invalid("global particle ledger sum overflowed");
        return result;
      }
    }
  }

  std::sort(result.globalObservations.begin(), result.globalObservations.end(),
            [](const Output::ParticleObservation& left,
               const Output::ParticleObservation& right) {
              return left.stableId < right.stableId;
            });
  for (std::size_t i = 0; i < result.globalObservations.size(); ++i) {
    if (result.globalObservations[i].stableId == 0 ||
        (i != 0 && result.globalObservations[i].stableId ==
                       result.globalObservations[i - 1].stableId)) {
      result.status = Invalid(
          "global observation merge found zero or duplicate stable identity");
      return result;
    }
  }

  for (const auto& item : ledgers) {
    const Adapters::LedgerRow& row = item.second;
    if (AddOverflow(row.activeStart, row.injected) ||
        AddOverflow(row.activeEnd, row.escaped) ||
        AddOverflow(row.activeEnd + row.escaped, row.absorbed) ||
        AddOverflow(row.activeEnd + row.escaped + row.absorbed, row.failed) ||
        row.activeStart + row.injected !=
            row.activeEnd + row.escaped + row.absorbed + row.failed) {
      result.status = Invalid("globally reduced particle ledger does not close exactly");
      return result;
    }
    result.globalLedgerRows.push_back(row);
  }

  const double meanParticles = static_cast<double>(
      result.globalObservations.size()) / static_cast<double>(partitions.size());
  result.particleMaxToMean = meanParticles == 0.0 ? 0.0 :
      static_cast<double>(maximumParticles) / meanParticles;
  if (result.particleMaxToMean > budget.maximumParticleMaxToMean) {
    std::ostringstream message;
    message << "particle max/mean=" << result.particleMaxToMean
            << " exceeds " << budget.maximumParticleMaxToMean;
    result.budgetViolations.push_back(message.str());
  }
  if (result.maximumRankWallSeconds > budget.maximumWallSeconds) {
    std::ostringstream message;
    message << "maximum rank wall seconds=" << result.maximumRankWallSeconds
            << " exceeds " << budget.maximumWallSeconds;
    result.budgetViolations.push_back(message.str());
  }
  if (result.totalPeakResidentBytes > budget.maximumTotalResidentBytes) {
    std::ostringstream message;
    message << "summed peak resident bytes=" << result.totalPeakResidentBytes
            << " exceeds " << budget.maximumTotalResidentBytes;
    result.budgetViolations.push_back(message.str());
  }
  result.withinBudget = result.budgetViolations.empty();
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Validation
}  // namespace SEP3D
