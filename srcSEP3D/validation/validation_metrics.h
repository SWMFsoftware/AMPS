// ============================================================================
// Phase-V integration and scientific-validation contracts.
//
// This layer is deliberately independent of AMPS and MPI.  A coupled host
// exports immutable rank-local observations and scientific curves; these
// routines perform the canonical merge and comparison.  Keeping the audit
// here gives a single numerical definition to standalone tests, MPI campaign
// post-processing, and observational validation reports.
//
// Two distinctions are fundamental:
//   * `status` describes whether evidence is structurally usable.
//   * `accepted`/`withinBudget` describe whether valid evidence satisfies the
//     declared scientific or integration thresholds.
// A malformed or incomplete record must never be converted into a failed
// physics comparison with invented zero-valued metrics.
// ============================================================================

#ifndef SEP3D_VALIDATION_VALIDATION_METRICS_H
#define SEP3D_VALIDATION_VALIDATION_METRICS_H

#include "../output/sampling.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace SEP3D {
namespace Validation {

// A unit-bearing meaning is supplied by the campaign manifest.  The numerical
// layer needs only a strictly increasing coordinate and a positive observable;
// examples are time/intensity, radius/intensity, or energy/differential flux.
struct SeriesPoint {
  double coordinate = 0.0;
  double value = 0.0;
};

enum class Normalization {
  Absolute,
  // One global multiplicative factor is fitted in log space.  This is useful
  // only where the source area/normalization is explicitly unavailable; it
  // cannot change timing, spectral shape, or individual channels.
  OneGlobalAmplitude,
  // Each complete curve is divided by its own peak.  This is a shape-only
  // diagnostic and must not be described as absolute-intensity validation.
  UnitPeak
};

struct ComparisonPolicy {
  Normalization normalization = Normalization::Absolute;
  double minimumCoverage = 1.0;
  double onsetFractionOfPeak = 0.01;
  double maximumLog10Rmse = std::numeric_limits<double>::infinity();
  double maximumMedianAbsoluteLog10Error =
      std::numeric_limits<double>::infinity();
  double minimumLog10Correlation = -1.0;
  double maximumOnsetCoordinateError =
      std::numeric_limits<double>::infinity();
  double maximumPeakCoordinateError =
      std::numeric_limits<double>::infinity();
  double maximumAbsoluteLog10PeakRatio =
      std::numeric_limits<double>::infinity();
  double maximumAbsoluteLog10FluenceRatio =
      std::numeric_limits<double>::infinity();
};

struct ComparisonMetrics {
  Core::Status status;
  bool accepted = false;
  std::size_t referencePoints = 0;
  std::size_t matchedPoints = 0;
  double coverage = 0.0;
  double amplitudeScale = 1.0;
  double log10Rmse = 0.0;
  double medianAbsoluteLog10Error = 0.0;
  double log10Correlation = 0.0;
  double onsetCoordinateError = 0.0;
  double peakCoordinateError = 0.0;
  double log10PeakRatio = 0.0;
  double log10FluenceRatio = 0.0;
  std::vector<std::string> violations;
};

// Interpolate the model in log(value) at every covered reference coordinate,
// apply exactly the requested normalization, and evaluate all metrics on that
// common grid.  Inputs are never extrapolated.  This prevents a short model
// curve from obtaining an artificially good score by silently omitting onset
// or decay observations.
ComparisonMetrics ComparePositiveSeries(
    const std::vector<SeriesPoint>& model,
    const std::vector<SeriesPoint>& reference,
    const ComparisonPolicy& policy);

// One host/MPI partition contributes immutable observations, already-closed
// conservation rows, and measured resource usage.  Rank is an evidence label,
// not an array position; partition order may therefore change freely.
struct RankPartition {
  unsigned rank = 0;
  std::vector<Output::ParticleObservation> observations;
  std::vector<Adapters::LedgerRow> ledgerRows;
  double elapsedSeconds = 0.0;
  std::size_t peakResidentBytes = 0;
};

struct IntegrationBudget {
  double maximumParticleMaxToMean = 2.0;
  double maximumWallSeconds = std::numeric_limits<double>::infinity();
  std::size_t maximumTotalResidentBytes =
      std::numeric_limits<std::size_t>::max();
};

struct IntegrationAudit {
  Core::Status status;
  bool withinBudget = false;
  std::vector<Output::ParticleObservation> globalObservations;
  std::vector<Adapters::LedgerRow> globalLedgerRows;
  double particleMaxToMean = 0.0;
  double maximumRankWallSeconds = 0.0;
  std::size_t totalPeakResidentBytes = 0;
  std::vector<std::string> budgetViolations;
};

// Canonically merge a coupled run without using floating MPI reductions:
// particles are sorted by stable physical ID and integer conservation rows are
// summed by (step,species).  Duplicate IDs/ranks, open ledgers, overflow, or a
// non-closing global ledger are structural errors.  Resource-budget misses are
// returned separately so their measured values remain reportable.
IntegrationAudit AuditIntegration(
    const std::vector<RankPartition>& partitions,
    const IntegrationBudget& budget);

}  // namespace Validation
}  // namespace SEP3D

#endif  // SEP3D_VALIDATION_VALIDATION_METRICS_H
