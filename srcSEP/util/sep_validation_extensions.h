#ifndef SEP_UTIL_SEP_VALIDATION_EXTENSIONS_H
#define SEP_UTIL_SEP_VALIDATION_EXTENSIONS_H

#include "sep_evidence.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace ValidationExtensions {

// WP57: lineage, rather than a split child, is the independent sampling unit.
// contribution is already expressed in the product's physical units.  Speed
// and Jacobian remain explicit so broad energy bins can be accumulated from
// particle-resolved relativistic kinematics instead of one representative
// speed.
enum class EstimatorKind { SnapshotDensity, CrossingFlux };

struct LineageSample {
  std::uint64_t rootLineage = 0;
  std::uint64_t samplingWindow = 0;
  std::size_t bin = 0;
  double statisticalWeight = 0.0;
  double physicalValue = 0.0;
  double speedMPerS = 0.0;
  double jacobian = 1.0;
};

struct LineageEstimate {
  Transport::Status status;
  EstimatorKind kind = EstimatorKind::SnapshotDensity;
  std::string independenceUnit = "root-lineage";
  std::vector<double> estimate;
  std::vector<double> variance;
  std::vector<double> covariance;
  std::vector<double> effectiveSampleSize;
  std::uint64_t independentRoots = 0;
};

LineageEstimate EstimateByLineage(const std::vector<LineageSample>& samples,
                                  std::size_t bins,
                                  EstimatorKind kind);

// WP58: detector response is a versioned matrix.  The flattened index is
// [channel][species][direction][true-energy].  Values are detection
// probabilities and can therefore express contamination and angular response.
enum class DeadTimeModel { None, Nonparalyzable, Paralyzable };

struct InstrumentResponseMatrix {
  std::string instrumentId;
  std::string calibrationVersion;
  std::string checksum;
  std::string dataLevel;
  double validFromS = 0.0;
  double validUntilS = 0.0;
  std::size_t channels = 0;
  std::size_t species = 0;
  std::size_t directions = 0;
  std::size_t trueEnergyBins = 0;
  std::vector<double> probability;
};

struct CountObservation {
  double counts = 0.0;
  double backgroundCounts = 0.0;
  bool saturated = false;
  double saturationThreshold = 0.0;
};

struct ResponseFoldResult {
  Transport::Status status;
  std::vector<double> incidentCounts;
  std::vector<double> recordedCounts;
  std::vector<double> poissonLogLikelihood;
  std::vector<bool> censored;
};

Transport::Status ValidateResponseMatrix(
    const InstrumentResponseMatrix& matrix);
ResponseFoldResult FoldInstrumentCounts(
    const InstrumentResponseMatrix& matrix,
    const std::vector<double>& trueCounts,
    const std::vector<CountObservation>& observations,
    double exposureS, double deadTimeS, DeadTimeModel deadTimeModel);

// WP59: a native row records callbacks observed in the linked application.
// Source-only tests may construct this structure, but ValidateNativeMatrixRow
// refuses to promote them unless observedLevel is NativeAmps or higher.
struct NativeMatrixTrace {
  Evidence::Level observedLevel = Evidence::Level::SourceIntegration;
  std::string configurationFingerprint;
  std::string executableChecksum;
  std::string compiler;
  std::string flags;
  std::string moverId;
  std::string coefficientProviderId;
  std::string backgroundOwnerId;
  std::string turbulenceOwnerId;
  bool moverEntered = false;
  bool coefficientEntered = false;
  bool turbulenceEntered = false;
  bool couplingEntered = false;
};

Transport::Status ValidateNativeMatrixRow(
    const NativeMatrixTrace& trace, bool couplingExpected);

// WP60: discrete identities and counters require exact equality across
// decompositions.  Floating ledgers use a declared relative bound when the
// configured global reduction is not bitwise reproducible.
struct DecompositionSignature {
  int ranks = 1;
  int threads = 1;
  std::string configurationFingerprint;
  std::string particleIdentityHash;
  std::string transactionOrderHash;
  std::string restartHash;
  std::uint64_t sourceEvents = 0;
  std::uint64_t transactions = 0;
  std::vector<double> ledger;
};

Transport::Status CompareDecompositions(
    const DecompositionSignature& reference,
    const DecompositionSignature& candidate,
    double ledgerRelativeTolerance);

// WP61: refinement order is fitted from all levels by least squares in
// log(error) versus log(resolution), rather than inferred from one favorable
// pair.  The manufactured case records which PDE terms are nonzero so an
// apparently converged but incomplete equation cannot pass.
struct ManufacturedCase {
  std::string id;
  std::vector<std::string> activeTerms;
  std::vector<double> resolution;
  std::vector<double> error;
};

struct RefinementFit {
  Transport::Status status;
  double observedOrder = 0.0;
  double intercept = 0.0;
  double rSquared = 0.0;
};

RefinementFit FitRefinementOrder(const ManufacturedCase& testCase);

// WP62: a frozen statistical specification makes power, significance, and
// multiplicity decisions before outcomes are known.  The normal approximation
// is appropriate for bounded smoke planning; event campaigns may replace it
// with exact or simulation-based power while retaining this metadata schema.
struct PowerSpecification {
  std::string observable;
  double standardizedEffect = 0.0;
  double alpha = 0.05;
  double power = 0.8;
  std::size_t comparisons = 1;
  std::string seedPanelVersion;
};

struct PowerResult {
  Transport::Status status;
  std::size_t requiredSamples = 0;
  double adjustedAlpha = 0.0;
};

PowerResult PlanNormalMeanPower(const PowerSpecification& specification);
Transport::ScalarResult KolmogorovSmirnovStatistic(
    const std::vector<double>& sortedSamples,
    const std::vector<double>& expectedCdfAtSamples);

// WP63: external evidence remains fail closed.  A manifest must identify input
// lineage, checksums, calibration/validation role, preprocessing, exclusions,
// configuration prior, and metrics before execution.
enum class EventRole { Calibration, Validation, HeldOut };

struct ExternalEventManifest {
  std::string eventId;
  EventRole role = EventRole::Calibration;
  std::string backgroundChecksum;
  std::string geometryChecksum;
  std::string shockChecksum;
  std::string calibrationChecksum;
  std::string preprocessing;
  std::string exclusions;
  std::string configurationPrior;
  std::vector<std::string> metrics;
  bool synthetic = true;
};

Transport::Status ValidateExternalEventManifest(
    const ExternalEventManifest& manifest,
    Evidence::Level claimedLevel);

struct EventMetrics {
  double onsetS = 0.0;
  double peak = 0.0;
  double fluence = 0.0;
  double anisotropy = 0.0;
  double spectralIndex = 0.0;
  double profileRmse = 0.0;
};

Transport::ScalarResult NormalizedMetricDistance(
    const EventMetrics& model, const EventMetrics& reference,
    const EventMetrics& oneSigma);

// WP64: scalability evidence separates deterministic algorithmic work from
// noisy wall time.  Communication and memory complexity are normalized by the
// declared problem size; a wall-time gate is valid only for matching
// environment fingerprints.
struct ScalingSample {
  std::string workloadId;
  std::string environmentFingerprint;
  int ranks = 1;
  int threadsPerRank = 1;
  std::uint64_t problemItems = 0;
  std::uint64_t operations = 0;
  std::uint64_t communicationBytes = 0;
  std::uint64_t peakResidentBytes = 0;
  std::uint64_t peakQueueRecords = 0;
  double wallSeconds = 0.0;
};

struct ScalingGate {
  Transport::Status status;
  double parallelEfficiency = 0.0;
  double operationsPerItem = 0.0;
  double bytesPerItem = 0.0;
  double memoryPerItem = 0.0;
  double queuePerItem = 0.0;
};

ScalingGate CompareScaling(const ScalingSample& reference,
                           const ScalingSample& candidate,
                           double minimumEfficiency,
                           double maximumWorkRatio,
                           double maximumMemoryRatio);

}  // namespace ValidationExtensions
}  // namespace SEP

#endif  // SEP_UTIL_SEP_VALIDATION_EXTENSIONS_H
