#ifndef SHIELDSIM_STAGE6_PHYSICS_HH
#define SHIELDSIM_STAGE6_PHYSICS_HH

/*
 * Deterministic Stage-6 LET and spectral-hardness mathematics.
 *
 * This module intentionally contains no Geant4 types.  It is the independently
 * testable boundary between transport tallies and the reported quantities:
 *
 *   - LogarithmicLETBinning defines explicit, non-clamping linear-LET bins.
 *   - IntegratePowerLawAbove performs exact partial-bin threshold integration.
 *   - EstimateRatioOfMeans propagates event-level numerator/denominator
 *     variance and covariance and constructs a bounded-domain Fieller interval.
 *
 * Production transport supplies step fluence/LET values and primary-proton
 * rear-face fluence tallies.  Keeping the numerical definitions here prevents
 * histogram centers or an uncorrelated ratio approximation from silently
 * replacing the approved Stage-6 equations.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace shield_sim {

class LogarithmicLETBinning {
public:
  LogarithmicLETBinning() = default;
  LogarithmicLETBinning(double minimumKeVPerUm,double maximumKeVPerUm,
                        int binsPerDecade);

  // Returns -1 for underflow and BinCount() for overflow.  A value exactly on
  // the final edge is included in the final bin; no value is ever clamped.
  long Locate(double linearLETKeVPerUm) const;
  std::size_t BinCount() const { return edges_.empty()?0:edges_.size()-1; }
  double Lower(std::size_t index) const;
  double Upper(std::size_t index) const;
  double Width(std::size_t index) const;
  double GeometricCenter(std::size_t index) const;
  const std::vector<double>& Edges() const { return edges_; }

private:
  std::vector<double> edges_;
};

// One differential-fluence interval represented by the exact local power law
// phi(E)=valueAtReference*(E/referenceEnergy)^spectralIndex.  Setting the index
// to zero gives a piecewise-constant interval.
struct PowerLawFluenceBin {
  double lowerMeV = 0.0;
  double upperMeV = 0.0;
  double referenceEnergyMeV = 0.0;
  double valueAtReferencePerCm2PerMeV = 0.0;
  double spectralIndex = 0.0;
};

double IntegratePowerLawAbove(const PowerLawFluenceBin& bin,double thresholdMeV);

struct HardnessFoldResult {
  bool defined = false;
  double numeratorFluencePerCm2 = 0.0;
  double denominatorFluencePerCm2 = 0.0;
  double value = 0.0;
  std::string reason;
};

HardnessFoldResult FoldH100Over10(const std::vector<PowerLawFluenceBin>& bins);

// Sufficient event-level moments for a ratio R=mean(N)/mean(D).  The sums of
// squares and products include every history, including implicit zero scores.
struct RatioMoments {
  long histories = 0;
  double numeratorSum = 0.0;
  double denominatorSum = 0.0;
  double numeratorSumSquares = 0.0;
  double denominatorSumSquares = 0.0;
  double numeratorDenominatorSumProducts = 0.0;
};

enum class RatioStatus {
  Available,
  ZeroDenominator,
  StatisticallyUnresolvedDenominator,
  InvalidMoments
};

std::string RatioStatusName(RatioStatus status);

struct RatioEstimate {
  RatioStatus status = RatioStatus::InvalidMoments;
  double numeratorMean = 0.0;
  double denominatorMean = 0.0;
  double numeratorVarianceOfMean = 0.0;
  double denominatorVarianceOfMean = 0.0;
  double covarianceOfMeans = 0.0;
  double value = 0.0;
  double standardError = 0.0;
  double interval95Low = 0.0;
  double interval95High = 0.0;
  std::string reason;
};

// z=1.959963984540054 is the approved two-sided 95% normal quantile.  Fieller
// is used because the denominator is estimated and correlated with the
// numerator.  If its confidence set includes zero, no finite ratio is emitted.
RatioEstimate EstimateRatioOfMeans(
    const RatioMoments& moments,double z=1.959963984540054);

} // namespace shield_sim

#endif
