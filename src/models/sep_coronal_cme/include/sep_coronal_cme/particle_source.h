#ifndef SEP_CORONAL_CME_PARTICLE_SOURCE_H
#define SEP_CORONAL_CME_PARTICLE_SOURCE_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

// MomentumSpectrum is a normalized dN/dp law.  For local test-particle DSA,
// f(p)~p^-q with q=3r/(r-1), so the exact spherical momentum Jacobian gives
// dN/dp~4*pi*p^2*f(p)~p^(2-q).
class MomentumSpectrum {
 public:
  static Core::Result<MomentumSpectrum> FixedPowerLaw(
      double minimumMomentumKgMPerS, double maximumMomentumKgMPerS,
      double numberExponent);
  static Core::Result<MomentumSpectrum> LocalCompressionDsa(
      double minimumMomentumKgMPerS, double maximumMomentumKgMPerS,
      double compressionRatio);

  double DensityPerMomentum(double momentumKgMPerS) const;
  double Cdf(double momentumKgMPerS) const;
  double InverseCdf(double unitUniform) const;
  double NumberExponent() const noexcept { return exponent_; }

 private:
  double minimum_ = 0.0;
  double maximum_ = 0.0;
  double exponent_ = 0.0;
  double normalization_ = 0.0;
};

double RadialSourceEnvelope(double radiusM, double fullStrengthRadiusM,
                            double zeroRadiusM, bool smoothQuintic);

// A counter-based stream derives every draw from the full semantic key.  No
// mutable generator state exists, so adding a draw to one stream cannot move
// momentum, pitch, position, or species draws in another stream.
double KeyedUniform01(std::uint64_t campaignSeed, std::uint64_t streamId,
                      std::uint64_t generation, std::uint64_t tick,
                      std::uint64_t species, std::uint64_t patch,
                      std::uint64_t sample);

Core::Result<std::vector<std::uint64_t>> AllocateLargestRemainder(
    std::uint64_t samples, const std::vector<double>& physicalRates);

struct CompiledSpeciesIdentity {
  std::string stableId;
  int compiledSlot = -1;
  std::string chemicalSymbol;
  double massKg = 0.0;
  double chargeC = 0.0;
  int nucleonCount = 0;
  bool sourceEnabled = true;
};
Core::Status ValidateSourceSpecies(
    const std::vector<CompiledSpeciesIdentity>& species,
    int compiledSpeciesCount);

struct SourceBudget {
  double requestedNumberRatePerS = 0.0;
  double availableUpstreamNumberRatePerS = 0.0;
  double maximumNumberFluxFraction = 0.0;
  double requestedKineticPowerW = 0.0;
  double availableKineticPowerW = 0.0;
  double maximumNonthermalEnergyFraction = 0.0;
};
Core::Status ValidateSourceBudget(const SourceBudget& budget);

enum class SourcePatchTopology {
  OpenEligible,
  ClosedDiagnoseOnly,
  TransitionClearance
};
double EligiblePhysicalSourceRate(double candidateRatePerS,
                                  SourcePatchTopology topology,
                                  bool sourceTerminated);

enum class SourcePopulationMeaning {
  NetFirstPassageAtReferenceSurface,
  GrossShockEmission,
  AcceleratorDerivedFreeEscape
};
Core::Status ValidateSourcePopulationMeaning(SourcePopulationMeaning meaning,
                                             bool productionIntent);

struct ReferenceSurfacePoint {
  std::uint64_t stablePatchId = 0;
  Vec3 frontPointM;
  Vec3 referencePointM;
  Vec3 outwardNormal;
};
Core::Status ValidateReferenceSurface(
    const std::vector<ReferenceSurfacePoint>& points,
    double referenceDistanceM, double absoluteToleranceM);

struct TransportDepthSample {
  double distanceM = 0.0;
  double inflowMPerS = 0.0;
  double normalDiffusionM2PerS = 0.0;
};
Core::Result<double> IntegrateTransportDepth(
    const std::vector<TransportDepthSample>& samples);

enum class FocusedEscapeStatus {
  Admissible,
  NoFocusedEscape,
  UnresolvedPositiveFluxNumerics
};
struct FocusedEscapeDistribution {
  FocusedEscapeStatus status = FocusedEscapeStatus::NoFocusedEscape;
  double normalizationMPerS = 0.0;
  bool diagnosticMuCApplicable = false;
  double diagnosticMuC = 0.0;
  std::vector<double> mu;
  std::vector<double> cdf;
};
Core::Result<FocusedEscapeDistribution> BuildFocusedEscapeDistribution(
    double particleSpeedMPerS, int magneticPolarity,
    double fieldNormalProjection, double advectiveNormalSpeedMPerS,
    int quadratureIntervals, double absoluteFluxToleranceMPerS);
Core::Result<double> SampleFocusedMu(
    const FocusedEscapeDistribution& distribution, double unitUniform);

struct DiffusionTensor {
  double xx = 0.0, xy = 0.0, xz = 0.0;
  double yy = 0.0, yz = 0.0, zz = 0.0;
};
Core::Result<Vec3> ParkerConormalDirection(const DiffusionTensor& tensor,
                                           Vec3 normal,
                                           double minimumNormalDiffusionM2PerS);

Core::Result<double> ShockAdjacentEscapeProbability(
    double inflowMPerS, double placementDistanceM,
    double referenceDistanceM, double normalDiffusionM2PerS);

struct FourMomentum {
  double totalEnergyJ = 0.0;
  Vec3 momentumKgMPerS;
};
Core::Result<FourMomentum> BoostFourMomentum(
    const FourMomentum& sourceFrame, Vec3 sourceVelocityMPerS);
double FourMomentumInvariant(const FourMomentum& state);

double IntegratePiecewiseLinearRate(
    const std::vector<std::pair<double, double>>& timeRate,
    double beginTimeS, double endTimeS);

class InjectionCommitRegistry {
 public:
  bool CommitOnce(std::uint64_t generation, std::uint64_t tick);
 private:
  std::set<std::pair<std::uint64_t, std::uint64_t>> committed_;
};

enum class ReleaseLedgerTerm {
  CandidateReferenceRelease,
  NoFocusedEscapeExcluded,
  CommittedFirstPassageRelease,
  ImmediateShockAdjacentReturn,
  DelayedFrontReturn,
  NetFirstPassageRelease,
  SurvivingUpstreamInventory,
  GrossShockAdjacentEmission,
  HitEscapeSurfaceBeforeShock,
  TransitionSheetContact
};

struct CohortKey {
  std::string speciesId;
  std::uint64_t patchLineage = 0;
  std::uint64_t generation = 0;
  std::uint64_t birthTick = 0;
  bool operator<(const CohortKey& other) const;
};
struct LedgerMeasure {
  double representedNumber = 0.0;
  double birthKineticEnergyJ = 0.0;
  double eventKineticEnergyJ = 0.0;
  FourMomentum eventFourMomentum;
};
class CohortLedger {
 public:
  Core::Status Add(const CohortKey& key, ReleaseLedgerTerm term,
                   const LedgerMeasure& measure);
  LedgerMeasure Get(const CohortKey& key, ReleaseLedgerTerm term) const;
 private:
  std::map<std::pair<CohortKey, ReleaseLedgerTerm>, LedgerMeasure> entries_;
};

enum class FractionValidity { Valid, InapplicableZeroDenominator };
struct FiniteHorizonReturn {
  FractionValidity validity = FractionValidity::Valid;
  double committedNumber = 0.0;
  double delayedReturnNumber = 0.0;
  double survivingInventoryNumber = 0.0;
  double noFrontReturnFraction = 0.0;
  double availableFollowupS = 0.0;
  bool rightCensored = false;
};
Core::Result<FiniteHorizonReturn> ReduceFiniteHorizonReturn(
    const CohortLedger& ledger, const CohortKey& key,
    double birthEndTimeS, double horizonTimeS, double requestedAgeS);

struct LossCapResult {
  FractionValidity numberValidity = FractionValidity::Valid;
  FractionValidity birthEnergyValidity = FractionValidity::Valid;
  double numberFraction = 0.0;
  double birthEnergyFraction = 0.0;
  bool withinCaps = false;
};
Core::Result<LossCapResult> EvaluateRepresentedLossCaps(
    double bornRepresentedNumber, double lostRepresentedNumber,
    double bornKineticEnergyJ, double lostBirthKineticEnergyJ,
    double maximumNumberFraction, double maximumBirthEnergyFraction,
    bool removingFrontBranch);

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_PARTICLE_SOURCE_H
