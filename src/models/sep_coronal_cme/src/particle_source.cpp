#include "sep_coronal_cme/particle_source.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>

namespace SEP { namespace CoronalCME { namespace {

bool Finite(double value) { return std::isfinite(value); }

double IntegralPower(double minimum, double maximum, double exponent) {
  if (std::abs(exponent + 1.0) < 1.0e-14)
    return std::log(maximum / minimum);
  return (std::pow(maximum, exponent + 1.0) -
          std::pow(minimum, exponent + 1.0)) / (exponent + 1.0);
}

std::uint64_t Mix(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

Vec3 Apply(const DiffusionTensor& tensor, Vec3 vector) {
  return {tensor.xx * vector.x + tensor.xy * vector.y + tensor.xz * vector.z,
          tensor.xy * vector.x + tensor.yy * vector.y + tensor.yz * vector.z,
          tensor.xz * vector.x + tensor.yz * vector.y + tensor.zz * vector.z};
}

}  // namespace

Core::Result<MomentumSpectrum> MomentumSpectrum::FixedPowerLaw(
    double minimum, double maximum, double exponent) {
  if (!(minimum > 0.0 && maximum > minimum && Finite(exponent)))
    return Core::Result<MomentumSpectrum>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "momentum spectrum requires finite 0<p_min<p_max and exponent");
  const double integral = IntegralPower(minimum, maximum, exponent);
  if (!(integral > 0.0) || !Finite(integral))
    return Core::Result<MomentumSpectrum>::Failure(
        Core::StatusCode::NumericalFailure,
        "momentum spectrum normalization is nonpositive or nonfinite");
  MomentumSpectrum result;
  result.minimum_ = minimum;
  result.maximum_ = maximum;
  result.exponent_ = exponent;
  result.normalization_ = 1.0 / integral;
  return Core::Result<MomentumSpectrum>::Success(result);
}

Core::Result<MomentumSpectrum> MomentumSpectrum::LocalCompressionDsa(
    double minimum, double maximum, double compression) {
  if (!(compression > 1.0))
    return Core::Result<MomentumSpectrum>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "DSA compression ratio must exceed one");
  const double phaseSpaceExponent = 3.0 * compression / (compression - 1.0);
  return FixedPowerLaw(minimum, maximum, 2.0 - phaseSpaceExponent);
}

double MomentumSpectrum::DensityPerMomentum(double momentum) const {
  if (momentum < minimum_ || momentum > maximum_) return 0.0;
  return normalization_ * std::pow(momentum, exponent_);
}

double MomentumSpectrum::Cdf(double momentum) const {
  if (momentum <= minimum_) return 0.0;
  if (momentum >= maximum_) return 1.0;
  return normalization_ * IntegralPower(minimum_, momentum, exponent_);
}

double MomentumSpectrum::InverseCdf(double uniform) const {
  const double clipped = std::max(0.0, std::min(1.0, uniform));
  if (std::abs(exponent_ + 1.0) < 1.0e-14)
    return minimum_ * std::pow(maximum_ / minimum_, clipped);
  const double power = exponent_ + 1.0;
  return std::pow(std::pow(minimum_, power) + clipped *
      (std::pow(maximum_, power) - std::pow(minimum_, power)), 1.0 / power);
}

double RadialSourceEnvelope(double radius, double full, double zero,
                            bool smooth) {
  if (!(Finite(radius) && full >= 0.0 && zero >= full)) return 0.0;
  if (radius <= full) return 1.0;
  if (radius >= zero || zero == full) return 0.0;
  if (!smooth) return 1.0;
  const double t = (radius - full) / (zero - full);
  return 1.0 - (10.0 * t * t * t - 15.0 * std::pow(t, 4) +
                6.0 * std::pow(t, 5));
}

double KeyedUniform01(std::uint64_t seed, std::uint64_t stream,
                      std::uint64_t generation, std::uint64_t tick,
                      std::uint64_t species, std::uint64_t patch,
                      std::uint64_t sample) {
  std::uint64_t state = Mix(seed);
  for (std::uint64_t value : {stream, generation, tick, species, patch, sample})
    state = Mix(state ^ Mix(value));
  // Use the leading 53 bits so conversion has the same exact resolution as a
  // double mantissa and never returns one.
  return static_cast<double>(state >> 11U) * 0x1.0p-53;
}

Core::Result<std::vector<std::uint64_t>> AllocateLargestRemainder(
    std::uint64_t samples, const std::vector<double>& rates) {
  if (samples == 0 || rates.empty())
    return Core::Result<std::vector<std::uint64_t>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "sample allocation requires samples and physical rates");
  double total = 0.0;
  for (double rate : rates) {
    if (rate < 0.0 || !Finite(rate))
      return Core::Result<std::vector<std::uint64_t>>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "physical allocation rates must be finite and nonnegative");
    total += rate;
  }
  if (!(total > 0.0))
    return Core::Result<std::vector<std::uint64_t>>::Failure(
        Core::StatusCode::InvalidState, "physical allocation rate is zero");
  std::vector<std::uint64_t> allocated(rates.size());
  std::vector<std::pair<double, std::size_t>> remainders;
  std::uint64_t assigned = 0;
  for (std::size_t i = 0; i < rates.size(); ++i) {
    const double exact = static_cast<double>(samples) * rates[i] / total;
    allocated[i] = static_cast<std::uint64_t>(std::floor(exact));
    assigned += allocated[i];
    remainders.push_back({exact - allocated[i], i});
  }
  std::sort(remainders.begin(), remainders.end(),
      [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
      });
  for (std::uint64_t i = assigned; i < samples; ++i)
    ++allocated[remainders[static_cast<std::size_t>(i - assigned)].second];
  return Core::Result<std::vector<std::uint64_t>>::Success(allocated);
}

Core::Status ValidateSourceSpecies(
    const std::vector<CompiledSpeciesIdentity>& species, int count) {
  if (count <= 0 || species.size() != static_cast<std::size_t>(count))
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "source species must cover every compiled slot");
  std::set<int> slots;
  std::set<std::string> ids;
  for (const auto& item : species) {
    if (item.stableId.empty() || item.chemicalSymbol.empty() ||
        item.compiledSlot < 0 || item.compiledSlot >= count ||
        !(item.massKg > 0.0) || item.nucleonCount < 0 ||
        !slots.insert(item.compiledSlot).second ||
        !ids.insert(item.stableId).second ||
        (item.sourceEnabled && item.chargeC == 0.0))
      return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
          "source species identity is incomplete, duplicate, or neutral");
  }
  return Core::Status::Success();
}

Core::Status ValidateSourceBudget(const SourceBudget& budget) {
  if (!(budget.requestedNumberRatePerS >= 0.0 &&
        budget.availableUpstreamNumberRatePerS >= 0.0 &&
        budget.maximumNumberFluxFraction >= 0.0 &&
        budget.maximumNumberFluxFraction <= 1.0 &&
        budget.requestedKineticPowerW >= 0.0 &&
        budget.availableKineticPowerW >= 0.0 &&
        budget.maximumNonthermalEnergyFraction >= 0.0 &&
        budget.maximumNonthermalEnergyFraction <= 1.0))
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "source budgets require nonnegative SI measures/unit caps");
  if (budget.requestedNumberRatePerS >
          budget.maximumNumberFluxFraction *
          budget.availableUpstreamNumberRatePerS ||
      budget.requestedKineticPowerW >
          budget.maximumNonthermalEnergyFraction *
          budget.availableKineticPowerW)
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "source exceeds number-flux or nonthermal-energy budget");
  return Core::Status::Success();
}

double EligiblePhysicalSourceRate(double candidateRate,
                                  SourcePatchTopology topology,
                                  bool terminated) {
  if (!(candidateRate >= 0.0) || !Finite(candidateRate) || terminated ||
      topology != SourcePatchTopology::OpenEligible) return 0.0;
  return candidateRate;
}

Core::Status ValidateSourcePopulationMeaning(SourcePopulationMeaning meaning,
                                             bool productionIntent) {
  if (meaning == SourcePopulationMeaning::NetFirstPassageAtReferenceSurface)
    return Core::Status::Success();
  if (productionIntent)
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
        "schema-5 production source is net first passage at finite L_ref");
  return Core::Status::Success();
}

Core::Status ValidateReferenceSurface(
    const std::vector<ReferenceSurfacePoint>& points, double distance,
    double tolerance) {
  if (points.empty() || !(distance > 0.0) || tolerance < 0.0)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "reference surface requires points and fixed L_ref>0");
  std::set<std::uint64_t> ids;
  for (const auto& point : points) {
    const Vec3 normal = Unit(point.outwardNormal);
    const Vec3 offset = point.referencePointM - point.frontPointM;
    if (point.stablePatchId == 0 || !ids.insert(point.stablePatchId).second ||
        Norm(normal) == 0.0 || Dot(offset, normal) <= 0.0 ||
        std::abs(Dot(offset, normal) - distance) > tolerance ||
        Norm(offset - distance * normal) > tolerance)
      return Core::Status::Failure(Core::StatusCode::InvalidState,
          "reference surface folds, overlaps, or violates fixed upstream offset");
  }
  return Core::Status::Success();
}

Core::Result<double> IntegrateTransportDepth(
    const std::vector<TransportDepthSample>& samples) {
  if (samples.size() < 2)
    return Core::Result<double>::Failure(Core::StatusCode::InvalidConfiguration,
                                        "transport-depth quadrature needs two samples");
  double result = 0.0;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (!(samples[i].inflowMPerS > 0.0) ||
        !(samples[i].normalDiffusionM2PerS > 0.0) ||
        !Finite(samples[i].inflowMPerS) ||
        !Finite(samples[i].normalDiffusionM2PerS) ||
        (i && !(samples[i].distanceM > samples[i - 1].distanceM)))
      return Core::Result<double>::Failure(Core::StatusCode::InvalidState,
          "transport depth requires ordered distance and positive u_in/kappa_nn");
    if (i) {
      const double left = samples[i - 1].inflowMPerS /
          samples[i - 1].normalDiffusionM2PerS;
      const double right = samples[i].inflowMPerS /
          samples[i].normalDiffusionM2PerS;
      result += 0.5 * (left + right) *
          (samples[i].distanceM - samples[i - 1].distanceM);
    }
  }
  return Core::Result<double>::Success(result);
}

Core::Result<FocusedEscapeDistribution> BuildFocusedEscapeDistribution(
    double speed, int polarity, double bDotN, double advection,
    int intervals, double absoluteTolerance) {
  if (!(speed > 0.0) || (polarity != -1 && polarity != 1) ||
      !Finite(bDotN) || !Finite(advection) || intervals < 4 ||
      intervals % 2 != 0 || absoluteTolerance < 0.0)
    return Core::Result<FocusedEscapeDistribution>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "focused escape requires speed, polarity, and even quadrature intervals");
  FocusedEscapeDistribution result;
  result.mu.resize(static_cast<std::size_t>(intervals + 1));
  result.cdf.resize(result.mu.size());
  const double h = 2.0 / intervals;
  double accumulated = 0.0;
  double previousFlux = 0.0;
  for (int i = 0; i <= intervals; ++i) {
    const double mu = -1.0 + i * h;
    const double normalSpeed = advection + speed * polarity * mu * bDotN;
    const double flux = 0.5 * std::max(0.0, normalSpeed);
    result.mu[static_cast<std::size_t>(i)] = mu;
    if (i) accumulated += 0.5 * h * (previousFlux + flux);
    result.cdf[static_cast<std::size_t>(i)] = accumulated;
    previousFlux = flux;
  }
  result.normalizationMPerS = accumulated;
  if (std::abs(bDotN) > 1.0e-12) {
    result.diagnosticMuCApplicable = true;
    result.diagnosticMuC = -advection / (speed * polarity * bDotN);
  }
  if (accumulated == 0.0) {
    result.status = FocusedEscapeStatus::NoFocusedEscape;
    return Core::Result<FocusedEscapeDistribution>::Success(result);
  }
  if (accumulated <= absoluteTolerance) {
    result.status = FocusedEscapeStatus::UnresolvedPositiveFluxNumerics;
    return Core::Result<FocusedEscapeDistribution>::Success(result);
  }
  result.status = FocusedEscapeStatus::Admissible;
  for (double& value : result.cdf) value /= accumulated;
  result.cdf.back() = 1.0;
  return Core::Result<FocusedEscapeDistribution>::Success(result);
}

Core::Result<double> SampleFocusedMu(
    const FocusedEscapeDistribution& distribution, double uniform) {
  if (distribution.status != FocusedEscapeStatus::Admissible ||
      distribution.mu.size() < 2 || distribution.mu.size() !=
          distribution.cdf.size() || uniform < 0.0 || uniform >= 1.0)
    return Core::Result<double>::Failure(Core::StatusCode::InvalidState,
                                        "focused distribution is not sampleable");
  auto upper = std::upper_bound(distribution.cdf.begin(),
                                distribution.cdf.end(), uniform);
  std::size_t right = static_cast<std::size_t>(
      std::distance(distribution.cdf.begin(), upper));
  if (right == 0) right = 1;
  if (right >= distribution.cdf.size()) right = distribution.cdf.size() - 1;
  const double c0 = distribution.cdf[right - 1];
  const double c1 = distribution.cdf[right];
  const double f = c1 > c0 ? (uniform - c0) / (c1 - c0) : 0.0;
  return Core::Result<double>::Success(
      distribution.mu[right - 1] + f *
      (distribution.mu[right] - distribution.mu[right - 1]));
}

Core::Result<Vec3> ParkerConormalDirection(const DiffusionTensor& tensor,
                                           Vec3 normal, double minimum) {
  const Vec3 unit = Unit(normal);
  const Vec3 applied = Apply(tensor, unit);
  const double normalDiffusion = Dot(unit, applied);
  if (Norm(unit) == 0.0 || !Finite(normalDiffusion) ||
      !(normalDiffusion >= minimum) || !(minimum > 0.0))
    return Core::Result<Vec3>::Failure(Core::StatusCode::InvalidState,
        "conormal requires finite n.kappa.n above the configured threshold");
  return Core::Result<Vec3>::Success(applied / normalDiffusion);
}

Core::Result<double> ShockAdjacentEscapeProbability(
    double inflow, double delta, double length, double diffusion) {
  if (!(inflow > 0.0 && delta >= 0.0 && length > delta && diffusion > 0.0))
    return Core::Result<double>::Failure(Core::StatusCode::InvalidConfiguration,
                                        "absorbing verification parameters are invalid");
  return Core::Result<double>::Success(
      std::expm1(inflow * delta / diffusion) /
      std::expm1(inflow * length / diffusion));
}

double FourMomentumInvariant(const FourMomentum& state) {
  const double c = Constants::kSpeedOfLightMPerS;
  return state.totalEnergyJ * state.totalEnergyJ -
      Dot(state.momentumKgMPerS, state.momentumKgMPerS) * c * c;
}

Core::Result<FourMomentum> BoostFourMomentum(
    const FourMomentum& source, Vec3 velocity) {
  const double c = Constants::kSpeedOfLightMPerS;
  const double speed = Norm(velocity);
  if (!(source.totalEnergyJ > 0.0) || speed >= c || !Finite(speed))
    return Core::Result<FourMomentum>::Failure(Core::StatusCode::InvalidState,
                                              "four-momentum boost is unphysical");
  if (speed == 0.0) return Core::Result<FourMomentum>::Success(source);
  const Vec3 direction = velocity / speed;
  const double beta = speed / c;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  const double parallel = Dot(source.momentumKgMPerS, direction);
  const Vec3 perpendicular = source.momentumKgMPerS - parallel * direction;
  FourMomentum result;
  result.totalEnergyJ = gamma * (source.totalEnergyJ +
      speed * parallel);
  const double boostedParallel = gamma *
      (parallel + speed * source.totalEnergyJ / (c * c));
  result.momentumKgMPerS = perpendicular + boostedParallel * direction;
  return Core::Result<FourMomentum>::Success(result);
}

double IntegratePiecewiseLinearRate(
    const std::vector<std::pair<double, double>>& points,
    double begin, double end) {
  if (points.size() < 2 || end <= begin) return 0.0;
  double result = 0.0;
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double left = std::max(begin, points[i - 1].first);
    const double right = std::min(end, points[i].first);
    if (right <= left) continue;
    const double dt = points[i].first - points[i - 1].first;
    if (!(dt > 0.0)) return 0.0;
    const auto value = [&](double time) {
      const double f = (time - points[i - 1].first) / dt;
      return points[i - 1].second + f *
          (points[i].second - points[i - 1].second);
    };
    result += 0.5 * (value(left) + value(right)) * (right - left);
  }
  return result;
}

bool InjectionCommitRegistry::CommitOnce(std::uint64_t generation,
                                         std::uint64_t tick) {
  return generation != 0 && committed_.insert({generation, tick}).second;
}

bool CohortKey::operator<(const CohortKey& other) const {
  return std::tie(speciesId, patchLineage, generation, birthTick) <
      std::tie(other.speciesId, other.patchLineage,
               other.generation, other.birthTick);
}

Core::Status CohortLedger::Add(const CohortKey& key, ReleaseLedgerTerm term,
                               const LedgerMeasure& measure) {
  if (key.speciesId.empty() || key.patchLineage == 0 || key.generation == 0 ||
      measure.representedNumber < 0.0 || measure.birthKineticEnergyJ < 0.0 ||
      measure.eventKineticEnergyJ < 0.0)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "cohort ledger key/measure is invalid");
  auto& destination = entries_[{key, term}];
  destination.representedNumber += measure.representedNumber;
  destination.birthKineticEnergyJ += measure.birthKineticEnergyJ;
  destination.eventKineticEnergyJ += measure.eventKineticEnergyJ;
  destination.eventFourMomentum.totalEnergyJ +=
      measure.eventFourMomentum.totalEnergyJ;
  destination.eventFourMomentum.momentumKgMPerS =
      destination.eventFourMomentum.momentumKgMPerS +
      measure.eventFourMomentum.momentumKgMPerS;
  return Core::Status::Success();
}

LedgerMeasure CohortLedger::Get(const CohortKey& key,
                                ReleaseLedgerTerm term) const {
  const auto found = entries_.find({key, term});
  return found == entries_.end() ? LedgerMeasure{} : found->second;
}

Core::Result<FiniteHorizonReturn> ReduceFiniteHorizonReturn(
    const CohortLedger& ledger, const CohortKey& key,
    double birthEnd, double horizon, double requestedAge) {
  if (!(horizon >= birthEnd && requestedAge >= 0.0))
    return Core::Result<FiniteHorizonReturn>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "finite-horizon reduction requires nonnegative follow-up");
  FiniteHorizonReturn result;
  result.committedNumber = ledger.Get(
      key, ReleaseLedgerTerm::CommittedFirstPassageRelease).representedNumber;
  result.delayedReturnNumber = ledger.Get(
      key, ReleaseLedgerTerm::DelayedFrontReturn).representedNumber;
  result.survivingInventoryNumber = ledger.Get(
      key, ReleaseLedgerTerm::SurvivingUpstreamInventory).representedNumber;
  result.availableFollowupS = horizon - birthEnd;
  result.rightCensored = result.availableFollowupS < requestedAge;
  if (result.delayedReturnNumber > result.committedNumber)
    return Core::Result<FiniteHorizonReturn>::Failure(
        Core::StatusCode::DataIntegrityFailure,
        "delayed return exceeds committed first-passage release");
  if (result.committedNumber == 0.0) {
    result.validity = FractionValidity::InapplicableZeroDenominator;
  } else {
    result.noFrontReturnFraction = 1.0 -
        result.delayedReturnNumber / result.committedNumber;
  }
  return Core::Result<FiniteHorizonReturn>::Success(result);
}

Core::Result<LossCapResult> EvaluateRepresentedLossCaps(
    double bornNumber, double lostNumber, double bornEnergy,
    double lostBirthEnergy, double numberCap, double energyCap,
    bool removingFront) {
  if (!removingFront)
    return Core::Result<LossCapResult>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "front-loss caps are inactive for a nonabsorbing branch");
  if (!(bornNumber >= 0.0 && lostNumber >= 0.0 && lostNumber <= bornNumber &&
        bornEnergy >= 0.0 && lostBirthEnergy >= 0.0 &&
        lostBirthEnergy <= bornEnergy && numberCap >= 0.0 && numberCap <= 1.0 &&
        energyCap >= 0.0 && energyCap <= 1.0))
    return Core::Result<LossCapResult>::Failure(
        Core::StatusCode::InvalidConfiguration, "loss-cap measure is invalid");
  LossCapResult result;
  if (bornNumber == 0.0)
    result.numberValidity = FractionValidity::InapplicableZeroDenominator;
  else result.numberFraction = lostNumber / bornNumber;
  if (bornEnergy == 0.0)
    result.birthEnergyValidity = FractionValidity::InapplicableZeroDenominator;
  else result.birthEnergyFraction = lostBirthEnergy / bornEnergy;
  result.withinCaps =
      (result.numberValidity != FractionValidity::Valid ||
       result.numberFraction <= numberCap) &&
      (result.birthEnergyValidity != FractionValidity::Valid ||
       result.birthEnergyFraction <= energyCap);
  return Core::Result<LossCapResult>::Success(result);
}

} }  // namespace SEP::CoronalCME
