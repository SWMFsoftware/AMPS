// ============================================================================
// Phase-V integration and scientific-validation acceptance tests.
//
// These callbacks use production Phase-P/A/O kernels and the Phase-V audit
// implementation.  They are not substitutes for a linked MPI campaign or an
// observational event: those evidence classes are registered by the Python
// campaign runner and explicitly SKIP when their external prerequisites are
// absent.  The tests here establish that the numerical comparison machinery,
// rank merge, analytical transport limits, and SWCME source statistics are
// correct before expensive evidence is interpreted.
// ============================================================================

#include "sep3d_test_registry.h"

#include "focused_transport.h"
#include "parker_transport.h"
#include "swcme_source_adapter.h"
#include "validation_metrics.h"

#include "sep_injection_spectrum.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace A = SEP3D::Adapters;
namespace C = SEP3D::Core;
namespace T = SEP3D::Transport;
namespace V = SEP3D::Validation;
using SEP3D::Testing::Metric;
using SEP3D::Testing::Result;

Result Pass(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Pass;
  result.message = message;
  return result;
}

Result Fail(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Fail;
  result.message = message;
  return result;
}

T::RandomKey Key(std::uint64_t particle, T::RandomPurpose purpose,
                 std::uint64_t step = 0) {
  T::RandomKey key;
  key.campaignSeed = UINT64_C(0x56414c4944415445);  // "VALIDATE"
  key.particleId = particle;
  key.step = step;
  key.purpose = purpose;
  return key;
}

SEP3D::Output::ParticleObservation Observation(std::uint64_t id) {
  SEP3D::Output::ParticleObservation observation;
  observation.stableId = id;
  observation.cellId = 1;
  observation.species = 0;
  observation.positionM = C::Vec3(static_cast<double>(id), 0.0, 0.0);
  observation.momentumKgMPerS = 1.0e-20;
  observation.restMassKg = C::Const::m_p;
  observation.statisticalWeight = 1.0;
  return observation;
}

A::LedgerRow ClosedRow(std::uint64_t step, std::uint64_t activeStart,
                       std::uint64_t injected, std::uint64_t activeEnd,
                       std::uint64_t escaped, std::uint64_t absorbed) {
  A::LedgerRow row;
  row.key = {step, 0};
  row.activeStart = activeStart;
  row.injected = injected;
  row.activeEnd = activeEnd;
  row.advanced = activeEnd;
  row.escaped = escaped;
  row.absorbed = absorbed;
  row.closed = true;
  return row;
}

Result RunINT3D01() {
  V::RankPartition first, second;
  first.rank = 7;
  second.rank = 2;
  for (std::uint64_t id : {6, 2, 4}) first.observations.push_back(Observation(id));
  for (std::uint64_t id : {5, 1, 3}) second.observations.push_back(Observation(id));
  const V::IntegrationBudget budget;
  const V::IntegrationAudit left = V::AuditIntegration({first, second}, budget);

  // The same physical set is repartitioned and rank records are reversed.
  // A correct global gather depends only on stable identity, not on sender or
  // receive order.
  V::RankPartition alternate0, alternate1;
  alternate0.rank = 2;
  alternate1.rank = 7;
  for (std::uint64_t id : {1, 4, 5})
    alternate0.observations.push_back(Observation(id));
  for (std::uint64_t id : {6, 3, 2})
    alternate1.observations.push_back(Observation(id));
  const V::IntegrationAudit right =
      V::AuditIntegration({alternate1, alternate0}, budget);
  if (!left.status.ok() || !right.status.ok() ||
      left.globalObservations.size() != 6 ||
      right.globalObservations.size() != 6)
    return Fail("valid rank partitions could not be merged");
  for (std::size_t i = 0; i < 6; ++i) {
    if (left.globalObservations[i].stableId != i + 1 ||
        right.globalObservations[i].stableId != i + 1)
      return Fail("global observation order depends on rank partitioning");
  }
  return Pass("stable-ID gather is bitwise ordered across rank order and repartitioning");
}

Result RunINT3D02() {
  V::RankPartition first, second;
  first.rank = 0;
  second.rank = 1;
  first.observations = {Observation(1), Observation(2)};
  second.observations = {Observation(3)};
  // Rank zero: 4+1 = 4+1. Rank one: 3+0 = 2+0+1.
  first.ledgerRows.push_back(ClosedRow(11, 4, 1, 4, 1, 0));
  second.ledgerRows.push_back(ClosedRow(11, 3, 0, 2, 0, 1));
  const V::IntegrationAudit audit =
      V::AuditIntegration({first, second}, V::IntegrationBudget{});
  if (!audit.status.ok() || audit.globalLedgerRows.size() != 1) {
    return Fail("closed local ledgers did not form one global row");
  }
  const A::LedgerRow& row = audit.globalLedgerRows.front();
  if (row.activeStart != 7 || row.injected != 1 || row.activeEnd != 6 ||
      row.escaped != 1 || row.absorbed != 1 || row.advanced != 6)
    return Fail("global integer conservation counters are incorrect");

  // Stable identities are globally unique.  Detecting this only within a rank
  // would allow a migration/restart bug to double-count a physical particle.
  second.observations.push_back(Observation(2));
  if (V::AuditIntegration({first, second}, V::IntegrationBudget{}).status.ok())
    return Fail("duplicate stable identity on another rank was accepted");
  return Pass("rank ledgers close exactly and cross-rank duplicate identity is rejected");
}

Result RunINT3D03() {
  V::RankPartition first, second;
  first.rank = 0;
  second.rank = 1;
  for (std::uint64_t id = 1; id <= 8; ++id)
    first.observations.push_back(Observation(id));
  for (std::uint64_t id = 9; id <= 12; ++id)
    second.observations.push_back(Observation(id));
  first.elapsedSeconds = 2.0;
  second.elapsedSeconds = 3.0;
  first.peakResidentBytes = 100;
  second.peakResidentBytes = 200;

  V::IntegrationBudget accepted;
  accepted.maximumParticleMaxToMean = 1.5;
  accepted.maximumWallSeconds = 4.0;
  accepted.maximumTotalResidentBytes = 400;
  const V::IntegrationAudit within =
      V::AuditIntegration({first, second}, accepted);
  if (!within.status.ok() || !within.withinBudget ||
      std::fabs(within.particleMaxToMean - 4.0 / 3.0) > 1.0e-14 ||
      within.maximumRankWallSeconds != 3.0 ||
      within.totalPeakResidentBytes != 300)
    return Fail("integration budget metrics are incorrect");

  V::IntegrationBudget strict = accepted;
  strict.maximumTotalResidentBytes = 250;
  const V::IntegrationAudit exceeded =
      V::AuditIntegration({first, second}, strict);
  if (!exceeded.status.ok() || exceeded.withinBudget ||
      exceeded.budgetViolations.size() != 1)
    return Fail("budget miss was hidden or misclassified as malformed evidence");
  return Pass("load balance, wall time, and memory budgets remain explicit reportable metrics");
}

Result RunVFY3D01() {
  const std::vector<V::SeriesPoint> reference = {
      {0.0, 1.0}, {1.0, 10.0}, {2.0, 100.0}, {3.0, 10.0}, {4.0, 1.0}};
  std::vector<V::SeriesPoint> model = reference;
  for (V::SeriesPoint& point : model) point.value *= 7.0;
  V::ComparisonPolicy policy;
  policy.normalization = V::Normalization::OneGlobalAmplitude;
  policy.maximumLog10Rmse = 1.0e-13;
  policy.maximumMedianAbsoluteLog10Error = 1.0e-13;
  policy.minimumLog10Correlation = 1.0 - 1.0e-13;
  policy.maximumOnsetCoordinateError = 0.0;
  policy.maximumPeakCoordinateError = 0.0;
  policy.maximumAbsoluteLog10PeakRatio = 1.0e-13;
  policy.maximumAbsoluteLog10FluenceRatio = 1.0e-13;
  const V::ComparisonMetrics metrics =
      V::ComparePositiveSeries(model, reference, policy);
  if (!metrics.status.ok() || !metrics.accepted || metrics.coverage != 1.0 ||
      std::fabs(metrics.amplitudeScale - 1.0 / 7.0) > 1.0e-14)
    return Fail("one global amplitude did not recover an identical curve shape");
  Result result = Pass(
      "log-space profile metrics recover exact shape under one declared amplitude");
  result.metrics.push_back(Metric{"log10_rmse", metrics.log10Rmse,
                                  policy.maximumLog10Rmse, "<=", "dex"});
  result.metrics.push_back(Metric{"coverage", metrics.coverage,
                                  policy.minimumCoverage, ">=", "fraction"});
  return result;
}

Result RunVFY3D02() {
  const std::vector<V::SeriesPoint> reference = {
      {0.0, 1.0}, {1.0, 2.0}, {2.0, 4.0}, {3.0, 2.0}, {4.0, 1.0}};
  const std::vector<V::SeriesPoint> shortModel = {
      {1.0, 2.0}, {2.0, 4.0}, {3.0, 2.0}};
  V::ComparisonPolicy policy;
  policy.minimumCoverage = 0.8;
  const V::ComparisonMetrics incomplete =
      V::ComparePositiveSeries(shortModel, reference, policy);
  if (!incomplete.status.ok() || incomplete.accepted ||
      std::fabs(incomplete.coverage - 0.6) > 1.0e-14 ||
      incomplete.violations.empty())
    return Fail("incomplete temporal coverage was accepted or hidden");

  std::vector<V::SeriesPoint> malformed = reference;
  malformed[2].coordinate = malformed[1].coordinate;
  if (V::ComparePositiveSeries(malformed, reference, policy).status.ok())
    return Fail("unordered scientific evidence received numerical metrics");
  return Pass("coverage loss is a scientific rejection and malformed curves are typed errors");
}

Result RunVFY3D03() {
  // For constant U and kappa, the field-aligned Green function is Gaussian:
  // mean=U*t and sigma=sqrt(2*kappa*t).  The 3-D kernel is projected onto
  // three unrelated field directions and checked with empirical CDF values;
  // this tests an entire distribution rather than only its second moment.
  constexpr std::size_t count = 40000;
  const double drift = 1.5, kappa = 2.0, time = 3.0;
  const double mean = drift * time;
  const double sigma = std::sqrt(2.0 * kappa * time);
  const std::array<C::Vec3, 3> directions = {
      C::Vec3(1.0, 0.0, 0.0),
      C::Vec3(1.0, 2.0, -1.0).Normalized(),
      C::Vec3(-2.0, 1.0, 3.0).Normalized()};
  double maximumCdfError = 0.0;
  for (const C::Vec3& direction : directions) {
    std::vector<double> projected;
    projected.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
      T::ParkerParticleState state;
      state.momentumKgMPerS = 1.0;
      T::ParkerLocalState local;
      local.bHat = direction;
      local.bulkVelocityMPerS = direction * drift;
      local.kappaParallelM2PerS = kappa;
      T::KeyedRandomStream random(
          Key(i, T::RandomPurpose::ParkerParallel));
      const T::ParkerStepResult moved =
          T::AdvanceParker(state, local, time, &random);
      if (!moved.status.ok()) return Fail("Parker validation trajectory failed");
      projected.push_back(moved.state.positionM.Dot(direction));
    }
    std::sort(projected.begin(), projected.end());
    for (double z : {-1.0, 0.0, 1.0}) {
      const double threshold = mean + z * sigma;
      const double empirical = static_cast<double>(
          std::upper_bound(projected.begin(), projected.end(), threshold) -
          projected.begin()) / count;
      const double analytic = 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
      maximumCdfError = std::max(maximumCdfError,
                                 std::fabs(empirical - analytic));
    }
  }
  if (maximumCdfError > 0.012)
    return Fail("projected Parker propagator differs from its Gaussian Green function");
  Result result = Pass(
      "3-D Parker projections match the independent 1-D Gaussian propagator");
  result.metrics.push_back(Metric{"maximum_cdf_error", maximumCdfError,
                                  0.012, "<=", "fraction"});
  return result;
}

double FocusedError(double dt) {
  const double mass = C::Const::m_p;
  const double momentum = mass;
  const double speed = T::RelativisticSpeed(momentum, mass);
  const double total = 2.0;
  T::FocusedParticleState state;
  state.momentumKgMPerS = momentum;
  state.mu = 0.2;
  T::FocusedLocalState local;
  local.bHat = C::Vec3(1.0, 0.0, 0.0);
  local.divBhatPerM = 0.2 / speed;
  const int steps = static_cast<int>(std::llround(total / dt));
  for (int i = 0; i < steps; ++i) {
    const T::FocusedStepResult moved =
        T::AdvanceFocused(state, local, mass, dt, nullptr);
    if (!moved.status.ok()) return std::numeric_limits<double>::infinity();
    state = moved.state;
  }
  // dmu/dt=0.1(1-mu^2) has the exact characteristic below.
  const double exact = std::tanh(std::atanh(0.2) + 0.1 * total);
  return std::fabs(state.mu - exact);
}

Result RunVFY3D04() {
  const double coarse = FocusedError(0.2);
  const double medium = FocusedError(0.1);
  const double fine = FocusedError(0.05);
  if (!std::isfinite(coarse) || !std::isfinite(medium) ||
      !std::isfinite(fine) || !(coarse > medium && medium > fine))
    return Fail("focused characteristic did not converge monotonically");
  const double order = std::log(medium / fine) / std::log(2.0);
  if (order < 1.8)
    return Fail("focused split failed its second-order characteristic limit");
  Result result = Pass(
      "focused transport converges at second order to the analytic focusing characteristic");
  result.metrics.push_back(Metric{"observed_order", order, 1.8, ">=", "order"});
  result.metrics.push_back(Metric{"fine_mu_error", fine, 1.0e-5, "<=", "absolute"});
  return result;
}

swcme::sep::SEPSourceState Source() {
  swcme::sep::SEPSourceState source;
  source.status = swcme::ModelStatus::success();
  source.active = true;
  source.source_id = 91;
  source.position_m = {{10.0, 0.0, 0.0}};
  source.normal = {{1.0, 0.0, 0.0}};
  source.relative_patch_weight = 0.75;
  source.compression = 4.0;
  source.normal_speed_m_s = 8.0e5;
  source.q_phase_space = 4.0;
  source.spectrum.particle_mass_kg = C::Const::m_p;
  source.spectrum.kinetic_energy_min_MeV = 1.0;
  source.spectrum.kinetic_energy_max_MeV = 100.0;
  source.spectrum.reference_energy_MeV = 10.0;
  return source;
}

Result RunVFY3D05() {
  constexpr std::uint64_t count = 20000;
  const A::ShockSourceRecord source =
      A::MakeShockSourceRecord(Source(), 17, 431, count, 0.2);
  if (!source.status.ok()) return Fail(source.status.message);
  std::vector<double> samples;
  samples.reserve(count);
  double representedWeight = 0.0;
  for (std::uint64_t i = 0; i < count; ++i) {
    const A::InjectedParticle particle =
        A::SampleInjectedParticle(source, i, 0);
    if (!particle.status.ok()) return Fail(particle.status.message);
    samples.push_back(particle.particle.momentumKgMPerS);
    representedWeight += particle.particle.statisticalWeight;
  }
  std::sort(samples.begin(), samples.end());
  double maximumCdfError = 0.0;
  for (double probability : {0.1, 0.25, 0.5, 0.75, 0.9}) {
    const SEP::Transport::ScalarResult quantile =
        SEP::Injection::InverseCdf(source.injection.spectrum, probability);
    if (!quantile.status.ok()) return Fail(quantile.status.message);
    const double empirical = static_cast<double>(
        std::upper_bound(samples.begin(), samples.end(), quantile.value) -
        samples.begin()) / count;
    maximumCdfError = std::max(maximumCdfError,
                               std::fabs(empirical - probability));
  }
  const double expectedWeight =
      source.relativePatchWeight * source.injection.injectionEfficiency;
  if (maximumCdfError > 0.015 ||
      std::fabs(representedWeight - expectedWeight) > 1.0e-12)
    return Fail("SWCME source spectrum or normalization misses its declared law");
  Result result = Pass(
      "SWCME injection follows the independent truncated DSA CDF and closes weight");
  result.metrics.push_back(Metric{"maximum_cdf_error", maximumCdfError,
                                  0.015, "<=", "fraction"});
  return result;
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterValidationTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* group, const char* name,
                 SEP3D::Testing::TestCallback callback) {
    D descriptor;
    descriptor.id = id;
    descriptor.name = name;
    descriptor.group = group;
    descriptor.description =
        "Phase-V controlled integration/scientific validation prerequisite";
    descriptor.initialization = I::None;
    descriptor.supportedBuildModes = "standalone-no-AMPS";
    descriptor.runtime = RC::Routine;
    descriptor.seedPolicy = "fixed campaign and particle-keyed streams";
    descriptor.stateIsolation =
        "fresh rank records, curves, and particle ensembles per callback";
    descriptor.callback = std::move(callback);
    return descriptor;
  };
  return {
      make("INT3D01", "INT3D", "Deterministic rank gather", RunINT3D01),
      make("INT3D02", "INT3D", "Global conservation audit", RunINT3D02),
      make("INT3D03", "INT3D", "Load and resource budgets", RunINT3D03),
      make("VFY3D01", "VFY3D", "Profile metric normalization", RunVFY3D01),
      make("VFY3D02", "VFY3D", "Coverage and malformed evidence", RunVFY3D02),
      make("VFY3D03", "VFY3D", "Parker Green-function validation", RunVFY3D03),
      make("VFY3D04", "VFY3D", "Focused characteristic convergence", RunVFY3D04),
      make("VFY3D05", "VFY3D", "SWCME DSA source distribution", RunVFY3D05),
  };
}
