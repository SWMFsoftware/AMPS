#include "sep_rng_validation.h"

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include "sep_parker_core.h"
#include "sep_focused_transport_core.h"
#include "sep_focused_transport_mfp_core.h"
#include "sep_source_lifecycle.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace SEP {
namespace Testing {
namespace {

namespace T = Transport;
const double kC = 299792458.0;
const double kMass = 1.67262192369e-27;

Result Complete(bool ok, const char* message, std::uint64_t seedDomain=0) {
  Result result;
  result.status = ok ? Status::Pass : Status::Fail;
  result.message = message;
  result.hasSeed = true;
  result.seed = TestSeed(seedDomain);
  result.configuration.push_back("evidence_level=analytical-core");
  result.metrics.push_back({"assertion_failures", ok ? 0.0 : 1.0,
                            0.0, "<=", "count"});
  return result;
}

bool Same(const T::KeyedRandomStreamState& a,
          const T::KeyedRandomStreamState& b) {
  return a.schema == b.schema && a.keyIdentity == b.keyIdentity &&
      a.integerState == b.integerState &&
      a.hasSpareNormal == b.hasSpareNormal && a.spareNormal == b.spareNormal;
}

T::OperatorRandomState Initial(std::uint64_t particle, std::uint64_t operation) {
  T::KeyedRandomStream random(TestSeed(0), particle, operation, 0);
  T::OperatorRandomState state;
  state.stream = random.SaveState();
  return state;
}

class Spatial final : public T::SpatialDiffusionProvider {
 public:
  T::SpatialDiffusionSample Evaluate(double, double) const override {
    T::SpatialDiffusionSample sample;
    sample.status = T::Status::Ok();
    sample.kappaParallelM2PerS = 1.0;
    sample.provenance = "rng-lifecycle:constant-kappa-v1";
    return sample;
  }
};

class Pitch final : public T::PitchAngleDiffusionProvider {
 public:
  T::PitchAngleDiffusionSample Evaluate(double, double, double mu) const override {
    T::PitchAngleDiffusionSample sample;
    sample.status = T::Status::Ok();
    sample.dMuMuPerS = 0.2 * (1.0 - mu * mu);
    sample.dDmuMuDmuPerS = -0.4 * mu;
    sample.provenance = "rng-lifecycle:Legendre-v1";
    sample.turbulenceStateIdentity = "rng-lifecycle:frozen:g1";
    return sample;
  }
};

class Mfp final : public T::MeanFreePathProvider {
 public:
  T::MeanFreePathSample Evaluate(double, double, double) const override {
    T::MeanFreePathSample sample;
    sample.status = T::Status::Ok();
    sample.lambdaParallelM = 1.0e7;
    sample.hasBranchResolvedRates = true;
    sample.nuPlusPerS = sample.nuMinusPerS = 0.5;
    sample.provenance = "rng-lifecycle:unit-rate-v1";
    sample.turbulenceStateIdentity = "rng-lifecycle:frozen:g1";
    return sample;
  }
};

Result RoundTrip() {
  bool ok = true;
  double firstNormal=0.0;
  for (unsigned prefix = 0; prefix < 8; ++prefix) {
    T::KeyedRandomStream uninterrupted(TestSeed(0), 17, 8, 0);
    for (unsigned i = 0; i < prefix; ++i) uninterrupted.Normal01();
    const T::KeyedRandomStreamState saved = uninterrupted.SaveState();
    ok = ok && saved.hasSpareNormal == prefix % 2;
    T::KeyedRandomStream resumed(TestSeed(0), 17, 8, 0);
    ok = ok && resumed.RestoreState(saved).ok();
    for (unsigned i = 0; i < 100; ++i) {
      const double actual=uninterrupted.Normal01(),continued=resumed.Normal01();
      if(prefix==0&&i==0)firstNormal=actual;
      ok = ok && actual == continued;
      if (i % 3 == 0)
        ok = ok && uninterrupted.UniformOpen01() == resumed.UniformOpen01();
    }
    ok = ok && Same(uninterrupted.SaveState(), resumed.SaveState());
  }
  Result result = Complete(ok,
      "mixed uniform/normal continuation preserves the complete Box-Muller state");
  result.configuration.push_back("normal_prefixes=0..7;continued_draws=100");
  result.metrics.push_back({"actual_seeded_observable",firstNormal,0.0,"observed","dimensionless"});
  return result;
}

Result ParkerCalls() {
  const unsigned n = 20000;
  Spatial provider;
  long double sum = 0, square = 0, negativeSum = 0, negativeSquare = 0;
  long double lagProduct = 0;
  bool ok = true;
  for (unsigned particle = 1; particle <= n; ++particle) {
    T::OperatorRandomState saved = Initial(particle, 7);
    T::ParkerState state(0.0, 1.0), negative(0.0, 1.0);
    double previousIncrement = 0.0;
    for (unsigned call = 0; call < 8; ++call) {
      T::KeyedRandomStream resumed(TestSeed(0), particle, 7, 0);
      std::uint64_t accepted = 0;
      ok = ok && T::ResumeOperatorRandomStream(saved, &resumed, &accepted).ok();
      const T::ParkerIncrement step = T::AdvanceParker(
          state, T::ParkerBackground(), 1.0, 1.0, provider, resumed);
      if (!step.status.ok()) return Complete(false, step.status.message.c_str());
      ok = ok && T::StageOperatorRandomStream(&saved, resumed, accepted + 1).ok();
      if (call == 1) lagProduct += previousIncrement * step.displacementM;
      previousIncrement = step.displacementM;
      state = step.state;
      // Deliberately reproduce the old adapter's per-call restart. This
      // negative control must violate the independent Brownian variance law.
      T::KeyedRandomStream restarted(TestSeed(0), particle, 7, 0);
      const T::ParkerIncrement wrong = T::AdvanceParker(
          negative, T::ParkerBackground(), 1.0, 1.0, provider, restarted);
      if (!wrong.status.ok()) return Complete(false, wrong.status.message.c_str());
      negative = wrong.state;
    }
    sum += state.arcLengthM;
    square += state.arcLengthM * state.arcLengthM;
    negativeSum += negative.arcLengthM;
    negativeSquare += negative.arcLengthM * negative.arcLengthM;
    ok = ok && saved.acceptedIntervals == 8;
  }
  const double variance = static_cast<double>((square - sum * sum / n) / (n - 1));
  const double wrongVariance = static_cast<double>(
      (negativeSquare - negativeSum * negativeSum / n) / (n - 1));
  // Gaussian sample variance has relative standard deviation sqrt(2/(N-1)).
  // The 5-sigma bound and lag-product 5-sigma bound are fixed analytically
  // before observing this seed; no tolerance is fitted to the repaired output.
  const double expected = 16.0;
  const double varianceLimit = 5.0 * expected * std::sqrt(2.0 / (n - 1));
  const double lag = static_cast<double>(lagProduct / n);
  const double lagLimit = 5.0 * 2.0 / std::sqrt(static_cast<double>(n));
  ok = ok && std::fabs(variance - expected) <= varianceLimit &&
      std::fabs(lag) <= lagLimit &&
      std::fabs(wrongVariance - expected) > 8.0 * varianceLimit;
  Result result = Complete(ok,
      "repeated-call Parker variance is 2*kappa*T and restart control is rejected");
  result.configuration.push_back("particles=20000;calls=8;dt_s=1;kappa_m2_per_s=1");
  result.metrics.push_back({"variance_absolute_error", std::fabs(variance - expected),
                            varianceLimit, "<=", "m^2"});
  result.metrics.push_back({"lag_product_absolute", std::fabs(lag), lagLimit,
                            "<=", "m^2"});
  result.metrics.push_back({"negative_control_variance_error",
      std::fabs(wrongVariance - expected), 8.0 * varianceLimit, ">", "m^2"});
  return result;
}

Result DmumuCalls() {
  const unsigned n = 12000;
  const unsigned steps = 100;
  const double dt = 0.005;
  const double amplitude = 0.8;
  const T::ScalarResult momentum = T::MomentumFromSpeed(1.0e7, kMass, kC);
  if (!momentum.status.ok()) return Complete(false, momentum.status.message.c_str());
  Pitch provider;
  long double sum = 0.0, square = 0.0;
  bool ok = true;
  for (unsigned particle = 1; particle <= n; ++particle) {
    T::KeyedRandomStream initial(TestSeed(0), particle, 108, 0);
    double mu = 0.0;
    for (;;) {
      mu = 2.0 * initial.UniformOpen01() - 1.0;
      const double p2 = 0.5 * (3.0 * mu * mu - 1.0);
      if (initial.UniformOpen01() * (1.0 + amplitude) <= 1.0 + amplitude * p2) break;
    }
    T::OperatorRandomState saved = Initial(particle, 8);
    T::FocusedTransportState state(0.0, momentum.value, mu), uninterrupted = state;
    T::KeyedRandomStream direct(TestSeed(0), particle, 8, 0);
    for (unsigned call = 0; call < steps; ++call) {
      T::KeyedRandomStream resumed(TestSeed(0), particle, 8, 0);
      std::uint64_t accepted = 0;
      ok = ok && T::ResumeOperatorRandomStream(saved, &resumed, &accepted).ok();
      const T::FocusedTransportIncrement step = T::AdvanceFocusedTransportDmumu(
          state, T::FocusedTransportBackground(), kMass, kC, dt,
          provider, resumed, NULL);
      const T::FocusedTransportIncrement directStep = T::AdvanceFocusedTransportDmumu(
          uninterrupted, T::FocusedTransportBackground(), kMass, kC, dt,
          provider, direct, NULL);
      if (!step.status.ok() || !directStep.status.ok())
        return Complete(false, "focused lifecycle step failed");
      ok = ok && T::StageOperatorRandomStream(&saved, resumed, accepted + 1).ok();
      state = step.state;
      uninterrupted = directStep.state;
      ok = ok && state.mu == uninterrupted.mu &&
          state.arcLengthM == uninterrupted.arcLengthM &&
          state.momentumKgMPerS == uninterrupted.momentumKgMPerS;
    }
    const double p2 = 0.5 * (3.0 * state.mu * state.mu - 1.0);
    sum += p2;
    square += p2 * p2;
    ok = ok && saved.acceptedIntervals == steps;
  }
  const double mean = static_cast<double>(sum / n);
  const double standardError = std::sqrt(static_cast<double>(
      (square - sum * sum / n) / (n - 1) / n));
  // Independent Legendre orthogonality/eigenvalue reference, matching the
  // declared weak-diffusion benchmark rather than a production-generated PDF.
  const double expected = amplitude / 5.0 * std::exp(-6.0 * 0.2 * steps * dt);
  ok = ok && standardError > 0.0 &&
      std::fabs(mean - expected) <= 5.0 * standardError;
  Result result = Complete(ok,
      "repeated-call Dmumu preserves exact stream continuation and Legendre decay");
  result.configuration.push_back("particles=12000;steps=100;dt_s=0.005;D0_per_s=0.2");
  result.metrics.push_back({"P2_moment_z_score", std::fabs(mean - expected) / standardError,
                            5.0, "<=", "sigma"});
  return result;
}

Result MfpCalls() {
  Mfp provider;
  const T::ScalarResult momentum = T::MomentumFromSpeed(1.0e7, kMass, kC);
  if (!momentum.status.ok()) return Complete(false, momentum.status.message.c_str());
  const T::FocusedTransportMfpState initial(0.0, momentum.value, 0.2);
  const T::FocusedTransportMfpBackground background;
  T::KeyedRandomStream direct(TestSeed(0), 123, 9, 0);
  T::ThreadLocalWaveAccumulator directEvents;
  const T::FocusedTransportMfpIncrement whole = T::AdvanceFocusedTransportMfp(
      initial, background, kMass, kC, 10.0, 1.0, provider, direct, &directEvents);
  if (!whole.status.ok()) return Complete(false, whole.status.message.c_str());
  T::FocusedTransportMfpState state = initial, wrong = initial;
  T::OperatorRandomState saved = Initial(123, 9);
  T::ThreadLocalWaveAccumulator partitionedEvents;
  bool ok = true;
  T::KeyedRandomStreamState restartedFinal;
  for (unsigned call = 0; call < 10; ++call) {
    T::KeyedRandomStream resumed(TestSeed(0), 123, 9, 0);
    std::uint64_t accepted = 0;
    ok = ok && T::ResumeOperatorRandomStream(saved, &resumed, &accepted).ok();
    const T::FocusedTransportMfpIncrement step = T::AdvanceFocusedTransportMfp(
        state, background, kMass, kC, 1.0, 1.0, provider, resumed, &partitionedEvents);
    if (!step.status.ok()) return Complete(false, step.status.message.c_str());
    ok = ok && T::StageOperatorRandomStream(&saved, resumed,
        accepted + step.diagnostics.deterministicIntervals).ok();
    state = step.state;
    T::KeyedRandomStream restarted(TestSeed(0), 123, 9, 0);
    const T::FocusedTransportMfpIncrement negative = T::AdvanceFocusedTransportMfp(
        wrong, background, kMass, kC, 1.0, 1.0, provider, restarted, NULL);
    if (!negative.status.ok()) return Complete(false, negative.status.message.c_str());
    wrong = negative.state;
    restartedFinal = restarted.SaveState();
  }
  std::vector<std::pair<std::uint64_t, int> > events, partitioned;
  for (const T::WaveContribution& event : directEvents.Contributions())
    if (event.scatteringEventAtEnd)
      events.push_back(std::make_pair(event.eventIndex, event.resonantBranch));
  for (const T::WaveContribution& event : partitionedEvents.Contributions())
    if (event.scatteringEventAtEnd)
      partitioned.push_back(std::make_pair(event.eventIndex, event.resonantBranch));
  const double epsilon = std::numeric_limits<double>::epsilon();
  const double positionLimit = 128.0 * epsilon * 1.0e8;
  const double hazardLimit = 128.0 * epsilon * 20.0;
  ok = ok && events == partitioned &&
      state.nextEventIndex == whole.state.nextEventIndex &&
      Same(saved.stream, direct.SaveState()) &&
      std::fabs(state.mu - whole.state.mu) <= 128.0 * epsilon &&
      std::fabs(state.arcLengthM - whole.state.arcLengthM) <= positionLimit &&
      std::fabs(state.remainingOpticalDepth - whole.state.remainingOpticalDepth) <= hazardLimit &&
      !Same(restartedFinal, direct.SaveState());
  Result result = Complete(ok,
      "MFP outer partitions preserve draws/events and reject restarted-stream history");
  result.configuration.push_back("particle=123;nu_per_s=1;duration_s=10;outer_partitions=10");
  result.metrics.push_back({"position_partition_error",
      std::fabs(state.arcLengthM - whole.state.arcLengthM), positionLimit, "<=", "m"});
  result.metrics.push_back({"optical_depth_partition_error",
      std::fabs(state.remainingOpticalDepth - whole.state.remainingOpticalDepth),
      hazardLimit, "<=", "dimensionless"});
  // Event counts can coincide for a different valid campaign seed. The
  // negative control is the lost draw history itself, not that accidental
  // observable. The full positive event/path/hazard checks above are retained.
  result.metrics.push_back({"negative_control_random_history_mismatch",
      Same(restartedFinal,direct.SaveState())?0.0:1.0,1.0,"==","count"});
  result.metrics.push_back({"actual_seeded_observable",whole.state.remainingOpticalDepth,
      0.0,"observed","dimensionless"});
  return result;
}

Result StateValidation() {
  T::KeyedRandomStream random(TestSeed(0), 19, 8, 0);
  random.Normal01();
  T::OperatorRandomState saved = Initial(19, 8);
  bool ok = T::StageOperatorRandomStream(&saved, random, 37).ok();
  // AMPS particle migration/checkpoint owns raw requested bytes. Exercise that
  // byte boundary without claiming this fixture is the native MPI transport.
  unsigned char bytes[sizeof(T::OperatorRandomState)];
  std::memcpy(bytes, &saved, sizeof(saved));
  T::OperatorRandomState restored;
  std::memcpy(&restored, bytes, sizeof(restored));
  T::KeyedRandomStream resumed(TestSeed(0), 19, 8, 0);
  std::uint64_t accepted = 0;
  ok = ok && T::ResumeOperatorRandomStream(restored, &resumed, &accepted).ok() && accepted == 37;
  for (unsigned i = 0; i < 20; ++i) ok = ok && random.Normal01() == resumed.Normal01();
  const T::KeyedRandomStreamState before = resumed.SaveState();
  T::KeyedRandomStreamState invalid = before;
  invalid.schema = 0;
  ok = ok && !resumed.RestoreState(invalid).ok() && Same(before, resumed.SaveState());
  invalid = before;
  invalid.keyIdentity ^= 1;
  ok = ok && !resumed.RestoreState(invalid).ok() && Same(before, resumed.SaveState());
  invalid = before;
  invalid.hasSpareNormal = 2;
  ok = ok && !resumed.RestoreState(invalid).ok() && Same(before, resumed.SaveState());
  invalid = before;
  invalid.spareNormal = std::numeric_limits<double>::quiet_NaN();
  ok = ok && !resumed.RestoreState(invalid).ok() && Same(before, resumed.SaveState());
  const T::OperatorRandomState beforeRollback = saved;
  ok = ok && !T::StageOperatorRandomStream(&saved, resumed, 36).ok() &&
      Same(beforeRollback.stream, saved.stream) && saved.acceptedIntervals == 37;
  return Complete(ok,
      "serialized cache resumes exactly and invalid/key-mismatched/rollback state is rejected atomically");
}

Result Scheduling() {
  const unsigned n = 64;
  std::vector<T::OperatorRandomState> states;
  std::vector<double> reference(n, 0.0), permuted(n, 0.0);
  for (unsigned i = 0; i < n; ++i) states.push_back(Initial(i + 1, 8));
  bool ok = true;
  for (unsigned i = 0; i < n; ++i) {
    T::KeyedRandomStream uninterrupted(TestSeed(0), i + 1, 8, 0);
    for (unsigned j = 0; j < 7; ++j) reference[i] += uninterrupted.Normal01();
  }
  for (unsigned call = 0; call < 7; ++call) {
    for (unsigned offset = 0; offset < n; ++offset) {
      const unsigned i = (offset * 17 + call * 11) % n;
      T::KeyedRandomStream resumed(TestSeed(0), i + 1, 8, 0);
      std::uint64_t accepted = 0;
      ok = ok && T::ResumeOperatorRandomStream(states[i], &resumed, &accepted).ok();
      permuted[i] += resumed.Normal01();
      ok = ok && T::StageOperatorRandomStream(&states[i], resumed, accepted + 1).ok();
    }
  }
  ok = ok && reference == permuted;
  T::KeyedRandomStream otherOperator(TestSeed(0), 1, 7, 0);
  const T::KeyedRandomStreamState before = otherOperator.SaveState();
  ok = ok && !otherOperator.RestoreState(states[0].stream).ok() &&
      Same(before, otherOperator.SaveState());
  Result result = Complete(ok,
      "particle-order permutation preserves stream identity and operator streams cannot cross-restore");
  result.configuration.push_back("particles=64;calls=7;permutation_multiplier=17");
  return result;
}

Result EventIdentityOverflow() {
  Mfp provider;
  const T::ScalarResult momentum = T::MomentumFromSpeed(1.0e7, kMass, kC);
  if (!momentum.status.ok()) return Complete(false, momentum.status.message.c_str());
  T::FocusedTransportMfpState initial(0.0, momentum.value, 0.2);
  initial.remainingOpticalDepth = 0.01;
  initial.nextEventIndex = std::numeric_limits<std::uint64_t>::max();
  T::KeyedRandomStream random(TestSeed(0), 123, 9, 0);
  const T::KeyedRandomStreamState before = random.SaveState();
  T::ThreadLocalWaveAccumulator pending;
  const T::FocusedTransportMfpIncrement step = T::AdvanceFocusedTransportMfp(
      initial, T::FocusedTransportMfpBackground(), kMass, kC,
      0.1, 0.1, provider, random, &pending);
  const bool ok = step.status.code == T::StatusCode::StepUnderflow &&
      step.state.nextEventIndex == initial.nextEventIndex &&
      step.state.arcLengthM == initial.arcLengthM &&
      step.state.momentumKgMPerS == initial.momentumKgMPerS &&
      step.state.mu == initial.mu &&
      step.state.remainingOpticalDepth == initial.remainingOpticalDepth &&
      Same(before, random.SaveState()) && pending.Contributions().empty();
  return Complete(ok,
      "exhausted scattering identity rejects before phase-space, RNG, or queue mutation");
}

Result OpenIntervalEndpoints() {
  double output = 0.375;
  bool ok = !T::UniformOpen01FromBits(0, &output) && output == 0.375;
  ok = ok && !T::UniformOpen01FromBits(UINT64_C(2047), &output) && output == 0.375;
  ok = ok && !T::UniformOpen01FromBits(UINT64_C(2048), NULL);
  ok = ok && T::UniformOpen01FromBits(UINT64_C(2048), &output) &&
      output == std::ldexp(1.0, -53) && output > 0.0;
  ok = ok && T::UniformOpen01FromBits(
      std::numeric_limits<std::uint64_t>::max(), &output) &&
      output == 1.0 - std::ldexp(1.0, -53) && output < 1.0;
  // Check each bin-edge bit position and its least-significant bits against an
  // independent power-of-two reference. Low 11 bits must not affect mapping.
  for (unsigned bit = 11; bit < 64; ++bit) {
    const std::uint64_t word = UINT64_C(1) << bit;
    ok = ok && T::UniformOpen01FromBits(word, &output) &&
        output == std::ldexp(1.0, static_cast<int>(bit) - 64);
    const double expected = output;
    ok = ok && T::UniformOpen01FromBits(word | UINT64_C(2047), &output) &&
        output == expected;
  }
  T::KeyedRandomStream random(TestSeed(0), 1, 9, 0);
  const T::KeyedRandomStreamState before = random.SaveState();
  T::KeyedRandomStreamState old = before;
  old.schema = 1;
  ok = ok && !random.RestoreState(old).ok() && Same(before, random.SaveState());
  Result result = Complete(ok,
      "bit-extreme uniform mapping is strictly open and old algorithm state is rejected");
  result.configuration.push_back("algorithm_schema=2;grid=integer_m*2^-53;reject_m=0");
  return result;
}

Result SourceEventIdentity() {
  namespace S=SourceLifecycle;
  S::Scheduler scheduler;
  bool ok=scheduler.Initialize(TestSeed(0),"ordered-fixed-geometry-v1").ok();
  std::set<std::uint64_t> ids;
  // The old line/macro/momentum key was identical for these independent,
  // constant-energy, field-aligned births. No phase-space value enters the new
  // key, and a zero-birth event still consumes its own physical event index.
  for(unsigned event=0;event<4;++event)for(unsigned kind=1;kind<=2;++kind)
    for(unsigned species=0;species<3;++species)for(unsigned line=0;line<3;++line) {
      S::Channel channel;channel.kind=static_cast<S::Kind>(kind);
      channel.line=line;channel.species=species;
      S::Event source;ok=ok&&scheduler.Begin(channel,&source).ok()&&source.index==event;
      for(unsigned macro=0;macro<16;++macro)
        ok=ok&&ids.insert(S::ParticleIdentity(source,macro)).second;
      auto energy=S::Random(source,0,2),direction=S::Random(source,0,3);
      ok=ok&&energy.UniformOpen01()!=direction.UniformOpen01();
    }
  Result result=Complete(ok,"source kind, physical line, species, event, and macro separate repeated fixed-state births");
  result.metrics.push_back({"unique_birth_keys",static_cast<double>(ids.size()),1152.0,"==","count"});
  return result;
}

Result SourceSchedulerCheckpoint() {
  namespace S=SourceLifecycle;
  S::Scheduler uninterrupted;
  bool ok=uninterrupted.Initialize(TestSeed(9),"fixed-layout-A").ok();
  S::Channel channel;channel.kind=S::Kind::FixedBeginning;channel.line=3;channel.species=2;
  S::Event event;
  S::Scheduler absent;event.index=321;
  ok=ok&&!absent.Begin(channel,&event).ok()&&event.index==321;
  for(unsigned i=0;i<7;++i)ok=ok&&uninterrupted.Begin(channel,&event).ok();
  FILE* file=std::tmpfile();
  if(!file)return Complete(false,"cannot allocate source checkpoint fixture");
  const S::State snapshot=uninterrupted.Capture();
  ok=ok&&S::WriteState(file,snapshot).ok();std::rewind(file);
  S::State loaded;ok=ok&&S::ReadState(file,&loaded).ok();std::fclose(file);
  S::Scheduler resumed;ok=ok&&resumed.Restore(loaded,TestSeed(9),"fixed-layout-A").ok();
  double firstDraw=0.0;
  for(unsigned i=0;i<11;++i) {
    S::Event a,b;ok=ok&&uninterrupted.Begin(channel,&a).ok()&&resumed.Begin(channel,&b).ok();
    ok=ok&&a.index==b.index&&S::ParticleIdentity(a,4)==S::ParticleIdentity(b,4);
    auto r0=S::Random(a,4,2),r1=S::Random(b,4,2);
    const double actual=r0.UniformOpen01(),continued=r1.UniformOpen01();
    if(i==0)firstDraw=actual;
    ok=ok&&actual==continued;
  }
  auto bad=loaded;bad.schema=0;
  ok=ok&&!resumed.Restore(bad,TestSeed(9),"fixed-layout-A").ok();
  ok=ok&&!resumed.Restore(loaded,TestSeed(9)^UINT64_C(1),"fixed-layout-A").ok();
  ok=ok&&!resumed.Restore(loaded,TestSeed(9),"reordered-layout-B").ok();
  bad=loaded;bad.counters.push_back(bad.counters.front());
  ok=ok&&!resumed.Restore(bad,TestSeed(9),"fixed-layout-A").ok();
  const auto before=resumed.Capture();
  ok=ok&&before.counters.front().next==18;
  // Overflow rejection preserves both scheduler and caller output.
  bad=loaded;bad.counters.front().next=std::numeric_limits<std::uint64_t>::max();
  ok=ok&&resumed.Restore(bad,TestSeed(9),"fixed-layout-A").ok();
  event.index=123;
  ok=ok&&!resumed.Begin(channel,&event).ok()&&event.index==123&&
      resumed.Capture().counters.front().next==std::numeric_limits<std::uint64_t>::max();
  // Failed/truncated reads must not publish a partial scheduler state.
  file=std::tmpfile();if(!file)return Complete(false,"cannot allocate truncated source fixture");
  const S::State untouched=loaded;
  ok=ok&&!S::ReadState(file,&loaded).ok()&&loaded.counters.front().next==untouched.counters.front().next;
  std::fclose(file);
  auto result=Complete(ok,"source scheduler serialization continues exact next events and rejects incompatible or missing history atomically",9);
  result.metrics.push_back({"actual_seeded_observable",firstDraw,0.0,"observed","dimensionless"});
  return result;
}

Descriptor DescriptorFor(const char* id, const char* name, TestCallback callback) {
  Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = "rng-lifecycle";
  descriptor.description = name;
  descriptor.initialization = InitializationLevel::None;
  descriptor.supportedBuildModes = "native CLI and dependency-light C++11 sanitizer runner";
  descriptor.runtime = RuntimeClass::Routine;
  descriptor.seedPolicy = "TestSeed schema1; defaults0 (source checkpoint9); stable particle/operator keys";
  const bool deterministic=descriptor.id=="RNG07"||descriptor.id=="RNG08";
  descriptor.seedCapability=deterministic?SeedCapability::Deterministic:SeedCapability::Context;
  descriptor.repeatSupported=true;
  descriptor.stateIsolation = "private operator state and pure production cores; no PIC/MPI claim";
  descriptor.callback = [callback,deterministic]() {
    auto result=callback();
    if(deterministic) {
      result.hasSeed=false;
      result.configuration.push_back("random_consumption=none;fixed-identity-validation-only");
    }
    return result;
  };
  return descriptor;
}

}  // namespace

std::vector<Descriptor> RngLifecycleDescriptors() {
  std::vector<Descriptor> descriptors;
  descriptors.push_back(DescriptorFor("RNG01", "Complete random-stream round trip", RoundTrip));
  descriptors.push_back(DescriptorFor("RNG02", "Repeated-call Parker diffusion", ParkerCalls));
  descriptors.push_back(DescriptorFor("RNG03", "Repeated-call Dmumu modal decay", DmumuCalls));
  descriptors.push_back(DescriptorFor("RNG04", "MFP event history across outer partitions", MfpCalls));
  descriptors.push_back(DescriptorFor("RNG05", "Particle RNG serialization and rejection", StateValidation));
  descriptors.push_back(DescriptorFor("RNG06", "Scheduler-independent operator continuation", Scheduling));
  descriptors.push_back(DescriptorFor("RNG07", "Scattering event identity overflow", EventIdentityOverflow));
  descriptors.push_back(DescriptorFor("RNG08", "Strict open-interval random endpoints", OpenIntervalEndpoints));
  descriptors.push_back(DescriptorFor("RNG09", "Persistent physical source event identities", SourceEventIdentity));
  descriptors.push_back(DescriptorFor("RNG10", "Source scheduler checkpoint lifecycle", SourceSchedulerCheckpoint));
  return descriptors;
}

}  // namespace Testing
}  // namespace SEP
