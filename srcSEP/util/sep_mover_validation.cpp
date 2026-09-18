#include "sep_mover_validation.h"

#include "sep_focused_transport_core.h"
#include "sep_focused_transport_mfp_core.h"
#include "sep_parker_core.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include "../QLT.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace SEP {
namespace Testing {
namespace {

const double kSpeedOfLightMPerS = 299792458.0;
const double kProtonMassKg = 1.67262192369e-27;

bool NearAbsolute(double actual, double expected, double tolerance) {
  return std::fabs(actual - expected) <= tolerance;
}

Result Complete(bool pass, const std::string& success,
                const std::string& failure) {
  Result result;
  result.status = pass ? Status::Pass : Status::Fail;
  result.message = pass ? success : failure;
  result.metrics.push_back(
      {"assertion_failures", pass ? 0.0 : 1.0, 0.0, "<=", "count"});
  return result;
}

Descriptor MakeDescriptor(const char* id, const char* name, const char* group,
                          const char* description, RuntimeClass runtime,
                          const char* seed_policy, TestCallback callback) {
  Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = group;
  descriptor.description = description;
  descriptor.initialization = InitializationLevel::None;
  descriptor.supportedBuildModes =
      "linked CLI and source-only C++11 ASan/UBSan runner";
  descriptor.runtime = runtime;
  descriptor.seedPolicy = seed_policy;
  descriptor.stateIsolation =
      "stack/vector-owned controlled state and injected coefficient providers";
  descriptor.callback = callback;
  return descriptor;
}

class ConstantSpatialProvider final
    : public Transport::SpatialDiffusionProvider {
 public:
  ConstantSpatialProvider(double kappa, double gradient)
      : kappa_(kappa), gradient_(gradient) {}

  Transport::SpatialDiffusionSample Evaluate(double, double) const override {
    Transport::SpatialDiffusionSample sample;
    sample.status = Transport::Status::Ok();
    sample.kappaParallelM2PerS = kappa_;
    sample.dKappaParallelDsMPerS = gradient_;
    sample.provenance = "controlled-test:constant-spatial-v1";
    return sample;
  }

 private:
  double kappa_;
  double gradient_;
};

class LinearDriftSpatialProvider final
    : public Transport::SpatialDiffusionProvider {
 public:
  explicit LinearDriftSpatialProvider(double rate_per_s)
      : rate_per_s_(rate_per_s) {}

  Transport::SpatialDiffusionSample Evaluate(double s_m, double) const override {
    Transport::SpatialDiffusionSample sample;
    sample.status = Transport::Status::Ok();
    sample.kappaParallelM2PerS = 0.0;
    sample.dKappaParallelDsMPerS = rate_per_s_ * s_m;
    sample.provenance = "controlled-test:linear-gradient-v1";
    return sample;
  }

 private:
  double rate_per_s_;
};

class ConstantPitchProvider final
    : public Transport::PitchAngleDiffusionProvider {
 public:
  ConstantPitchProvider(double diffusion, double derivative,
                        const std::string& identity =
                            "controlled-test:turbulence:g1")
      : diffusion_(diffusion), derivative_(derivative), identity_(identity) {}

  Transport::PitchAngleDiffusionSample Evaluate(double, double,
                                                  double) const override {
    Transport::PitchAngleDiffusionSample sample;
    sample.status = Transport::Status::Ok();
    sample.dMuMuPerS = diffusion_;
    sample.dDmuMuDmuPerS = derivative_;
    sample.provenance = "controlled-test:constant-dmumu-v1";
    sample.turbulenceStateIdentity = identity_;
    return sample;
  }

 private:
  double diffusion_;
  double derivative_;
  std::string identity_;
};

class LegendrePitchProvider final
    : public Transport::PitchAngleDiffusionProvider {
 public:
  explicit LegendrePitchProvider(double d0_per_s) : d0_per_s_(d0_per_s) {}

  Transport::PitchAngleDiffusionSample Evaluate(double, double,
                                                  double mu) const override {
    Transport::PitchAngleDiffusionSample sample;
    sample.status = Transport::Status::Ok();
    // D_mumu=D0(1-mu^2) makes each Legendre polynomial an eigenfunction
    // of the continuum pitch-angle operator.  The derivative is supplied
    // independently and explicitly because it is the Ito drift.
    sample.dMuMuPerS = d0_per_s_ * std::max(0.0, 1.0 - mu * mu);
    sample.dDmuMuDmuPerS = -2.0 * d0_per_s_ * mu;
    sample.provenance = "controlled-test:legendre-dmumu-v1";
    sample.turbulenceStateIdentity = "controlled-test:frozen-legendre:g1";
    return sample;
  }

 private:
  double d0_per_s_;
};

class ConstantMfpProvider final : public Transport::MeanFreePathProvider {
 public:
  explicit ConstantMfpProvider(double lambda_m) : lambda_m_(lambda_m) {}

  Transport::MeanFreePathSample Evaluate(double, double,
                                          double) const override {
    Transport::MeanFreePathSample sample;
    sample.status = Transport::Status::Ok();
    sample.lambdaParallelM = lambda_m_;
    sample.provenance = "controlled-test:constant-lambda-v1";
    sample.turbulenceStateIdentity = "controlled-test:turbulence:g1";
    return sample;
  }

 private:
  double lambda_m_;
};

Transport::FocusedTransportIncrement FocusedStep(
    const Transport::FocusedTransportState& state,
    const Transport::FocusedTransportBackground& background,
    double dt_s, const Transport::PitchAngleDiffusionProvider& provider,
    std::uint64_t particle,
    Transport::ThreadLocalWaveAccumulator* accumulator = NULL) {
  Transport::KeyedRandomStream random(31, particle, 8, 0);
  return Transport::AdvanceFocusedTransportDmumu(
      state, background, kProtonMassKg, kSpeedOfLightMPerS, dt_s,
      provider, random, accumulator);
}

Transport::FocusedTransportMfpIncrement AdvanceMfp(
    double speed_m_per_s, double mu, double dt_s, double maximum_interval_s,
    double lambda_m,
    const Transport::FocusedTransportMfpBackground& background,
    std::uint64_t particle_id) {
  const Transport::ScalarResult momentum = Transport::MomentumFromSpeed(
      speed_m_per_s, kProtonMassKg, kSpeedOfLightMPerS);
  ConstantMfpProvider provider(lambda_m);
  Transport::KeyedRandomStream random(91, particle_id, 9, 0);
  return Transport::AdvanceFocusedTransportMfp(
      {0.0, momentum.value, mu}, background, kProtonMassKg,
      kSpeedOfLightMPerS, dt_s, maximum_interval_s, provider, random, NULL);
}

Result RunPark01() {
  ConstantSpatialProvider provider(0.0, 0.0);
  Transport::KeyedRandomStream random(1, 1, 7, 0);
  const Transport::ParkerIncrement step = Transport::AdvanceParker(
      {2.0, 3.0e-19}, {400000.0, 0.0}, 1.0e7, 5.0, provider, random);
  const double expected_displacement_m = 2.0e6;
  const double error_m = std::fabs(step.displacementM - expected_displacement_m);
  Result result = Complete(step.status.ok() && error_m <= 1.0e-9,
      "zero-diffusion Parker state follows the exact convection characteristic",
      "Parker convection differs from U*dt");
  result.hasSeed = true;
  result.seed = 1;
  result.configuration = {"U_parallel_m_per_s=400000", "dt_s=5", "kappa=0"};
  result.metrics.push_back(
      {"displacement_absolute_error", error_m, 1.0e-9, "<=", "m"});
  return result;
}

Result RunPark02() {
  const std::uint64_t seed = 2;
  const int samples = 120000;
  const double kappa = 2.5e10;
  const double dt_s = 3.0;
  ConstantSpatialProvider provider(kappa, 0.0);
  long double sum = 0.0L;
  long double sum2 = 0.0L;
  bool status_ok = true;
  for (int i = 0; i < samples; ++i) {
    Transport::KeyedRandomStream random(seed, static_cast<unsigned>(i), 7, 0);
    const Transport::ParkerIncrement step = Transport::AdvanceParker(
        {0.0, 1.0e-19}, {0.0, 0.0}, 1.0e6, dt_s, provider, random);
    status_ok = status_ok && step.status.ok();
    sum += step.displacementM;
    sum2 += static_cast<long double>(step.displacementM) * step.displacementM;
  }
  const double mean = static_cast<double>(sum / samples);
  const double variance = static_cast<double>(sum2 / samples) - mean * mean;
  const double expected_variance = 2.0 * kappa * dt_s;
  const double mean_scaled = std::fabs(mean) / std::sqrt(expected_variance);
  const double variance_error =
      std::fabs(variance / expected_variance - 1.0);
  Result result = Complete(status_ok && mean_scaled < 0.015 &&
      variance_error < 0.015,
      "Gaussian Parker displacement matches its analytical moments",
      "Parker displacement moments exceed their statistical tolerances");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration = {"samples=120000", "kappa_m2_per_s=2.5e10", "dt_s=3"};
  result.metrics.push_back(
      {"scaled_absolute_mean", mean_scaled, 0.015, "<", "sigma"});
  result.metrics.push_back(
      {"variance_relative_error", variance_error, 0.015, "<", "dimensionless"});
  return result;
}

Result RunPark03() {
  ConstantSpatialProvider provider(0.0, 12.5);
  Transport::KeyedRandomStream random(3, 1, 7, 0);
  const Transport::ParkerIncrement step = Transport::AdvanceParker(
      {0.0, 1.0e-19}, {0.0, 0.0}, 1.0e6, 4.0, provider, random);
  const double error_m = std::fabs(step.displacementM - 50.0);
  Result result = Complete(step.status.ok() && error_m <= 1.0e-14,
      "Ito variable-kappa drift has the analytical positive gradient sign",
      "Ito variable-kappa drift has the wrong sign or magnitude");
  result.hasSeed = true;
  result.seed = 3;
  result.configuration = {"d_kappa_ds_m_per_s=12.5", "dt_s=4", "kappa=0"};
  result.metrics.push_back(
      {"drift_displacement_error", error_m, 1.0e-14, "<=", "m"});
  return result;
}

Result RunPark04() {
  ConstantSpatialProvider provider(0.0, 0.0);
  Transport::KeyedRandomStream random(4, 1, 7, 0);
  const double p0 = 2.0e-19;
  const Transport::ParkerIncrement step = Transport::AdvanceParker(
      {0.0, p0}, {0.0, 3.0e-4}, 1.0e6, 5.0, provider, random);
  const double expected = p0 * std::exp(-5.0e-4);
  const double error = std::fabs(step.state.momentumKgMPerS - expected);
  Result result = Complete(step.status.ok() && error <= 1.0e-32,
      "Parker momentum follows the exact adiabatic characteristic",
      "Parker momentum differs from the analytical cooling solution");
  result.hasSeed = true;
  result.seed = 4;
  result.configuration = {"p0_kg_m_per_s=2e-19", "divU_s^-1=3e-4", "dt_s=5"};
  result.metrics.push_back(
      {"momentum_absolute_error", error, 1.0e-32, "<=", "kg m s^-1"});
  return result;
}

Result RunPark05() {
  const Transport::CoordinateAdvance out = Transport::AdvanceCoordinate(
      0.9, 0.2, 1.0, 0.0, 1.0, Transport::BoundaryPolicy::Absorb);
  const bool pass = out.status.code == Transport::StatusCode::OutOfDomain &&
      out.crossedBoundary;
  Result result = Complete(pass,
      "outward Parker crossing reports the controlled absorbing boundary",
      "absorbing boundary failed to report the crossing");
  result.configuration = {"domain_m=[0,1]", "initial_s_m=0.9", "displacement_m=0.2"};
  result.metrics.push_back(
      {"crossing_detected", out.crossedBoundary ? 1.0 : 0.0, 1.0, "==", "boolean"});
  return result;
}

double IntegrateLinearParker(double dt_s) {
  LinearDriftSpatialProvider provider(0.4);
  Transport::KeyedRandomStream random(6, 1, 7, 0);
  Transport::ParkerState state = {1.0, 1.0e-19};
  for (double time_s = 0.0; time_s < 1.0 - 0.5 * dt_s; time_s += dt_s)
    state = Transport::AdvanceParker(
        state, {0.0, 0.0}, 1.0, dt_s, provider, random).state;
  return state.arcLengthM;
}

Result RunPark06() {
  const double exact = std::exp(0.4);
  const double coarse = std::fabs(IntegrateLinearParker(0.1) - exact);
  const double fine = std::fabs(IntegrateLinearParker(0.05) - exact);
  const RefinementOrderEstimate order =
      EstimateRefinementOrder(coarse, fine, 0.1, 0.05);
  // Forward Euler is first order for this manufactured drift problem.  A
  // modest lower bound permits pre-asymptotic error while expressing the
  // numerical requirement in the physically meaningful observed order.
  const double minimum_order = 0.85;
  Result result = Complete(order.valid && order.observedOrder >= minimum_order,
      "variable-kappa manufactured drift converges under refinement",
      "variable-kappa drift does not show the required observed order");
  result.hasSeed = true;
  result.seed = 6;
  result.configuration = {"ds_dt_s^-1=0.4*s", "duration_s=1", "dt_s=0.1,0.05"};
  result.metrics.push_back({"observed_refinement_order", order.observedOrder,
      minimum_order, ">=", "dimensionless"});
  result.metrics.push_back({"coarse_to_fine_error_ratio",
      order.coarseToFineErrorRatio, 1.0, ">", "dimensionless"});
  return result;
}

Result RunPark07() {
  const std::uint64_t seed = 407;
  const std::size_t particles = 20000;
  const double minimum_m = -1.0;
  const double maximum_m = 1.0;
  const double initial_m = 0.2;
  const double kappa_m2_per_s = 0.25;
  const double dt_s = 0.002;
  const double maximum_time_s = 30.0;
  ConstantSpatialProvider provider(kappa_m2_per_s, 0.0);

  std::size_t upper_exits = 0;
  std::size_t completed = 0;
  long double sum_time = 0.0L;
  bool status_ok = true;
  for (std::size_t particle = 0; particle < particles && status_ok; ++particle) {
    Transport::KeyedRandomStream random(seed, particle, 7, 4070);
    Transport::ParkerState state = {initial_m, 1.0e-19};
    for (double time_s = 0.0; time_s < maximum_time_s; time_s += dt_s) {
      const Transport::ParkerIncrement step = Transport::AdvanceParker(
          state, {0.0, 0.0}, 1.0, dt_s, provider, random);
      status_ok = status_ok && step.status.ok();
      if (!status_ok) break;

      // Boundary detection is performed by the same shared coordinate kernel
      // used by the PIC adapters.  In an Euler path simulation the crossing
      // is observable only at the end of a step, so dt is recorded and the
      // mean-time tolerance includes this controlled discretization error.
      const Transport::CoordinateAdvance boundary =
          Transport::AdvanceCoordinate(
              state.arcLengthM, step.displacementM,
              step.displacementM / dt_s, minimum_m, maximum_m,
              Transport::BoundaryPolicy::Absorb);
      if (boundary.status.code == Transport::StatusCode::OutOfDomain) {
        const double exit_time_s = time_s + dt_s;
        const double attempted_m = state.arcLengthM + step.displacementM;
        if (attempted_m > maximum_m) ++upper_exits;
        ++completed;
        sum_time += exit_time_s;
        break;
      }
      status_ok = status_ok && boundary.status.ok();
      state = step.state;
    }
  }

  // For drift-free Brownian motion dX=sqrt(2*kappa)dW on [a,b],
  // P(exit at b)=(x-a)/(b-a) and E[T]=(x-a)(b-x)/(2*kappa).
  // These references are evaluated independently of AdvanceParker.
  const double expected_upper =
      (initial_m - minimum_m) / (maximum_m - minimum_m);
  const double expected_mean_s =
      (initial_m - minimum_m) * (maximum_m - initial_m) /
      (2.0 * kappa_m2_per_s);
  const double actual_upper =
      static_cast<double>(upper_exits) / static_cast<double>(particles);
  const double actual_mean_s = completed == 0 ? 0.0
      : static_cast<double>(sum_time / completed);
  // For the symmetric [-L,L] interval the exact first-passage-time variance is
  // (L^4-x^4)/(6*kappa^2).  Using it rather than the simulated sample variance
  // keeps the statistical acceptance limit independent of production output.
  const double half_width_m = 0.5 * (maximum_m - minimum_m);
  const double centered_initial_m =
      initial_m - 0.5 * (minimum_m + maximum_m);
  const double expected_variance_s2 =
      (std::pow(half_width_m, 4) - std::pow(centered_initial_m, 4)) /
      (6.0 * kappa_m2_per_s * kappa_m2_per_s);
  const double probability_standard_error = std::sqrt(
      expected_upper * (1.0 - expected_upper) / particles);
  const double mean_standard_error_s =
      std::sqrt(expected_variance_s2 / particles);
  const double probability_tolerance =
      5.0 * probability_standard_error + 0.005;
  const double mean_tolerance_s = 5.0 * mean_standard_error_s + 0.06;
  const double probability_error = std::fabs(actual_upper - expected_upper);
  const double mean_error_s = std::fabs(actual_mean_s - expected_mean_s);
  const bool pass = status_ok && completed == particles &&
      probability_error <= probability_tolerance &&
      mean_error_s <= mean_tolerance_s;
  Result result = Complete(pass,
      "Parker absorbing first passage matches exact exit statistics",
      "Parker first-passage probability, mean time, or completion is outside tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration = {"domain_m=[-1,1]", "initial_s_m=0.2",
      "kappa_m2_per_s=0.25", "dt_s=0.002", "particles=20000",
      "maximum_time_s=30"};
  result.metrics.push_back({"upper_exit_probability_absolute_error",
      probability_error, probability_tolerance, "<=", "probability"});
  result.metrics.push_back({"mean_first_passage_time_absolute_error",
      mean_error_s, mean_tolerance_s, "<=", "s"});
  result.metrics.push_back({"completed_fraction",
      static_cast<double>(completed) / particles, 1.0, "==", "fraction"});
  return result;
}

Result RunFted01() {
  const std::uint64_t seed = 31;
  const int samples = 120000;
  const double diffusion = 0.2;
  const double dt_s = 0.01;
  ConstantPitchProvider provider(diffusion, 0.0);
  long double sum = 0.0L;
  long double sum2 = 0.0L;
  bool status_ok = true;
  for (int i = 0; i < samples; ++i) {
    const Transport::FocusedTransportIncrement step = FocusedStep(
        {0.0, 1.0e-19, 0.0}, {0.0, 0.0, 0.0, 0.0}, dt_s,
        provider, static_cast<unsigned>(i));
    status_ok = status_ok && step.status.ok();
    sum += step.state.mu;
    sum2 += step.state.mu * step.state.mu;
  }
  const double mean = static_cast<double>(sum / samples);
  const double variance = static_cast<double>(sum2 / samples) - mean * mean;
  const double expected_variance = 2.0 * diffusion * dt_s;
  const double mean_scaled = std::fabs(mean) / std::sqrt(expected_variance);
  const double variance_error = std::fabs(variance / expected_variance - 1.0);
  Result result = Complete(status_ok && mean_scaled < 0.015 &&
      variance_error < 0.015,
      "constant-Dmumu ensemble reproduces analytical Ito moments",
      "constant-Dmumu moments exceed their statistical tolerances");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration = {"samples=120000", "Dmumu_s^-1=0.2", "dt_s=0.01"};
  result.metrics.push_back({"scaled_absolute_mean", mean_scaled, 0.015, "<", "sigma"});
  result.metrics.push_back({"variance_relative_error", variance_error, 0.015, "<", "dimensionless"});
  return result;
}

Result RunFted02() {
  ConstantPitchProvider provider(0.0, 0.3);
  const Transport::FocusedTransportIncrement step = FocusedStep(
      {0.0, 1.0e-19, 0.2}, {0.0, 0.0, 0.0, 0.0}, 0.1, provider, 1);
  const double error = std::fabs(step.state.mu - 0.23);
  Result result = Complete(step.status.ok() && error <= 1.0e-13,
      "Ito pitch-angle drift includes the supplied dDmumu/dmu",
      "pitch-angle drift differs from mu0+dDmumu_dmu*dt");
  result.configuration = {"mu0=0.2", "dDmumu_dmu_s^-1=0.3", "dt_s=0.1"};
  result.metrics.push_back({"mu_absolute_error", error, 1.0e-13, "<=", "dimensionless"});
  return result;
}

Result RunFted03() {
  ConstantPitchProvider provider(0.0, 0.0);
  const double momentum = 1.0e-19;
  const Transport::ScalarResult speed = Transport::SpeedFromMomentum(
      momentum, kProtonMassKg, kSpeedOfLightMPerS);
  const double gradient = -1.0e-9;
  const double dt_s = 0.02;
  const double mu0 = 0.3;
  const Transport::FocusedTransportIncrement step = FocusedStep(
      {0.0, momentum, mu0}, {gradient, 0.0, 0.0, 0.0}, dt_s, provider, 1);
  const auto half_midpoint = [&](double start) {
    const double interval = 0.5 * dt_s;
    const double rate0 = -0.5 * (1.0 - start * start) * speed.value * gradient;
    const double midpoint = start + 0.5 * interval * rate0;
    const double rate_mid =
        -0.5 * (1.0 - midpoint * midpoint) * speed.value * gradient;
    return start + interval * rate_mid;
  };
  const double expected = half_midpoint(half_midpoint(mu0));
  const double error = std::fabs(step.state.mu - expected);
  Result result = Complete(step.status.ok() && error <= 1.0e-13,
      "focusing follows the independently evaluated split midpoint update",
      "focused mover uses an inconsistent focusing convention");
  result.configuration = {"mu0=0.3", "dlnB_ds_m^-1=-1e-9", "dt_s=0.02"};
  result.metrics.push_back({"mu_absolute_error", error, 1.0e-13, "<=", "dimensionless"});
  return result;
}

Result RunFted04() {
  ConstantPitchProvider provider(0.0, 0.0);
  const Transport::FocusedTransportIncrement step = FocusedStep(
      {10.0, 1.0e-19, 0.4}, {-2.0e-10, 400000.0, 3.0e-5, 2.0e-5},
      0.5, provider, 1);
  const double expected_momentum = 1.0e-19 * std::exp(-2.0e-5 * 0.5 / 3.0);
  const double momentum_error =
      std::fabs(step.state.momentumKgMPerS - expected_momentum);
  const bool pass = step.status.ok() && momentum_error <= 1.0e-32 &&
      std::isfinite(step.displacementM) && step.state.arcLengthM > 10.0;
  Result result = Complete(pass,
      "combined controlled focused step preserves its exact cooling characteristic",
      "combined focusing/gradient/streaming/cooling step is invalid");
  result.configuration = {"mu0=0.4", "dlnB_ds_m^-1=-2e-10",
      "U_m_per_s=400000", "dU_ds_s^-1=3e-5", "divU_s^-1=2e-5", "dt_s=0.5"};
  result.metrics.push_back({"momentum_absolute_error", momentum_error, 1.0e-32, "<=", "kg m s^-1"});
  result.metrics.push_back({"positive_displacement", step.state.arcLengthM > 10.0 ? 1.0 : 0.0, 1.0, "==", "boolean"});
  return result;
}

Result RunFted05() {
  unsigned reflections_a = 0;
  unsigned reflections_b = 0;
  const double a = Transport::ReflectPitchAngle(1.2, &reflections_a);
  const double b = Transport::ReflectPitchAngle(-5.7, &reflections_b);
  const double zero = Transport::ReflectPitchAngle(0.0, NULL);
  const bool pass = NearAbsolute(a, 0.8, 1.0e-15) && b >= -1.0 && b <= 1.0 &&
      reflections_a > 0 && reflections_b > 0 && zero == 0.0;
  Result result = Complete(pass,
      "reflective pitch boundaries handle zero and arbitrary overshoot",
      "pitch-angle reflection left the physical interval");
  result.configuration = {"samples_mu=1.2,-5.7,0"};
  result.metrics.push_back({"maximum_mu", std::max(std::fabs(a), std::fabs(b)), 1.0, "<=", "dimensionless"});
  return result;
}

Result RunFted06() {
  ConstantPitchProvider provider(0.0, 0.0, "turbulence:controlled:g17");
  Transport::ThreadLocalWaveAccumulator accumulator;
  const Transport::FocusedTransportIncrement step = FocusedStep(
      {0.0, 1.0e-19, 0.5}, {0.0, 0.0, 0.0, 0.0}, 1.0,
      provider, 1, &accumulator);
  const std::vector<Transport::WaveContribution>& values =
      accumulator.Contributions();
  const double magnetic_field_t = 5.0e-9;
  const double fluctuation_t = 2.0e-9;
  const double speed_m_per_s = 1.0e5;
  const double pitch_mu = 0.5;
  const double omega = QLT::proton_charge * magnetic_field_t / QLT::proton_mass;
  const double wave_number = omega / (speed_m_per_s * std::fabs(pitch_mu));
  const double bandwidth = std::pow(QLT::k_min_1AU, -2.0 / 3.0) -
                           std::pow(QLT::k_max_1AU, -2.0 / 3.0);
  const double normalization =
      (2.0 / 3.0) * fluctuation_t * fluctuation_t / bandwidth;
  const double spectrum = normalization * std::pow(wave_number, -5.0 / 3.0);
  const double expected_dmumu = 1.57079632679489661923 * omega * omega *
      (1.0 - pitch_mu * pitch_mu) * spectrum /
      (magnetic_field_t * magnetic_field_t * speed_m_per_s * std::fabs(pitch_mu));
  const double actual_dmumu = QLT::calculateDmuMu(
      magnetic_field_t, fluctuation_t, speed_m_per_s, pitch_mu, QLT::r0);
  const double relative_error = std::fabs(actual_dmumu - expected_dmumu) /
      std::max(1.0, std::fabs(expected_dmumu));
  const bool pass = step.status.ok() && values.size() == 1 &&
      values[0].turbulenceStateIdentity == "turbulence:controlled:g17" &&
      step.coefficientProvenance == "controlled-test:constant-dmumu-v1" &&
      relative_error <= 2.0e-14;
  Result result = Complete(pass,
      "QLT normalization and turbulence identity close independently",
      "QLT coefficient or wave-state identity differs from its controlled reference");
  result.configuration = {"B_T=5e-9", "deltaB_T=2e-9", "speed_m_per_s=1e5", "mu=0.5"};
  result.metrics.push_back({"Dmumu_scaled_error", relative_error, 2.0e-14, "<=", "dimensionless"});
  result.metrics.push_back({"wave_contribution_count", static_cast<double>(values.size()), 1.0, "==", "count"});
  return result;
}

double IntegrateFocused(double dt_s) {
  ConstantPitchProvider provider(0.0, 0.0);
  Transport::FocusedTransportState state = {0.0, 1.0e-19, 0.2};
  for (double time_s = 0.0; time_s < 1.0 - 0.5 * dt_s; time_s += dt_s)
    state = FocusedStep(state, {-2.0e-9, 0.0, 0.0, 0.0}, dt_s,
                        provider, static_cast<std::uint64_t>(time_s / dt_s)).state;
  return state.mu;
}

Result RunFted07() {
  const double reference = IntegrateFocused(0.0005);
  const double coarse = std::fabs(IntegrateFocused(0.02) - reference);
  const double fine = std::fabs(IntegrateFocused(0.01) - reference);
  const RefinementOrderEstimate order =
      EstimateRefinementOrder(coarse, fine, 0.02, 0.01);
  const double minimum_order = 1.80;
  Result result = Complete(order.valid && order.observedOrder >= minimum_order,
      "symmetric deterministic focusing split converges under refinement",
      "focused deterministic split misses its observed-order target");
  result.hasSeed = true;
  result.seed = 31;
  result.configuration = {"dlnB_ds_m^-1=-2e-9", "duration_s=1", "dt_s=0.02,0.01,0.0005"};
  result.metrics.push_back({"observed_refinement_order", order.observedOrder,
      minimum_order, ">=", "dimensionless"});
  result.metrics.push_back({"coarse_to_fine_error_ratio",
      order.coarseToFineErrorRatio, 1.0, ">", "dimensionless"});
  return result;
}

double LegendreP2(double mu) { return 0.5 * (3.0 * mu * mu - 1.0); }

double SampleLegendreP2(double amplitude, Transport::RandomStream& random) {
  // Rejection sampling is independent of the mover.  The envelope is exact
  // because max(P2)=1 on [-1,1], and amplitude is chosen in [0,1] so the
  // prescribed density 0.5*(1+a*P2) is everywhere non-negative.
  for (;;) {
    const double mu = 2.0 * random.UniformOpen01() - 1.0;
    const double acceptance =
        (1.0 + amplitude * LegendreP2(mu)) / (1.0 + amplitude);
    if (random.UniformOpen01() <= acceptance) return mu;
  }
}

Result RunFted08() {
  const std::uint64_t seed = 308;
  const std::size_t particles = 40000;
  const std::size_t steps = 100;
  const double dt_s = 0.005;
  const double duration_s = steps * dt_s;
  const double d0_per_s = 0.2;
  const double amplitude = 0.8;
  LegendrePitchProvider provider(d0_per_s);
  const Transport::ScalarResult momentum = Transport::MomentumFromSpeed(
      1.0e6, kProtonMassKg, kSpeedOfLightMPerS);
  long double sum = 0.0L;
  long double sum2 = 0.0L;
  bool status_ok = momentum.status.ok();
  for (std::size_t particle = 0; particle < particles; ++particle) {
    Transport::KeyedRandomStream initial_random(seed, particle, 0, 3080);
    double mu = SampleLegendreP2(amplitude, initial_random);
    for (std::size_t step_index = 0; step_index < steps; ++step_index) {
      Transport::KeyedRandomStream mover_random(
          seed, particle, step_index, 3081);
      const Transport::FocusedTransportIncrement step =
          Transport::AdvanceFocusedTransportDmumu(
              {0.0, momentum.value, mu}, {0.0, 0.0, 0.0, 0.0},
              kProtonMassKg, kSpeedOfLightMPerS, dt_s, provider,
              mover_random, NULL);
      if (!step.status.ok()) {
        status_ok = false;
        break;
      }
      mu = step.state.mu;
    }
    const double mode = LegendreP2(mu);
    sum += mode;
    sum2 += mode * mode;
  }
  const double measured = static_cast<double>(sum / particles);
  const double variance = std::max(0.0,
      static_cast<double>(sum2 / particles) - measured * measured);
  // Orthogonality gives <P_l>(0)=a/(2l+1).  For
  // D_mumu=D0(1-mu^2), the l-th mode decays as
  // exp[-l(l+1)D0*t]; here l=2.
  const double expected = amplitude / 5.0 *
      std::exp(-6.0 * d0_per_s * duration_s);
  const double standard_error = std::sqrt(variance / particles);
  const double z_score = std::fabs(measured - expected) / standard_error;
  const bool pass = status_ok && std::isfinite(z_score) && z_score <= 5.0;
  Result result = Complete(pass,
      "P2 pitch-angle eigenmode follows its analytical Legendre decay",
      "P2 pitch-angle mode exceeds its preregistered statistical tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration = {"particles=40000", "steps=100", "dt_s=0.005",
      "Dmumu_s^-1=0.2*(1-mu^2)", "initial_f=0.5*(1+0.8*P2(mu))"};
  result.metrics.push_back({"P2_moment_z_score", z_score, 5.0, "<=", "sigma"});
  result.metrics.push_back({"P2_moment_absolute_error",
      std::fabs(measured - expected), 5.0 * standard_error, "<=", "dimensionless"});
  return result;
}

Result RunFtem01() {
  const std::uint64_t seed = 101;
  const int samples = 120000;
  const double rate_per_s = 2.5;
  long double sum = 0.0L;
  bool status_ok = true;
  for (int i = 0; i < samples; ++i) {
    Transport::KeyedRandomStream random(seed, i, 9, 0);
    const Transport::ScalarResult wait =
        Transport::SampleExponentialWaitingTime(rate_per_s, random);
    status_ok = status_ok && wait.status.ok();
    sum += wait.value;
  }
  const double error_s =
      std::fabs(static_cast<double>(sum / samples) - 1.0 / rate_per_s);
  Result result = Complete(status_ok && error_s < 0.004,
      "exponential waiting-time mean matches the analytical inverse rate",
      "waiting-time ensemble mean exceeds tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration = {"samples=120000", "nu_s^-1=2.5"};
  result.metrics.push_back({"mean_wait_absolute_error", error_s, 0.004, "<", "s"});
  return result;
}

Result RunFtem02() {
  const int histories = 30000;
  const double speed_m_per_s = 1.0e6;
  const double rate_per_s = 2.0;
  long long events = 0;
  bool status_ok = true;
  for (int i = 0; i < histories; ++i) {
    const Transport::FocusedTransportMfpIncrement step = AdvanceMfp(
        speed_m_per_s, 0.4, 1.0, 1.0, speed_m_per_s / rate_per_s,
        {0.0, 0.0, 0.0, 0.0, 0.0}, static_cast<std::uint64_t>(i));
    status_ok = status_ok && step.status.ok();
    events += static_cast<long long>(step.diagnostics.scatteringEvents);
  }
  const double mean = static_cast<double>(events) / histories;
  const double error = std::fabs(mean - rate_per_s);
  Result result = Complete(status_ok && error < 0.025,
      "event counts reproduce the analytical Poisson mean",
      "event-count mean exceeds its statistical tolerance");
  result.hasSeed = true;
  result.seed = 91;
  result.configuration = {"histories=30000", "nu_s^-1=2", "duration_s=1"};
  result.metrics.push_back({"event_mean_absolute_error", error, 0.025, "<", "count"});
  return result;
}

Result RunFtem03() {
  const double speed = 2.0e6;
  const double mu = 0.25;
  const double advection = 4.0e5;
  const double dt_s = 3.0;
  const Transport::FocusedTransportMfpIncrement step = AdvanceMfp(
      speed, mu, dt_s, dt_s, std::numeric_limits<double>::infinity(),
      {0.0, advection, 0.0, 0.0, 0.0}, 1);
  const double expected = (advection + speed * mu) * dt_s;
  const double relative_error = std::fabs(step.displacementM - expected) /
      std::max(1.0, std::fabs(expected));
  const bool pass = step.status.ok() && relative_error <= 2.0e-14 &&
      step.diagnostics.scatteringEvents == 0 &&
      step.diagnostics.ballisticIntervals == 1;
  Result result = Complete(pass,
      "infinite mean free path follows the exact ballistic characteristic",
      "ballistic MFP limit differs from (U+v*mu)*dt");
  result.configuration = {"speed_m_per_s=2e6", "mu=0.25", "U_m_per_s=4e5", "dt_s=3", "lambda=inf"};
  result.metrics.push_back({"displacement_relative_error", relative_error, 2.0e-14, "<=", "dimensionless"});
  return result;
}

Result RunFtem04() {
  Transport::KeyedRandomStream random(104, 1, 9, 0);
  const Transport::WaveFrameScatterResult scatter =
      Transport::ScatterIsotropicallyInWaveFrame(
          0.7 * kSpeedOfLightMPerS, 0.6, -4.0e5,
          kSpeedOfLightMPerS, random);
  const double relative_error = std::fabs(
      scatter.waveFrameSpeedAfterMPerS - scatter.waveFrameSpeedBeforeMPerS) /
      std::max(1.0, std::fabs(scatter.waveFrameSpeedBeforeMPerS));
  const bool pass = scatter.status.ok() && scatter.speedMPerS > 0.0 &&
      scatter.speedMPerS < kSpeedOfLightMPerS && scatter.mu >= -1.0 &&
      scatter.mu <= 1.0 && relative_error <= 3.0e-14;
  Result result = Complete(pass,
      "Lorentz scattering preserves speed in the selected Alfven-wave frame",
      "wave-frame elastic-scattering invariant was violated");
  result.hasSeed = true;
  result.seed = 104;
  result.configuration = {"speed=0.7c", "mu=0.6", "wave_speed_m_per_s=-4e5"};
  result.metrics.push_back({"wave_frame_speed_relative_error", relative_error, 3.0e-14, "<=", "dimensionless"});
  return result;
}

double FocusedMfpMu(double maximum_interval_s) {
  return AdvanceMfp(1.0e6, 0.2, 2.0, maximum_interval_s,
      std::numeric_limits<double>::infinity(),
      {-8.0e-8, 0.0, 0.0, 0.0, 0.0}, 1).state.mu;
}

Result RunFtem05() {
  const double reference = FocusedMfpMu(0.0005);
  const double coarse = std::fabs(FocusedMfpMu(0.1) - reference);
  const double fine = std::fabs(FocusedMfpMu(0.05) - reference);
  const RefinementOrderEstimate order =
      EstimateRefinementOrder(coarse, fine, 0.1, 0.05);
  const double minimum_order = 1.80;
  Result result = Complete(order.valid && order.observedOrder >= minimum_order,
      "MFP focusing characteristic converges under interval refinement",
      "MFP focusing split misses its observed-order target");
  result.configuration = {"dlnB_ds_m^-1=-8e-8", "duration_s=2", "max_interval_s=0.1,0.05,0.0005"};
  result.metrics.push_back({"observed_refinement_order", order.observedOrder,
      minimum_order, ">=", "dimensionless"});
  result.metrics.push_back({"coarse_to_fine_error_ratio",
      order.coarseToFineErrorRatio, 1.0, ">", "dimensionless"});
  return result;
}

Result RunFtem06() {
  const double speed = 3.0e6;
  const double divergence = 4.0e-4;
  const double dt_s = 5.0;
  const Transport::ScalarResult p0 = Transport::MomentumFromSpeed(
      speed, kProtonMassKg, kSpeedOfLightMPerS);
  const Transport::FocusedTransportMfpIncrement step = AdvanceMfp(
      speed, 0.3, dt_s, 0.2, std::numeric_limits<double>::infinity(),
      {0.0, 0.0, 0.0, divergence, 0.0}, 1);
  const double expected = p0.value * std::exp(-divergence * dt_s / 3.0);
  const double relative_error = std::fabs(step.state.momentumKgMPerS - expected) /
      std::max(std::fabs(expected), std::numeric_limits<double>::min());
  Result result = Complete(step.status.ok() && relative_error <= 5.0e-14,
      "MFP mover applies the exact plasma-frame cooling characteristic",
      "MFP momentum differs from its analytical cooling solution");
  result.configuration = {"speed_m_per_s=3e6", "divU_s^-1=4e-4", "dt_s=5", "lambda=inf"};
  result.metrics.push_back({"momentum_relative_error", relative_error, 5.0e-14, "<=", "dimensionless"});
  return result;
}

double CombinedMfpState(double maximum_interval_s) {
  const Transport::FocusedTransportMfpIncrement step = AdvanceMfp(
      1.5e6, -0.35, 1.0, maximum_interval_s,
      std::numeric_limits<double>::infinity(),
      {-2.0e-8, 3.0e5, 1.0e-3, 2.0e-4, 0.0}, 1);
  return step.state.arcLengthM + 1.0e6 * step.state.mu;
}

Result RunFtem07() {
  const double reference = CombinedMfpState(0.0005);
  const double coarse = std::fabs(CombinedMfpState(0.05) - reference);
  const double fine = std::fabs(CombinedMfpState(0.025) - reference);
  const RefinementOrderEstimate order =
      EstimateRefinementOrder(coarse, fine, 0.05, 0.025);
  const double minimum_order = 0.85;
  Result result = Complete(order.valid && order.observedOrder >= minimum_order,
      "combined MFP characteristic converges under interval refinement",
      "combined MFP split misses its observed-order target");
  result.configuration = {"duration_s=1", "lambda=inf", "max_interval_s=0.05,0.025,0.0005"};
  result.metrics.push_back({"observed_refinement_order", order.observedOrder,
      minimum_order, ">=", "dimensionless"});
  result.metrics.push_back({"coarse_to_fine_error_ratio",
      order.coarseToFineErrorRatio, 1.0, ">", "dimensionless"});
  return result;
}

struct PersistentFlightMoment {
  bool ok = false;
  double meanSquareDisplacementM2 = 0.0;
};

PersistentFlightMoment MeasurePersistentFlightMsd(
    double duration_s, std::uint64_t stream_purpose) {
  const std::uint64_t seed = 908;
  const std::size_t particles = 40000;
  const double speed_m_per_s = 1.0e6;
  const double lambda_m = 1.0e7;
  const Transport::ScalarResult momentum = Transport::MomentumFromSpeed(
      speed_m_per_s, kProtonMassKg, kSpeedOfLightMPerS);
  ConstantMfpProvider provider(lambda_m);
  long double sum_square = 0.0L;
  bool status_ok = momentum.status.ok();
  for (std::size_t particle = 0; particle < particles && status_ok; ++particle) {
    // Initial pitch and event history use distinct purpose keys.  This avoids
    // introducing an artificial correlation between the isotropic initial
    // ensemble and the first exponentially distributed scattering time.
    Transport::KeyedRandomStream initial_random(
        seed, particle, 9, stream_purpose);
    const double initial_mu =
        2.0 * initial_random.UniformOpen01() - 1.0;
    Transport::KeyedRandomStream mover_random(
        seed, particle, 9, stream_purpose + 1);
    const Transport::FocusedTransportMfpIncrement step =
        Transport::AdvanceFocusedTransportMfp(
            {0.0, momentum.value, initial_mu},
            {0.0, 0.0, 0.0, 0.0, 0.0}, kProtonMassKg,
            kSpeedOfLightMPerS, duration_s, duration_s,
            provider, mover_random, NULL);
    status_ok = status_ok && step.status.ok();
    sum_square += static_cast<long double>(step.state.arcLengthM) *
                  step.state.arcLengthM;
  }
  PersistentFlightMoment result;
  result.ok = status_ok;
  result.meanSquareDisplacementM2 =
      static_cast<double>(sum_square / particles);
  return result;
}

double PersistentFlightMsdExact(double speed_m_per_s,
                                double rate_per_s,
                                double duration_s) {
  // Isotropic pitch resets give the stationary correlation
  // <v_parallel(0)v_parallel(t)>=v^2 exp(-nu*t)/3.  Integrating this
  // correlation twice yields the exact finite-time MSD below.  At nu*t>>1 it
  // approaches 2*kappa_parallel*t with kappa_parallel=v^2/(3*nu)=v*lambda/3.
  const double q = rate_per_s * duration_s;
  return 2.0 * speed_m_per_s * speed_m_per_s /
      (3.0 * rate_per_s * rate_per_s) * (q - 1.0 + std::exp(-q));
}

Result RunFtem08() {
  const double speed_m_per_s = 1.0e6;
  const double lambda_m = 1.0e7;
  const double rate_per_s = speed_m_per_s / lambda_m;
  const double short_time_s = 2.0;   // nu*t=0.2: persistent/ballistic regime.
  const double long_time_s = 500.0;  // nu*t=50: diffusion asymptote.
  const PersistentFlightMoment short_time =
      MeasurePersistentFlightMsd(short_time_s, 9080);
  const PersistentFlightMoment long_time =
      MeasurePersistentFlightMsd(long_time_s, 9090);
  const double exact_short = PersistentFlightMsdExact(
      speed_m_per_s, rate_per_s, short_time_s);
  const double exact_long = PersistentFlightMsdExact(
      speed_m_per_s, rate_per_s, long_time_s);
  const double diffusion_coefficient = speed_m_per_s * lambda_m / 3.0;
  const double short_relative_error = std::fabs(
      short_time.meanSquareDisplacementM2 / exact_short - 1.0);
  const double long_relative_error = std::fabs(
      long_time.meanSquareDisplacementM2 / exact_long - 1.0);
  const double diffusion_limit_error = std::fabs(
      long_time.meanSquareDisplacementM2 /
          (2.0 * diffusion_coefficient * long_time_s) - 1.0);
  const double moment_tolerance = 0.04;
  const double diffusion_limit_tolerance = 0.07;
  const bool pass = short_time.ok && long_time.ok &&
      short_relative_error <= moment_tolerance &&
      long_relative_error <= moment_tolerance &&
      diffusion_limit_error <= diffusion_limit_tolerance;
  Result result = Complete(pass,
      "event-driven MFP transport matches persistent-flight and diffusion-limit MSD",
      "event-driven MFP mean-square displacement misses its analytical limits");
  result.hasSeed = true;
  result.seed = 908;
  result.configuration = {"particles=40000", "speed_m_per_s=1e6",
      "lambda_parallel_m=1e7", "nu_s^-1=0.1",
      "short_duration_s=2", "long_duration_s=500",
      "initial_mu=isotropic", "alfven_speed_m_per_s=0"};
  result.metrics.push_back({"persistent_flight_MSD_relative_error",
      short_relative_error, moment_tolerance, "<=", "dimensionless"});
  result.metrics.push_back({"long_time_MSD_relative_error",
      long_relative_error, moment_tolerance, "<=", "dimensionless"});
  result.metrics.push_back({"diffusion_limit_relative_error",
      diffusion_limit_error, diffusion_limit_tolerance, "<=", "dimensionless"});
  return result;
}

}  // namespace

std::vector<Descriptor> ControlledMoverDescriptors() {
  std::vector<Descriptor> descriptors;
  descriptors.push_back(MakeDescriptor("PARK01", "Exact Parker convection", "parker",
      "Compare zero-diffusion transport with ds=U dt.", RuntimeClass::Routine,
      "fixed keyed seed 1; stochastic term disabled", RunPark01));
  descriptors.push_back(MakeDescriptor("PARK02", "Parker Gaussian moments", "parker",
      "Compare constant-diffusion ensemble moments with variance 2*kappa*dt.", RuntimeClass::Extended,
      "fixed keyed seed 2", RunPark02));
  descriptors.push_back(MakeDescriptor("PARK03", "Parker Ito gradient drift", "parker",
      "Check the sign and magnitude of the analytical d(kappa)/ds drift.", RuntimeClass::Routine,
      "fixed keyed seed 3; stochastic term disabled", RunPark03));
  descriptors.push_back(MakeDescriptor("PARK04", "Parker adiabatic cooling", "parker",
      "Compare momentum with p0*exp[-div(U)t/3].", RuntimeClass::Routine,
      "fixed keyed seed 4; stochastic term disabled", RunPark04));
  descriptors.push_back(MakeDescriptor("PARK05", "Parker absorbing boundary", "parker",
      "Verify a controlled outward crossing reports absorbing flux.", RuntimeClass::Routine,
      "deterministic; no RNG", RunPark05));
  descriptors.push_back(MakeDescriptor("PARK06", "Parker drift refinement", "parker",
      "Verify refinement convergence for ds/dt=0.4s.", RuntimeClass::Routine,
      "fixed keyed seed 6; stochastic term disabled", RunPark06));
  descriptors.push_back(MakeDescriptor("PARK07", "Parker first passage", "parker",
      "Compare absorbing-boundary exit side and mean time with exact Brownian first passage.", RuntimeClass::Extended,
      "fixed keyed seed 407", RunPark07));

  descriptors.push_back(MakeDescriptor("FTED01", "Dmumu Ito moments", "fte-dmumu",
      "Compare a constant-Dmumu ensemble with analytical Ito moments.", RuntimeClass::Extended,
      "fixed keyed seed 31", RunFted01));
  descriptors.push_back(MakeDescriptor("FTED02", "Dmumu Ito drift", "fte-dmumu",
      "Check the supplied dDmumu/dmu drift independently.", RuntimeClass::Routine,
      "fixed keyed seed 31; stochastic term disabled", RunFted02));
  descriptors.push_back(MakeDescriptor("FTED03", "Dmumu focusing convention", "fte-dmumu",
      "Compare focusing with an independent split-midpoint calculation.", RuntimeClass::Routine,
      "fixed keyed seed 31; stochastic term disabled", RunFted03));
  descriptors.push_back(MakeDescriptor("FTED04", "Dmumu combined controlled step", "fte-dmumu",
      "Exercise focusing, flow gradient, streaming, and exact cooling together.", RuntimeClass::Routine,
      "fixed keyed seed 31; stochastic term disabled", RunFted04));
  descriptors.push_back(MakeDescriptor("FTED05", "Dmumu reflecting boundary", "fte-dmumu",
      "Verify regular reflection for zero and arbitrary pitch overshoot.", RuntimeClass::Routine,
      "deterministic; no RNG", RunFted05));
  descriptors.push_back(MakeDescriptor("FTED06", "Dmumu QLT and wave identity", "fte-dmumu",
      "Reconstruct the QLT coefficient and verify wave-state identity.", RuntimeClass::Routine,
      "fixed keyed seed 31; stochastic term disabled", RunFted06));
  descriptors.push_back(MakeDescriptor("FTED07", "Dmumu split refinement", "fte-dmumu",
      "Verify deterministic symmetric-split convergence.", RuntimeClass::Routine,
      "fixed keyed seed 31; stochastic term disabled", RunFted07));
  descriptors.push_back(MakeDescriptor("FTED08", "Dmumu Legendre eigenmode", "fte-dmumu",
      "Compare the P2 pitch-angle moment with exp[-6 D0 t] decay.", RuntimeClass::Extended,
      "fixed keyed seed 308", RunFted08));

  descriptors.push_back(MakeDescriptor("FTEM01", "MFP exponential waiting time", "fte-mfp",
      "Compare the sampled waiting-time mean with 1/nu.", RuntimeClass::Extended,
      "fixed keyed seed 101", RunFtem01));
  descriptors.push_back(MakeDescriptor("FTEM02", "MFP Poisson event count", "fte-mfp",
      "Compare event-count mean with nu*t.", RuntimeClass::Extended,
      "fixed keyed seed 91", RunFtem02));
  descriptors.push_back(MakeDescriptor("FTEM03", "MFP ballistic limit", "fte-mfp",
      "Compare infinite-lambda transport with (U+v*mu)dt.", RuntimeClass::Routine,
      "fixed keyed seed 91; no event occurs", RunFtem03));
  descriptors.push_back(MakeDescriptor("FTEM04", "MFP wave-frame elasticity", "fte-mfp",
      "Verify exact speed conservation in the scattering wave frame.", RuntimeClass::Routine,
      "fixed keyed seed 104", RunFtem04));
  descriptors.push_back(MakeDescriptor("FTEM05", "MFP focusing refinement", "fte-mfp",
      "Verify focusing convergence between scattering events.", RuntimeClass::Routine,
      "fixed keyed seed 91; ballistic closure", RunFtem05));
  descriptors.push_back(MakeDescriptor("FTEM06", "MFP adiabatic cooling", "fte-mfp",
      "Compare momentum with the exact plasma-frame cooling solution.", RuntimeClass::Routine,
      "fixed keyed seed 91; ballistic closure", RunFtem06));
  descriptors.push_back(MakeDescriptor("FTEM07", "MFP combined refinement", "fte-mfp",
      "Verify convergence of focusing, streaming, and cooling splitting.", RuntimeClass::Routine,
      "fixed keyed seed 91; ballistic closure", RunFtem07));
  descriptors.push_back(MakeDescriptor("FTEM08", "MFP persistent-flight limit", "fte-mfp",
      "Compare finite-time telegraph-like MSD and its diffusion asymptote with exact correlation theory.", RuntimeClass::Extended,
      "fixed keyed seed 908", RunFtem08));
  return descriptors;
}

}  // namespace Testing
}  // namespace SEP
