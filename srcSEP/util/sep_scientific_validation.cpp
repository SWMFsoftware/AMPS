#include "sep_scientific_validation.h"

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_registry.h)
#include "sep_focused_transport_core.h"
#include "sep_focused_transport_mfp_core.h"
#include "sep_parker_core.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace SEP {
namespace Testing {
namespace {

const double kProtonMassKg = 1.67262192369e-27;
const double kSpeedOfLightMPerS = 299792458.0;

Result ValidationResult(bool pass, const std::string& success,
                        const std::string& failure) {
  Result result;
  result.status = pass ? Status::Pass : Status::Fail;
  result.message = pass ? success : failure;
  result.metrics.push_back(
      {"assertion_failures", pass ? 0.0 : 1.0, 0.0, "<=", "count"});
  return result;
}

Descriptor MakeDescriptor(const char* id, const char* name,
                          const char* description, RuntimeClass runtime,
                          const char* seed_policy, TestCallback callback) {
  Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = "validation";
  descriptor.description = description;
  descriptor.initialization = InitializationLevel::None;
  descriptor.supportedBuildModes =
      "serial/MPI linked CLI and source-only sanitizer runner";
  descriptor.runtime = runtime;
  descriptor.seedPolicy = seed_policy;
  descriptor.stateIsolation =
      "stack/vector-owned ensemble state; no production globals or files mutated";
  descriptor.callback = callback;
  return descriptor;
}

double RelativeError(double actual, double expected) {
  return std::fabs(actual - expected) /
      std::max(std::fabs(expected), std::numeric_limits<double>::min());
}

class ConstantSpatialDiffusion final
    : public Transport::SpatialDiffusionProvider {
 public:
  explicit ConstantSpatialDiffusion(double kappa_m2_per_s)
      : kappa_(kappa_m2_per_s) {}

  Transport::SpatialDiffusionSample Evaluate(double, double) const override {
    Transport::SpatialDiffusionSample sample;
    sample.status = Transport::Status::Ok();
    sample.kappaParallelM2PerS = kappa_;
    sample.dKappaParallelDsMPerS = 0.0;
    sample.provenance = "validation:constant-kappa-v1";
    return sample;
  }

 private:
  double kappa_;
};

class IsotropicPitchAngleDiffusion final
    : public Transport::PitchAngleDiffusionProvider {
 public:
  explicit IsotropicPitchAngleDiffusion(double scattering_rate_per_s)
      : rate_(scattering_rate_per_s) {}

  Transport::PitchAngleDiffusionSample Evaluate(double, double,
                                                  double mu) const override {
    Transport::PitchAngleDiffusionSample sample;
    sample.status = Transport::Status::Ok();
    // D_mumu=(nu/2)(1-mu^2) is the isotropic closure whose l=1
    // autocorrelation decays as exp(-nu*t).  Its derivative is mandatory for
    // the Ito drift and must be supplied with the same coefficient identity.
    sample.dMuMuPerS = 0.5 * rate_ * std::max(0.0, 1.0 - mu * mu);
    sample.dDmuMuDmuPerS = -rate_ * mu;
    sample.provenance = "validation:isotropic-dmumu-v1";
    sample.turbulenceStateIdentity = "validation:frozen-coefficients:g1";
    return sample;
  }

 private:
  double rate_;
};

class ConstantMeanFreePath final : public Transport::MeanFreePathProvider {
 public:
  explicit ConstantMeanFreePath(double lambda_parallel_m)
      : lambda_(lambda_parallel_m) {}

  Transport::MeanFreePathSample Evaluate(double, double,
                                          double) const override {
    Transport::MeanFreePathSample sample;
    sample.status = Transport::Status::Ok();
    sample.lambdaParallelM = lambda_;
    sample.provenance = "validation:constant-lambda-v1";
    sample.turbulenceStateIdentity = "validation:frozen-coefficients:g1";
    return sample;
  }

 private:
  double lambda_;
};

Result RunVal01ManufacturedParker() {
  // The reference is evaluated independently from the stochastic kernel:
  // X(t) is Gaussian with mean U*t and variance 2*kappa*t, while constant
  // plasma divergence gives p(t)=p0*exp[-div(U)t/3].  The ensemble tolerance
  // is stated in standard errors so increasing the sample count cannot hide a
  // systematic displacement bias behind an arbitrary percentage threshold.
  const std::uint64_t seed = 1501001;
  const std::size_t samples = 80000;
  const double kappa = 2.5e13;
  const double plasma_speed = 4.0e5;
  const double divergence = 1.5e-5;
  const double duration = 40.0;
  const double initial_momentum = 3.0e-19;
  ConstantSpatialDiffusion provider(kappa);

  long double displacement_sum = 0.0L;
  long double displacement_square_sum = 0.0L;
  double maximum_momentum_relative_error = 0.0;
  bool kernel_ok = true;
  const double expected_momentum = initial_momentum *
      std::exp(-divergence * duration / 3.0);
  for (std::size_t particle = 0; particle < samples; ++particle) {
    Transport::KeyedRandomStream random(
        seed, static_cast<std::uint64_t>(particle), 0, 1501);
    const Transport::ParkerIncrement step = Transport::AdvanceParker(
        Transport::ParkerState(0.0, initial_momentum),
        Transport::ParkerBackground(plasma_speed, divergence),
        1.2e7, duration, provider, random);
    if (!step.status.ok()) {
      kernel_ok = false;
      break;
    }
    displacement_sum += step.displacementM;
    displacement_square_sum +=
        static_cast<long double>(step.displacementM) * step.displacementM;
    maximum_momentum_relative_error = std::max(
        maximum_momentum_relative_error,
        RelativeError(step.state.momentumKgMPerS, expected_momentum));
  }

  const double mean = static_cast<double>(displacement_sum / samples);
  const double second =
      static_cast<double>(displacement_square_sum / samples);
  const double variance = std::max(0.0, second - mean * mean);
  const double expected_mean = plasma_speed * duration;
  const double expected_variance = 2.0 * kappa * duration;
  const double mean_standard_error =
      std::sqrt(expected_variance / static_cast<double>(samples));
  const double mean_z_score =
      std::fabs(mean - expected_mean) / mean_standard_error;
  const double variance_relative_error =
      RelativeError(variance, expected_variance);
  const double momentum_tolerance =
      64.0 * std::numeric_limits<double>::epsilon();
  const bool pass = kernel_ok && mean_z_score <= 4.5 &&
      variance_relative_error <= 0.02 &&
      maximum_momentum_relative_error <= momentum_tolerance;

  Result result = ValidationResult(
      pass,
      "Parker ensemble matches independent diffusion and cooling solutions",
      "Parker manufactured benchmark exceeded its preregistered tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("samples=80000");
  result.configuration.push_back("kappa_parallel_m2_per_s=2.5e13");
  result.configuration.push_back("plasma_speed_m_per_s=4e5");
  result.configuration.push_back("div_velocity_s^-1=1.5e-5");
  result.configuration.push_back("duration_s=40");
  result.metrics.push_back(
      {"mean_displacement_z_score", mean_z_score, 4.5, "<=", "sigma"});
  result.metrics.push_back({"variance_relative_error",
      variance_relative_error, 0.02, "<=", "dimensionless"});
  result.metrics.push_back({"momentum_relative_error",
      maximum_momentum_relative_error, momentum_tolerance, "<=",
      "dimensionless"});
  return result;
}

struct ProfileSummary {
  std::vector<double> normalized_intensity;
  double onset_s = std::numeric_limits<double>::quiet_NaN();
  double peak_s = std::numeric_limits<double>::quiet_NaN();
  double fluence = 0.0;
};

ProfileSummary SummarizeProfile(const std::vector<double>& counts,
                                double sampling_interval_s) {
  ProfileSummary summary;
  if (counts.empty()) return summary;
  std::vector<double> intensity(counts.size(), 0.0);
  for (std::size_t i = 0; i < counts.size(); ++i) {
    // A fixed three-cadence moving average is declared before either mover is
    // inspected.  Peak time from a raw 4000-particle occupancy series is an
    // unstable order statistic; this forward operator represents the finite
    // cadence/averaging that must also be applied to observational profiles.
    const std::size_t first = i == 0 ? 0 : i - 1;
    const std::size_t last = std::min(counts.size() - 1, i + 1);
    for (std::size_t j = first; j <= last; ++j) intensity[i] += counts[j];
    intensity[i] /= static_cast<double>(last - first + 1);
  }
  const double peak = *std::max_element(intensity.begin(), intensity.end());
  summary.normalized_intensity.resize(counts.size(), 0.0);
  if (!(peak > 0.0)) return summary;

  // Onset is defined before looking at either mover as the first sample at ten
  // percent of that mover's peak.  Fluence is the trapezoidal time integral of
  // the normalized detector occupancy, so absolute macroparticle weight and
  // bin width cancel in the cross-mover comparison.
  const double onset_threshold = 0.10 * peak;
  for (std::size_t i = 0; i < counts.size(); ++i) {
    summary.normalized_intensity[i] = intensity[i] / peak;
    if (!std::isfinite(summary.onset_s) && intensity[i] >= onset_threshold)
      summary.onset_s = (i + 1) * sampling_interval_s;
  }
  const std::size_t peak_index = static_cast<std::size_t>(
      std::max_element(intensity.begin(), intensity.end()) - intensity.begin());
  summary.peak_s = (peak_index + 1) * sampling_interval_s;
  for (std::size_t i = 1; i < summary.normalized_intensity.size(); ++i) {
    summary.fluence += 0.5 * sampling_interval_s *
        (summary.normalized_intensity[i - 1] +
         summary.normalized_intensity[i]);
  }
  return summary;
}

Result RunVal02CrossMoverCampaign() {
  const std::uint64_t seed = 1502001;
  const std::size_t particles = 4000;
  const std::size_t steps = 1600;
  const std::size_t sampling_stride = 20;
  const double speed = 1.0e7;
  const double lambda = 1.0e8;
  const double rate = speed / lambda;
  const double dt = 0.25;
  const double observer_m = 5.0e8;
  const double detector_half_width_m = 5.0e7;
  const Transport::ScalarResult momentum = Transport::MomentumFromSpeed(
      speed, kProtonMassKg, kSpeedOfLightMPerS);
  IsotropicPitchAngleDiffusion dmumu_provider(rate);
  ConstantMeanFreePath mfp_provider(lambda);

  std::vector<Transport::FocusedTransportState> dmumu_states;
  std::vector<Transport::FocusedTransportMfpState> mfp_states;
  dmumu_states.reserve(particles);
  mfp_states.reserve(particles);
  for (std::size_t particle = 0; particle < particles; ++particle) {
    Transport::KeyedRandomStream initial_random(
        seed, static_cast<std::uint64_t>(particle), 0, 15020);
    const double mu = -1.0 + 2.0 * initial_random.UniformOpen01();
    dmumu_states.push_back(
        Transport::FocusedTransportState(0.0, momentum.value, mu));
    mfp_states.push_back(
        Transport::FocusedTransportMfpState(0.0, momentum.value, mu));
  }

  std::vector<double> dmumu_counts;
  std::vector<double> mfp_counts;
  bool kernels_ok = momentum.status.ok();
  for (std::size_t step_index = 0; step_index < steps && kernels_ok;
       ++step_index) {
    for (std::size_t particle = 0; particle < particles; ++particle) {
      Transport::KeyedRandomStream dmumu_random(
          seed, static_cast<std::uint64_t>(particle), step_index, 15021);
      const Transport::FocusedTransportIncrement dmumu =
          Transport::AdvanceFocusedTransportDmumu(
              dmumu_states[particle],
              Transport::FocusedTransportBackground(0.0, 0.0, 0.0, 0.0),
              kProtonMassKg, kSpeedOfLightMPerS, dt, dmumu_provider,
              dmumu_random, NULL);
      Transport::KeyedRandomStream mfp_random(
          seed, static_cast<std::uint64_t>(particle), step_index, 15022);
      const Transport::FocusedTransportMfpIncrement mfp =
          Transport::AdvanceFocusedTransportMfp(
              mfp_states[particle],
              Transport::FocusedTransportMfpBackground(
                  0.0, 0.0, 0.0, 0.0, 0.0),
              kProtonMassKg, kSpeedOfLightMPerS, dt, dt, mfp_provider,
              mfp_random, NULL);
      if (!dmumu.status.ok() || !mfp.status.ok()) {
        kernels_ok = false;
        break;
      }
      dmumu_states[particle] = dmumu.state;
      mfp_states[particle] = mfp.state;
    }

    if ((step_index + 1) % sampling_stride == 0) {
      double dmumu_count = 0.0;
      double mfp_count = 0.0;
      for (std::size_t particle = 0; particle < particles; ++particle) {
        if (std::fabs(dmumu_states[particle].arcLengthM - observer_m) <=
            detector_half_width_m) ++dmumu_count;
        if (std::fabs(mfp_states[particle].arcLengthM - observer_m) <=
            detector_half_width_m) ++mfp_count;
      }
      dmumu_counts.push_back(dmumu_count);
      mfp_counts.push_back(mfp_count);
    }
  }

  long double dmumu_mu = 0.0L, mfp_mu = 0.0L;
  long double dmumu_mu2 = 0.0L, mfp_mu2 = 0.0L;
  long double dmumu_s2 = 0.0L, mfp_s2 = 0.0L;
  for (std::size_t particle = 0; particle < particles; ++particle) {
    dmumu_mu += dmumu_states[particle].mu;
    mfp_mu += mfp_states[particle].mu;
    dmumu_mu2 += dmumu_states[particle].mu * dmumu_states[particle].mu;
    mfp_mu2 += mfp_states[particle].mu * mfp_states[particle].mu;
    dmumu_s2 += dmumu_states[particle].arcLengthM *
                 dmumu_states[particle].arcLengthM;
    mfp_s2 += mfp_states[particle].arcLengthM *
              mfp_states[particle].arcLengthM;
  }
  const double inv_particles = 1.0 / static_cast<double>(particles);
  const double mean_mu_dmumu = static_cast<double>(dmumu_mu) * inv_particles;
  const double mean_mu_mfp = static_cast<double>(mfp_mu) * inv_particles;
  const double mean_mu2_dmumu = static_cast<double>(dmumu_mu2) * inv_particles;
  const double mean_mu2_mfp = static_cast<double>(mfp_mu2) * inv_particles;
  const double duration = steps * dt;
  const double kappa = speed * lambda / 3.0;
  const double expected_msd = 2.0 * kappa *
      (duration - (1.0 - std::exp(-rate * duration)) / rate);
  const double dmumu_msd_error = RelativeError(
      static_cast<double>(dmumu_s2) * inv_particles, expected_msd);
  const double mfp_msd_error = RelativeError(
      static_cast<double>(mfp_s2) * inv_particles, expected_msd);

  const double sampling_interval = sampling_stride * dt;
  const ProfileSummary dmumu_profile =
      SummarizeProfile(dmumu_counts, sampling_interval);
  const ProfileSummary mfp_profile =
      SummarizeProfile(mfp_counts, sampling_interval);
  const double onset_relative_difference = RelativeError(
      dmumu_profile.onset_s, mfp_profile.onset_s);
  const double peak_relative_difference = RelativeError(
      dmumu_profile.peak_s, mfp_profile.peak_s);
  const double fluence_relative_difference = RelativeError(
      dmumu_profile.fluence, mfp_profile.fluence);
  const double mu_mean_difference =
      std::fabs(mean_mu_dmumu - mean_mu_mfp);
  const double mu2_difference =
      std::fabs(mean_mu2_dmumu - mean_mu2_mfp);
  const bool profiles_finite = std::isfinite(onset_relative_difference) &&
      std::isfinite(peak_relative_difference) &&
      std::isfinite(fluence_relative_difference);

  // This is a diffusion-limit comparison, not an assertion that two distinct
  // collision operators have identical short-time distributions.  The 15%
  // MSD band includes finite-ensemble error and the explicitly measured
  // finite-dt bias of the reflected Ito path; it is paired with much tighter
  // pitch-moment checks so an incorrect scattering rate cannot pass merely by
  // broadening the spatial tolerance.  Onset permits a wider band because the
  // event-reset and continuous-diffusion operators are not equivalent before
  // the diffusive limit, while peak and fluence use the same fixed cadence
  // forward operator.
  const bool pass = kernels_ok && profiles_finite &&
      std::fabs(mean_mu_dmumu) <= 0.03 && std::fabs(mean_mu_mfp) <= 0.03 &&
      std::fabs(mean_mu2_dmumu - 1.0 / 3.0) <= 0.03 &&
      std::fabs(mean_mu2_mfp - 1.0 / 3.0) <= 0.03 &&
      mu_mean_difference <= 0.04 && mu2_difference <= 0.04 &&
      dmumu_msd_error <= 0.15 && mfp_msd_error <= 0.15 &&
      onset_relative_difference <= 0.50 &&
      peak_relative_difference <= 0.35 &&
      fluence_relative_difference <= 0.30;
  Result result = ValidationResult(
      pass,
      "matched Dmumu and MFP campaigns agree in the diffusion-limit observables",
      "cross-mover pitch-angle or time-profile metric exceeded tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("particles=4000");
  result.configuration.push_back("speed_m_per_s=1e7");
  result.configuration.push_back("lambda_parallel_m=1e8");
  result.configuration.push_back("nu_s^-1=0.1");
  result.configuration.push_back("duration_s=400;dt_s=0.25");
  result.configuration.push_back("detector_center_m=5e8;half_width_m=5e7");
  result.metrics.push_back(
      {"absolute_mean_mu_difference", mu_mean_difference, 0.04, "<=", "1"});
  result.metrics.push_back(
      {"absolute_mean_mu2_difference", mu2_difference, 0.04, "<=", "1"});
  result.metrics.push_back({"dmumu_msd_relative_error",
      dmumu_msd_error, 0.15, "<=", "dimensionless"});
  result.metrics.push_back({"mfp_msd_relative_error",
      mfp_msd_error, 0.15, "<=", "dimensionless"});
  result.metrics.push_back({"onset_relative_difference",
      onset_relative_difference, 0.50, "<=", "dimensionless"});
  result.metrics.push_back({"peak_time_relative_difference",
      peak_relative_difference, 0.35, "<=", "dimensionless"});
  result.metrics.push_back({"fluence_relative_difference",
      fluence_relative_difference, 0.30, "<=", "dimensionless"});
  return result;
}

double InitialPitchDensity(double mu, double slope) {
  return 0.5 * (1.0 + slope * mu);
}

double SampleInitialPitch(double slope, Transport::RandomStream& random) {
  // Rejection sampling is independent of the production transport solver.  A
  // linear, positive initial density has a known normalization and exercises a
  // nonzero l=1 anisotropy without singular endpoint data.
  for (;;) {
    const double mu = -1.0 + 2.0 * random.UniformOpen01();
    if (random.UniformOpen01() * (1.0 + slope) <= 1.0 + slope * mu)
      return mu;
  }
}

std::vector<double> IndependentPitchAngleFiniteVolume(
    std::size_t cells, double d0_per_s, double duration_s,
    double initial_slope) {
  const double dmu = 2.0 / static_cast<double>(cells);
  std::vector<double> density(cells, 0.0);
  for (std::size_t i = 0; i < cells; ++i) {
    const double mu = -1.0 + (i + 0.5) * dmu;
    density[i] = InitialPitchDensity(mu, initial_slope);
  }

  // This explicit conservative finite-volume solver shares no production
  // flux, index, reflection, or random-stream code.  Face diffusion vanishes
  // analytically at mu=+-1, giving the required zero-flux boundary.  The CFL
  // factor is fixed before comparison and the final fractional step lands
  // exactly on the requested physical time.
  const double stable_dt = 0.20 * dmu * dmu / d0_per_s;
  double elapsed = 0.0;
  std::vector<double> next(cells, 0.0);
  std::vector<double> flux(cells + 1, 0.0);
  while (elapsed < duration_s) {
    const double dt = std::min(stable_dt, duration_s - elapsed);
    flux.front() = 0.0;
    flux.back() = 0.0;
    for (std::size_t face = 1; face < cells; ++face) {
      const double mu_face = -1.0 + face * dmu;
      const double diffusion = d0_per_s * (1.0 - mu_face * mu_face);
      flux[face] = -diffusion *
          (density[face] - density[face - 1]) / dmu;
    }
    for (std::size_t i = 0; i < cells; ++i)
      next[i] = density[i] - dt * (flux[i + 1] - flux[i]) / dmu;
    density.swap(next);
    elapsed += dt;
  }
  return density;
}

Result RunVal03IndependentSolver() {
  const std::uint64_t seed = 1503001;
  const std::size_t particles = 30000;
  const std::size_t production_steps = 100;
  const std::size_t comparison_bins = 48;
  const std::size_t reference_cells = 192;
  const double d0 = 0.4;
  const double duration = 0.4;
  const double dt = duration / production_steps;
  const double initial_slope = 0.8;
  const double dmu_bin = 2.0 / comparison_bins;
  const Transport::ScalarResult momentum = Transport::MomentumFromSpeed(
      1.0e7, kProtonMassKg, kSpeedOfLightMPerS);
  IsotropicPitchAngleDiffusion provider(2.0 * d0);

  std::vector<double> pitch(particles, 0.0);
  bool kernel_ok = momentum.status.ok();
  for (std::size_t particle = 0; particle < particles; ++particle) {
    Transport::KeyedRandomStream initial_random(
        seed, static_cast<std::uint64_t>(particle), 0, 15030);
    pitch[particle] = SampleInitialPitch(initial_slope, initial_random);
  }
  for (std::size_t step_index = 0;
       step_index < production_steps && kernel_ok; ++step_index) {
    for (std::size_t particle = 0; particle < particles; ++particle) {
      Transport::KeyedRandomStream random(
          seed, static_cast<std::uint64_t>(particle), step_index, 15031);
      const Transport::FocusedTransportIncrement step =
          Transport::AdvanceFocusedTransportDmumu(
              Transport::FocusedTransportState(0.0, momentum.value,
                                                pitch[particle]),
              Transport::FocusedTransportBackground(0.0, 0.0, 0.0, 0.0),
              kProtonMassKg, kSpeedOfLightMPerS, dt, provider, random, NULL);
      if (!step.status.ok()) {
        kernel_ok = false;
        break;
      }
      pitch[particle] = step.state.mu;
    }
  }

  std::vector<double> model_density(comparison_bins, 0.0);
  for (std::size_t particle = 0; particle < particles; ++particle) {
    const double scaled = (pitch[particle] + 1.0) / 2.0 * comparison_bins;
    const std::size_t bin = std::min(
        comparison_bins - 1, static_cast<std::size_t>(scaled));
    model_density[bin] += 1.0 /
        (static_cast<double>(particles) * dmu_bin);
  }
  const std::vector<double> fine_reference =
      IndependentPitchAngleFiniteVolume(
          reference_cells, d0, duration, initial_slope);
  std::vector<double> reference_density(comparison_bins, 0.0);
  const std::size_t fine_per_bin = reference_cells / comparison_bins;
  for (std::size_t bin = 0; bin < comparison_bins; ++bin) {
    for (std::size_t fine = 0; fine < fine_per_bin; ++fine)
      reference_density[bin] +=
          fine_reference[bin * fine_per_bin + fine] / fine_per_bin;
  }

  double l1_error = 0.0;
  double model_mean = 0.0, reference_mean = 0.0;
  double model_second = 0.0, reference_second = 0.0;
  for (std::size_t bin = 0; bin < comparison_bins; ++bin) {
    const double mu = -1.0 + (bin + 0.5) * dmu_bin;
    l1_error += std::fabs(model_density[bin] - reference_density[bin]) *
                dmu_bin;
    model_mean += mu * model_density[bin] * dmu_bin;
    reference_mean += mu * reference_density[bin] * dmu_bin;
    model_second += mu * mu * model_density[bin] * dmu_bin;
    reference_second += mu * mu * reference_density[bin] * dmu_bin;
  }
  const double mean_error = std::fabs(model_mean - reference_mean);
  const double second_moment_error =
      std::fabs(model_second - reference_second);
  const bool pass = kernel_ok && l1_error <= 0.08 &&
      mean_error <= 0.015 && second_moment_error <= 0.015;
  Result result = ValidationResult(
      pass,
      "srcSEP Dmumu evolution agrees with an independent finite-volume solver",
      "independent focused-transport comparison exceeded tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("particles=30000");
  result.configuration.push_back("production_steps=100;duration_s=0.4");
  result.configuration.push_back("reference_cells=192;comparison_bins=48");
  result.configuration.push_back("Dmumu_s^-1=0.4*(1-mu^2)");
  result.configuration.push_back("initial_density=0.5*(1+0.8*mu)");
  result.metrics.push_back(
      {"pitch_density_L1_error", l1_error, 0.08, "<=", "dimensionless"});
  result.metrics.push_back(
      {"mean_mu_absolute_error", mean_error, 0.015, "<=", "dimensionless"});
  result.metrics.push_back({"mean_mu2_absolute_error",
      second_moment_error, 0.015, "<=", "dimensionless"});
  return result;
}

}  // namespace

std::vector<Descriptor> ScientificValidationDescriptors() {
  std::vector<Descriptor> descriptors;
  descriptors.push_back(MakeDescriptor(
      "VAL01", "Manufactured Parker benchmark",
      "Compare Parker diffusion moments and adiabatic cooling with independent analytical solutions.",
      RuntimeClass::Routine, "fixed keyed campaign seed 1501001",
      RunVal01ManufacturedParker));
  descriptors.push_back(MakeDescriptor(
      "VAL02", "Matched focused-mover campaign",
      "Compare Dmumu and MFP pitch moments, diffusion, onset, peak, and fluence under one lambda closure.",
      RuntimeClass::Extended, "fixed keyed campaign seed 1502001",
      RunVal02CrossMoverCampaign));
  descriptors.push_back(MakeDescriptor(
      "VAL03", "Independent focused-transport solver",
      "Compare the Dmumu ensemble with a separately implemented conservative finite-volume equation solver.",
      RuntimeClass::Extended, "fixed keyed campaign seed 1503001",
      RunVal03IndependentSolver));
  return descriptors;
}

}  // namespace Testing
}  // namespace SEP
