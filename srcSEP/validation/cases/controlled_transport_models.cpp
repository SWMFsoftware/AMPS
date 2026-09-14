#include "controlled_transport_models.h"

// Repository-relative headers are required because the enclosing AMPS build
// does not add srcSEP/util as a flat include directory for nested validation
// sources. Keeping these paths local also makes the focused compile gate match
// the production build that previously exposed this integration error.
#include "../../util/sep_focused_transport_core.h"
#include "../../util/sep_parker_core.h"
#include "../../util/sep_transport_common.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::Transport::AdvanceFocusedTransportDmumu;
using SEP::Transport::AdvanceParker;
using SEP::Transport::FocusedEquationMode;
using SEP::Transport::FocusedTransportBackground;
using SEP::Transport::FocusedTransportIncrement;
using SEP::Transport::FocusedTransportState;
using SEP::Transport::KeyedRandomStream;
using SEP::Transport::ParkerBackground;
using SEP::Transport::ParkerIncrement;
using SEP::Transport::ParkerState;
using SEP::Transport::PitchAngleDiffusionProvider;
using SEP::Transport::PitchAngleDiffusionSample;
using SEP::Transport::ScalarResult;
using SEP::Transport::SpatialDiffusionProvider;
using SEP::Transport::SpatialDiffusionSample;
using SEP::Transport::Status;

const double Pi = 3.141592653589793238462643383279502884;

std::map<std::string, std::string> ParseArguments(
    const std::vector<std::string>& arguments) {
  if (arguments.empty() || arguments.size() % 2 != 0)
    throw std::runtime_error("controlled model arguments must be name/value pairs");
  std::map<std::string, std::string> values;
  for (std::size_t i = 0; i < arguments.size(); i += 2) {
    const std::string& option = arguments[i];
    if (option.size() < 3 || option.substr(0, 2) != "--")
      throw std::runtime_error("controlled model option must begin with --");
    if (values.count(option.substr(2)) != 0)
      throw std::runtime_error("duplicate controlled model option: " + option);
    values[option.substr(2)] = arguments[i + 1];
  }
  return values;
}

std::string Require(const std::map<std::string, std::string>& values,
                    const std::string& name) {
  const std::map<std::string, std::string>::const_iterator found =
      values.find(name);
  if (found == values.end() || found->second.empty())
    throw std::runtime_error("missing controlled model option --" + name);
  return found->second;
}

double ParseDouble(const std::map<std::string, std::string>& values,
                   const std::string& name) {
  const std::string text = Require(values, name);
  char* end = NULL;
  errno = 0;
  const double value = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || !end || *end != '\0' || !std::isfinite(value))
    throw std::runtime_error("invalid finite number for --" + name);
  return value;
}

std::uint64_t ParseUnsigned(const std::map<std::string, std::string>& values,
                            const std::string& name) {
  const std::string text = Require(values, name);
  if (text.empty() || text[0] == '-')
    throw std::runtime_error("invalid unsigned integer for --" + name);
  char* end = NULL;
  errno = 0;
  const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || !end || *end != '\0' ||
      value > std::numeric_limits<std::uint64_t>::max())
    throw std::runtime_error("invalid unsigned integer for --" + name);
  return static_cast<std::uint64_t>(value);
}

std::vector<double> ParseDoubleList(
    const std::map<std::string, std::string>& values,
    const std::string& name) {
  std::vector<double> parsed;
  std::istringstream input(Require(values, name));
  std::string field;
  while (std::getline(input, field, ',')) {
    char* end = NULL;
    errno = 0;
    const double value = std::strtod(field.c_str(), &end);
    if (field.empty() || errno == ERANGE || !end || *end != '\0' ||
        !std::isfinite(value))
      throw std::runtime_error("invalid list element for --" + name);
    parsed.push_back(value);
  }
  if (parsed.empty()) throw std::runtime_error("empty list for --" + name);
  return parsed;
}

std::vector<std::uint64_t> ParseUnsignedList(
    const std::map<std::string, std::string>& values,
    const std::string& name) {
  std::vector<std::uint64_t> parsed;
  std::istringstream input(Require(values, name));
  std::string field;
  while (std::getline(input, field, ',')) {
    if (field.empty() || field[0] == '-')
      throw std::runtime_error("invalid list element for --" + name);
    char* end = NULL;
    errno = 0;
    const unsigned long long value = std::strtoull(field.c_str(), &end, 10);
    if (errno == ERANGE || !end || *end != '\0' || value == 0)
      throw std::runtime_error("invalid positive list element for --" + name);
    parsed.push_back(static_cast<std::uint64_t>(value));
  }
  if (parsed.empty()) throw std::runtime_error("empty list for --" + name);
  return parsed;
}

void Commit(std::ofstream* output, const std::string& temporaryPath,
            const std::string& outputPath) {
  output->close();
  if (output->fail())
    throw std::runtime_error("failed while flushing controlled model CSV");
  if (std::rename(temporaryPath.c_str(), outputPath.c_str()) != 0)
    throw std::runtime_error("cannot publish controlled model CSV: " + outputPath);
}

double WrapPeriodic(double position, double minimum, double maximum) {
  const double width = maximum - minimum;
  double wrapped = std::fmod(position - minimum, width);
  if (wrapped < 0.0) wrapped += width;
  return minimum + wrapped;
}

struct Moments {
  double mean = 0.0;
  double variance = 0.0;
  double skewness = 0.0;
  double kurtosis = 0.0;
};

Moments CalculateMoments(const std::vector<double>& values) {
  if (values.empty()) throw std::runtime_error("cannot measure an empty ensemble");
  Moments result;
  for (std::size_t i = 0; i < values.size(); ++i) result.mean += values[i];
  result.mean /= static_cast<double>(values.size());
  double third = 0.0;
  double fourth = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    const double delta = values[i] - result.mean;
    const double square = delta * delta;
    result.variance += square;
    third += square * delta;
    fourth += square * square;
  }
  result.variance /= static_cast<double>(values.size());
  if (result.variance > 0.0) {
    const double sigma = std::sqrt(result.variance);
    result.skewness = third / static_cast<double>(values.size()) /
                      (sigma * sigma * sigma);
    result.kurtosis = fourth / static_cast<double>(values.size()) /
                      (result.variance * result.variance);
  }
  return result;
}

class ConstantSpatialDiffusion final : public SpatialDiffusionProvider {
 public:
  explicit ConstantSpatialDiffusion(double kappa) : kappa_(kappa) {}
  SpatialDiffusionSample Evaluate(double, double) const override {
    SpatialDiffusionSample sample;
    sample.status = Status::Ok();
    sample.kappaParallelM2PerS = kappa_;
    sample.dKappaParallelDsMPerS = 0.0;
    sample.provenance = "validation:constant-kappa-linked-v1";
    return sample;
  }
 private:
  double kappa_;
};

class SinusoidalSpatialDiffusion final : public SpatialDiffusionProvider {
 public:
  SinusoidalSpatialDiffusion(double kappa0, double amplitude, double length,
                             bool numericalDerivative, double derivativeStep,
                             double derivativeSign)
      : kappa0_(kappa0), amplitude_(amplitude), length_(length),
        numericalDerivative_(numericalDerivative),
        derivativeStep_(derivativeStep), derivativeSign_(derivativeSign) {}

  SpatialDiffusionSample Evaluate(double position, double) const override {
    SpatialDiffusionSample sample;
    sample.status = Status::Ok();
    const double wrapped = WrapPeriodic(position, 0.0, length_);
    sample.kappaParallelM2PerS = Kappa(wrapped);
    if (numericalDerivative_) {
      // The centered periodic difference is intentionally independent of the
      // analytic derivative path while sampling the same smooth coefficient.
      const double plus = WrapPeriodic(wrapped + derivativeStep_, 0.0, length_);
      const double minus = WrapPeriodic(wrapped - derivativeStep_, 0.0, length_);
      sample.dKappaParallelDsMPerS = derivativeSign_ *
          (Kappa(plus) - Kappa(minus)) / (2.0 * derivativeStep_);
      sample.provenance = "validation:sinusoidal-kappa-numerical-gradient-v1";
    } else {
      sample.dKappaParallelDsMPerS = derivativeSign_ * kappa0_ * amplitude_ *
          (2.0 * Pi / length_) * std::cos(2.0 * Pi * wrapped / length_);
      sample.provenance = "validation:sinusoidal-kappa-analytic-gradient-v1";
    }
    return sample;
  }

 private:
  double Kappa(double position) const {
    return kappa0_ * (1.0 + amplitude_ *
        std::sin(2.0 * Pi * position / length_));
  }
  double kappa0_;
  double amplitude_;
  double length_;
  bool numericalDerivative_;
  double derivativeStep_;
  double derivativeSign_;
};

class ZeroSpatialDiffusion final : public SpatialDiffusionProvider {
 public:
  SpatialDiffusionSample Evaluate(double, double) const override {
    SpatialDiffusionSample sample;
    sample.status = Status::Ok();
    sample.kappaParallelM2PerS = 0.0;
    sample.dKappaParallelDsMPerS = 0.0;
    sample.provenance = "validation:zero-kappa-linked-v1";
    return sample;
  }
};

class ZeroPitchDiffusion final : public PitchAngleDiffusionProvider {
 public:
  PitchAngleDiffusionSample Evaluate(double, double, double) const override {
    PitchAngleDiffusionSample sample;
    sample.status = Status::Ok();
    sample.dMuMuPerS = 0.0;
    sample.dDmuMuDmuPerS = 0.0;
    sample.provenance = "validation:zero-dmumu-linked-v1";
    sample.turbulenceStateIdentity = "validation:frozen-controlled-background-v1";
    return sample;
  }
};

void RunCV02(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double kappa = ParseDouble(values, "kappa-m2-per-s");
  const double finalTime = ParseDouble(values, "final-time-s");
  const double minimum = ParseDouble(values, "minimum-m");
  const double maximum = ParseDouble(values, "maximum-m");
  const double center = ParseDouble(values, "center-m");
  const double profileHalfWidth = ParseDouble(values, "profile-half-width-m");
  const std::uint64_t baseSeed = ParseUnsigned(values, "campaign-seed");
  const std::uint64_t seedCount = ParseUnsigned(values, "seed-count");
  const std::uint64_t binCount = ParseUnsigned(values, "profile-bins");
  const std::uint64_t sampleCount = ParseUnsigned(values, "sample-count");
  const std::vector<double> timeSteps = ParseDoubleList(values, "time-steps-s");
  const std::vector<std::uint64_t> particleCounts =
      ParseUnsignedList(values, "particle-counts");
  if (!(kappa > 0.0) || !(minimum < center && center < maximum) ||
      !(profileHalfWidth > 0.0) || !(center - profileHalfWidth > minimum) ||
      !(center + profileHalfWidth < maximum) || binCount < 16 || seedCount < 1 ||
      sampleCount < 1)
    throw std::runtime_error("CV02 configuration violates its physical domain");

  const std::string temporaryPath = outputPath + ".tmp";
  std::ofstream output(temporaryPath.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good()) throw std::runtime_error("cannot create CV02 model CSV");
  output << std::setprecision(17)
         << "row_type,seed,particle_count,dt_s,time_s,bin_index,bin_left_m,"
            "bin_right_m,probability,mean_m,variance_m2,skewness,kurtosis,"
            "escaped_weight\n";
  ConstantSpatialDiffusion provider(kappa);
  const ParkerBackground background(0.0, 0.0);
  const double binWidth = 2.0 * profileHalfWidth /
                          static_cast<double>(binCount);

  for (std::size_t nIndex = 0; nIndex < particleCounts.size(); ++nIndex) {
    const std::uint64_t particleCount = particleCounts[nIndex];
    for (std::size_t dtIndex = 0; dtIndex < timeSteps.size(); ++dtIndex) {
      const double dt = timeSteps[dtIndex];
      const double stepsReal = finalTime / dt;
      const std::uint64_t steps = static_cast<std::uint64_t>(std::llround(stepsReal));
      if (!(dt > 0.0) || std::fabs(stepsReal - steps) > 1.0e-10 ||
          steps % sampleCount != 0)
        throw std::runtime_error("CV02 timestep/sample count is not integral");
      const std::uint64_t sampleEvery = steps / sampleCount;
      for (std::uint64_t seedOffset = 0; seedOffset < seedCount; ++seedOffset) {
        std::vector<double> positions(static_cast<std::size_t>(particleCount), center);
        for (std::uint64_t step = 0; step < steps; ++step) {
          for (std::uint64_t particle = 0; particle < particleCount; ++particle) {
            ParkerState state(positions[static_cast<std::size_t>(particle)], 1.0);
            KeyedRandomStream random(baseSeed + seedOffset, particle + 1, 202, step);
            const ParkerIncrement increment = AdvanceParker(
                state, background, 1.0, dt, provider, random);
            if (!increment.status.ok())
              throw std::runtime_error("CV02 production Parker step failed: " +
                                       increment.status.message);
            positions[static_cast<std::size_t>(particle)] = increment.state.arcLengthM;
          }
          if ((step + 1) % sampleEvery == 0) {
            const double time = static_cast<double>(step + 1) * dt;
            const Moments moments = CalculateMoments(positions);
            std::size_t escaped = 0;
            for (std::size_t i = 0; i < positions.size(); ++i)
              if (positions[i] <= minimum || positions[i] >= maximum) ++escaped;
            output << "moment," << baseSeed + seedOffset << ',' << particleCount
                   << ',' << dt << ',' << time << ",,,,," << moments.mean << ','
                   << moments.variance << ',' << moments.skewness << ','
                   << moments.kurtosis << ','
                   << static_cast<double>(escaped) / particleCount << '\n';
          }
        }

        std::vector<std::uint64_t> bins(static_cast<std::size_t>(binCount), 0);
        for (std::size_t i = 0; i < positions.size(); ++i) {
          const long index = static_cast<long>(
              std::floor((positions[i] - (center - profileHalfWidth)) / binWidth));
          if (index >= 0 && index < static_cast<long>(binCount))
            ++bins[static_cast<std::size_t>(index)];
        }
        for (std::uint64_t bin = 0; bin < binCount; ++bin) {
          const double left = center - profileHalfWidth + bin * binWidth;
          output << "profile," << baseSeed + seedOffset << ',' << particleCount
                 << ',' << dt << ',' << finalTime << ',' << bin << ',' << left
                 << ',' << left + binWidth << ','
                 << static_cast<double>(bins[static_cast<std::size_t>(bin)]) /
                    particleCount << ",,,,,\n";
        }
      }
    }
  }
  Commit(&output, temporaryPath, outputPath);
}

void WriteCV03Profile(std::ofstream* output, const std::string& mode,
                      const std::string& derivativePath, double driftSign,
                      std::uint64_t seed, std::uint64_t particleCount,
                      double dt, double finalTime, double length,
                      const std::vector<double>& positions,
                      std::uint64_t binCount) {
  std::vector<std::uint64_t> bins(static_cast<std::size_t>(binCount), 0);
  for (std::size_t i = 0; i < positions.size(); ++i) {
    long bin = static_cast<long>(std::floor(positions[i] / length * binCount));
    if (bin == static_cast<long>(binCount)) bin = 0;
    if (bin < 0 || bin >= static_cast<long>(binCount))
      throw std::runtime_error("CV03 periodic coordinate escaped its domain");
    ++bins[static_cast<std::size_t>(bin)];
  }
  const Moments moments = CalculateMoments(positions);
  for (std::uint64_t bin = 0; bin < binCount; ++bin) {
    const double left = length * bin / static_cast<double>(binCount);
    *output << mode << ',' << derivativePath << ',' << driftSign << ',' << seed
            << ',' << particleCount << ',' << dt << ',' << finalTime << ','
            << bin << ',' << left << ','
            << length * (bin + 1) / static_cast<double>(binCount) << ','
            << static_cast<double>(bins[static_cast<std::size_t>(bin)]) /
               particleCount << ',' << moments.mean << ',' << moments.variance
            << ",0\n";
  }
}

void RunCV03(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double length = ParseDouble(values, "line-length-m");
  const double kappa0 = ParseDouble(values, "kappa0-m2-per-s");
  const double amplitude = ParseDouble(values, "amplitude");
  const double derivativeStep = ParseDouble(values, "derivative-step-m");
  const double finalTime = ParseDouble(values, "final-time-s");
  const double center = ParseDouble(values, "center-m");
  const std::uint64_t particleCount = ParseUnsigned(values, "particle-count");
  const std::uint64_t seedCount = ParseUnsigned(values, "seed-count");
  const std::uint64_t baseSeed = ParseUnsigned(values, "campaign-seed");
  const std::uint64_t binCount = ParseUnsigned(values, "profile-bins");
  const std::vector<double> timeSteps = ParseDoubleList(values, "time-steps-s");
  if (!(length > 0.0) || !(kappa0 > 0.0) || !(amplitude > 0.0 && amplitude <= 0.5) ||
      !(derivativeStep > 0.0 && derivativeStep < 0.1 * length) ||
      !(center >= 0.0 && center < length) || !(finalTime > 0.0) ||
      particleCount < 100 || seedCount < 2 || binCount < 16)
    throw std::runtime_error("CV03 configuration violates its physical domain");

  const std::string temporaryPath = outputPath + ".tmp";
  std::ofstream output(temporaryPath.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good()) throw std::runtime_error("cannot create CV03 model CSV");
  output << std::setprecision(17)
         << "mode,derivative_path,drift_sign,seed,particle_count,dt_s,time_s,"
            "bin_index,bin_left_m,bin_right_m,probability,mean_m,variance_m2,"
            "escaped_weight\n";
  const ParkerBackground background(0.0, 0.0);

  for (int derivativeMode = 0; derivativeMode < 2; ++derivativeMode) {
    const std::string derivativeName = derivativeMode == 0 ? "analytic" : "numerical";
    for (std::size_t dtIndex = 0; dtIndex < timeSteps.size(); ++dtIndex) {
      const double dt = timeSteps[dtIndex];
      const std::uint64_t steps = static_cast<std::uint64_t>(std::llround(finalTime / dt));
      if (!(dt > 0.0) || std::fabs(finalTime / dt - steps) > 1.0e-10)
        throw std::runtime_error("CV03 final time must be divisible by every dt");
      for (std::uint64_t seedOffset = 0; seedOffset < seedCount; ++seedOffset) {
        std::vector<double> positions(static_cast<std::size_t>(particleCount), center);
        SinusoidalSpatialDiffusion provider(kappa0, amplitude, length,
            derivativeMode != 0, derivativeStep, 1.0);
        for (std::uint64_t step = 0; step < steps; ++step) {
          for (std::uint64_t particle = 0; particle < particleCount; ++particle) {
            ParkerState state(positions[static_cast<std::size_t>(particle)], 1.0);
            // Use common random numbers for the analytic- versus numerical-
            // derivative comparison.  Only d(kappa)/ds should differ between
            // these paths; changing the random stream would add avoidable
            // Monte Carlo noise and weaken that focused diagnostic.
            KeyedRandomStream random(baseSeed + seedOffset, particle + 1,
                                     303, step);
            const ParkerIncrement increment = AdvanceParker(
                state, background, 1.0, dt, provider, random);
            if (!increment.status.ok())
              throw std::runtime_error("CV03 transient Parker step failed: " +
                                       increment.status.message);
            positions[static_cast<std::size_t>(particle)] = WrapPeriodic(
                increment.state.arcLengthM, 0.0, length);
          }
        }
        WriteCV03Profile(&output, "transient", derivativeName, 1.0,
            baseSeed + seedOffset, particleCount, dt, finalTime, length,
            positions, binCount);
      }
    }
  }

  // A uniform initial density is the exact zero-flux equilibrium of
  // partial_t f = partial_s(kappa partial_s f). Stratified initial positions
  // remove avoidable t=0 bin noise while all subsequent motion remains the
  // production stochastic Parker update.
  const double finestDt = *std::min_element(timeSteps.begin(), timeSteps.end());
  const std::uint64_t finestSteps =
      static_cast<std::uint64_t>(std::llround(finalTime / finestDt));
  for (std::uint64_t seedOffset = 0; seedOffset < seedCount; ++seedOffset) {
    std::vector<double> positions(static_cast<std::size_t>(particleCount));
    for (std::uint64_t particle = 0; particle < particleCount; ++particle)
      positions[static_cast<std::size_t>(particle)] =
          length * (particle + 0.5) / static_cast<double>(particleCount);
    SinusoidalSpatialDiffusion provider(kappa0, amplitude, length, false,
                                       derivativeStep, 1.0);
    for (std::uint64_t step = 0; step < finestSteps; ++step) {
      for (std::uint64_t particle = 0; particle < particleCount; ++particle) {
        ParkerState state(positions[static_cast<std::size_t>(particle)], 1.0);
        KeyedRandomStream random(baseSeed + 1000 + seedOffset, particle + 1,
                                 305, step);
        const ParkerIncrement increment = AdvanceParker(
            state, background, 1.0, finestDt, provider, random);
        if (!increment.status.ok())
          throw std::runtime_error("CV03 equilibrium Parker step failed: " +
                                   increment.status.message);
        positions[static_cast<std::size_t>(particle)] = WrapPeriodic(
            increment.state.arcLengthM, 0.0, length);
      }
    }
    WriteCV03Profile(&output, "equilibrium", "analytic", 1.0,
        baseSeed + 1000 + seedOffset, particleCount, finestDt, finalTime,
        length, positions, binCount);
  }

  // The sign-reversed Ito drift is saved as a negative control. The scorer
  // must demonstrate that it disagrees with the conservative finite-volume
  // reference by substantially more than the nominal solution.
  std::vector<double> negativePositions(static_cast<std::size_t>(particleCount), center);
  SinusoidalSpatialDiffusion negativeProvider(kappa0, amplitude, length, false,
                                              derivativeStep, -1.0);
  for (std::uint64_t step = 0; step < finestSteps; ++step) {
    for (std::uint64_t particle = 0; particle < particleCount; ++particle) {
      ParkerState state(negativePositions[static_cast<std::size_t>(particle)], 1.0);
      KeyedRandomStream random(baseSeed, particle + 1, 306, step);
      const ParkerIncrement increment = AdvanceParker(
          state, background, 1.0, finestDt, negativeProvider, random);
      if (!increment.status.ok())
        throw std::runtime_error("CV03 negative-control Parker step failed: " +
                                 increment.status.message);
      negativePositions[static_cast<std::size_t>(particle)] = WrapPeriodic(
          increment.state.arcLengthM, 0.0, length);
    }
  }
  WriteCV03Profile(&output, "negative-control", "analytic", -1.0, baseSeed,
      particleCount, finestDt, finalTime, length, negativePositions, binCount);
  Commit(&output, temporaryPath, outputPath);
}

double SpeedFromEnergyPerNucleon(double energyMeVPerNucleon,
                                 double nucleonMassKg, double lightSpeed) {
  const double electronVoltJ = 1.602176634e-19;
  const double kineticJ = energyMeVPerNucleon * 1.0e6 * electronVoltJ;
  const double gamma = 1.0 + kineticJ / (nucleonMassKg * lightSpeed * lightSpeed);
  return lightSpeed * std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

double EnergyPerNucleonFromMomentum(double momentum, double totalMass,
                                    double atomicMassNumber,
                                    double lightSpeed) {
  const double electronVoltJ = 1.602176634e-19;
  const double ratio = momentum / (totalMass * lightSpeed);
  // Rationalize sqrt(1+x^2)-1.  The direct subtraction loses significant
  // digits for the 0.1 MeV/nucleon case even though momentum itself is valid;
  // x^2/(sqrt(1+x^2)+1) is algebraically identical and remains accurate from
  // nonrelativistic through relativistic energies.
  const double kineticJ = ratio * ratio /
      (std::sqrt(1.0 + ratio * ratio) + 1.0) *
      totalMass * lightSpeed * lightSpeed;
  return kineticJ / atomicMassNumber / (1.0e6 * electronVoltJ);
}

void RunCV04(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double lightSpeed = ParseDouble(values, "light-speed-m-per-s");
  const double protonMass = ParseDouble(values, "proton-mass-kg");
  const double constantDivergence = ParseDouble(values, "constant-divergence-per-s");
  const double constantFinalTime = ParseDouble(values, "constant-final-time-s");
  const double sphericalFinalTime = ParseDouble(values, "spherical-final-time-s");
  const double radialSpeed = ParseDouble(values, "radial-wind-speed-m-per-s");
  const double initialRadius = ParseDouble(values, "initial-radius-m");
  const double spectralIndex = ParseDouble(values, "spectral-index");
  const std::uint64_t campaignSeed = ParseUnsigned(values, "campaign-seed");
  const std::vector<double> constantSteps =
      ParseDoubleList(values, "constant-time-steps-s");
  const std::vector<double> sphericalSteps =
      ParseDoubleList(values, "spherical-time-steps-s");
  const std::vector<double> energies =
      ParseDoubleList(values, "energies-mev-per-nucleon");
  const std::vector<double> massNumbers =
      ParseDoubleList(values, "atomic-mass-numbers");
  const std::uint64_t sampleCount = ParseUnsigned(values, "sample-count");
  if (!(lightSpeed > 0.0) || !(protonMass > 0.0) ||
      !(constantDivergence > 0.0) || !(radialSpeed > 0.0) ||
      !(initialRadius > 0.0) || !(spectralIndex > 0.0) || sampleCount < 2)
    throw std::runtime_error("CV04 configuration violates its physical domain");

  const std::string temporaryPath = outputPath + ".tmp";
  std::ofstream output(temporaryPath.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good()) throw std::runtime_error("cannot create CV04 model CSV");
  output << std::setprecision(17)
         << "scenario,atomic_mass_number,energy0_mev_per_nucleon,dt_s,time_s,"
            "radius_m,momentum_kg_m_per_s,kinetic_energy_mev_per_nucleon,weight\n";
  ZeroSpatialDiffusion provider;

  for (std::size_t species = 0; species < massNumbers.size(); ++species) {
    const double atomicMass = massNumbers[species];
    if (!(atomicMass >= 1.0)) throw std::runtime_error("CV04 mass number must be >=1");
    const double totalMass = atomicMass * protonMass;
    for (std::size_t energyIndex = 0; energyIndex < energies.size(); ++energyIndex) {
      const double speed = SpeedFromEnergyPerNucleon(
          energies[energyIndex], protonMass, lightSpeed);
      const ScalarResult momentum = SEP::Transport::MomentumFromSpeed(
          speed, totalMass, lightSpeed);
      if (!momentum.status.ok()) throw std::runtime_error(momentum.status.message);
      // A power law in momentum remains a power law under the multiplicative
      // adiabatic characteristic. Scaling by 1e-19 kg m/s keeps weights near
      // unity without changing the fitted logarithmic slope.
      const double weight = std::pow(momentum.value / 1.0e-19, -spectralIndex);

      for (std::size_t dtIndex = 0; dtIndex < constantSteps.size(); ++dtIndex) {
        const double dt = constantSteps[dtIndex];
        const std::uint64_t steps =
            static_cast<std::uint64_t>(std::llround(constantFinalTime / dt));
        if (!(dt > 0.0) || std::fabs(constantFinalTime / dt - steps) > 1.0e-10 ||
            steps % sampleCount != 0)
          throw std::runtime_error("CV04 constant-divergence cadence is not integral");
        const std::uint64_t sampleEvery = steps / sampleCount;
        ParkerState state(initialRadius, momentum.value);
        output << "constant," << atomicMass << ',' << energies[energyIndex] << ','
               << dt << ",0," << state.arcLengthM << ',' << state.momentumKgMPerS
               << ',' << energies[energyIndex] << ',' << weight << '\n';
        for (std::uint64_t step = 0; step < steps; ++step) {
          KeyedRandomStream random(campaignSeed, energyIndex + 1,
                                   404 + species, step);
          const ParkerIncrement increment = AdvanceParker(
              state, ParkerBackground(0.0, constantDivergence), speed, dt,
              provider, random);
          if (!increment.status.ok())
            throw std::runtime_error("CV04 constant Parker step failed: " +
                                     increment.status.message);
          state = increment.state;
          if ((step + 1) % sampleEvery == 0)
            output << "constant," << atomicMass << ',' << energies[energyIndex]
                   << ',' << dt << ',' << (step + 1) * dt << ','
                   << state.arcLengthM << ',' << state.momentumKgMPerS << ','
                   << EnergyPerNucleonFromMomentum(state.momentumKgMPerS,
                                                   totalMass, atomicMass,
                                                   lightSpeed) << ',' << weight << '\n';
        }
      }

      for (std::size_t dtIndex = 0; dtIndex < sphericalSteps.size(); ++dtIndex) {
        const double dt = sphericalSteps[dtIndex];
        const std::uint64_t steps =
            static_cast<std::uint64_t>(std::llround(sphericalFinalTime / dt));
        if (!(dt > 0.0) || std::fabs(sphericalFinalTime / dt - steps) > 1.0e-10 ||
            steps % sampleCount != 0)
          throw std::runtime_error("CV04 spherical-wind cadence is not integral");
        const std::uint64_t sampleEvery = steps / sampleCount;
        ParkerState state(initialRadius, momentum.value);
        output << "spherical," << atomicMass << ',' << energies[energyIndex] << ','
               << dt << ",0," << state.arcLengthM << ',' << state.momentumKgMPerS
               << ',' << energies[energyIndex] << ',' << weight << '\n';
        for (std::uint64_t step = 0; step < steps; ++step) {
          // Evaluate div(U)=2U/r at the deterministic radial midpoint. This is
          // the declared second-order coefficient location; AdvanceParker then
          // applies its exact frozen-divergence exponential over the step.
          const double midpointRadius = state.arcLengthM + 0.5 * radialSpeed * dt;
          KeyedRandomStream random(campaignSeed, energyIndex + 1,
                                   405 + species, step);
          const ParkerIncrement increment = AdvanceParker(
              state, ParkerBackground(radialSpeed,
                                      2.0 * radialSpeed / midpointRadius),
              speed, dt, provider, random);
          if (!increment.status.ok())
            throw std::runtime_error("CV04 spherical Parker step failed: " +
                                     increment.status.message);
          state = increment.state;
          if ((step + 1) % sampleEvery == 0)
            output << "spherical," << atomicMass << ',' << energies[energyIndex]
                   << ',' << dt << ',' << (step + 1) * dt << ','
                   << state.arcLengthM << ',' << state.momentumKgMPerS << ','
                   << EnergyPerNucleonFromMomentum(state.momentumKgMPerS,
                                                   totalMass, atomicMass,
                                                   lightSpeed) << ',' << weight << '\n';
        }
      }
    }
  }
  Commit(&output, temporaryPath, outputPath);
}

void RunCV05(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double speed = ParseDouble(values, "speed-m-per-s");
  const double protonMass = ParseDouble(values, "proton-mass-kg");
  const double lightSpeed = ParseDouble(values, "light-speed-m-per-s");
  const double focusingLength = ParseDouble(values, "focusing-length-m");
  const double initialPosition = ParseDouble(values, "initial-position-m");
  const double finalTime = ParseDouble(values, "final-time-s");
  const std::vector<double> timeSteps = ParseDoubleList(values, "time-steps-s");
  const std::vector<double> pitchAngles = ParseDoubleList(values, "initial-mus");
  const std::vector<double> gradientSigns = ParseDoubleList(values, "gradient-signs");
  const std::uint64_t sampleCount = ParseUnsigned(values, "sample-count");
  const std::uint64_t seed = ParseUnsigned(values, "campaign-seed");
  if (!(speed > 0.0 && speed < lightSpeed) || !(protonMass > 0.0) ||
      !(focusingLength > 0.0) || !(finalTime > 0.0) || sampleCount < 2)
    throw std::runtime_error("CV05 configuration violates its physical domain");
  const ScalarResult momentum = SEP::Transport::MomentumFromSpeed(
      speed, protonMass, lightSpeed);
  if (!momentum.status.ok()) throw std::runtime_error(momentum.status.message);

  const std::string temporaryPath = outputPath + ".tmp";
  std::ofstream output(temporaryPath.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good()) throw std::runtime_error("cannot create CV05 model CSV");
  output << std::setprecision(17)
         << "gradient_sign,dt_s,particle_id,initial_mu,weight,time_s,position_m,"
            "mu,momentum_kg_m_per_s,magnetic_moment_invariant\n";
  ZeroPitchDiffusion provider;

  for (std::size_t signIndex = 0; signIndex < gradientSigns.size(); ++signIndex) {
    const double sign = gradientSigns[signIndex];
    if (sign != -1.0 && sign != 1.0)
      throw std::runtime_error("CV05 gradient signs must be exactly -1 or +1");
    const double dLnBds = sign / focusingLength;
    for (std::size_t dtIndex = 0; dtIndex < timeSteps.size(); ++dtIndex) {
      const double dt = timeSteps[dtIndex];
      const std::uint64_t steps =
          static_cast<std::uint64_t>(std::llround(finalTime / dt));
      if (!(dt > 0.0) || std::fabs(finalTime / dt - steps) > 1.0e-10 ||
          steps % sampleCount != 0)
        throw std::runtime_error("CV05 cadence is not integral");
      const std::uint64_t sampleEvery = steps / sampleCount;
      for (std::size_t muIndex = 0; muIndex < pitchAngles.size(); ++muIndex) {
        const double initialMu = pitchAngles[muIndex];
        if (initialMu < -1.0 || initialMu > 1.0)
          throw std::runtime_error("CV05 initial pitch angle is outside [-1,1]");
        // Smooth positive quadrature weights represent f(mu)=1+0.4mu. They
        // are retained exactly, allowing the external scorer to compare angular
        // moments and weight conservation without stochastic sampling noise.
        const double weight = 1.0 + 0.4 * initialMu;
        FocusedTransportState state(initialPosition, momentum.value, initialMu);
        const double initialInvariant =
            (1.0 - initialMu * initialMu) /
            std::exp(dLnBds * initialPosition);
        output << sign << ',' << dt << ',' << muIndex + 1 << ',' << initialMu
               << ',' << weight << ",0," << state.arcLengthM << ',' << state.mu
               << ',' << state.momentumKgMPerS << ',' << initialInvariant << '\n';
        for (std::uint64_t step = 0; step < steps; ++step) {
          FocusedTransportBackground background;
          background.dLnAbsBdsPerM = dLnBds;
          background.plasmaAdvectionMPerS = 0.0;
          background.parallelVelocityGradientPerS = 0.0;
          background.velocityDivergencePerS = 0.0;
          background.fieldAlignedStrainPerS = 0.0;
          background.equationMode = FocusedEquationMode::FullGyrotropic;
          KeyedRandomStream random(seed, muIndex + 1, 505 + signIndex, step);
          const FocusedTransportIncrement increment = AdvanceFocusedTransportDmumu(
              state, background, protonMass, lightSpeed, dt, provider, random, NULL);
          if (!increment.status.ok())
            throw std::runtime_error("CV05 production focused step failed: " +
                                     increment.status.message);
          state = increment.state;
          if ((step + 1) % sampleEvery == 0) {
            const double invariant = (1.0 - state.mu * state.mu) /
                std::exp(dLnBds * state.arcLengthM);
            output << sign << ',' << dt << ',' << muIndex + 1 << ',' << initialMu
                   << ',' << weight << ',' << (step + 1) * dt << ','
                   << state.arcLengthM << ',' << state.mu << ','
                   << state.momentumKgMPerS << ',' << invariant << '\n';
          }
        }
      }
    }
  }
  Commit(&output, temporaryPath, outputPath);
}

}  // namespace

namespace SEP {
namespace Validation {

bool RunControlledTransportModel(const std::string& caseId,
                                 const std::vector<std::string>& arguments,
                                 const std::string& outputPath,
                                 std::string* error) {
  try {
    const std::map<std::string, std::string> values = ParseArguments(arguments);
    if (caseId == "CV02") RunCV02(values, outputPath);
    else if (caseId == "CV03") RunCV03(values, outputPath);
    else if (caseId == "CV04") RunCV04(values, outputPath);
    else if (caseId == "CV05") RunCV05(values, outputPath);
    else throw std::runtime_error("unsupported controlled transport case: " + caseId);
    if (error) error->clear();
    return true;
  } catch (const std::exception& exception) {
    if (error) *error = exception.what();
    return false;
  }
}

}  // namespace Validation
}  // namespace SEP

#ifdef SRCSEP_CONTROLLED_MODELS_STANDALONE_TEST_HARNESS
// Dependency-light harness for strict compile/sanitizer testing only. Physics
// validation never executes this binary; it requires the linked srcSEP/AMPS
// registry and passes the same argument vector through its execution context.
int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: controlled_models CASE OUTPUT --name value ...\n";
    return 2;
  }
  std::vector<std::string> arguments;
  for (int i = 3; i < argc; ++i) arguments.push_back(argv[i]);
  std::string error;
  if (!SEP::Validation::RunControlledTransportModel(
          argv[1], arguments, argv[2], &error)) {
    std::cerr << argv[1] << " model error: " << error << '\n';
    return 2;
  }
  return 0;
}
#endif
