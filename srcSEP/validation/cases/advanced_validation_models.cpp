#include "advanced_validation_models.h"

// Repository-relative includes deliberately reproduce the enclosing AMPS
// build contract. Nested validation sources must not depend on an unconfigured
// srcSEP/util include path.
#include "../../util/sep_focused_transport_core.h"
#include "../../util/sep_focused_transport_mfp_core.h"
#include "../../util/sep_parker_core.h"
#include "../../util/sep_shock_source_core.h"
#include "../../util/sep_transport_common.h"
#include "../../util/sep_turbulence_core.h"

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

using SEP::Transport::KeyedRandomStream;
using SEP::Transport::Status;

const double Pi = 3.141592653589793238462643383279502884;
const double ProtonMassKg = 1.67262192369e-27;
const double LightSpeedMPerS = 299792458.0;

std::map<std::string, std::string> ParseArguments(
    const std::vector<std::string>& arguments) {
  if (arguments.empty() || arguments.size() % 2 != 0)
    throw std::runtime_error("advanced validation arguments must be name/value pairs");
  std::map<std::string, std::string> values;
  for (std::size_t i = 0; i < arguments.size(); i += 2) {
    if (arguments[i].size() < 3 || arguments[i].substr(0, 2) != "--")
      throw std::runtime_error("advanced validation option must begin with --");
    const std::string name = arguments[i].substr(2);
    if (values.count(name))
      throw std::runtime_error("duplicate advanced validation option --" + name);
    values[name] = arguments[i + 1];
  }
  return values;
}

std::string Require(const std::map<std::string, std::string>& values,
                    const std::string& name) {
  const std::map<std::string, std::string>::const_iterator found = values.find(name);
  if (found == values.end() || found->second.empty())
    throw std::runtime_error("missing advanced validation option --" + name);
  return found->second;
}

double Number(const std::map<std::string, std::string>& values,
              const std::string& name) {
  const std::string text = Require(values, name);
  char* end = NULL;
  errno = 0;
  const double result = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || !end || *end != '\0' || !std::isfinite(result))
    throw std::runtime_error("invalid finite number for --" + name);
  return result;
}

std::uint64_t Unsigned(const std::map<std::string, std::string>& values,
                       const std::string& name) {
  const std::string text = Require(values, name);
  if (text.empty() || text[0] == '-')
    throw std::runtime_error("invalid unsigned integer for --" + name);
  char* end = NULL;
  errno = 0;
  const unsigned long long result = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || !end || *end != '\0')
    throw std::runtime_error("invalid unsigned integer for --" + name);
  return static_cast<std::uint64_t>(result);
}

std::vector<double> NumberList(const std::map<std::string, std::string>& values,
                               const std::string& name) {
  std::vector<double> result;
  std::istringstream input(Require(values, name));
  std::string item;
  while (std::getline(input, item, ',')) {
    char* end = NULL;
    errno = 0;
    const double value = std::strtod(item.c_str(), &end);
    if (item.empty() || errno == ERANGE || !end || *end != '\0' ||
        !std::isfinite(value))
      throw std::runtime_error("invalid list element for --" + name);
    result.push_back(value);
  }
  if (result.empty()) throw std::runtime_error("empty list for --" + name);
  return result;
}

std::vector<std::uint64_t> UnsignedList(
    const std::map<std::string, std::string>& values, const std::string& name) {
  std::vector<std::uint64_t> result;
  std::istringstream input(Require(values, name));
  std::string item;
  while (std::getline(input, item, ',')) {
    if (item.empty() || item[0] == '-')
      throw std::runtime_error("invalid list element for --" + name);
    char* end = NULL;
    errno = 0;
    const unsigned long long value = std::strtoull(item.c_str(), &end, 10);
    if (errno == ERANGE || !end || *end != '\0')
      throw std::runtime_error("invalid list element for --" + name);
    result.push_back(static_cast<std::uint64_t>(value));
  }
  if (result.empty()) throw std::runtime_error("empty list for --" + name);
  return result;
}

void Commit(std::ofstream* output, const std::string& temporary,
            const std::string& destination) {
  output->close();
  if (output->fail())
    throw std::runtime_error("failed while flushing advanced validation CSV");
  if (std::rename(temporary.c_str(), destination.c_str()) != 0)
    throw std::runtime_error("cannot publish advanced validation CSV");
}

double Legendre(unsigned mode, double mu) {
  if (mode == 0) return 1.0;
  if (mode == 1) return mu;
  double previous = 1.0;
  double current = mu;
  for (unsigned ell = 2; ell <= mode; ++ell) {
    const double next = ((2.0 * ell - 1.0) * mu * current -
                         (ell - 1.0) * previous) / ell;
    previous = current;
    current = next;
  }
  return current;
}

class LegendreDiffusionProvider final
    : public SEP::Transport::PitchAngleDiffusionProvider {
 public:
  explicit LegendreDiffusionProvider(double d0) : d0_(d0) {}
  SEP::Transport::PitchAngleDiffusionSample Evaluate(
      double, double, double mu) const override {
    SEP::Transport::PitchAngleDiffusionSample sample;
    sample.status = Status::Ok();
    sample.dMuMuPerS = d0_ * std::max(0.0, 1.0 - mu * mu);
    sample.dDmuMuDmuPerS = -2.0 * d0_ * mu;
    sample.provenance = "validation:CV06:legendre-diffusion-v1";
    sample.turbulenceStateIdentity = "validation:CV06:frozen-spectrum-v1";
    return sample;
  }
 private:
  double d0_;
};

double SampleLegendre(unsigned mode, double epsilon,
                      KeyedRandomStream* random) {
  // |P_l(mu)|<=1, hence 1+epsilon*P_l is nonnegative and bounded by
  // 1+epsilon for 0<epsilon<1. Rejection sampling is therefore exact and
  // independent of the production scattering update.
  for (;;) {
    const double mu = 2.0 * random->UniformOpen01() - 1.0;
    if (random->UniformOpen01() <=
        (1.0 + epsilon * Legendre(mode, mu)) / (1.0 + epsilon)) return mu;
  }
}

void RunCV06(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double d0 = Number(values, "d0-per-s");
  const double epsilon = Number(values, "epsilon");
  const double finalTime = Number(values, "final-time-s");
  const std::vector<double> timeSteps = NumberList(values, "time-steps-s");
  const std::vector<std::uint64_t> modes = UnsignedList(values, "modes");
  const std::uint64_t samples = Unsigned(values, "sample-count");
  const std::uint64_t particles = Unsigned(values, "particle-count");
  const std::uint64_t seeds = Unsigned(values, "seed-count");
  const std::uint64_t baseSeed = Unsigned(values, "campaign-seed");
  if (!(d0 > 0.0) || !(epsilon > 0.0 && epsilon < 1.0) ||
      !(finalTime > 0.0) || particles < 100 || seeds < 1 || samples < 1)
    throw std::runtime_error("CV06 configuration violates its physical domain");
  const SEP::Transport::ScalarResult momentum = SEP::Transport::MomentumFromSpeed(
      1.0e6, ProtonMassKg, LightSpeedMPerS);
  if (!momentum.status.ok()) throw std::runtime_error(momentum.status.message);
  LegendreDiffusionProvider provider(d0);
  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  if (!output.good()) throw std::runtime_error("cannot create CV06 model CSV");
  output << std::setprecision(17)
         << "initial_mode,measured_mode,seed,dt_s,time_s,coefficient,"
            "boundary_reflections,particle_count\n";
  for (std::size_t modeIndex = 0; modeIndex < modes.size(); ++modeIndex) {
    const unsigned initialMode = static_cast<unsigned>(modes[modeIndex]);
    if (initialMode < 1 || initialMode > 6)
      throw std::runtime_error("CV06 modes must lie in [1,6]");
    for (std::size_t dtIndex = 0; dtIndex < timeSteps.size(); ++dtIndex) {
      const double dt = timeSteps[dtIndex];
      const std::uint64_t steps = static_cast<std::uint64_t>(std::llround(finalTime / dt));
      if (!(dt > 0.0) || std::fabs(finalTime / dt - steps) > 1.0e-10 ||
          steps % samples != 0)
        throw std::runtime_error("CV06 cadence must be integral for every dt");
      const std::uint64_t sampleEvery = steps / samples;
      for (std::uint64_t seedOffset = 0; seedOffset < seeds; ++seedOffset) {
        std::vector<double> mu(static_cast<std::size_t>(particles));
        for (std::uint64_t particle = 0; particle < particles; ++particle) {
          KeyedRandomStream initial(baseSeed + seedOffset, particle + 1,
                                    6060 + initialMode, 0);
          mu[static_cast<std::size_t>(particle)] =
              SampleLegendre(initialMode, epsilon, &initial);
        }
        std::uint64_t reflections = 0;
        for (std::uint64_t step = 0; step < steps; ++step) {
          for (std::uint64_t particle = 0; particle < particles; ++particle) {
            KeyedRandomStream random(baseSeed + seedOffset, particle + 1,
                                     6061 + initialMode, step);
            const SEP::Transport::FocusedTransportIncrement increment =
                SEP::Transport::AdvanceFocusedTransportDmumu(
                    {0.0, momentum.value, mu[static_cast<std::size_t>(particle)]},
                    {0.0, 0.0, 0.0, 0.0}, ProtonMassKg, LightSpeedMPerS,
                    dt, provider, random, NULL);
            if (!increment.status.ok())
              throw std::runtime_error("CV06 production Dmumu step failed: " +
                                       increment.status.message);
            mu[static_cast<std::size_t>(particle)] = increment.state.mu;
            reflections += increment.pitchAngleReflections;
          }
          if ((step + 1) % sampleEvery == 0) {
            const double time = (step + 1) * dt;
            for (unsigned measured = 0; measured <= 6; ++measured) {
              long double sum = 0.0L;
              for (std::size_t particle = 0; particle < mu.size(); ++particle)
                sum += Legendre(measured, mu[particle]);
              output << initialMode << ',' << measured << ','
                     << baseSeed + seedOffset << ',' << dt << ',' << time << ','
                     << static_cast<double>(sum / particles) << ','
                     << reflections << ',' << particles << '\n';
            }
          }
        }
      }
    }
  }
  Commit(&output, temporary, outputPath);
}

void RunCV07(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double speed = Number(values, "speed-m-per-s");
  const std::vector<double> rates = NumberList(values, "switching-rates-per-s");
  const std::vector<double> normalizedTimes = NumberList(values, "normalized-times");
  const std::uint64_t particles = Unsigned(values, "particle-count");
  const std::uint64_t seeds = Unsigned(values, "seed-count");
  const std::uint64_t bins = Unsigned(values, "profile-bins");
  const std::uint64_t baseSeed = Unsigned(values, "campaign-seed");
  if (!(speed > 0.0) || particles < 100 || seeds < 1 || bins < 20)
    throw std::runtime_error("CV07 configuration violates its physical domain");
  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  if (!output.good()) throw std::runtime_error("cannot create CV07 model CSV");
  output << std::setprecision(17)
         << "rate_per_s,normalized_time,time_s,seed,row_type,bin_index,"
            "bin_left_m,bin_right_m,probability,front_mass,msd_m2,kurtosis,"
            "outside_support,event_mean\n";
  for (std::size_t rateIndex = 0; rateIndex < rates.size(); ++rateIndex) {
    const double rate = rates[rateIndex];
    if (!(rate > 0.0)) throw std::runtime_error("CV07 switching rate must be positive");
    for (std::size_t timeIndex = 0; timeIndex < normalizedTimes.size(); ++timeIndex) {
      const double q = normalizedTimes[timeIndex];
      const double duration = q / rate;
      const double extent = speed * duration;
      for (std::uint64_t seedOffset = 0; seedOffset < seeds; ++seedOffset) {
        std::vector<std::uint64_t> histogram(static_cast<std::size_t>(bins), 0);
        long double second = 0.0L, fourth = 0.0L;
        std::uint64_t fronts = 0, outside = 0, events = 0;
        for (std::uint64_t particle = 0; particle < particles; ++particle) {
          KeyedRandomStream random(baseSeed + seedOffset, particle + 1,
                                   707 + rateIndex, timeIndex);
          double direction = random.UniformOpen01() < 0.5 ? -1.0 : 1.0;
          double elapsed = 0.0, position = 0.0;
          std::uint64_t particleEvents = 0;
          for (;;) {
            const SEP::Transport::ScalarResult wait =
                SEP::Transport::SampleExponentialWaitingTime(rate, random);
            if (!wait.status.ok())
              throw std::runtime_error("CV07 exponential wait failed");
            const double interval = std::min(wait.value, duration - elapsed);
            position += direction * speed * interval;
            elapsed += interval;
            if (elapsed >= duration) break;
            direction = -direction;
            ++particleEvents;
          }
          events += particleEvents;
          if (particleEvents == 0) ++fronts;
          if (std::fabs(position) > extent * (1.0 + 8.0e-15)) ++outside;
          const long bin = static_cast<long>(std::floor(
              (position + extent) / (2.0 * extent) * bins));
          const std::size_t index = static_cast<std::size_t>(
              std::max<long>(0, std::min<long>(bins - 1, bin)));
          ++histogram[index];
          second += position * position;
          fourth += position * position * position * position;
        }
        const double msd = static_cast<double>(second / particles);
        const double kurtosis = static_cast<double>(fourth / particles) /
            (msd * msd);
        output << rate << ',' << q << ',' << duration << ','
               // Four empty profile fields (bin index, left, right, and
               // probability) keep moment rows aligned with the CSV schema.
               << baseSeed + seedOffset << ",moment,,,,,"
               << static_cast<double>(fronts) / particles << ',' << msd << ','
               << kurtosis << ',' << static_cast<double>(outside) / particles
               << ',' << static_cast<double>(events) / particles << '\n';
        const double width = 2.0 * extent / bins;
        for (std::uint64_t binIndex = 0; binIndex < bins; ++binIndex) {
          const double left = -extent + binIndex * width;
          output << rate << ',' << q << ',' << duration << ','
                 << baseSeed + seedOffset << ",profile," << binIndex << ','
                 << left << ',' << left + width << ','
                 << static_cast<double>(histogram[binIndex]) / particles
                 << ",,,,,\n";
        }
      }
    }
  }
  Commit(&output, temporary, outputPath);
}

class ConstantSpatialProvider final : public SEP::Transport::SpatialDiffusionProvider {
 public:
  explicit ConstantSpatialProvider(double kappa) : kappa_(kappa) {}
  SEP::Transport::SpatialDiffusionSample Evaluate(double, double) const override {
    SEP::Transport::SpatialDiffusionSample sample;
    sample.status = Status::Ok();
    sample.kappaParallelM2PerS = kappa_;
    sample.dKappaParallelDsMPerS = 0.0;
    sample.provenance = "validation:CV08:constant-drift-diffusion-v1";
    return sample;
  }
 private:
  double kappa_;
};

void RunCV08(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double boundary = Number(values, "boundary-m");
  const double kappa = Number(values, "kappa-m2-per-s");
  const std::vector<double> drifts = NumberList(values, "drifts-m-per-s");
  const std::vector<double> timeSteps = NumberList(values, "time-steps-s");
  const double maximumTime = Number(values, "maximum-time-s");
  const std::uint64_t particles = Unsigned(values, "particle-count");
  const std::uint64_t seeds = Unsigned(values, "seed-count");
  const std::uint64_t baseSeed = Unsigned(values, "campaign-seed");
  if (!(boundary > 0.0) || !(kappa > 0.0) || !(maximumTime > 0.0) ||
      particles < 100 || seeds < 10)
    throw std::runtime_error("CV08 configuration violates its physical domain");
  ConstantSpatialProvider provider(kappa);
  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  if (!output.good()) throw std::runtime_error("cannot create CV08 model CSV");
  output << std::setprecision(17)
         << "drift_m_per_s,dt_s,seed,particle_id,arrived,arrival_time_s,"
            "censor_time_s,overshoot_m\n";
  for (std::size_t driftIndex = 0; driftIndex < drifts.size(); ++driftIndex) {
    const double drift = drifts[driftIndex];
    if (!(drift > 0.0)) throw std::runtime_error("CV08 drift must be positive");
    for (std::size_t dtIndex = 0; dtIndex < timeSteps.size(); ++dtIndex) {
      const double dt = timeSteps[dtIndex];
      if (!(dt > 0.0)) throw std::runtime_error("CV08 dt must be positive");
      const std::uint64_t steps = static_cast<std::uint64_t>(
          std::ceil(maximumTime / dt));
      for (std::uint64_t seedOffset = 0; seedOffset < seeds; ++seedOffset) {
        for (std::uint64_t particle = 0; particle < particles; ++particle) {
          double position = 0.0;
          bool arrived = false;
          double arrival = maximumTime;
          double overshoot = 0.0;
          for (std::uint64_t step = 0; step < steps; ++step) {
            const double startTime = step * dt;
            if (startTime >= maximumTime) break;
            const double localDt = std::min(dt, maximumTime - startTime);
            KeyedRandomStream random(baseSeed + seedOffset, particle + 1,
                                     808 + driftIndex, step);
            const SEP::Transport::ParkerIncrement increment =
                SEP::Transport::AdvanceParker(
                    {position, 1.0}, {drift, 0.0}, 1.0, localDt,
                    provider, random);
            if (!increment.status.ok())
              throw std::runtime_error("CV08 production Parker step failed");
            const double candidate = increment.state.arcLengthM;
            if (candidate >= boundary) {
              overshoot = candidate - boundary;
              // This is the production straight-segment crossing convention:
              // interpolate the fraction of the accepted particle step whose
              // spatial chord reaches the boundary, then stop the history.
              const double fraction = (boundary - position) /
                  (candidate - position);
              arrival = startTime + std::max(0.0, std::min(1.0, fraction)) * localDt;
              arrived = true;
              break;
            }
            position = candidate;
          }
          output << drift << ',' << dt << ',' << baseSeed + seedOffset << ','
                 << particle + 1 << ',' << (arrived ? 1 : 0) << ',';
          if (arrived) output << arrival;
          output << ',' << maximumTime << ',' << overshoot << '\n';
        }
      }
    }
  }
  Commit(&output, temporary, outputPath);
}

void RunCV09(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const std::vector<double> ratios = NumberList(values, "compression-ratios");
  const double upstreamSpeed = Number(values, "upstream-speed-m-per-s");
  const double particleSpeed = Number(values, "particle-speed-m-per-s");
  const double kappaUp = Number(values, "kappa-upstream-m2-per-s");
  const double kappaDown = Number(values, "kappa-downstream-m2-per-s");
  const double injectionMomentum = Number(values, "injection-momentum-kg-m-per-s");
  const std::uint64_t particles = Unsigned(values, "particle-count");
  const std::uint64_t baseSeed = Unsigned(values, "campaign-seed");
  if (!(upstreamSpeed > 0.0) || !(particleSpeed > 20.0 * upstreamSpeed) ||
      !(kappaUp > 0.0) || !(kappaDown > 0.0) ||
      !(injectionMomentum > 0.0) || particles < 1000)
    throw std::runtime_error("CV09 configuration violates its resolved DSA domain");
  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  if (!output.good()) throw std::runtime_error("cannot create CV09 model CSV");
  output << std::setprecision(17)
         << "compression_ratio,particle_id,momentum_kg_m_per_s,"
            "acceleration_time_s,cycles,escaped_weight,energy_proxy_j,"
            "upstream_diffusion_length_m,downstream_diffusion_length_m\n";
  for (std::size_t ratioIndex = 0; ratioIndex < ratios.size(); ++ratioIndex) {
    const double ratio = ratios[ratioIndex];
    if (!(ratio > 1.0)) throw std::runtime_error("CV09 compression ratio must exceed one");
    const double downstreamSpeed = upstreamSpeed / ratio;
    const double gain = 4.0 * (upstreamSpeed - downstreamSpeed) /
                        (3.0 * particleSpeed);
    const double escapeProbability = 4.0 * downstreamSpeed / particleSpeed;
    const double accelerationTime = 3.0 / (upstreamSpeed - downstreamSpeed) *
        (kappaUp / upstreamSpeed + kappaDown / downstreamSpeed);
    const double cycleTime = accelerationTime * std::log1p(gain);
    for (std::uint64_t particle = 0; particle < particles; ++particle) {
      KeyedRandomStream random(baseSeed, particle + 1, 909, ratioIndex);
      // Invert the geometric survival law rather than looping over thousands
      // of shock crossings. This is exactly equivalent for constant return
      // probability and preserves deterministic keyed randomness.
      const std::uint64_t cycles = static_cast<std::uint64_t>(std::floor(
          std::log(random.UniformOpen01()) /
          std::log(1.0 - escapeProbability)));
      const double momentum = injectionMomentum *
          std::pow(1.0 + gain, static_cast<double>(cycles));
      const double time = cycles * cycleTime;
      // The energy proxy p*v is sufficient for a closed accounting check in
      // this constant-speed analytical DSA cycle model; it is explicitly not
      // presented as a relativistic kinetic energy.
      output << ratio << ',' << particle + 1 << ',' << momentum << ',' << time
             << ',' << cycles << ",1," << momentum * particleSpeed << ','
             << kappaUp / upstreamSpeed << ','
             << kappaDown / downstreamSpeed << '\n';
    }
  }
  Commit(&output, temporary, outputPath);
}

double SineAverage(std::size_t cell, std::size_t cells, double shift) {
  const double dx = 1.0 / cells;
  const double left = cell * dx - shift;
  const double right = (cell + 1) * dx - shift;
  return 1.0 + 0.25 * (std::cos(2.0 * Pi * left) -
      std::cos(2.0 * Pi * right)) / (2.0 * Pi * dx);
}

SEP::Turbulence::CellState WaveCell(double length, double area,
                                    double plusJ, double minusJ) {
  SEP::Turbulence::CellState cell;
  cell.lengthM = length;
  cell.volumeM3 = length * area;
  cell.plasmaSpeedMPerS = 1.0;
  cell.alfvenSpeedMPerS = 0.0;
  cell.magneticFieldT = 5.0e-9;
  cell.massDensityKgPerM3 = 1.0;
  cell.ePlusJ = plusJ;
  cell.eMinusJ = minusJ;
  return cell;
}

SEP::Turbulence::State WaveAdvectionState(std::size_t cells, bool variableArea,
                                          bool plusBranch, std::size_t bins) {
  SEP::Turbulence::State state;
  state.configuration.source = bins > 1
      ? SEP::Turbulence::Source::SelfConsistentSpectral
      : SEP::Turbulence::Source::SelfConsistentIntegrated;
  state.configuration.representation = bins > 1
      ? SEP::Turbulence::Representation::Spectral
      : SEP::Turbulence::Representation::Integrated;
  state.configuration.spectralBins = bins;
  state.configuration.advectionEnabled = true;
  state.configuration.reflectionEnabled = false;
  state.configuration.cascadeEnabled = false;
  state.configuration.shockInjectionEnabled = false;
  state.configuration.periodicBoundaries = true;
  state.configuration.coupling = SEP::Turbulence::CouplingPolicy::Disabled;
  state.provenance = "validation:CV10:manufactured-wave-advection-v1";
  const double dx = 1.0 / cells;
  for (std::size_t i = 0; i < cells; ++i) {
    const double center = (i + 0.5) * dx;
    const double area = variableArea ? 1.0 + 0.3 * std::sin(2.0 * Pi * center) : 1.0;
    const double conserved = SineAverage(i, cells, 0.0) * dx;
    SEP::Turbulence::CellState cell = WaveCell(
        dx, area, plusBranch ? conserved : 0.0,
        plusBranch ? 0.0 : conserved);
    // U=0 and VA=1 m/s make the two production characteristic velocities
    // U+VA and U-VA equal and opposite. This explicitly tests propagation
    // sense while keeping the manufactured period and SI reference simple.
    cell.plasmaSpeedMPerS = 0.0;
    cell.alfvenSpeedMPerS = 1.0;
    if (bins > 1) {
      cell.spectralEnergyJ.assign(2 * bins, 0.0);
      for (std::size_t k = 0; k < bins; ++k) {
        const double fraction = (k + 1.0) / (0.5 * bins * (bins + 1.0));
        cell.spectralEnergyJ[(plusBranch ? 0 : bins) + k] = conserved * fraction;
      }
    }
    state.cells.push_back(cell);
  }
  return state;
}

void RunCV10(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const std::vector<std::uint64_t> resolutions = UnsignedList(values, "resolutions");
  const std::uint64_t spectralBins = Unsigned(values, "spectral-bins");
  const double duration = Number(values, "duration-s");
  const double cfl = Number(values, "cfl");
  if (!(duration > 0.0) || !(cfl > 0.0 && cfl < 1.0) ||
      spectralBins < 2 || resolutions.size() < 3)
    throw std::runtime_error("CV10 configuration violates its numerical domain");
  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  if (!output.good()) throw std::runtime_error("cannot create CV10 model CSV");
  output << std::setprecision(17)
         << "scenario,branch,resolution,spectral_bin,cell,left_m,right_m,"
            "wave_action_j,total_branch_energy_j,negative_cells\n";
  const char* scenarios[] = {"fixed", "expanding-area", "moving-grid"};
  for (int scenarioIndex = 0; scenarioIndex < 3; ++scenarioIndex) {
    for (int branchIndex = 0; branchIndex < 2; ++branchIndex) {
      const bool plus = branchIndex == 0;
      for (std::size_t resolutionIndex = 0;
           resolutionIndex < resolutions.size(); ++resolutionIndex) {
        const std::size_t cells = static_cast<std::size_t>(resolutions[resolutionIndex]);
        if (cells < 8) throw std::runtime_error("CV10 resolution must be >=8");
        const bool moving = scenarioIndex == 2;
        const bool variableArea = scenarioIndex == 1;
        const std::size_t initialCells = moving ? cells / 2 : cells;
        SEP::Turbulence::State state = WaveAdvectionState(
            initialCells, variableArea, plus, spectralBins);
        if (!SEP::Turbulence::InitializeState(&state).ok())
          throw std::runtime_error("CV10 turbulence initialization failed");
        const double dt = cfl / cells;
        const std::uint64_t steps = static_cast<std::uint64_t>(std::llround(duration / dt));
        if (std::fabs(duration / dt - steps) > 1.0e-10 || steps % 2 != 0)
          throw std::runtime_error("CV10 duration must contain an even integral step count");
        for (std::uint64_t step = 0; step < steps; ++step) {
          if (moving && step == steps / 2) {
            std::vector<SEP::Turbulence::CellState> geometry;
            const double dx = 1.0 / cells;
            for (std::size_t i = 0; i < cells; ++i) {
              SEP::Turbulence::CellState target =
                  WaveCell(dx, 1.0, 0.0, 0.0);
              // Remapping replaces geometry as well as conserved values. Keep
              // the manufactured U=0, VA=1 m/s characteristic speeds on the
              // refined mesh or the minus branch would silently reverse.
              target.plasmaSpeedMPerS = 0.0;
              target.alfvenSpeedMPerS = 1.0;
              geometry.push_back(target);
            }
            SEP::Turbulence::State remapped;
            SEP::Turbulence::EnergyLedger ledger;
            const Status status = SEP::Turbulence::RemapConservatively(
                state, geometry, &remapped, &ledger);
            if (!status.ok()) throw std::runtime_error("CV10 moving-grid remap failed");
            state = remapped;
          }
          const SEP::Turbulence::StepResult advanced =
              SEP::Turbulence::Advance(&state, dt);
          if (!advanced.status.ok())
            throw std::runtime_error("CV10 production turbulence advance failed: " +
                                     advanced.status.message);
        }
        const std::size_t finalCells = state.cells.size();
        for (std::size_t k = 0; k < spectralBins; ++k) {
          double total = 0.0;
          std::uint64_t negative = 0;
          for (std::size_t i = 0; i < finalCells; ++i) {
            const double energy = state.cells[i].spectralEnergyJ[
                (plus ? 0 : spectralBins) + k];
            total += energy;
            if (energy < 0.0) ++negative;
          }
          for (std::size_t i = 0; i < finalCells; ++i) {
            const double dx = 1.0 / finalCells;
            const double energy = state.cells[i].spectralEnergyJ[
                (plus ? 0 : spectralBins) + k];
            output << scenarios[scenarioIndex] << ',' << (plus ? "plus" : "minus")
                   << ',' << cells << ',' << k << ',' << i << ',' << i * dx << ','
                   << (i + 1) * dx << ',' << energy << ',' << total << ','
                   << negative << '\n';
          }
        }
      }
    }
  }
  Commit(&output, temporary, outputPath);
}

SEP::Turbulence::State OneHotSpectrum(std::size_t bins, std::size_t activeBin,
                                      double initialEnergy) {
  SEP::Turbulence::State state;
  state.configuration.source = SEP::Turbulence::Source::SelfConsistentSpectral;
  state.configuration.representation = SEP::Turbulence::Representation::Spectral;
  state.configuration.spectralBins = bins;
  state.configuration.advectionEnabled = false;
  state.configuration.reflectionEnabled = false;
  state.configuration.cascadeEnabled = false;
  state.configuration.shockInjectionEnabled = false;
  state.configuration.coupling = SEP::Turbulence::CouplingPolicy::StreamingEnergyExchange;
  state.provenance = "validation:CV11:one-hot-resonant-spectrum-v1";
  SEP::Turbulence::CellState cell = WaveCell(1.0, 1.0, initialEnergy, 0.0);
  cell.spectralEnergyJ.assign(2 * bins, 0.0);
  cell.spectralEnergyJ[activeBin] = initialEnergy;
  state.cells.push_back(cell);
  return state;
}

void RunCV11(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const std::uint64_t bins = Unsigned(values, "spectral-bins");
  const std::uint64_t activeBin = Unsigned(values, "active-bin");
  const double initial = Number(values, "initial-energy-j");
  const double growth = Number(values, "growth-rate-per-s");
  const double damping = Number(values, "damping-rate-per-s");
  const double sinusoidal = Number(values, "sinusoidal-rate-per-s");
  const double omega = Number(values, "angular-frequency-per-s");
  const double finalTime = Number(values, "final-time-s");
  const std::vector<double> timeSteps = NumberList(values, "time-steps-s");
  if (bins < 4 || activeBin >= bins || !(initial > 0.0) || !(omega > 0.0) ||
      !(finalTime > 0.0))
    throw std::runtime_error("CV11 configuration violates its physical domain");
  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  if (!output.good()) throw std::runtime_error("cannot create CV11 model CSV");
  output << std::setprecision(17)
         << "scenario,dt_s,time_s,spectral_bin,wave_energy_j,integrated_energy_j,"
            "ledger_residual_j,positive\n";
  const char* scenarios[] = {"growth", "damping", "cancellation", "sign-change"};
  for (int scenario = 0; scenario < 4; ++scenario) {
    for (std::size_t dtIndex = 0; dtIndex < timeSteps.size(); ++dtIndex) {
      const double dt = timeSteps[dtIndex];
      const std::uint64_t steps = static_cast<std::uint64_t>(std::llround(finalTime / dt));
      if (!(dt > 0.0) || std::fabs(finalTime / dt - steps) > 1.0e-10)
        throw std::runtime_error("CV11 final time must be divisible by dt");
      SEP::Turbulence::State state = OneHotSpectrum(
          bins, activeBin, initial);
      if (!SEP::Turbulence::InitializeState(&state).ok())
        throw std::runtime_error("CV11 turbulence initialization failed");
      for (std::uint64_t step = 0; step <= steps; ++step) {
        const double time = step * dt;
        double total = 0.0;
        for (std::size_t k = 0; k < bins; ++k)
          total += state.cells[0].spectralEnergyJ[k];
        for (std::size_t k = 0; k < bins; ++k)
          output << scenarios[scenario] << ',' << dt << ',' << time << ',' << k
                 << ',' << state.cells[0].spectralEnergyJ[k] << ',' << total
                 << ',' << state.accumulatedLedger.closureResidualJ << ','
                 << (state.cells[0].spectralEnergyJ[k] >= 0.0 ? 1 : 0) << '\n';
        if (step == steps) break;
        const double t0 = step * dt;
        const double t1 = (step + 1) * dt;
        double exponent = 0.0;
        if (scenario == 0) exponent = 2.0 * growth * dt;
        if (scenario == 1) exponent = -2.0 * damping * dt;
        if (scenario == 2) exponent = 2.0 * (growth - growth) * dt;
        if (scenario == 3) {
          // Integrate the time-varying rate with the midpoint rule, then send
          // the resulting energy transfer through the production ledger. The
          // independent reference integrates the sine exactly, turning this
          // scenario into a real second-order temporal-refinement check.
          const double midpoint = 0.5 * (t0 + t1);
          exponent = 2.0 *
              (sinusoidal * std::sin(omega * midpoint) - damping) * dt;
        }
        const double current = state.cells[0].ePlusJ;
        state.cells[0].pendingParticlePlusJ = current * std::expm1(exponent);
        const SEP::Turbulence::StepResult result = SEP::Turbulence::Advance(&state, dt);
        if (!result.status.ok())
          throw std::runtime_error("CV11 production turbulence source update failed");
      }
    }
  }
  Commit(&output, temporary, outputPath);
}

double RelativisticEnergy(double speed) {
  const double beta = speed / LightSpeedMPerS;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  return ProtonMassKg * LightSpeedMPerS * LightSpeedMPerS *
      gamma * gamma * beta * beta / (gamma + 1.0);
}

void RunCV12(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double initialSpeed = Number(values, "particle-speed-m-per-s");
  const double alfvenSpeed = Number(values, "alfven-speed-m-per-s");
  const double initialWave = Number(values, "initial-wave-energy-j");
  const double targetExchange = Number(values, "target-exchange-j");
  const std::vector<double> timeSteps = NumberList(values, "time-steps-s");
  const std::vector<std::uint64_t> particleCounts = UnsignedList(values, "particle-counts");
  const std::uint64_t steps = Unsigned(values, "steps");
  const std::uint64_t baseSeed = Unsigned(values, "campaign-seed");
  if (!(initialSpeed > 0.0 && initialSpeed < LightSpeedMPerS) ||
      !(std::fabs(alfvenSpeed) < initialSpeed) || !(initialWave > 0.0) ||
      !(targetExchange > 0.0 && targetExchange < initialWave) || steps < 1)
    throw std::runtime_error("CV12 configuration violates its closed-system domain");
  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  if (!output.good()) throw std::runtime_error("cannot create CV12 model CSV");
  output << std::setprecision(17)
         << "scenario,particle_count,dt_s,step,time_s,particle_energy_j,"
            "wave_energy_j,total_energy_j,exchange_j,ledger_exchange_j,"
            "step_residual_j,resonant_branch\n";
  const char* scenarios[] = {"particle-only", "wave-only", "coupled", "negative-control"};
  for (std::size_t countIndex = 0; countIndex < particleCounts.size(); ++countIndex) {
    const std::uint64_t count = particleCounts[countIndex];
    for (std::size_t dtIndex = 0; dtIndex < timeSteps.size(); ++dtIndex) {
      const double dt = timeSteps[dtIndex];
      for (int scenario = 0; scenario < 4; ++scenario) {
        double particleEnergy = count * RelativisticEnergy(initialSpeed);
        double waveEnergy = initialWave;
        double speed = initialSpeed;
        const double initialTotal = particleEnergy + waveEnergy;
        for (std::uint64_t step = 0; step <= steps; ++step) {
          if (step == 0) {
            output << scenarios[scenario] << ',' << count << ',' << dt
                   << ",0,0," << particleEnergy << ',' << waveEnergy << ','
                   << particleEnergy + waveEnergy << ",0,0,0,none\n";
            continue;
          }
          double exchange = 0.0, ledgerExchange = 0.0, residual = 0.0;
          std::string branch = "none";
          if (scenario >= 2) {
            KeyedRandomStream random(baseSeed, count, 1212, step);
            const double initialMu = step % 2 ? 0.65 : -0.55;
            const double signedWave = initialMu > 0.0 ? -alfvenSpeed : alfvenSpeed;
            branch = initialMu > 0.0 ? "minus" : "plus";
            const SEP::Transport::WaveFrameScatterResult scattered =
                SEP::Transport::ScatterIsotropicallyInWaveFrame(
                    speed, initialMu, signedWave, LightSpeedMPerS, random);
            if (!scattered.status.ok())
              throw std::runtime_error("CV12 wave-frame scattering failed");
            const double singleChange = RelativisticEnergy(scattered.speedMPerS) -
                                        RelativisticEnergy(speed);
            if (singleChange == 0.0)
              throw std::runtime_error("CV12 controlled scattering produced zero exchange");
            const double macroWeight = targetExchange / std::fabs(singleChange);
            exchange = macroWeight * singleChange;
            particleEnergy += exchange;
            SEP::Turbulence::State wave;
            wave.configuration.advectionEnabled = false;
            wave.configuration.reflectionEnabled = false;
            wave.configuration.cascadeEnabled = false;
            wave.configuration.shockInjectionEnabled = false;
            wave.configuration.coupling =
                SEP::Turbulence::CouplingPolicy::StreamingEnergyExchange;
            wave.provenance = "validation:CV12:closed-exchange-v1";
            wave.cells.push_back(WaveCell(1.0, 1.0,
                branch == "plus" ? waveEnergy : 0.0,
                branch == "minus" ? waveEnergy : 0.0));
            if (!SEP::Turbulence::InitializeState(&wave).ok())
              throw std::runtime_error("CV12 wave initialization failed");
            if (scenario == 2) {
              if (branch == "plus") wave.cells[0].pendingParticlePlusJ = -exchange;
              else wave.cells[0].pendingParticleMinusJ = -exchange;
            }
            const SEP::Turbulence::StepResult advanced =
                SEP::Turbulence::Advance(&wave, dt);
            if (!advanced.status.ok())
              throw std::runtime_error("CV12 wave ledger update failed");
            ledgerExchange = advanced.ledger.particleExchangeJ;
            waveEnergy = wave.cells[0].ePlusJ + wave.cells[0].eMinusJ;
            residual = particleEnergy + waveEnergy - initialTotal;
            speed = scattered.speedMPerS;
          }
          output << scenarios[scenario] << ',' << count << ',' << dt << ','
                 << step << ',' << step * dt << ',' << particleEnergy << ','
                 << waveEnergy << ',' << particleEnergy + waveEnergy << ','
                 << exchange << ',' << ledgerExchange << ',' << residual << ','
                 << branch << '\n';
        }
      }
    }
  }
  Commit(&output, temporary, outputPath);
}

}  // namespace

namespace SEP {
namespace Validation {

bool RunAdvancedValidationModel(const std::string& caseId,
                                const std::vector<std::string>& arguments,
                                const std::string& outputPath,
                                std::string* error) {
  try {
    const std::map<std::string, std::string> values = ParseArguments(arguments);
    if (caseId == "CV06") RunCV06(values, outputPath);
    else if (caseId == "CV07") RunCV07(values, outputPath);
    else if (caseId == "CV08") RunCV08(values, outputPath);
    else if (caseId == "CV09") RunCV09(values, outputPath);
    else if (caseId == "CV10") RunCV10(values, outputPath);
    else if (caseId == "CV11") RunCV11(values, outputPath);
    else if (caseId == "CV12") RunCV12(values, outputPath);
    else throw std::runtime_error("unsupported advanced validation case: " + caseId);
    if (error) error->clear();
    return true;
  } catch (const std::exception& exception) {
    if (error) *error = exception.what();
    return false;
  }
}

}  // namespace Validation
}  // namespace SEP

#ifdef SRCSEP_ADVANCED_MODELS_STANDALONE_TEST_HARNESS
// This executable exists only for strict source and orchestration checks. A
// scientific run must enter through the linked srcSEP/AMPS test registry.
int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: advanced_models CASE OUTPUT --name value ...\n";
    return 2;
  }
  std::vector<std::string> arguments;
  for (int i = 3; i < argc; ++i) arguments.push_back(argv[i]);
  std::string error;
  if (!SEP::Validation::RunAdvancedValidationModel(
          argv[1], arguments, argv[2], &error)) {
    std::cerr << argv[1] << " model error: " << error << '\n';
    return 2;
  }
  return 0;
}
#endif
