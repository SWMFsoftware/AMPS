// CV01 linked-application model: ballistic focused transport on a uniform line.
//
// The implementation is compiled into srcSEP's main library and reached from
// the production component-test registry. Every physical step calls
// AdvanceFocusedTransportDmumu(), the same kernel used by the PIC adapter. The
// code below supplies only the controlled background, immutable particle
// fixture, boundary policy, and serialization. The analytical reference stays
// in a separate Python process and never enters the linked application.

#include "cv01_model.h"

// Use paths relative to this translation unit rather than relying on an
// application-specific ``-I srcSEP/util`` flag. The enclosing AMPS build
// compiles sources in nested srcSEP directories with its global include list,
// which intentionally does not expose srcSEP/util as a flat include root.
// Repository-relative includes therefore work both in the linked production
// build and in the dependency-light CV01 compile gate.
#include "../../../util/sep_common_header_path.h"
#include "../../../util/sep_focused_transport_core.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cmath>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cstdio>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::Transport::AdvanceCoordinate;
using SEP::Transport::AdvanceFocusedTransportDmumu;
using SEP::Transport::BoundaryPolicy;
using SEP::Transport::FocusedEquationMode;
using SEP::Transport::FocusedTransportBackground;
using SEP::Transport::FocusedTransportIncrement;
using SEP::Transport::FocusedTransportState;
using SEP::Transport::KeyedRandomStream;
using SEP::Transport::PitchAngleDiffusionProvider;
using SEP::Transport::PitchAngleDiffusionSample;
using SEP::Transport::ScalarResult;
using SEP::Transport::Status;
using SEP::Transport::StatusCode;

struct InitialParticle {
  std::uint64_t id = 0;
  double initialPositionM = 0.0;
  double mu = 0.0;
  double weight = 0.0;
};

struct Particle {
  InitialParticle initial;
  FocusedTransportState state;
  double unwrappedPositionM = 0.0;
  double boundaryPositionM = std::numeric_limits<double>::quiet_NaN();
  double crossingTimeS = std::numeric_limits<double>::quiet_NaN();
  bool active = true;
};

// A valid zero-scattering provider must still publish coefficient provenance
// and turbulence-state identity.  This ensures CV01 exercises the same strict
// provider contract as a physical D_mumu calculation while turning only the
// scattering operator off.
class ZeroPitchAngleDiffusion final : public PitchAngleDiffusionProvider {
 public:
  PitchAngleDiffusionSample Evaluate(double, double, double) const override {
    PitchAngleDiffusionSample sample;
    sample.status = Status::Ok();
    sample.dMuMuPerS = 0.0;
    sample.dDmuMuDmuPerS = 0.0;
    sample.provenance = "CV01:prescribed-zero-dmumu-v1";
    sample.turbulenceStateIdentity = "CV01:frozen-uniform-background-v1";
    return sample;
  }
};

double ParseDouble(const std::map<std::string, std::string>& options,
                   const std::string& name) {
  const std::map<std::string, std::string>::const_iterator found =
      options.find(name);
  if (found == options.end()) throw std::runtime_error("missing option --" + name);
  char* end = NULL;
  const double value = std::strtod(found->second.c_str(), &end);
  if (!end || *end != '\0' || !std::isfinite(value))
    throw std::runtime_error("invalid finite number for --" + name);
  return value;
}

std::string Require(const std::map<std::string, std::string>& options,
                    const std::string& name) {
  const std::map<std::string, std::string>::const_iterator found =
      options.find(name);
  if (found == options.end() || found->second.empty())
    throw std::runtime_error("missing option --" + name);
  return found->second;
}

std::uint64_t ParseUnsigned64(
    const std::map<std::string, std::string>& options,
    const std::string& name) {
  const std::string text = Require(options, name);
  if (text.empty() || text[0] == '-')
    throw std::runtime_error("invalid unsigned integer for --" + name);
  char* end = NULL;
  errno = 0;
  const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || !end || *end != '\0' ||
      value > std::numeric_limits<std::uint64_t>::max())
    throw std::runtime_error("invalid unsigned integer for --" + name);
  // Parsing the seed as an integer, rather than through double, preserves all
  // 64 bits of the deterministic keyed-stream identity in reviewed variants.
  return static_cast<std::uint64_t>(value);
}

std::map<std::string, std::string> ParseOptions(int argc, char** argv) {
  std::map<std::string, std::string> values;
  for (int i = 1; i < argc; i += 2) {
    const std::string token = argv[i];
    if (token.size() < 3 || token.substr(0, 2) != "--" || i + 1 >= argc)
      throw std::runtime_error("options require --name value pairs");
    values[token.substr(2)] = argv[i + 1];
  }
  return values;
}

std::vector<InitialParticle> ReadParticles(const std::string& path) {
  std::ifstream input(path.c_str());
  if (!input.good()) throw std::runtime_error("cannot open initial particle CSV: " + path);
  std::string line;
  if (!std::getline(input, line))
    throw std::runtime_error("initial particle CSV has no header");
  // Python's standards-compliant csv.writer may emit CRLF even on Linux.
  // std::getline removes LF only, so discard a terminal CR before comparing
  // the schema rather than making the cross-language input platform-specific.
  if (!line.empty() && line[line.size() - 1] == '\r') line.erase(line.size() - 1);
  if (line != "particle_id,initial_s_m,mu,weight")
    throw std::runtime_error("unexpected initial particle CSV header");
  std::vector<InitialParticle> particles;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::istringstream row(line);
    std::string field;
    InitialParticle particle;
    if (!std::getline(row, field, ',')) throw std::runtime_error("missing particle ID");
    particle.id = static_cast<std::uint64_t>(std::strtoull(field.c_str(), NULL, 10));
    if (!std::getline(row, field, ',')) throw std::runtime_error("missing initial position");
    particle.initialPositionM = std::strtod(field.c_str(), NULL);
    if (!std::getline(row, field, ',')) throw std::runtime_error("missing pitch angle");
    particle.mu = std::strtod(field.c_str(), NULL);
    if (!std::getline(row, field, ',')) throw std::runtime_error("missing weight");
    particle.weight = std::strtod(field.c_str(), NULL);
    if (particle.id == 0 || !std::isfinite(particle.initialPositionM) ||
        !std::isfinite(particle.mu) || particle.mu < -1.0 || particle.mu > 1.0 ||
        !std::isfinite(particle.weight) || particle.weight <= 0.0)
      throw std::runtime_error("initial particle row is outside the CV01 physical domain");
    particles.push_back(particle);
  }
  if (particles.empty()) throw std::runtime_error("initial particle CSV is empty");
  return particles;
}

double WrapPeriodic(double positionM, double minimumM, double maximumM) {
  const double width = maximumM - minimumM;
  double wrapped = std::fmod(positionM - minimumM, width);
  if (wrapped < 0.0) wrapped += width;
  return minimumM + wrapped;
}

void WriteParticle(std::ostream& output, const Particle& particle,
                   double timeS, double dtS, const std::string& boundary) {
  const double activeWeight = particle.active ? particle.initial.weight : 0.0;
  const double escapedWeight = particle.active ? 0.0 : particle.initial.weight;
  output << timeS << ',' << dtS << ',' << boundary << ','
         << particle.initial.id << ',' << particle.initial.mu << ','
         << (particle.active ? 1 : 0) << ',' << particle.state.arcLengthM << ','
         << particle.unwrappedPositionM << ',' << particle.state.mu << ','
         << particle.state.momentumKgMPerS << ',' << activeWeight << ','
         << escapedWeight << ',';
  if (std::isfinite(particle.crossingTimeS)) output << particle.crossingTimeS;
  output << '\n';
}

int Run(const std::vector<std::string>& arguments,
        const std::string& outputPath) {
  // Convert the registry-owned immutable argument vector to the historical
  // parser representation. Keeping this small adapter preserves the already
  // reviewed option/domain checks while removing the standalone executable
  // from the scientific-validation path.
  std::vector<char*> argv;
  argv.push_back(const_cast<char*>("CV01"));
  for (std::size_t i = 0; i < arguments.size(); ++i)
    argv.push_back(const_cast<char*>(arguments[i].c_str()));
  const std::map<std::string, std::string> options =
      ParseOptions(static_cast<int>(argv.size()), argv.data());
  const std::string boundary = Require(options, "boundary");
  if (boundary != "periodic" && boundary != "open")
    throw std::runtime_error("--boundary must be periodic or open");

  const double minimumM = ParseDouble(options, "minimum-m");
  const double maximumM = ParseDouble(options, "maximum-m");
  const double speedMPerS = ParseDouble(options, "speed-m-per-s");
  const double massKg = ParseDouble(options, "mass-kg");
  const double lightSpeedMPerS = ParseDouble(options, "light-speed-m-per-s");
  const double plasmaSpeedMPerS = ParseDouble(options, "plasma-speed-m-per-s");
  const double dtS = ParseDouble(options, "dt-s");
  const double finalTimeS = ParseDouble(options, "final-time-s");
  const double sampleIntervalS = ParseDouble(options, "sample-interval-s");
  const std::uint64_t seed = ParseUnsigned64(options, "campaign-seed");
  if (!(minimumM < maximumM) || !(speedMPerS >= 0.0) ||
      !(speedMPerS < lightSpeedMPerS) || !(massKg > 0.0) || !(dtS > 0.0) ||
      !(finalTimeS > 0.0) || !(sampleIntervalS > 0.0))
    throw std::runtime_error("CV01 dimensions, speeds, mass, and times must be physical");

  const double totalStepsReal = finalTimeS / dtS;
  const double sampleStepsReal = sampleIntervalS / dtS;
  const std::uint64_t totalSteps = static_cast<std::uint64_t>(std::llround(totalStepsReal));
  const std::uint64_t sampleSteps = static_cast<std::uint64_t>(std::llround(sampleStepsReal));
  if (sampleSteps == 0 || std::fabs(totalStepsReal - totalSteps) > 1.0e-10 ||
      std::fabs(sampleStepsReal - sampleSteps) > 1.0e-10)
    throw std::runtime_error("final and sample times must be integer multiples of dt");

  const ScalarResult initialMomentum = SEP::Transport::MomentumFromSpeed(
      speedMPerS, massKg, lightSpeedMPerS);
  if (!initialMomentum.status.ok())
    throw std::runtime_error(initialMomentum.status.message);

  std::vector<InitialParticle> initial = ReadParticles(
      Require(options, "initial-particles"));
  std::vector<Particle> particles;
  for (std::size_t i = 0; i < initial.size(); ++i) {
    if (initial[i].initialPositionM < minimumM ||
        initial[i].initialPositionM > maximumM)
      throw std::runtime_error("initial particle is outside the field-line domain");
    Particle particle;
    particle.initial = initial[i];
    particle.state = FocusedTransportState(
        initial[i].initialPositionM, initialMomentum.value, initial[i].mu);
    particle.unwrappedPositionM = initial[i].initialPositionM;
    particles.push_back(particle);
  }

  // Write to a same-directory staging path and publish with rename only after
  // every sampled row is flushed successfully. A killed validation process can
  // therefore leave a .tmp file, but never a truncated CSV that looks valid.
  const std::string temporaryPath = outputPath + ".tmp";
  std::ofstream output(temporaryPath.c_str(), std::ios::out | std::ios::trunc);
  if (!output.good())
    throw std::runtime_error("cannot create staged model CSV: " + temporaryPath);
  output << std::setprecision(17);
  output << "time_s,dt_s,boundary,particle_id,initial_mu,active,position_m,"
            "unwrapped_position_m,mu,momentum_kg_m_per_s,active_weight,"
            "escaped_weight,crossing_time_s\n";

  FocusedTransportBackground background;
  background.dLnAbsBdsPerM = 0.0;
  background.plasmaAdvectionMPerS = plasmaSpeedMPerS;
  background.parallelVelocityGradientPerS = 0.0;
  background.velocityDivergencePerS = 0.0;
  background.fieldAlignedStrainPerS = 0.0;
  background.equationMode = FocusedEquationMode::FullGyrotropic;
  ZeroPitchAngleDiffusion coefficient;

  for (std::size_t i = 0; i < particles.size(); ++i)
    WriteParticle(output, particles[i], 0.0, dtS, boundary);

  for (std::uint64_t stepIndex = 0; stepIndex < totalSteps; ++stepIndex) {
    const double startTimeS = static_cast<double>(stepIndex) * dtS;
    for (std::size_t i = 0; i < particles.size(); ++i) {
      Particle& particle = particles[i];
      if (!particle.active) continue;

      KeyedRandomStream random(seed, particle.initial.id, 8, stepIndex);
      const double startPositionM = particle.state.arcLengthM;
      const FocusedTransportIncrement increment = AdvanceFocusedTransportDmumu(
          particle.state, background, massKg, lightSpeedMPerS, dtS,
          coefficient, random, NULL);
      if (!increment.status.ok())
        throw std::runtime_error("production focused-transport kernel failed: " +
                                 increment.status.message);
      particle.unwrappedPositionM += increment.displacementM;

      if (boundary == "periodic") {
        particle.state = increment.state;
        particle.state.arcLengthM = WrapPeriodic(
            increment.state.arcLengthM, minimumM, maximumM);
      } else {
        const double directedSpeedMPerS = increment.displacementM / dtS;
        const SEP::Transport::CoordinateAdvance advance = AdvanceCoordinate(
            startPositionM, increment.displacementM, directedSpeedMPerS,
            minimumM, maximumM, BoundaryPolicy::Absorb);
        if (advance.status.code == StatusCode::OutOfDomain) {
          particle.active = false;
          const double boundaryPositionM = directedSpeedMPerS > 0.0
              ? maximumM : minimumM;
          const double inDomainDistanceM = boundaryPositionM - startPositionM;
          particle.crossingTimeS = startTimeS +
              inDomainDistanceM / directedSpeedMPerS;
          particle.boundaryPositionM = boundaryPositionM;
          particle.state = increment.state;
          particle.state.arcLengthM = boundaryPositionM;
        } else if (!advance.status.ok()) {
          throw std::runtime_error("open-boundary coordinate update failed: " +
                                   advance.status.message);
        } else {
          particle.state = increment.state;
          particle.state.arcLengthM = advance.positionM;
        }
      }
    }

    if (((stepIndex + 1) % sampleSteps) == 0) {
      const double sampleTimeS = static_cast<double>(stepIndex + 1) * dtS;
      for (std::size_t i = 0; i < particles.size(); ++i)
        WriteParticle(output, particles[i], sampleTimeS, dtS, boundary);
    }
  }
  output.close();
  if (output.fail()) throw std::runtime_error("failed while writing model CSV");
  if (std::rename(temporaryPath.c_str(), outputPath.c_str()) != 0)
    throw std::runtime_error("cannot publish staged model CSV: " + outputPath);
  return 0;
}

}  // namespace

namespace SEP {
namespace Validation {
namespace CV01 {

bool RunModel(const std::vector<std::string>& arguments,
              const std::string& outputPath,
              std::string* error) {
  try {
    Run(arguments, outputPath);
    if (error) error->clear();
    return true;
  } catch (const std::exception& exception) {
    if (error) *error = exception.what();
    return false;
  }
}

}  // namespace CV01
}  // namespace Validation
}  // namespace SEP

#ifdef SRCSEP_CV01_STANDALONE_TEST_HARNESS
// This main exists only for dependency-light compilation/sanitizer tests of the
// linked model implementation. The scientific validation runner never defines
// this macro and never executes this binary; it requires the actual srcSEP/AMPS
// executable and reaches RunModel through the production registry callback.
int main(int argc, char** argv) {
  std::vector<std::string> arguments;
  std::string outputPath;
  for (int i = 1; i < argc; ++i) {
    const std::string token = argv[i] ? argv[i] : "";
    if (token == "--output") {
      if (++i >= argc) {
        std::cerr << "CV01 test harness error: --output requires a path\n";
        return 2;
      }
      outputPath = argv[i];
    } else {
      arguments.push_back(token);
    }
  }
  if (outputPath.empty()) {
    std::cerr << "CV01 test harness error: --output is required\n";
    return 2;
  }
  std::string error;
  if (!SEP::Validation::CV01::RunModel(arguments, outputPath, &error)) {
    std::cerr << "CV01 test harness error: " << error << '\n';
    return 2;
  }
  return 0;
}
#endif
