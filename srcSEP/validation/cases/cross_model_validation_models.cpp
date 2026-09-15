#include "cross_model_validation_models.h"

// Repository-relative includes preserve the enclosing AMPS build contract:
// this nested source must compile without adding srcSEP/util to the include
// search path.
#include "../../util/sep_focused_transport_core.h"
#include "../../util/sep_transport_common.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::Transport::FocusedTransportBackground;
using SEP::Transport::FocusedTransportIncrement;
using SEP::Transport::FocusedTransportState;
using SEP::Transport::KeyedRandomStream;
using SEP::Transport::PitchAngleDiffusionProvider;
using SEP::Transport::PitchAngleDiffusionSample;
using SEP::Transport::Status;

const double ProtonMassKg = 1.67262192369e-27;
const double LightSpeedMPerS = 299792458.0;
const double TwoPi = 6.283185307179586476925286766559;

std::map<std::string, std::string> Parse(
    const std::vector<std::string>& arguments) {
  if (arguments.empty() || arguments.size() % 2 != 0)
    throw std::runtime_error("XM arguments must be non-empty name/value pairs");
  std::map<std::string, std::string> values;
  for (std::size_t i = 0; i < arguments.size(); i += 2) {
    if (arguments[i].size() < 3 || arguments[i].substr(0, 2) != "--")
      throw std::runtime_error("XM option must begin with --");
    if (!values.insert(std::make_pair(arguments[i].substr(2),
                                      arguments[i + 1])).second)
      throw std::runtime_error("duplicate XM option " + arguments[i]);
  }
  return values;
}

std::string Require(const std::map<std::string, std::string>& values,
                    const std::string& name) {
  const std::map<std::string, std::string>::const_iterator found =
      values.find(name);
  if (found == values.end() || found->second.empty())
    throw std::runtime_error("missing XM option --" + name);
  return found->second;
}

double Number(const std::map<std::string, std::string>& values,
              const std::string& name) {
  const std::string text = Require(values, name);
  char* end = NULL;
  errno = 0;
  const double value = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || !end || *end != '\0' || !std::isfinite(value))
    throw std::runtime_error("invalid finite XM number --" + name);
  return value;
}

std::uint64_t Unsigned(const std::map<std::string, std::string>& values,
                       const std::string& name) {
  const std::string text = Require(values, name);
  char* end = NULL;
  errno = 0;
  const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
  if (text.empty() || text[0] == '-' || errno == ERANGE || !end || *end != '\0')
    throw std::runtime_error("invalid unsigned XM integer --" + name);
  return static_cast<std::uint64_t>(value);
}

void Commit(std::ofstream* output, const std::string& temporary,
            const std::string& destination) {
  output->close();
  if (output->fail()) throw std::runtime_error("cannot flush XM output CSV");
  if (std::rename(temporary.c_str(), destination.c_str()) != 0)
    throw std::runtime_error("cannot publish XM output CSV");
}

class IsotropicPitchDiffusion final : public PitchAngleDiffusionProvider {
 public:
  explicit IsotropicPitchDiffusion(double d0PerS) : d0PerS_(d0PerS) {}

  PitchAngleDiffusionSample Evaluate(double, double, double mu) const override {
    // D_mumu=D0(1-mu^2) makes the pitch-angle boundaries natural because the
    // stochastic amplitude vanishes at |mu|=1. Both the coefficient and its
    // Ito drift derivative are expressed in s^-1 as required by the production
    // focused-transport core.
    PitchAngleDiffusionSample sample;
    sample.status = Status::Ok();
    sample.dMuMuPerS = d0PerS_ * std::max(0.0, 1.0 - mu * mu);
    sample.dDmuMuDmuPerS = -2.0 * d0PerS_ * mu;
    sample.provenance = "validation:XM01:D0(1-mu2)";
    sample.turbulenceStateIdentity = "validation:XM01:frozen";
    return sample;
  }

 private:
  double d0PerS_;
};

double Gaussian(KeyedRandomStream* random) {
  // Box-Muller is used only to create the reviewed initial packet. The
  // production mover receives a separate keyed stream, so initialization and
  // transport random draws cannot change one another when the timestep does.
  const double u1 = random->UniformOpen01();
  const double u2 = random->UniformOpen01();
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(TwoPi * u2);
}

double SampleInitialMu(KeyedRandomStream* random) {
  // The normalized angular density is (1+0.6*mu)/2. Rejection sampling from a
  // uniform proposal is exact and leaves a nonzero first moment for scattering
  // and focusing to evolve.
  for (;;) {
    const double mu = 2.0 * random->UniformOpen01() - 1.0;
    if (random->UniformOpen01() <= (1.0 + 0.6 * mu) / 1.6) return mu;
  }
}

struct XM01Scenario {
  const char* name;
  double d0PerS;
  double dLnBdsPerM;
  double divergencePerS;
};

void RunXM01(const std::map<std::string, std::string>& values,
             const std::string& outputPath) {
  const double lengthM = Number(values, "length-m");
  const double speedMPerS = Number(values, "speed-m-per-s");
  const double durationS = Number(values, "duration-s");
  const double timeStepS = Number(values, "time-step-s");
  const double d0PerS = Number(values, "d0-per-s");
  const double dLnBdsPerM = Number(values, "dlnb-ds-per-m");
  const double divergencePerS = Number(values, "divergence-per-s");
  const unsigned sBins = static_cast<unsigned>(Unsigned(values, "s-bins"));
  const unsigned muBins = static_cast<unsigned>(Unsigned(values, "mu-bins"));
  const std::uint64_t baseParticles = Unsigned(values, "particles");
  const std::uint64_t campaignSeed = Unsigned(values, "campaign-seed");
  if (!(lengthM > 0.0 && speedMPerS > 0.0 && durationS > 0.0 &&
        timeStepS > 0.0 && d0PerS >= 0.0) || sBins < 8 || muBins < 8 ||
      baseParticles < 1000)
    throw std::runtime_error("XM01 has an invalid physical/numerical domain");

  const SEP::Transport::ScalarResult momentum =
      SEP::Transport::MomentumFromSpeed(speedMPerS, ProtonMassKg,
                                        LightSpeedMPerS);
  if (!momentum.status.ok()) throw std::runtime_error(momentum.status.message);
  const XM01Scenario scenarios[] = {
      {"streaming", 0.0, 0.0, 0.0},
      {"scattering", d0PerS, 0.0, 0.0},
      {"focusing", 0.0, dLnBdsPerM, 0.0},
      {"adiabatic", 0.0, 0.0, divergencePerS},
      {"combined", d0PerS, dLnBdsPerM, divergencePerS},
  };

  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  output << std::setprecision(17)
         << "scenario,refinement,s_bin,mu_bin,s_left_m,s_right_m,"
            "mu_left,mu_right,probability,intensity,anisotropy,mean_log_p\n";
  for (std::size_t scenarioIndex = 0;
       scenarioIndex < sizeof(scenarios) / sizeof(scenarios[0]);
       ++scenarioIndex) {
    const XM01Scenario& scenario = scenarios[scenarioIndex];
    for (unsigned refinement = 1; refinement <= 2; ++refinement) {
      // Four times as many samples on the second level gives the expected
      // Monte-Carlo standard-error reduction of two. Both levels retain the
      // same physical mesh so refinement measures sampling convergence, not a
      // change in the observable definition.
      const std::uint64_t particleCount = baseParticles *
          (refinement == 1 ? 1U : 4U);
      std::vector<std::uint64_t> histogram(sBins * muBins, 0);
      std::vector<double> logMomentumSum(sBins * muBins, 0.0);
      double muSum = 0.0;
      IsotropicPitchDiffusion diffusion(scenario.d0PerS);
      for (std::uint64_t particle = 0; particle < particleCount; ++particle) {
        KeyedRandomStream initial(campaignSeed, particle + 1,
                                  1000 + scenarioIndex, 0);
        double s = 0.35 * lengthM + 0.06 * lengthM * Gaussian(&initial);
        s -= std::floor(s / lengthM) * lengthM;
        FocusedTransportState state(s, momentum.value, SampleInitialMu(&initial));
        double elapsed = 0.0;
        std::uint64_t step = 0;
        while (elapsed < durationS) {
          const double dt = std::min(timeStepS, durationS - elapsed);
          const FocusedTransportBackground background(
              scenario.dLnBdsPerM, 0.0, 0.0, scenario.divergencePerS);
          KeyedRandomStream transport(campaignSeed, particle + 1,
                                      2000 + scenarioIndex, step);
          const FocusedTransportIncrement increment =
              SEP::Transport::AdvanceFocusedTransportDmumu(
                  state, background, ProtonMassKg, LightSpeedMPerS, dt,
                  diffusion, transport, NULL);
          if (!increment.status.ok())
            throw std::runtime_error(increment.status.message);
          state = increment.state;
          state.arcLengthM -=
              std::floor(state.arcLengthM / lengthM) * lengthM;
          elapsed += dt;
          ++step;
        }
        const unsigned i = std::min(sBins - 1,
            static_cast<unsigned>(state.arcLengthM / lengthM * sBins));
        const unsigned j = std::min(muBins - 1,
            static_cast<unsigned>((state.mu + 1.0) * 0.5 * muBins));
        const std::size_t bin = i * muBins + j;
        ++histogram[bin];
        logMomentumSum[bin] += std::log(state.momentumKgMPerS);
        muSum += state.mu;
      }
      const double anisotropy = 3.0 * muSum / particleCount;
      for (unsigned i = 0; i < sBins; ++i) {
        double intensity = 0.0;
        for (unsigned j = 0; j < muBins; ++j)
          intensity += static_cast<double>(histogram[i * muBins + j]) /
                       particleCount;
        for (unsigned j = 0; j < muBins; ++j) {
          const std::size_t bin = i * muBins + j;
          const double probability =
              static_cast<double>(histogram[bin]) / particleCount;
          const double meanLogP = histogram[bin]
              ? logMomentumSum[bin] / histogram[bin]
              : std::numeric_limits<double>::quiet_NaN();
          output << scenario.name << ',' << refinement << ',' << i << ',' << j
                 << ',' << lengthM * i / sBins << ','
                 << lengthM * (i + 1) / sBins << ','
                 << -1.0 + 2.0 * j / muBins << ','
                 << -1.0 + 2.0 * (j + 1) / muBins << ',' << probability
                 << ',' << intensity << ',' << anisotropy << ',' << meanLogP
                 << '\n';
        }
      }
    }
  }
  Commit(&output, temporary, outputPath);
}

void NormalizeExternalCSV(const std::map<std::string, std::string>& values,
                          const std::string& outputPath,
                          const std::string& caseId) {
  const std::string sourcePath = Require(values, "source-csv");
  std::ifstream source(sourcePath.c_str());
  std::string header;
  if (!source.good() || !std::getline(source, header))
    throw std::runtime_error(caseId + " cannot read --source-csv");
  if (!header.empty() && header[header.size() - 1] == '\r') header.erase(header.size() - 1);
  const std::string required = caseId == "XM02"
      ? "elapsed_hours,series,intensity"
      : "observable,coordinate,value";
  if (header != required)
    throw std::runtime_error(caseId + " source CSV header must be exactly: " + required);

  const std::string temporary = outputPath + ".tmp";
  std::ofstream output(temporary.c_str());
  output << header << '\n';
  std::string row;
  std::size_t rowCount = 0;
  while (std::getline(source, row)) {
    if (!row.empty() && row[row.size() - 1] == '\r') row.erase(row.size() - 1);
    if (row.empty() || std::count(row.begin(), row.end(), ',') != 2)
      throw std::runtime_error(caseId + " source CSV contains a malformed row");
    // All XM02/XM03 comparison tables use a finite numeric coordinate and
    // value. Parse both here so NaN, infinity, locale-dependent text, and
    // truncated rows fail inside the linked application evidence boundary.
    const std::size_t first = row.find(',');
    const std::size_t second = row.find(',', first + 1);
    const std::string coordinate = caseId == "XM02"
        ? row.substr(0, first) : row.substr(first + 1, second - first - 1);
    const std::string series = caseId == "XM03" ? row.substr(0, first) : "";
    const std::string value = row.substr(second + 1);
    char* end = NULL;
    const double x = std::strtod(coordinate.c_str(), &end);
    if (!end || *end != '\0' || !std::isfinite(x))
      throw std::runtime_error(caseId + " source CSV has a non-finite coordinate");
    end = NULL;
    const double y = std::strtod(value.c_str(), &end);
    const bool signedObservable =
        caseId == "XM03" && series == "earth_fluence_spectral_index";
    if (!end || *end != '\0' || !std::isfinite(y) ||
        (!signedObservable && y < 0.0))
      throw std::runtime_error(caseId + " source CSV has an invalid nonnegative value");
    output << row << '\n';
    ++rowCount;
  }
  if (!source.eof() || rowCount < 2)
    throw std::runtime_error(caseId + " source CSV is truncated or has fewer than two rows");
  Commit(&output, temporary, outputPath);
}

}  // namespace

namespace SEP { namespace Validation {

bool RunCrossModelValidationModel(const std::string& caseId,
    const std::vector<std::string>& arguments, const std::string& outputPath,
    std::string* error) {
  try {
    const std::map<std::string, std::string> values = Parse(arguments);
    if (caseId == "XM01") RunXM01(values, outputPath);
    else if (caseId == "XM02" || caseId == "XM03")
      NormalizeExternalCSV(values, outputPath, caseId);
    else throw std::runtime_error("unsupported XM validation case");
    if (error) error->clear();
    return true;
  } catch (const std::exception& exception) {
    if (error) *error = exception.what();
    return false;
  }
}

}}

#ifdef SRCSEP_CROSS_MODEL_STANDALONE_TEST_HARNESS
#include <iostream>
int main(int argc, char** argv) {
  if (argc < 4) return 2;
  std::vector<std::string> arguments;
  for (int i = 3; i < argc; ++i) arguments.push_back(argv[i]);
  std::string error;
  if (!SEP::Validation::RunCrossModelValidationModel(
          argv[1], arguments, argv[2], &error)) {
    std::cerr << error << '\n';
    return 2;
  }
  return 0;
}
#endif
