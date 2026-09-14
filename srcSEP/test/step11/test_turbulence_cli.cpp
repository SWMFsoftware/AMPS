#include "sep_cli.h"

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool Parse(std::vector<std::string> arguments, SEP::Util::CLI::Options* options,
           std::string* error) {
  std::vector<char*> argv;
  for (std::size_t i = 0; i < arguments.size(); ++i)
    argv.push_back(&arguments[i][0]);
  std::ostringstream out, err;
  const bool result = SEP::Util::CLI::ParseCommandLine(
      static_cast<int>(argv.size()), argv.data(), *options, out, err);
  *error = err.str();
  return result;
}

}  // namespace

int main() {
  SEP::Util::CLI::Options options;
  std::string error;
  const bool valid = Parse({"amps", "--turbulence-source",
      "self-consistent-spectral", "--turbulence-coupling-policy=disabled",
      "--turbulence-inner-boundary", "specified-incoming-flux",
      "--turbulence-inner-value", "2.5", "--spectral-k-min", "1e-9",
      "--spectral-k-max", "1e-3", "--spectral-bins", "64",
      "--reflection-coefficient", "0.2", "--cascade-coefficient", "0.4",
      "--turbulence-correlation-length", "5e6",
      "--turbulence-conservation-tolerance", "1e-10"}, &options, &error);
  if (!valid || options.turbulence.source !=
                    SEP::Turbulence::Source::SelfConsistentSpectral ||
      options.turbulence.representation != SEP::Turbulence::Representation::Spectral ||
      options.turbulence.coupling != SEP::Turbulence::CouplingPolicy::Disabled ||
      options.turbulence.spectralBins != 64 ||
      options.turbulence.innerBoundary.value != 2.5)
    return EXIT_FAILURE;

  SEP::Util::CLI::Options invalid;
  if (Parse({"amps", "--spectral-k-min", "2", "--spectral-k-max", "1"},
            &invalid, &error))
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
