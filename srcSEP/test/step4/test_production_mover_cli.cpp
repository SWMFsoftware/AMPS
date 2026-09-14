#include "../../util/sep_cli.h"
#include "../../util/sep_production_mover.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void Check(bool condition, const char* id, const std::string& detail) {
  if (!condition) {
    std::cerr << "FAIL " << id << ": " << detail << '\n';
    ++failures;
  } else {
    std::cout << "PASS " << id << ": " << detail << '\n';
  }
}

bool Parse(const std::vector<std::string>& arguments,
           SEP::Util::CLI::Options& options,
           std::string& output,
           std::string& error) {
  std::vector<std::string> storage(arguments);
  std::vector<char*> argv;
  for (std::size_t i = 0; i < storage.size(); ++i) {
    argv.push_back(&storage[i][0]);
  }
  std::ostringstream out;
  std::ostringstream err;
  const bool ok = SEP::Util::CLI::ParseCommandLine(
      static_cast<int>(argv.size()), argv.data(), options, out, err);
  output = out.str();
  error = err.str();
  return ok;
}

}  // namespace

int main() {
  using SEP::Mover::CoefficientContract;
  using SEP::Mover::ProductionMover;

  // MOVCLI01: every canonical spelling parses in both supported value forms,
  // and discovery exposes exactly those same three stable names.
  SEP::Util::CLI::Options options;
  std::string output, error;
  bool canonical_ok =
      Parse({"sep", "--particle-mover", "parker"}, options, output, error) &&
      options.particleMover == ProductionMover::Parker;
  options = SEP::Util::CLI::Options();
  canonical_ok = canonical_ok &&
      Parse({"sep", "--particle-mover=fte-dmumu"}, options, output, error) &&
      options.particleMover == ProductionMover::FocusedTransportDiffusion;
  options = SEP::Util::CLI::Options();
  canonical_ok = canonical_ok &&
      Parse({"sep", "--mover", "fte-mfp"}, options, output, error) &&
      options.particleMover == ProductionMover::FocusedTransportMeanFreePath;
  options = SEP::Util::CLI::Options();
  canonical_ok = canonical_ok &&
      Parse({"sep", "--list-movers"}, options, output, error) &&
      options.listMovers &&
      !SEP::Util::CLI::IsComponentTestExecutionRequested(options);
  Check(canonical_ok && SEP::Mover::Registry().size() == 3,
        "MOVCLI01", "three canonical movers parse and form the complete registry");

  // MOVCLI02: the Step 14 migration window is closed.  Formerly accepted
  // transition aliases must fail just like every other non-canonical spelling;
  // the migration manifest is now the sole source of replacement guidance.
  const char* former_aliases[] = {
      "parker-dxx", "parker-diffusion", "dxx", "fte", "default",
      "focused-transport", "focused-transport-equation", "legacy-fte",
      "focused-transport-event-driven", "event-driven", "event-driven-fte",
      "fte-event-driven"};
  bool aliases_rejected = true;
  for (std::size_t i = 0;
       i < sizeof(former_aliases) / sizeof(former_aliases[0]); ++i) {
    options = SEP::Util::CLI::Options();
    aliases_rejected = aliases_rejected &&
        !Parse({"sep", "--particle-mover", former_aliases[i]},
               options, output, error);
  }
  Check(aliases_rejected, "MOVCLI02",
        "all retired aliases require an explicit canonical replacement");

  // MOVCLI03: consumers can choose algorithms from coefficient and state
  // capabilities without comparing implementation function addresses.
  const SEP::Mover::MoverCapabilities& parker =
      SEP::Mover::Describe(ProductionMover::Parker).capabilities;
  const SEP::Mover::MoverCapabilities& dmumu =
      SEP::Mover::Describe(ProductionMover::FocusedTransportDiffusion).capabilities;
  const SEP::Mover::MoverCapabilities& mfp =
      SEP::Mover::Describe(ProductionMover::FocusedTransportMeanFreePath).capabilities;
  Check(parker.requiresFieldLineAttachment && !parker.usesPitchAngleState &&
            parker.coefficientContract == CoefficientContract::SpatialDiffusion &&
            dmumu.usesPitchAngleState &&
            dmumu.coefficientContract == CoefficientContract::PitchAngleDiffusion &&
            mfp.usesPitchAngleState &&
            mfp.coefficientContract == CoefficientContract::MeanFreePath &&
            !parker.evolvesWaveStateDirectly &&
            !dmumu.evolvesWaveStateDirectly &&
            !mfp.evolvesWaveStateDirectly,
        "MOVCLI03", "capabilities distinguish representation and coefficient contracts");

  // MOVCLI04: historical movers with different or ambiguous physics are not
  // silently redirected to one of the production algorithms.
  const char* rejected[] = {
      "focused-transport-wave-scattering", "coupled-fte",
      "parker-mean-free-path", "mean-free-path-scattering",
      "tenishev-2005-fl", "boris", "parker3d", "drift"};
  bool all_rejected = true;
  for (std::size_t i = 0; i < sizeof(rejected) / sizeof(rejected[0]); ++i) {
    options = SEP::Util::CLI::Options();
    all_rejected = all_rejected &&
        !Parse({"sep", "--particle-mover", rejected[i]},
               options, output, error);
  }
  Check(all_rejected, "MOVCLI04", "unsupported and ambiguous historical movers are rejected");

  std::ostringstream help;
  SEP::Util::CLI::PrintHelp("sep", help);
  std::ostringstream listing;
  SEP::Mover::PrintProductionMovers(listing);
  Check(help.str().find("Choices: parker, fte-dmumu, fte-mfp") != std::string::npos &&
            listing.str().find("parker") != std::string::npos &&
            listing.str().find("fte-dmumu") != std::string::npos &&
            listing.str().find("fte-mfp") != std::string::npos,
        "MOVCLI01", "help and --list-movers share canonical names");

  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
