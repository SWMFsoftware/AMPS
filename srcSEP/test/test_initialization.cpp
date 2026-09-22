#include "../util/sep_initialization_validation.h"
#include "../util/sep_initialization.h"

#include <iostream>

int main(int argc, char** argv) {
  SEP::Testing::Registry registry(SEP::Testing::InitializationDescriptors());
  const SEP::Testing::Summary summary = registry.Run(
      registry.Select(std::vector<std::string>(),
                      std::vector<std::string>(), true), std::cout);
  if (summary.ExitCode() != 0) return summary.ExitCode();
  if (argc != 2) {
    std::cerr << "expected shipped initialization example path\n";
    return 2;
  }
  SEP::Initialization::Configuration example;
  const SEP::Transport::Status loaded =
      SEP::Initialization::LoadFile(argv[1], &example);
  if (!loaded.ok()) {
    std::cerr << "shipped initialization example failed: "
              << loaded.message << '\n';
    return 1;
  }
  // The command-line directory override must retain both reviewed leaf names;
  // otherwise a preview run could silently change the declared artifact set.
  const std::string meshLeaf = example.meshTecplotFile;
  const std::string lineLeaf = example.fieldLineTecplotFile;
  const std::string dataLeaf = example.dataTecplotFile;
  const SEP::Transport::Status redirected =
      SEP::Initialization::ApplyOutputDirectoryOverride("preview", &example);
  if (!redirected.ok() ||
      example.meshTecplotFile != "preview/" + meshLeaf ||
      example.fieldLineTecplotFile != "preview/" + lineLeaf ||
      example.dataTecplotFile != "preview/" + dataLeaf ||
      example.observers.size() != 2 ||
      example.observers[0].energySpacing !=
          SEP::Initialization::EnergyChannelSpacing::Logarithmic ||
      example.observers[1].energySpacing !=
          SEP::Initialization::EnergyChannelSpacing::Linear) {
    std::cerr << "initialization output-directory override failed: "
              << redirected.message << '\n';
    return 1;
  }
  std::cout << "INIT-EXAMPLE PASS fingerprint="
            << SEP::Initialization::Fingerprint(example) << '\n';
  return 0;
}
