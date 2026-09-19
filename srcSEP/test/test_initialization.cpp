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
  std::cout << "INIT-EXAMPLE PASS fingerprint="
            << SEP::Initialization::Fingerprint(example) << '\n';
  return 0;
}
