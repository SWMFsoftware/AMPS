#include "test_framework.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

void test_shk17(swcme_test::Context& context) {
  std::cout << "SHK17 shock-reference regeneration and independence\n";

  // Validation is normally launched from test/, while a developer may invoke
  // the same registry binary from the source root.  Resolve only these two
  // documented layouts and fail visibly if neither contains the audited
  // script; never search PATH for a similarly named unreviewed generator.
  std::filesystem::path verifier =
      std::filesystem::path("reference") / "verify_shock_references.py";
  if (!std::filesystem::is_regular_file(verifier)) {
    verifier = std::filesystem::path("test") / "reference" /
               "verify_shock_references.py";
  }
  context.expect_true(std::filesystem::is_regular_file(verifier),
                      "SHK17 finds the checked-in reference verifier");
  if (!std::filesystem::is_regular_file(verifier)) return;

  // The path is selected from fixed source-tree literals rather than user
  // input.  The Python verifier owns temporary isolation, byte/hash comparison,
  // higher-precision repeats, and dependency analysis; the C++ registry owns
  // propagation of its nonzero exit status into the unified validation report.
  const std::string command = "python3 \"" + verifier.string() + "\"";
  const int status = std::system(command.c_str());
  context.expect_true(status == 0,
                      "SHK17 canonical regeneration and precision audit pass");
}
