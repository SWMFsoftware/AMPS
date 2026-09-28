#include "test_framework.h"

#include <fstream>
#include <sstream>

namespace SCCMTest {
namespace {
std::string gFixturePath;
}

std::string ReadText(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  Require(stream.good(), "cannot open test fixture: " + path);
  std::ostringstream bytes;
  bytes << stream.rdbuf();
  return bytes.str();
}

std::string Fixture() {
  if (gFixturePath.empty()) gFixturePath = "test/data/schema5_stage0.in";
  return ReadText(gFixturePath);
}

std::string ReplaceOnce(std::string input, const std::string& oldText,
                        const std::string& newText) {
  const std::size_t position = input.find(oldText);
  Require(position != std::string::npos, "mutation token absent: " + oldText);
  input.replace(position, oldText.size(), newText);
  return input;
}

SEP::CoronalCME::VersionedConfiguration ParseGood(const std::string& input) {
  const auto result = SEP::CoronalCME::ParseConfiguration(input);
  Require(result.ok(), "expected valid input: " + result.status.message);
  return result.value;
}

void RequireRejected(const std::string& input, const std::string& context) {
  const auto result = SEP::CoronalCME::ParseConfiguration(input);
  Require(!result.ok(), "invalid input accepted: " + context);
}

}  // namespace SCCMTest
