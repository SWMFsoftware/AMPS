#ifndef SEP_CORONAL_CME_TEST_FRAMEWORK_H
#define SEP_CORONAL_CME_TEST_FRAMEWORK_H

#include "sep_coronal_cme/configuration_parser.h"

#include <functional>
#include <map>
#include <stdexcept>
#include <string>

namespace SCCMTest {

// A failed assertion throws so the runner can report the canonical test ID
// without terminating the complete cumulative suite.  Production kernels do
// not throw; this exception is deliberately confined to test infrastructure.
class Failure : public std::runtime_error {
 public:
  explicit Failure(const std::string& message) : std::runtime_error(message) {}
};

inline void Require(bool condition, const std::string& message) {
  if (!condition) throw Failure(message);
}

using TestFunction = std::function<void()>;
using Registry = std::map<std::string, TestFunction>;

std::string ReadText(const std::string& path);
std::string Fixture();
std::string ReplaceOnce(std::string input, const std::string& oldText,
                        const std::string& newText);
SEP::CoronalCME::VersionedConfiguration ParseGood(const std::string& input);
void RequireRejected(const std::string& input, const std::string& context);
void RegisterStage0(Registry* tests);
void RegisterStage1(Registry* tests);
void RegisterStage2(Registry* tests);

}  // namespace SCCMTest

#endif
