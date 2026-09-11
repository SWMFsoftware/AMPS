#ifndef SWCME_TEST_FRAMEWORK_HPP
#define SWCME_TEST_FRAMEWORK_HPP

#include <cmath>
#include <iostream>
#include <string>

namespace swcme_test {

class Context {
public:
  void expect_true(bool condition, const std::string& message) {
    if (!condition) {
      ++failures_;
      std::cerr << "    assertion failed: " << message << '\n';
    }
  }

  void expect_near(double actual, double expected, double tolerance,
                   const std::string& message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
      ++failures_;
      std::cerr << "    assertion failed: " << message
                << " (actual=" << actual << ", expected=" << expected
                << ", tolerance=" << tolerance << ")\n";
    }
  }

  int failures() const { return failures_; }

private:
  int failures_ = 0;
};

using TestFunction = void (*)(Context&);

struct TestCase {
  const char* name;
  TestFunction function;
};

}  // namespace swcme_test

#endif
