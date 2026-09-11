#ifndef SWCME_TEST_FRAMEWORK_HPP
#define SWCME_TEST_FRAMEWORK_HPP

#include <cmath>
#include <iostream>
#include <string>

namespace swcme_test {

class Context {
public:
  // Record every subcheck so individual tests can publish structured counts
  // without requiring a JSON/CSV dependency. A skipped subcheck is counted
  // separately and never turns an otherwise valid mandatory test into FAIL.
  void record_result(bool passed) {
    ++checks_;
    if (!passed) {
      ++failures_;
    } else {
      ++passes_;
    }
  }

  void record_skip() {
    ++checks_;
    ++skips_;
  }

  void expect_true(bool condition, const std::string& message) {
    record_result(condition);
    if (!condition) {
      std::cerr << "    assertion failed: " << message << '\n';
    }
  }

  void expect_near(double actual, double expected, double tolerance,
                   const std::string& message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
      record_result(false);
      std::cerr << "    assertion failed: " << message
                << " (actual=" << actual << ", expected=" << expected
                << ", tolerance=" << tolerance << ")\n";
    } else {
      record_result(true);
    }
  }

  int checks() const { return checks_; }
  int passes() const { return passes_; }
  int failures() const { return failures_; }
  int skips() const { return skips_; }

private:
  int checks_ = 0;
  int passes_ = 0;
  int failures_ = 0;
  int skips_ = 0;
};

using TestFunction = void (*)(Context&);

struct TestCase {
  const char* id;
  const char* classification;
  const char* name;
  TestFunction function;
};

}  // namespace swcme_test

#endif
