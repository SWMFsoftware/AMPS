#ifndef AMPS_CORE_TEST_HARNESS_H
#define AMPS_CORE_TEST_HARNESS_H

// A deliberately small C++ test harness for AMPS core code.  It has no
// GoogleTest, PIC, MPI, or application dependency, so dependency-light core
// utilities and configured AMPS applications can use the same assertions.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace AMPS {
namespace Test {

class Failure : public std::runtime_error {
 public:
  explicit Failure(const std::string& message) : std::runtime_error(message) {}
};

struct TestCase {
  std::string name;
  std::string description;
  std::function<void()> function;
};

class Registry {
 public:
  // A function-local static avoids cross-translation-unit initialization
  // order dependencies: every AutoRegister first constructs this registry,
  // regardless of which test object the linker loads first.
  static Registry& Instance() {
    static Registry registry;
    return registry;
  }

  void Add(const char* name, const char* description,
           const std::function<void()>& function) {
    // Stable IDs are evidence keys.  Treat duplicates as harness corruption
    // instead of executing an ambiguous case and publishing a false receipt.
    for (std::size_t i = 0; i < cases_.size(); ++i) {
      if (cases_[i].name == name) {
        throw Failure(std::string("duplicate test name: ") + name);
      }
    }
    TestCase item;
    item.name = name;
    item.description = description;
    item.function = function;
    cases_.push_back(item);
  }

  const std::vector<TestCase>& Cases() const { return cases_; }

 private:
  std::vector<TestCase> cases_;
};

class AutoRegister {
 public:
  AutoRegister(const char* name, const char* description,
               const std::function<void()>& function) {
    Registry::Instance().Add(name, description, function);
  }
};

inline std::string Location(const char* file, int line) {
  std::ostringstream out;
  out << file << ':' << line << ": ";
  return out.str();
}

template <typename Left, typename Right>
inline void ExpectEqual(const Left& left, const Right& right,
                        const char* left_text, const char* right_text,
                        const char* file, int line) {
  if (!(left == right)) {
    std::ostringstream out;
    out << Location(file, line) << "expected " << left_text << " == " << right_text
        << ", got " << left << " and " << right;
    throw Failure(out.str());
  }
}

inline void ExpectNear(double value, double reference, double tolerance,
                       const char* value_text, const char* reference_text,
                       const char* file, int line) {
  if (!std::isfinite(value) || !std::isfinite(reference) || tolerance < 0.0 ||
      std::fabs(value - reference) > tolerance) {
    std::ostringstream out;
    out << Location(file, line) << "expected |" << value_text << " - "
        << reference_text << "| <= " << tolerance << ", got " << value
        << " and " << reference;
    throw Failure(out.str());
  }
}

enum CaseOrder {
  CASE_ORDER_NORMAL,
  CASE_ORDER_REVERSE,
  CASE_ORDER_SHUFFLE
};

inline bool ParseCaseOrder(const std::string& text, CaseOrder& order) {
  if (text == "normal") {
    order = CASE_ORDER_NORMAL;
  } else if (text == "reverse") {
    order = CASE_ORDER_REVERSE;
  } else if (text == "shuffle") {
    order = CASE_ORDER_SHUFFLE;
  } else {
    return false;
  }
  return true;
}

inline int RunAll(const std::string& filter = std::string(),
                  bool list_only = false,
                  const std::vector<std::string>& exact_cases =
                      std::vector<std::string>(),
                  CaseOrder order = CASE_ORDER_NORMAL,
                  unsigned int shuffle_seed = 0) {
  const std::vector<TestCase>& cases = Registry::Instance().Cases();
  if (list_only) {
    for (std::size_t i = 0; i < cases.size(); ++i) {
      std::printf("%s\t%s\n", cases[i].name.c_str(), cases[i].description.c_str());
    }
    return 0;
  }

  if (!filter.empty() && !exact_cases.empty()) {
    std::fprintf(stderr, "--filter and --case are mutually exclusive\n");
    return 2;
  }

  // Store pointers rather than TestCase copies so std::function ownership and
  // the static registry remain unchanged.  Exact selection starts in CLI
  // order; filter selection starts in registry order.  Reverse and seeded
  // shuffle then operate on that explicit base sequence, making the isolation
  // gate deterministic on every machine and independent of linker order when
  // the caller provides the complete repeated --case list.
  std::vector<const TestCase*> selected_cases;
  if (!exact_cases.empty()) {
    std::set<std::string> requested;
    for (std::size_t requested_index = 0;
         requested_index < exact_cases.size(); ++requested_index) {
      if (!requested.insert(exact_cases[requested_index]).second) {
        std::fprintf(stderr, "duplicate --case '%s'\n",
                     exact_cases[requested_index].c_str());
        return 2;
      }
      const TestCase* match = NULL;
      for (std::size_t registered = 0; registered < cases.size(); ++registered)
        if (cases[registered].name == exact_cases[requested_index]) {
          match = &cases[registered];
          break;
        }
      if (match == NULL) {
        std::fprintf(stderr, "no C++ test has exact ID '%s'\n",
                     exact_cases[requested_index].c_str());
        return 2;
      }
      selected_cases.push_back(match);
    }
  } else {
    for (std::size_t i = 0; i < cases.size(); ++i)
      if (filter.empty() || cases[i].name.find(filter) != std::string::npos)
        selected_cases.push_back(&cases[i]);
  }

  if (order == CASE_ORDER_REVERSE) {
    std::reverse(selected_cases.begin(), selected_cases.end());
  } else if (order == CASE_ORDER_SHUFFLE) {
    std::mt19937 generator(shuffle_seed);
    std::shuffle(selected_cases.begin(), selected_cases.end(), generator);
  }

  const char* order_name = order == CASE_ORDER_REVERSE ? "reverse" :
      (order == CASE_ORDER_SHUFFLE ? "shuffle" : "normal");
  std::printf("[HARNESS] order=%s seed=%u selected=%lu\n", order_name,
              shuffle_seed,
              static_cast<unsigned long>(selected_cases.size()));
  std::fflush(stdout);

  int selected = 0;
  int failed = 0;
  for (std::size_t i = 0; i < selected_cases.size(); ++i) {
    const TestCase& test = *selected_cases[i];
    ++selected;
    // This stable record is intentionally emitted and flushed before calling
    // the test function.  Both suite adapters and the global runner stream it,
    // so long statistical/native cases identify themselves while running
    // rather than only publishing a family-level banner or final result.
    std::printf("[CASE] RUN  %s %s\n", test.name.c_str(),
                test.description.c_str());
    std::fflush(stdout);
    // Each case is isolated at the exception boundary so one failed assertion
    // does not suppress evidence from the remaining cases.  Fatal signals and
    // memory corruption are intentionally not converted into PASS/FAIL; the
    // outer runner reports an abnormal/missing receipt instead.
    try {
      test.function();
      std::printf("[%s] PASS ok\n", test.name.c_str());
    } catch (const std::exception& error) {
      ++failed;
      std::printf("[%s] FAIL %s\n", test.name.c_str(), error.what());
    } catch (...) {
      ++failed;
      std::printf("[%s] FAIL unknown exception\n", test.name.c_str());
    }
  }
  if (selected == 0) {
    std::fprintf(stderr, "no C++ tests matched the requested selection\n");
    return 2;
  }
  std::printf("AMPS core C++ tests: total=%d pass=%d fail=%d\n",
              selected, selected - failed, failed);
  return failed == 0 ? 0 : 1;
}

}  // namespace Test
}  // namespace AMPS

#define AMPS_TEST_JOIN_INNER(a, b) a##b
#define AMPS_TEST_JOIN(a, b) AMPS_TEST_JOIN_INNER(a, b)
#define AMPS_TEST_CASE(NAME, DESCRIPTION)                                      \
  static void AMPS_TEST_JOIN(amps_test_function_, __LINE__)();                 \
  static AMPS::Test::AutoRegister AMPS_TEST_JOIN(amps_test_registration_,      \
                                                  __LINE__)(                    \
      NAME, DESCRIPTION, AMPS_TEST_JOIN(amps_test_function_, __LINE__));       \
  static void AMPS_TEST_JOIN(amps_test_function_, __LINE__)()

#define AMPS_EXPECT_TRUE(EXPRESSION)                                           \
  do {                                                                         \
    if (!(EXPRESSION)) {                                                       \
      throw AMPS::Test::Failure(AMPS::Test::Location(__FILE__, __LINE__) +     \
                                "expected true: " #EXPRESSION);                \
    }                                                                          \
  } while (false)

#define AMPS_EXPECT_FALSE(EXPRESSION) AMPS_EXPECT_TRUE(!(EXPRESSION))

#define AMPS_EXPECT_EQ(LEFT, RIGHT)                                            \
  AMPS::Test::ExpectEqual((LEFT), (RIGHT), #LEFT, #RIGHT, __FILE__, __LINE__)

#define AMPS_EXPECT_NEAR(VALUE, REFERENCE, TOLERANCE)                          \
  AMPS::Test::ExpectNear((VALUE), (REFERENCE), (TOLERANCE), #VALUE,            \
                         #REFERENCE, __FILE__, __LINE__)

#endif
