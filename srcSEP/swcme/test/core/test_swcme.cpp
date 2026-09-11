#include "test_framework.hpp"

#include <exception>
#include <iostream>

void test_1d_ambient_at_one_au(swcme_test::Context& context);

int main() {
  const swcme_test::TestCase tests[] = {
      {"1d.ambient_at_one_au", test_1d_ambient_at_one_au},
  };

  int failed_tests = 0;
  for (const auto& test : tests) {
    std::cout << "[ RUN      ] " << test.name << '\n';
    swcme_test::Context context;

    try {
      test.function(context);
    } catch (const std::exception& error) {
      context.expect_true(false, std::string("unexpected exception: ") + error.what());
    } catch (...) {
      context.expect_true(false, "unexpected non-standard exception");
    }

    if (context.failures() == 0) {
      std::cout << "[       OK ] " << test.name << '\n';
    } else {
      ++failed_tests;
      std::cout << "[  FAILED  ] " << test.name << '\n';
    }
  }

  const int test_count = static_cast<int>(sizeof(tests) / sizeof(tests[0]));
  std::cout << "[==========] " << test_count << " test(s) ran; "
            << failed_tests << " failed.\n";
  return failed_tests == 0 ? 0 : 1;
}
