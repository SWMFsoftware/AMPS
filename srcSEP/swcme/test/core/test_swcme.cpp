#include "test_framework.hpp"

#include <cstring>
#include <exception>
#include <iostream>

void test_cfg01(swcme_test::Context& context);
void test_1d_ambient_at_one_au(swcme_test::Context& context);

int main(int argc, char** argv) {
  const swcme_test::TestCase tests[] = {
      {"CFG01", "COMMON", "Physical constants and unit-conversion consistency",
       test_cfg01},
      {"1D_AMBIENT_01", "1D", "Ambient values at 1 AU",
       test_1d_ambient_at_one_au},
  };

  const char* requested_id = nullptr;
  if (argc == 3 && std::strcmp(argv[1], "--test") == 0) {
    requested_id = argv[2];
  } else if (argc != 1) {
    std::cerr << "usage: " << argv[0] << " [--test TEST_ID]\n";
    return 2;
  }

  int failed_tests = 0;
  int selected_tests = 0;
  for (const auto& test : tests) {
    if (requested_id != nullptr && std::strcmp(requested_id, test.id) != 0) {
      continue;
    }

    ++selected_tests;
    std::cout << "[ RUN      ] " << test.id << " [" << test.classification
              << "] " << test.name << '\n';
    swcme_test::Context context;

    try {
      test.function(context);
    } catch (const std::exception& error) {
      context.expect_true(false, std::string("unexpected exception: ") + error.what());
    } catch (...) {
      context.expect_true(false, "unexpected non-standard exception");
    }

    if (context.failures() == 0) {
      std::cout << "[       OK ] " << test.id << '\n';
    } else {
      ++failed_tests;
      std::cout << "[  FAILED  ] " << test.id << '\n';
    }
  }

  if (selected_tests == 0) {
    std::cerr << "unknown test id: " << requested_id << '\n';
    return 2;
  }

  std::cout << "[==========] " << selected_tests << " test(s) ran; "
            << failed_tests << " failed.\n";
  return failed_tests == 0 ? 0 : 1;
}
