#include "test_framework.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  SCCMTest::Registry tests;
  SCCMTest::RegisterStage0(&tests);
  SCCMTest::RegisterStage1(&tests);
  SCCMTest::RegisterStage2(&tests);
  SCCMTest::RegisterStage3(&tests);
  SCCMTest::RegisterStage4(&tests);
  SCCMTest::RegisterStage5(&tests);
  SCCMTest::RegisterStage6(&tests);
  SCCMTest::RegisterStage7(&tests);
  SCCMTest::RegisterStage8(&tests);
  SCCMTest::RegisterStage9(&tests);
  SCCMTest::RegisterStage10(&tests);
  SCCMTest::RegisterStage11(&tests);
  SCCMTest::RegisterStage14(&tests);
  SCCMTest::RegisterStage14Drift(&tests);

  std::string selected;
  bool list = false;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--list") {
      list = true;
    } else if (argument == "--test" && index + 1 < argc) {
      selected = argv[++index];
    } else {
      std::cerr << "unknown test-runner argument: " << argument << '\n';
      return 2;
    }
  }
  if (list) {
    for (const auto& test : tests) std::cout << test.first << '\n';
    return 0;
  }
  if (selected.empty()) {
    std::cerr << "--test ID is required\n";
    return 2;
  }
  const auto found = tests.find(selected);
  if (found == tests.end()) {
    std::cerr << "unknown test ID: " << selected << '\n';
    return 2;
  }

  const auto start = std::chrono::steady_clock::now();
  try {
    found->second();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    std::cout << '[' << selected << "] PASS (" << std::fixed
              << std::setprecision(3) << elapsed << "s)\n";
    return 0;
  } catch (const std::exception& error) {
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    std::cerr << '[' << selected << "] FAIL (" << std::fixed
              << std::setprecision(3) << elapsed << "s) " << error.what()
              << '\n';
    return 1;
  }
}
