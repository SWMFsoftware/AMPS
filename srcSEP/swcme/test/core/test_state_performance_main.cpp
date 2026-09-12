#include "test_framework.hpp"

#include <iostream>

void test_pst08(swcme_test::Context& context);

// Minimal runner for the pinned PST08 optimized build.  The ordinary registry
// still executes the same test, while this target avoids compiling unrelated
// validation units with performance-specific flags and prints an unambiguous
// pass/fail result suitable for CI performance logs.
int main() {
  swcme_test::Context context;
  test_pst08(context);
  if (context.failures()!=0) {
    std::cerr << "PST08 RESULT: FAIL (" << context.failures()
              << " failed checks)\n";
    return 1;
  }
  std::cout << "PST08 RESULT: PASS (" << context.passes()
            << " checks)\n";
  return 0;
}
