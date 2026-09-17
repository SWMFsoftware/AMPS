#include "test_framework.hpp"

#include <iostream>

void test_thr01(swcme_test::Context& context);

// A minimal entry point lets TSan-capable CI run the scheduler matrix without
// recursively entering the complete sanitizer and coverage orchestrators.
int main() {
  swcme_test::Context context;
  test_thr01(context);
  std::cout << "THR01-TSAN checks=" << context.checks()
            << " failures=" << context.failures() << '\n';
  return context.failures()==0 ? 0 : 1;
}
