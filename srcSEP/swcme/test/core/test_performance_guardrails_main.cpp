#include "test_framework.hpp"

#include <iostream>

void test_perf01(swcme_test::Context& context);

int main() {
  swcme_test::Context context;
  test_perf01(context);
  std::cout << "PERF01 checks=" << context.checks()
            << " failures=" << context.failures() << '\n';
  return context.failures()==0 ? 0 : 1;
}
