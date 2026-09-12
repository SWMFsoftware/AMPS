#include "test_framework.hpp"

#include <iostream>

void test_out08(swcme_test::Context& context);

// The ordinary registry runs OUT08 with every other validation.  This compact
// entry point lets the strict-warning target compile and execute the identical
// writer assertion without pulling unrelated test translation units into the
// deliberately narrow writer/status/demo compiler gate.
int main() {
  swcme_test::Context context;
  test_out08(context);
  if (context.failures()!=0) {
    std::cerr << "OUT08 RESULT: FAIL (" << context.failures()
              << " failed checks)\n";
    return 1;
  }
  std::cout << "OUT08 RESULT: PASS (" << context.passes()
            << " checks)\n";
  return 0;
}
