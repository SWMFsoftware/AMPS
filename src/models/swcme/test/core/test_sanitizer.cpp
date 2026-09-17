#include "test_framework.hpp"

#include <cstdlib>
#include <iostream>

// SAN01 is both a registered validation and the entry point to a clean,
// separately instrumented build.  The ordinary executable invokes the
// Makefile campaign below; the sanitizer-built child sets SWCME_SAN01_CHILD so
// its own SAN01 entry verifies instrumentation and returns instead of recursing
// indefinitely.  All other registered tests continue normally in that child.
void test_san01(swcme_test::Context& context) {
  std::cout << "SAN01 address/undefined-behavior validation\n";
  const bool child=std::getenv("SWCME_SAN01_CHILD")!=nullptr;
  if (child) {
#if defined(SWCME_SAN01_INSTRUMENTED) && defined(__SANITIZE_ADDRESS__)
    context.expect_true(true,
                        "SAN01 child confirms ASan/UBSan instrumented build");
#else
    context.expect_true(false,
                        "SAN01 child was not compiled with required instrumentation");
#endif
    return;
  }

  // Use make rather than reimplementing compiler discovery or source lists in
  // C++.  The target compiles fresh objects with both sanitizers, runs the full
  // deterministic registry (including malformed-input tests), and executes all
  // three demonstration programs under the same runtime.  system() returns a
  // nonzero status for a compiler error, assertion failure, sanitizer report,
  // or demo failure, making every such finding fail this registered test.
  const int status=std::system("make --no-print-directory san01-sanitize");
  context.expect_true(status==0,
                      "SAN01 complete instrumented suite and demos pass");
}
