#ifndef STENCIL_TESTS_REGISTER_H
#define STENCIL_TESTS_REGISTER_H

#include <vector>
#include <string>
#include "test_harness.h"

/**
 * REGISTER_STENCIL_TEST(NS, NAME, BLURB)
 *
 * This registers NS::Run(...) and DECLARES NS::ForceLinkAllTests() so
 * the central hook can call it. Each test TU must provide a DEFINITION
 * of NS::ForceLinkAllTests() at file scope (non-static).
 */
#define REGISTER_STENCIL_TEST(NS, NAME, BLURB)                                   \
  namespace NS {                                                                  \
    int Run(const std::vector<std::string>&);                                     \
    void ForceLinkAllTests();   /* declaration only: definition in this TU */     \
  }                                                                               \
  static TestHarness::AutoRegister _auto_reg_##NS(                                \
      NAME,                                                                       \
      NS::Run,                                                                    \
      BLURB)

#endif // STENCIL_TESTS_REGISTER_H

